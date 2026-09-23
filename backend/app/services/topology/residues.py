"""
Canonical residue frame shared by every TM and SS provider.

Why this module exists
----------------------
Each block of the topology pipeline used to index residues its own way:

  * TM providers returned author residue numbers (int only: no chain, no insertion
    code) and each provider parsed the file itself with its own residue filter
    (the UniProt provider even kept HETATM ions whose residue name is "CA");
  * DSSP / STRIDE returned whatever numbering the tool printed - author numbers,
    sequential numbers, label_seq_id for mmCIF, "100A"-style strings, all chains;
  * the orchestrator rebuilt residue lists with ``range(start, end + 1)``, inventing
    residues that are not in the model and grouping by list adjacency.

That mismatch is NOT a constant offset (gaps, insertion codes, chains, tags and
different residue filters all break any single shift), so "add an offset before
merging" cannot fix it. Instead the structure is parsed ONCE into a ResidueFrame -
the ordered list of observed amino-acid residues of one chain - and every block
reports on frame POSITIONS 0..n-1:

  * TM providers return one label per position,
  * SS output is mapped onto positions by residue key (chain, resseq, icode) and
    VERIFIED against the residue identity; if the numbering does not agree the
    mapping falls back to a sequence alignment (see ``map_to_reference``),
  * author numbers reappear only when regions are emitted at the very end.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

import numpy as np

from app.core.constants import KYTE_DOOLITTLE

# CA-CA distance above which two consecutive residues are not covalently bonded
# (trans 3.8 A, cis ~2.9 A). Used to split SS elements / smoothing windows at gaps.
CA_BREAK_DISTANCE = 4.2

# Numbering is trusted only when residue identity confirms it.
MIN_NUMBERING_IDENTITY = 0.95
MIN_NUMBERING_COVERAGE = 0.80
# An alignment below this identity means "wrong sequence / wrong chain".
MIN_ALIGNMENT_IDENTITY = 0.70
# Above this many DP cells the numpy aligner hands over to Bio.Align.
MAX_ALIGN_CELLS = 20_000_000   # 3 uint8 traceback matrices -> ~60 MB

ONE_LETTER = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Modified / force-field residue names -> parent amino acid. MSE (selenomethionine)
# alone is in a large share of crystal structures; HSD/HIE/CYX come from MD output.
PARENT_RESIDUE = {
    "MSE": "MET", "SEP": "SER", "TPO": "THR", "PTR": "TYR", "HYP": "PRO", "MLY": "LYS",
    "M3L": "LYS", "ALY": "LYS", "KCX": "LYS", "LLP": "LYS", "CSO": "CYS", "CSD": "CYS",
    "OCS": "CYS", "CME": "CYS", "SEC": "CYS", "PYL": "LYS", "HIC": "HIS", "MLE": "LEU",
    "NLE": "LEU", "CGU": "GLU", "PCA": "GLU", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "CYX": "CYS", "CYM": "CYS", "ASH": "ASP",
    "GLH": "GLU", "LYN": "LYS",
}


def canonical_resname(raw: str) -> str:
    name = (raw or "").strip().upper()
    return PARENT_RESIDUE.get(name, name)


def one_letter(resname: str) -> str:
    return ONE_LETTER.get(canonical_resname(resname), "X")


# =============================================================================
# The frame
# =============================================================================
@dataclass(frozen=True)
class Residue:
    chain: str
    resseq: int
    icode: str          # "" when blank
    name: str           # canonical 3-letter name (MSE -> MET)
    one: str            # one-letter code

    @property
    def key(self) -> tuple[int, str]:
        return self.resseq, self.icode

    @property
    def label(self) -> str:
        return f"{self.resseq}{self.icode}"


class ResidueFrame:
    """Ordered observed residues of ONE chain; positions 0..n-1 are the shared index."""

    def __init__(self, chain_id: str, residues: Sequence[Residue],
                 coords: Optional[np.ndarray] = None):
        self.chain_id = chain_id
        self.residues: list[Residue] = list(residues)
        n = len(self.residues)
        if coords is None:
            coords = np.full((n, 3), np.nan)
        coords = np.asarray(coords, dtype=float).reshape(-1, 3) if n else np.empty((0, 3))
        if len(coords) != n:
            raise ValueError(f"coords has {len(coords)} rows for {n} residues")
        self.coords = coords
        self._pos: dict[tuple[int, str], int] = {}
        for i, r in enumerate(self.residues):
            self._pos.setdefault(r.key, i)

    # -- construction helpers -------------------------------------------------
    @classmethod
    def from_items(cls, chain_id: str, items: Iterable[tuple[int, str, str]],
                   coords: Optional[np.ndarray] = None) -> "ResidueFrame":
        """items: (resseq, icode, resname3). Handy for tests and non-PDB inputs."""
        residues = [Residue(chain_id, int(rs), (ic or "").strip(), canonical_resname(nm),
                            one_letter(nm)) for rs, ic, nm in items]
        return cls(chain_id, residues, coords)

    # -- basic accessors ------------------------------------------------------
    def __len__(self) -> int:
        return len(self.residues)

    @property
    def sequence(self) -> str:
        return "".join(r.one for r in self.residues)

    @property
    def names(self) -> list[str]:
        return [r.name for r in self.residues]

    def position(self, resseq: int, icode: str = "") -> Optional[int]:
        return self._pos.get((int(resseq), (icode or "").strip()))

    def residues_data(self) -> list[dict]:
        """Legacy [{'id': resseq, 'name': resname}] view used by older helpers."""
        return [{"id": r.resseq, "name": r.name, "icode": r.icode} for r in self.residues]

    def span_to_positions(self, start: int, end: int, start_icode: Optional[str] = None,
                          end_icode: Optional[str] = None) -> Optional[tuple[int, int]]:
        """Author residue range -> inclusive position range, using ONLY observed residues.

        Without an insertion code the range covers every insertion of the boundary
        residue (100, 100A, 100B ...). Returns None when no observed residue lies in
        the range (the region is outside the model and must not be painted)."""
        lo = hi = None
        for i, r in enumerate(self.residues):
            if not (start <= r.resseq <= end):
                continue
            if start_icode is not None and r.resseq == start and r.icode < start_icode:
                continue
            if end_icode is not None and r.resseq == end and r.icode > end_icode:
                continue
            lo = i if lo is None else lo
            hi = i
        return None if lo is None else (lo, hi)

    def chain_breaks(self) -> np.ndarray:
        """breaks[i] is True when residue i is not bonded to residue i-1 (breaks[0] True)."""
        n = len(self)
        breaks = np.ones(n, dtype=bool)
        if n < 2:
            return breaks
        if np.isfinite(self.coords).all():
            dist = np.linalg.norm(np.diff(self.coords, axis=0), axis=1)
            breaks[1:] = dist > CA_BREAK_DISTANCE
        else:  # no coordinates: fall back to numbering (insertions are contiguous)
            rs = np.array([r.resseq for r in self.residues])
            breaks[1:] = np.diff(rs) > 1
        return breaks

    def fragments(self) -> list[tuple[int, int]]:
        """Inclusive position ranges of covalently continuous stretches."""
        starts = list(np.flatnonzero(self.chain_breaks())) + [len(self)]
        return [(int(a), int(b) - 1) for a, b in zip(starts[:-1], starts[1:])]


def _pick_chain(model, chain_id: Optional[str]):
    if chain_id:
        try:
            return model[chain_id]
        except KeyError:
            available = [c.id for c in model]
            raise ValueError(f"Chain '{chain_id}' not found; available chains: {available}") from None
    for chain in model:  # first chain that actually contains amino acids
        for residue in chain:
            if "CA" in residue and canonical_resname(residue.get_resname()) in KYTE_DOOLITTLE:
                return chain
    chain = next(iter(model), None)
    if chain is None:
        raise ValueError("structure has no chains")
    return chain


def load_residue_frame(file_path: Path, chain_id: Optional[str] = None) -> ResidueFrame:
    """Parse the structure ONCE. First model; the requested chain or the first chain
    that contains amino acids. Keeps standard residues plus known modified residues
    (MSE...) that have a CA atom; drops waters, ligands and ions (e.g. Ca2+, whose
    residue and atom are both called "CA")."""
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning

    file_path = Path(file_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        is_cif = file_path.suffix.lower() in (".cif", ".mmcif")
        parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(file_path))

    model = next(iter(structure), None)
    if model is None:
        raise ValueError(f"{file_path.name}: no model/coordinates found")
    chain = _pick_chain(model, chain_id)

    residues, coords = [], []
    for residue in chain:
        hetflag, resseq, icode = residue.get_id()
        raw = residue.get_resname().strip().upper()
        name = canonical_resname(raw)
        if hetflag == "W" or name not in KYTE_DOOLITTLE or "CA" not in residue:
            continue
        if hetflag.strip() and raw not in PARENT_RESIDUE:
            continue  # free amino acid ligand (HETATM) - not part of the chain
        residues.append(Residue(chain.id, int(resseq), (icode or "").strip(), name, one_letter(name)))
        coords.append(residue["CA"].get_coord())
    return ResidueFrame(chain.id, residues, np.array(coords, dtype=float) if coords else None)


# =============================================================================
# Sequence alignment (semi-global, linear gaps) - numpy, no Biopython needed
# =============================================================================
_NEG = -(10 ** 8)


def align_sequences(a: str, b: str, match: int = 2, mismatch: int = -1,
                    gap_open: int = -6, gap_extend: int = -1) -> list[tuple[int, int]]:
    """Semi-global alignment with AFFINE gaps (Gotoh), vectorised per row.

    End gaps are free on both sequences (expression tags, truncated constructs);
    an internal gap costs ``gap_open`` for its first residue and ``gap_extend`` for
    each further one (unresolved loops). Affine costs matter: with a linear cost a
    gap can be split in two for free, which silently moves a residue next to an
    unresolved loop onto an identical residue a few positions away. Bio.Align's
    defaults (all gap scores 0) are worse still - a longest common subsequence.

    Returns aligned index pairs (i in a, j in b), mismatches included; 'X' scores 0."""
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return []
    if n * m > MAX_ALIGN_CELLS:
        return _align_with_biopython(a, b, match, mismatch, gap_open, gap_extend)

    A = np.frombuffer(a.encode("ascii", "replace"), dtype=np.uint8)
    B = np.frombuffer(b.encode("ascii", "replace"), dtype=np.uint8)
    XCH = ord("X")
    b_is_x = B == XCH
    idx = np.arange(m, dtype=np.int64)

    # states: 0 = M (a_i ~ b_j), 1 = X (a_i ~ gap, vertical), 2 = Y (b_j ~ gap, horizontal)
    pM = np.full((n + 1, m + 1), 0, dtype=np.uint8)
    pX = np.full((n + 1, m + 1), 0, dtype=np.uint8)
    pY = np.full((n + 1, m + 1), 0, dtype=np.uint8)

    M = np.full(m + 1, _NEG, dtype=np.int64)
    X = np.full(m + 1, _NEG, dtype=np.int64)
    Y = np.zeros(m + 1, dtype=np.int64)      # free leading b residues
    M[0], Y[0] = 0, _NEG
    Y[1:] = 0
    last_col = [(max(M[m], X[m], Y[m]), int(np.argmax([M[m], X[m], Y[m]])))]

    for i in range(1, n + 1):
        s = np.where(B == A[i - 1], match, mismatch).astype(np.int64)
        if A[i - 1] == XCH:
            s[:] = 0
        s[b_is_x] = 0

        stack = np.stack([M[:-1], X[:-1], Y[:-1]])
        pM[i, 1:] = np.argmax(stack, axis=0)
        newM = np.full(m + 1, _NEG, dtype=np.int64)
        newM[1:] = stack.max(axis=0) + s

        vstack = np.stack([M + gap_open, X + gap_extend, Y + gap_open])
        pX[i] = np.argmax(vstack, axis=0)
        newX = vstack.max(axis=0)
        newX[0] = 0                                  # free leading a residues
        pX[i, 0] = 1

        # horizontal: Y[j] = max_{k<j} (max(M[k], X[k]) + open + (j-1-k)*extend)
        opener = np.maximum(newM[:-1], newX[:-1]) + gap_open          # k = 0..m-1
        best_start = np.maximum.accumulate(opener - idx * gap_extend)
        newY = np.full(m + 1, _NEG, dtype=np.int64)
        newY[1:] = best_start + idx * gap_extend
        # pointer: extend if Y[j-1] + extend == Y[j] (j >= 2), else open from M/X at j-1
        ext = np.zeros(m + 1, dtype=bool)
        ext[2:] = newY[1:-1] + gap_extend == newY[2:]
        from_x = np.zeros(m + 1, dtype=bool)
        from_x[1:] = newX[:-1] > newM[:-1]
        pY[i] = np.where(ext, 2, np.where(from_x, 1, 0))

        M, X, Y = newM, newX, newY
        trio = [M[m], X[m], Y[m]]
        last_col.append((max(trio), int(np.argmax(trio))))

    # free trailing gaps: best cell on the last row or last column
    row_best = np.stack([M, X, Y])
    j_best = int(np.argmax(row_best.max(axis=0)))
    row_score = int(row_best[:, j_best].max())
    i_best = int(np.argmax([v for v, _ in last_col]))
    if row_score >= last_col[i_best][0]:
        i, j, state = n, j_best, int(np.argmax(row_best[:, j_best]))
    else:
        i, j, state = i_best, m, last_col[i_best][1]

    pairs = []
    while i > 0 and j > 0:
        if state == 0:
            pairs.append((i - 1, j - 1))
            state = int(pM[i, j])
            i, j = i - 1, j - 1
        elif state == 1:
            state = int(pX[i, j])
            i -= 1
        else:
            state = int(pY[i, j])
            j -= 1
    pairs.reverse()
    return pairs


def _align_with_biopython(a, b, match, mismatch, gap_open, gap_extend):
    from Bio import Align

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = match, mismatch
    aligner.open_gap_score, aligner.extend_gap_score = gap_open, gap_extend
    aligner.end_gap_score = 0
    alignment = next(iter(aligner.align(a, b)), None)   # never len()/bool(): can overflow
    if alignment is None:
        return []
    pairs = []
    for (a0, a1), (b0, b1) in zip(*alignment.aligned):
        pairs.extend((a0 + k, b0 + k) for k in range(a1 - a0))
    return pairs


# =============================================================================
# Mapping any per-residue reference onto the frame
# =============================================================================
@dataclass
class ReferenceMapping:
    """ref_index[p] = index into the reference list for frame position p (or None)."""
    ref_index: list[Optional[int]]
    method: str
    matched: int
    identity: Optional[float]
    warnings: list[str] = field(default_factory=list)

    @property
    def coverage(self) -> float:
        return self.matched / len(self.ref_index) if self.ref_index else 0.0


def _identity(frame: ResidueFrame, ref_aa: Sequence[Optional[str]], ref_index) -> Optional[float]:
    same = total = 0
    for p, j in enumerate(ref_index):
        if j is None:
            continue
        x, y = frame.residues[p].one, ref_aa[j]
        if not y or x == "X" or y == "X":
            continue
        total += 1
        same += x == y
    return same / total if total else None


def map_to_reference(frame: ResidueFrame, ref_keys: Sequence[tuple[int, str]],
                     ref_aa: Optional[Sequence[Optional[str]]] = None,
                     offsets: Sequence[int] = (0,), source: str = "reference",
                     allow_order_fallback: bool = True, note_remap: bool = True) -> ReferenceMapping:
    """Map frame positions to a reference list of (resseq, icode) keys.

    1. Numbering first: key (resseq + k, icode) for each candidate offset k (0 first).
       Accepted only when residue identity confirms it (if names are available).
    2. Otherwise a semi-global sequence alignment (needs residue names).
    3. Without residue names: numbering if it covers the chain, else by order when
       the counts agree (reported as a warning - it cannot be verified)."""
    n, m = len(frame), len(ref_keys)
    have_aa = ref_aa is not None and any(ref_aa)
    warns: list[str] = []
    if n == 0 or m == 0:
        return ReferenceMapping([None] * n, "none", 0, None, [f"{source}: nothing to map"])

    first_idx: dict[tuple[int, str], int] = {}
    duplicates = 0
    for j, key in enumerate(ref_keys):
        if key in first_idx:
            duplicates += 1
        else:
            first_idx[key] = j
    if duplicates:
        warns.append(f"{source}: {duplicates} duplicated residue numbers (several chains?) - "
                     "first occurrence used")

    need = min(n, m)
    best: Optional[ReferenceMapping] = None
    for k in dict.fromkeys(offsets):          # unique, order preserved
        idx = [first_idx.get((r.resseq + k, r.icode)) for r in frame.residues]
        matched = sum(j is not None for j in idx)
        ident = _identity(frame, ref_aa, idx) if have_aa else None
        method = "residue numbering" if k == 0 else f"residue numbering, offset {k:+d}"
        cand = ReferenceMapping(idx, method, matched, ident, list(warns))
        if have_aa and ident is not None:
            # identity-verified numbering (duplicates -> let the alignment pick the chain)
            if (matched >= MIN_NUMBERING_COVERAGE * need and ident >= MIN_NUMBERING_IDENTITY
                    and not duplicates):
                return cand
        elif matched == need:
            # unverifiable: accept only a complete match of the smaller side
            return cand
        if best is None or matched > best.matched:
            best = cand

    if have_aa:
        ref_seq = "".join((ref_aa[j] or "X") for j in range(m))
        pairs = align_sequences(frame.sequence, ref_seq)
        idx = [None] * n
        for p, j in pairs:
            idx[p] = j
        matched = len(pairs)
        ident = _identity(frame, ref_aa, idx)
        if ident is not None and ident >= MIN_ALIGNMENT_IDENTITY:
            deltas = {ref_keys[j][0] - frame.residues[p].resseq for p, j in pairs}
            if len(deltas) == 1:
                d = deltas.pop()
                method = "sequence alignment (numbering offset {:+d})".format(d) if d else "sequence alignment"
            else:
                method = "sequence alignment"
            if note_remap:
                if best is not None and best.matched:
                    warns.append(f"{source}: residue numbers did not match the structure "
                                 f"(identity {best.identity or 0:.0%}); re-mapped by {method}")
                else:
                    warns.append(f"{source}: residue numbers absent or unmatched; mapped by {method}")
            return ReferenceMapping(idx, method, matched, ident, warns)
        warns.append(f"{source}: sequence does not match chain {frame.chain_id} "
                     f"(best identity {ident or 0:.0%})")
        return ReferenceMapping([None] * n, "unmapped", 0, ident, warns)

    # No residue names: numbering cannot be verified.
    if allow_order_fallback and m == n:
        warns.append(f"{source}: residue numbers do not match; mapped by ORDER (same count, "
                     "unverified - expose residue names in the tool output to verify)")
        return ReferenceMapping(list(range(n)), "residue order", n, None, warns)
    if best is not None:
        warns.append(f"{source}: only {best.matched}/{n} residues matched by number")
        best.warnings = warns
        return best
    return ReferenceMapping([None] * n, "unmapped", 0, None, warns)


# =============================================================================
# Secondary-structure records
# =============================================================================
_RESNUM_RE = re.compile(r"^\s*(-?\d+)\s*([A-Za-z]?)\s*$")
_CHAIN_KEYS = ("chain", "chain_id", "chain_name", "auth_asym_id", "chainID", "asym_id")
_ICODE_KEYS = ("insertion_code", "icode", "ins_code", "pdbx_PDB_ins_code", "insertion")
_AA_KEYS = ("amino_acid", "aa", "residue_name", "resname", "res_name", "residue_type",
            "one_letter_code")
_NUM_KEYS = ("residue_number", "resnum", "residue_id", "res_id", "number", "auth_seq_id",
             "seq_id")
_CODE_KEYS = ("code", "ss", "secondary_structure", "sse", "ss_code")
_SS_WORDS = {"helix": "H", "alpha": "H", "alpha helix": "H", "3-10 helix": "G", "pi helix": "I",
             "strand": "E", "sheet": "E", "beta": "E", "beta strand": "E", "bridge": "B",
             "turn": "T", "bend": "S", "coil": "C", "loop": "C"}


def parse_resnum(value: Any) -> Optional[tuple[int, str]]:
    """100 -> (100, ''), '100A' -> (100, 'A'), ' -3 ' -> (-3, ''), (' ', 100, 'A') -> (100, 'A')."""
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, (int, np.integer)):
        return int(value), ""
    if isinstance(value, float):
        return (int(value), "") if value.is_integer() else None
    if isinstance(value, (tuple, list)):
        for pos, v in enumerate(value):
            if isinstance(v, (int, np.integer)) and not isinstance(v, bool):
                rest = [x for x in value[pos + 1:] if isinstance(x, str)]
                return int(v), (rest[0].strip() if rest else "")
        return None
    match = _RESNUM_RE.match(str(value))
    return (int(match.group(1)), match.group(2)) if match else None


def normalize_ss_code(code: Any) -> str:
    """Single upper-case DSSP-style code; blank / '-' -> 'C'; STRIDE 'b' -> 'B'."""
    if code is None:
        return "C"
    text = str(code).strip()
    if not text or text[0] in "-_.":
        return "C"
    if len(text) > 1 and text.lower() in _SS_WORDS:
        return _SS_WORDS[text.lower()]
    return text[0].upper()


def _normalize_aa(value: Any) -> tuple[Optional[str], bool]:
    """-> (one-letter or None, is_break_marker)."""
    if value is None:
        return None, False
    text = str(value).strip()
    if not text:
        return None, False
    if len(text) == 1:
        if text in "!*":
            return None, True                   # DSSP chain-break pseudo residue
        return ("C" if text.islower() else text.upper()), False   # DSSP: a..z = bonded Cys
    if len(text) == 3:
        return one_letter(text), False
    return None, False


def _first(record: Mapping, keys) -> Any:
    for key in keys:
        if key in record and record[key] not in (None, ""):
            return record[key]
    return None


@dataclass(frozen=True)
class SSRecord:
    resseq: int
    icode: str
    chain: Optional[str]
    code: str
    aa: Optional[str]


def parse_ss_records(raw: Iterable[Mapping]) -> tuple[list[SSRecord], int]:
    """Tolerant parser for DSSP/STRIDE ``to_dict()['residues']`` entries.
    Returns (records, n_skipped)."""
    out, skipped = [], 0
    for rec in raw or []:
        if not isinstance(rec, Mapping):
            skipped += 1
            continue
        aa, is_break = _normalize_aa(_first(rec, _AA_KEYS))
        num = parse_resnum(_first(rec, _NUM_KEYS))
        if is_break or num is None:
            skipped += 1
            continue
        resseq, icode = num
        extra_icode = _first(rec, _ICODE_KEYS)
        if not icode and extra_icode is not None:
            icode = str(extra_icode).strip()
        chain = _first(rec, _CHAIN_KEYS)
        out.append(SSRecord(resseq, icode, None if chain is None else str(chain).strip(),
                            normalize_ss_code(_first(rec, _CODE_KEYS)), aa))
    return out, skipped


@dataclass
class SSAssignment:
    """Per-position SS codes on a ResidueFrame (None = the tool said nothing)."""
    codes: list[Optional[str]]
    labeler: str
    method: str
    matched: int
    identity: Optional[float] = None
    warnings: list[str] = field(default_factory=list)

    @property
    def coverage(self) -> float:
        return self.matched / len(self.codes) if self.codes else 0.0


def map_ss_records(frame: ResidueFrame, records: Sequence[SSRecord],
                   labeler: str = "SS") -> SSAssignment:
    """Put DSSP/STRIDE output on the frame: filter to the frame's chain, then map by
    residue key verified by residue identity, falling back to sequence alignment."""
    warns: list[str] = []
    pool = list(records)
    if any(r.chain is not None for r in pool):
        same_chain = [r for r in pool if r.chain == frame.chain_id]
        if same_chain:
            pool = same_chain
        else:
            chains = sorted({r.chain for r in pool if r.chain is not None})
            warns.append(f"{labeler}: output has chains {chains}, not '{frame.chain_id}' - "
                         "matching by sequence")
    mapping = map_to_reference(
        frame, [(r.resseq, r.icode) for r in pool], [r.aa for r in pool],
        source=labeler)
    codes = [pool[j].code if j is not None else None for j in mapping.ref_index]
    return SSAssignment(codes, labeler, mapping.method, mapping.matched, mapping.identity,
                        warns + mapping.warnings)
