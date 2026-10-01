"""
TMHMM 2.0 transmembrane-helix provider (sequence-based hidden Markov model).

Model : TMHMM 2.0 - Sonnhammer, von Heijne & Krogh, ISMB 1998; Krogh, Larsson,
        von Heijne & Sonnhammer, J Mol Biol 305:567-580 (2001).
Engine: the vendored package ``_tmhmm`` - model parser taken verbatim from pyTMHMM
        1.3.6 (MIT) and a numpy port of its Viterbi / forward-backward code, so no C
        compiler is needed (``pip install pyTMHMM`` has no Windows wheel and does not
        build against NumPy 2). The port reproduces pyTMHMM's Viterbi path exactly and
        fixes two upstream bugs; see ``_tmhmm/hmm.py`` and scripts/check_tmhmm_port.py.

Which sequence is scored
------------------------
An HMM needs the CONTIGUOUS chain: gluing the resolved fragments of a structure
together would create hydrophobic stretches that do not exist and destroy the loop
lengths the model relies on. So the provider runs TMHMM on the full chain sequence
the file declares and maps the result onto the resolved residues:

  1. mmCIF ``_pdbx_poly_seq_scheme`` (full sequence + author numbering per residue),
  2. mmCIF ``_entity_poly.pdbx_seq_one_letter_code_can`` / PDB ``SEQRES``
     (positions mapped by sequence alignment),
  3. otherwise the resolved residues, every gap in the residue numbering filled with
     the unknown residue X (emission marginalised), so the loop keeps its length.

Output: per-position labels TM / Cytoplasmic (TMHMM "i") / Extracellular (TMHMM
"o", i.e. non-cytoplasmic) - TMHMM predicts the topology, not only the helices.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Optional

from fastapi import HTTPException

from app.schemas.topology import TMParams
from app.services.topology.labels import CYTO, EXTRA, TM, UNASSIGNED, build_regions
from app.services.topology.providers.tm import _tmhmm
from app.services.topology.providers.tm.base import TMPrediction, TMProvider
from app.services.topology.residues import (
    ResidueFrame, load_residue_frame, map_to_reference, one_letter,
)

MAX_GAP_FILL = 2000          # residue-number gaps above this are not filled (renumbered files)
_LABEL = {"M": TM, "i": CYTO, "o": EXTRA}


def _as_list(value):
    if value is None:
        return []
    return list(value) if isinstance(value, (list, tuple)) else [value]


def _blank(value) -> str:
    v = (value or "").strip()
    return "" if v in (".", "?") else v


def _int_or_none(value):
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


# =============================================================================
# Full chain sequence declared by the file
# =============================================================================
def _cif_sequence(file_path: Path, chain_id: str):
    """-> (sequence, keys or None, source) from mmCIF, or None."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    d = MMCIF2Dict(str(file_path))
    strands = _as_list(d.get("_pdbx_poly_seq_scheme.pdb_strand_id"))
    if strands:
        mons = _as_list(d.get("_pdbx_poly_seq_scheme.mon_id"))
        nums = _as_list(d.get("_pdbx_poly_seq_scheme.pdb_seq_num"))
        icodes = _as_list(d.get("_pdbx_poly_seq_scheme.pdb_ins_code")) or [""] * len(strands)
        seq, keys = [], []
        for strand, mon, num, ic in zip(strands, mons, nums, icodes):
            if strand != chain_id:
                continue
            seq.append(one_letter(mon))
            n = _int_or_none(num)
            keys.append((n, _blank(ic)) if n is not None else (None, ""))
        if seq:
            return "".join(seq), keys, "_pdbx_poly_seq_scheme"
    owners = _as_list(d.get("_entity_poly.pdbx_strand_id"))
    codes = _as_list(d.get("_entity_poly.pdbx_seq_one_letter_code_can"))
    for owner, code in zip(owners, codes):
        if chain_id in [c.strip() for c in owner.split(",")]:
            seq = "".join(code.split()).upper()
            if seq and seq not in ("?", "."):
                return seq, None, "_entity_poly"
    return None


def _pdb_sequence(file_path: Path, chain_id: str):
    names = []
    with open(file_path, errors="replace") as fh:
        for line in fh:
            rec = line[:6]
            if rec == "SEQRES" and line[11:12].strip() == (chain_id or "").strip():
                names.extend(line[19:].split())
            elif rec in ("ATOM  ", "HETATM"):
                break                       # header is over
    if not names:
        return None
    return "".join(one_letter(n) for n in names), None, "SEQRES"


def declared_sequence(file_path: Optional[Path], chain_id: str):
    """Full sequence of ``chain_id`` from the file header, or None."""
    if file_path is None:
        return None
    file_path = Path(file_path)
    try:
        if file_path.suffix.lower() in (".cif", ".mmcif"):
            return _cif_sequence(file_path, chain_id)
        return _pdb_sequence(file_path, chain_id)
    except Exception as err:                   # noqa: BLE001 - optional header data
        print(f"[tmhmm] could not read the chain sequence from {file_path.name}: {err}")
        return None


def gap_filled_sequence(frame: ResidueFrame) -> tuple[str, list[int]]:
    """Resolved residues with every residue-number gap filled with X.
    -> (sequence, index of each frame position in that sequence)."""
    seq, index = [], []
    prev = None
    for r in frame.residues:
        if prev is not None:
            gap = r.resseq - prev - 1
            if 0 < gap <= MAX_GAP_FILL:
                seq.extend("X" * gap)
        index.append(len(seq))
        seq.append(r.one)
        prev = r.resseq
    return "".join(seq), index


@lru_cache(maxsize=64)
def _run_tmhmm(sequence: str) -> "_tmhmm.TMHMMResult":
    try:
        model = _tmhmm.load_model()
    except FileNotFoundError as err:            # server set-up problem, not a bad request
        raise HTTPException(status_code=503, detail=str(err)) from err
    return _tmhmm.predict(sequence, model)


# =============================================================================
# Provider
# =============================================================================
class TMHMMProvider(TMProvider):
    LABELER = "TMHMM_2.0"

    def predict_tm(self, file_path: Path, params: Optional[TMParams] = None,
                   frame: Optional[ResidueFrame] = None, **kwargs) -> TMPrediction:
        frame = frame if frame is not None else load_residue_frame(file_path, kwargs.get("chain_id"))
        n = len(frame)
        if n == 0:
            return TMPrediction(labeler=self.LABELER, warnings=["no amino-acid residues in chain"])
        warns: list[str] = []

        # ---- sequence to score + frame position -> sequence index -------------
        ref_index = None
        seq = ""
        declared = declared_sequence(file_path or frame.source_path, frame.chain_id)
        if declared is not None:
            seq, keys, source = declared
            if keys is None or any(k[0] is None for k in keys):
                keys = [(j + 1, "") for j in range(len(seq))]
            mapping = map_to_reference(frame, keys, list(seq), offsets=(0,), source=source,
                                       allow_order_fallback=False, note_remap=False)
            if mapping.matched >= 0.9 * n:
                ref_index = mapping.ref_index
                warns.extend(mapping.warnings)
            else:
                warns.append(f"TMHMM: the {source} sequence does not match the resolved residues "
                             f"of chain {frame.chain_id} (identity {mapping.identity or 0:.0%}); "
                             "scored the resolved residues instead")
        if ref_index is None:
            seq, ref_index = gap_filled_sequence(frame)
            if declared is None:
                warns.append("TMHMM: the file has no SEQRES / entity sequence - scored the "
                             "resolved residues, gaps in the numbering as unknown residues (X)")

        result = _run_tmhmm(seq)
        helices = result.helices
        if not helices:
            return TMPrediction(
                labels=[UNASSIGNED] * n, segments=[], labeler=self.LABELER,
                warnings=warns + [f"TMHMM predicts no transmembrane helix (expected residues "
                                  f"in TM helices {result.exp_aa_tmh:.1f})"])

        # ---- labels and segments on frame positions ---------------------------
        labels = [_LABEL.get(result.path[j], UNASSIGNED) if j is not None else UNASSIGNED
                  for j in ref_index]
        pos_in_seq: dict[int, int] = {j: p for p, j in enumerate(ref_index) if j is not None}
        segments: list[tuple[int, int]] = []
        unresolved = 0
        for a, b in helices:
            hits = [pos_in_seq[j] for j in range(a, b + 1) if j in pos_in_seq]
            if hits:
                segments.append((min(hits), max(hits)))
            else:
                unresolved += 1
        groups = [-1] * n
        for g, (s, e) in enumerate(segments):
            for k in range(s, e + 1):
                labels[k] = TM          # a residue the alignment left unpaired inside a helix
                groups[k] = g

        if unresolved:
            warns.append(f"TMHMM: {unresolved} of {len(helices)} predicted TM helices lie in "
                         "residues not resolved in this structure")
        if result.possible_signal_peptide:
            warns.append(f"TMHMM: possible N-terminal signal peptide (expected TM residues in "
                         f"the first 60 = {result.exp_first60:.1f} > 10) - an N-terminal "
                         "'helix' may be a signal peptide, check before reading it as a TM helix")
        return TMPrediction(
            regions=build_regions(frame, labels, groups, skip=(UNASSIGNED,)),
            labels=labels,
            segments=segments,
            labeler=self.LABELER,
            warnings=warns,
        )
