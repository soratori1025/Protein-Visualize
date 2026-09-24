"""
Per-position label vocabulary and helpers shared by the orchestrator, the TM
providers and the legacy predictor functions.

All helpers work on lists indexed by ResidueFrame POSITION, never on residue IDs.
"""
from __future__ import annotations

from typing import Optional, Sequence

from app.core.constants import HELIX_CODES, MIN_EXTRA_SS_LEN, POSITIVE_RESIDUES, STRAND_CODES
from app.schemas.topology import TopologyRegion

TM = "Transmembrane"
INTRA = "Intramembrane"
SIGNAL = "Signal"
CYTO = "Cytoplasmic"
EXTRA = "Extracellular"
SIDE_A = "Side_A"
SIDE_B = "Side_B"
UNASSIGNED = "Unassigned"
KNOWN_SIDES = (CYTO, EXTRA)
GEOMETRIC_SIDES = (SIDE_A, SIDE_B)

TM_CLASS_LABEL = {"H": "Transmembrane Alpha Helix", "E": "Transmembrane Beta Strand",
                  # in the membrane, but not one regular helix/strand (mixed elements,
                  # no SS assignment, or SS and membrane evidence disagree)
                  "I": "Transmembrane Irregular",
                  "L": "Transmembrane Loop", None: TM}
INTERFACIAL = "Interfacial Helix"
# non-helical stretch INSIDE the bilayer between the two halves of a broken crossing
UNWOUND = "Transmembrane Unwound"

# UniProt "Topological domain" vocabulary -> which face it is topologically equivalent to.
# (Lumenal / periplasmic / intermembrane spaces are non-cytoplasmic = "outside".)
_INSIDE_WORDS = ("cytoplasm", "cytosol", "mitochondrial matrix", "stromal", "intravirion")
_OUTSIDE_WORDS = ("extracellular", "lumenal", "luminal", "periplasm", "intermembrane",
                  "vacuolar", "vesicular", "perinuclear", "exoplasmic", "virion surface",
                  "peroxisomal", "thylakoid")


def is_tm(label: Optional[str]) -> bool:
    return bool(label) and label.startswith(TM)


def side_from_text(text: Optional[str]) -> Optional[str]:
    t = (text or "").lower()
    if not t:
        return None
    if any(w in t for w in _OUTSIDE_WORDS):
        return EXTRA
    if any(w in t for w in _INSIDE_WORDS) or t.strip() == "nuclear":
        return CYTO
    return None


def normalize_base_label(label: Optional[str]) -> str:
    """Collapse any provider label to the base vocabulary used by the flows."""
    if not label:
        return UNASSIGNED
    if TM in label:
        return TM
    if INTRA in label:
        return INTRA
    if label.startswith(SIGNAL):
        return SIGNAL
    if label in (CYTO, EXTRA, SIDE_A, SIDE_B, UNASSIGNED):
        return label
    return side_from_text(label) or UNASSIGNED


def coarse_ss(code: Optional[str]) -> Optional[str]:
    """DSSP/STRIDE code -> 'H' | 'E' | 'C'; None stays None (no information)."""
    if code is None:
        return None
    if code in HELIX_CODES:
        return "H"
    if code in STRAND_CODES:
        return "E"
    return "C"


def ss_word(coarse: Optional[str]) -> str:
    return {"H": "Helix", "E": "Strand"}.get(coarse, "Coil")


def apply_positive_inside_rule(labels: Sequence[str], names: Sequence[str]) -> list[str]:
    """Rename Side_A/Side_B to Cytoplasmic/Extracellular (Lys/Arg-richer face = inside)."""
    labels = list(labels)
    if SIDE_A not in labels and SIDE_B not in labels:
        return labels

    def density(side):
        idx = [i for i, c in enumerate(labels) if c == side]
        if not idx:
            return 0.0
        return sum(1 for i in idx if names[i] in POSITIVE_RESIDUES) / len(idx)

    cyto = SIDE_A if density(SIDE_A) >= density(SIDE_B) else SIDE_B
    rename = {cyto: CYTO, (SIDE_B if cyto == SIDE_A else SIDE_A): EXTRA}
    return [rename.get(c, c) for c in labels]


def label_extramembrane_ss(labels: Sequence[str], coarse: Sequence[Optional[str]],
                           min_len: int = MIN_EXTRA_SS_LEN,
                           breaks: Optional[Sequence[bool]] = None) -> list[str]:
    """'Cytoplasmic' -> 'Cytoplasmic Helix' etc. (also Intramembrane), so loops can draw
    their own elements. Helix/strand runs shorter than ``min_len`` become Coil. Runs
    never continue across a side change or a chain break."""
    n = len(labels)
    out = list(labels)
    raw: list[Optional[str]] = [None] * n
    for i, lab in enumerate(labels):
        if lab in KNOWN_SIDES or lab == INTRA:
            raw[i] = ss_word(coarse[i])
    i = 0
    while i < n:
        if raw[i] in ("Helix", "Strand"):
            j = i + 1
            while (j < n and raw[j] == raw[i] and labels[j] == labels[i]
                   and not (breaks is not None and breaks[j])):
                j += 1
            if j - i < min_len:
                for k in range(i, j):
                    raw[k] = "Coil"
            i = j
        else:
            i += 1
    for i in range(n):
        if raw[i] is not None:
            out[i] = f"{out[i]} {raw[i]}"
    return out


def region_fields(desc: str) -> tuple[str, Optional[str], Optional[str]]:
    """description -> (type, side, ss) so the frontend needs no string parsing."""
    if desc.startswith(TM):
        rtype, side = TM, "membrane"
    elif desc.startswith(INTRA):
        rtype, side = INTRA, "membrane"
    elif desc.startswith(SIGNAL):
        rtype, side = SIGNAL, None
    else:
        rtype = "Topological domain"
        side = CYTO if desc.startswith(CYTO) else EXTRA if desc.startswith(EXTRA) else None
    ss = None
    for word in ("Irregular", "Helix", "Strand", "Loop", "Coil"):
        if word in desc:
            ss = word
            break
    return rtype, side, ss


def build_regions(frame, labels: Sequence[str], group_ids: Optional[Sequence[int]] = None,
                  skip: Sequence[str] = ()) -> list[TopologyRegion]:
    """Group consecutive POSITIONS into regions. A new region starts whenever the label
    changes OR the group id changes, so two adjacent TM spans never fuse. Start/end
    are author residue numbers (+ insertion codes when present)."""
    n = len(labels)
    if n == 0:
        return []
    if len(frame) != n:
        raise ValueError(f"{n} labels for a frame of {len(frame)} residues")
    groups = list(group_ids) if group_ids is not None else [-1] * n
    regions = []
    start = 0
    for i in range(1, n + 1):
        if i == n or labels[i] != labels[start] or groups[i] != groups[start]:
            desc = labels[start]
            if desc not in skip:
                first, last = frame.residues[start], frame.residues[i - 1]
                rtype, side, ss = region_fields(desc)
                regions.append(TopologyRegion(
                    type=rtype, start=first.resseq, end=last.resseq, description=desc,
                    side=side, ss=ss, start_icode=first.icode or None,
                    end_icode=last.icode or None))
            start = i
    return regions
