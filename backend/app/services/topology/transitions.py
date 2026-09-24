"""
Membrane transitions: what does the chain do between two consecutive SS fragments?

    A (helix/strand) -> gap -> B (helix/strand)

decided from membrane GEOMETRY first (signed depth along the membrane normal),
DSSP/STRIDE and the gap length only refine. State machine:

  1. structural fragments?    A and B >= MIN_FRAGMENT_LEN residues, no chain break
                              in the gap                              else SEPARATE
  2. gap inside the membrane? every gap residue within the envelope   else SEPARATE
  3. A and B compatible?      same sign of travel along the normal and both tilted
                              less than MAX_TILT_DEG from it
        no  -> same side + deep excursion of two PARTIAL fragments -> REENTRANT
               otherwise (hairpin / interfacial B)  -> SEPARATE
  4. A+B = ONE traversal?     ends on opposite sides of the mid-plane, net travel
                              >= min_cross_span, no reversal > REVERSAL_TOL
        no  -> same side + deep -> REENTRANT, else SEPARATE
  5. does A / B cross alone?  span >= full_cross_frac x thickness
        both   -> TWO_TM      (two independent crossings, never merged)
        one    -> AMBIGUOUS   (kept as one continuous crossing, confidence lowered)
        none   -> gap of >= MIN_UNWOUND non-regular residues centred in the CORE
                      -> BROKEN_TM  (one crossing, drawn as parts a / b)
                  otherwise (1-2 residue kink, or break at the interface)
                      -> CONTINUOUS (one crossing, drawn as one helix)

Only CONTINUOUS, AMBIGUOUS and BROKEN_TM merge A and B into one crossing.
Thresholds are algorithmic choices (state them as such in a paper), not biology.
Every decision returns its evidence so the result can be explained.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional

import numpy as np

from app.services.topology.membrane import EDGE_WIDTH, MembraneFrame

MIN_FRAGMENT_LEN = 5        # residues; shorter H/E bits are treated as part of a gap
MAX_TILT_DEG = 60.0         # fragment axis vs membrane normal; flatter = interfacial
REVERSAL_TOL = 5.0          # A: back-tracking allowed along the normal (helix wobble)
MIN_UNWOUND = 3             # non-regular residues for a real break (1-2 = kink)

MERGING = {"CONTINUOUS", "AMBIGUOUS", "BROKEN_TM"}


@dataclass
class TransitionEvidence:
    classification: str
    gap_residues: int
    gap_in_membrane: bool
    gap_in_core: bool
    same_orientation: bool
    ss_class_a: Optional[str]
    ss_class_b: Optional[str]
    same_ss_class: bool
    a_tilt_deg: Optional[float]
    b_tilt_deg: Optional[float]
    single_crossing: bool
    a_full_cross: bool
    b_full_cross: bool
    reason: str

    def as_dict(self) -> dict:
        d = asdict(self)
        for key in ("a_tilt_deg", "b_tilt_deg"):
            if d[key] is not None:
                d[key] = round(d[key], 1)
        return d


def _tilt_deg(coords, membrane: MembraneFrame, a: int, b: int) -> Optional[float]:
    if coords is None or b - a < 3:
        return None
    v = coords[b] - coords[a]
    norm = float(np.linalg.norm(v))
    if not np.isfinite(norm) or norm == 0:
        return None
    cos = abs(float(v @ membrane.normal)) / norm
    return float(np.degrees(np.arccos(min(1.0, cos))))


def _max_reversal(depths: np.ndarray, direction: float) -> float:
    """Largest back-tracking against the overall direction of travel (A)."""
    d = depths * np.sign(direction)
    return float(np.max(np.maximum.accumulate(d) - d)) if d.size else 0.0


def classify_membrane_transition(a: tuple[int, int], b: tuple[int, int], membrane: MembraneFrame,
                                 coords, breaks, min_cross: float,
                                 full_cross: float, ss_class_a: Optional[str] = None,
                                 ss_class_b: Optional[str] = None, ) -> TransitionEvidence:
    """Classify A -> gap -> B (inclusive frame positions, A before B)."""
    same_ss_class = (
        ss_class_a is None
        or ss_class_b is None
        or ss_class_a == ss_class_b
    )
    a0, a1 = a
    b0, b1 = b  
    depth = membrane.depth
    half = membrane.half_thickness
    gap = list(range(a1 + 1, b0))
    core_limit = half - EDGE_WIDTH

    gap_in_membrane = all(membrane.in_envelope(k) for k in gap)
    if gap:
        gap_center = float(np.nanmean(depth[gap]))
    else:
        gap_center = float((depth[a1] + depth[b0]) / 2.0)
    gap_in_core = abs(gap_center) <= core_limit

    du, dv = float(depth[a1] - depth[a0]), float(depth[b1] - depth[b0])
    ta, tb = _tilt_deg(coords, membrane, a0, a1), _tilt_deg(coords, membrane, b0, b1)
    same_sign = du * dv > 0
    tilt_ok = all(t is None or t <= MAX_TILT_DEG for t in (ta, tb))
    same_orientation = same_sign and tilt_ok

    path = depth[a0:b1 + 1]
    start, end = float(depth[a0]), float(depth[b1])
    opposite_sides = start * end < 0
    net = abs(end - start)
    reversal = _max_reversal(path, end - start)
    single_crossing = bool(opposite_sides and net >= min_cross and reversal <= REVERSAL_TOL)
    deep = bool(np.nanmin(np.abs(path)) <= core_limit)
    same_side = start * end > 0
    a_full = membrane.span(a0, a1) >= full_cross
    b_full = membrane.span(b0, b1) >= full_cross

    def result(cls, reason):
        return TransitionEvidence(
            cls,
            len(gap),
            gap_in_membrane,
            gap_in_core,
            same_orientation,
            ss_class_a,
            ss_class_b,
            same_ss_class,
            ta,
            tb,
            single_crossing,
            a_full,
            b_full,
            reason,
        )
    if min(a1 - a0, b1 - b0) + 1 < MIN_FRAGMENT_LEN:
        return result("SEPARATE", "fragment shorter than MIN_FRAGMENT_LEN")
    if breaks is not None and any(breaks[k] for k in range(a1 + 1, b0 + 1)):
        return result("SEPARATE", "chain break (unresolved residues) between the fragments")
    if not gap_in_membrane:
        return result("SEPARATE", "the chain leaves the membrane between the fragments")
    partial_pair = not a_full and not b_full
    if not same_orientation:
        if same_side and deep and partial_pair:
            return result("REENTRANT", "in and back out on the same side, reaching the core")
        if not same_sign:
            return result("SEPARATE", "antiparallel across the membrane (hairpin)")
        return result("SEPARATE", "a fragment lies almost parallel to the membrane (interfacial)")
    if not same_ss_class:
        return result(
            "SEPARATE",
            f"different SS classes: {ss_class_a} -> {ss_class_b}"
        )
    if not single_crossing:
        if same_side and deep and partial_pair:
            return result("REENTRANT", "A+B enter and leave on the same side")
        return result("SEPARATE", "A+B do not make exactly one traversal of the bilayer")
    if a_full and b_full:
        return result("TWO_TM", "A and B each cross the bilayer on their own")
    if a_full or b_full:
        return result("AMBIGUOUS", "one fragment crosses alone, the other only partly")
    if len(gap) >= MIN_UNWOUND and gap_in_core:
        return result("BROKEN_TM", "same-class partial fragments form one membrane traversal with an unwound gap in the core")
    if len(gap) >= MIN_UNWOUND:
        return result("CONTINUOUS", "break lies at the interface, not in the core")
    return result("CONTINUOUS", "1-2 residue kink inside one traversal")
