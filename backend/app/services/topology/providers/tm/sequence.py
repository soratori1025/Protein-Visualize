from pathlib import Path
from typing import Optional

import numpy as np

from app.core.constants import (
    KYTE_DOOLITTLE, SMOOTH_WINDOW, TM_HYDRO_THRESHOLD, TM_MERGE_GAP, TM_MIN_LENGTH,
)
from app.schemas.topology import TMParams
from app.services.topology.labels import TM, UNASSIGNED, build_regions
from app.services.topology.providers.tm.base import TMPrediction, TMProvider
from app.services.topology.residues import ResidueFrame, load_residue_frame


def _smooth_hydrophobicity(sequence, window: int = SMOOTH_WINDOW, fragments=None) -> np.ndarray:
    """Windowed Kyte-Doolittle mean. With ``fragments`` (inclusive position ranges of
    covalently continuous stretches) the window never reaches across a chain break,
    so residues on both sides of an unresolved loop are not averaged together."""
    scores = np.array([KYTE_DOOLITTLE.get(res, 0.0) for res in sequence], dtype=float)
    n = len(scores)
    if n == 0:
        return scores
    fragments = fragments or [(0, n - 1)]
    half = window // 2
    out = np.empty(n)
    for f0, f1 in fragments:
        for i in range(f0, f1 + 1):
            lo, hi = max(f0, i - half), min(f1, i + half)
            out[i] = scores[lo:hi + 1].mean()
    return out


def _find_hydrophobic_segments(hydro, threshold, min_length, merge_gap=TM_MERGE_GAP,
                               window: int = SMOOTH_WINDOW):
    """Positions (inclusive) of windows whose mean exceeds ``threshold``, widened by
    half a window (the window centre marks the middle of the segment)."""
    raw_segments = []
    in_seg, start = False, 0
    for i, h in enumerate(hydro):
        if h > threshold and not in_seg:
            in_seg, start = True, i
        elif h <= threshold and in_seg:
            in_seg = False
            raw_segments.append((start, i - 1))
    if in_seg:
        raw_segments.append((start, len(hydro) - 1))

    merged = []
    half = window // 2
    for s0, e0 in raw_segments:
        s = max(0, s0 - half)
        e = min(len(hydro) - 1, e0 + half)
        if merged and s - merged[-1][1] - 1 <= merge_gap:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return [(s, e) for s, e in merged if (e - s + 1) >= min_length]


class SequenceTMProvider(TMProvider):
    """Kyte-Doolittle hydropathy scan on the OBSERVED residues of the chain.

    Reports TM segments only; loop sides are left 'Unassigned' and are inferred by
    the orchestrator (alternation + positive-inside rule)."""

    LABELER = "Kyte-Doolittle_Sequence"

    def predict_tm(self, file_path: Path, params: Optional[TMParams] = None,
                   frame: Optional[ResidueFrame] = None, **kwargs) -> TMPrediction:
        frame = frame if frame is not None else load_residue_frame(file_path, kwargs.get("chain_id"))
        n = len(frame)
        if n == 0:
            return TMPrediction(labeler=self.LABELER, warnings=["no amino-acid residues in chain"])

        hydro = _smooth_hydrophobicity(frame.names, fragments=frame.fragments())
        segments = _find_hydrophobic_segments(hydro, TM_HYDRO_THRESHOLD, TM_MIN_LENGTH)

        labels = [UNASSIGNED] * n
        for s, e in segments:
            for k in range(s, e + 1):
                labels[k] = TM
        groups = [-1] * n
        for g, (s, e) in enumerate(segments):
            for k in range(s, e + 1):
                groups[k] = g
        return TMPrediction(
            regions=build_regions(frame, labels, groups, skip=(UNASSIGNED,)),
            labels=labels,
            segments=segments,
            labeler=self.LABELER,
        )
