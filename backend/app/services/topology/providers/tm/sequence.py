from pathlib import Path
from typing import Optional

import numpy as np

from app.core.constants import (
    KYTE_DOOLITTLE, SMOOTH_WINDOW, TM_EDGE_THRESHOLD, TM_EDGE_WINDOW, TM_HYDRO_THRESHOLD,
    TM_MAX_LENGTH, TM_MERGE_GAP, TM_SEGMENT_MIN_LENGTH,
)
from app.schemas.topology import TMParams
from app.services.topology.labels import TM, UNASSIGNED, build_regions
from app.services.topology.providers.tm.base import TMPrediction, TMProvider
from app.services.topology.residues import ResidueFrame, load_residue_frame

# modified residues seen in PDB files -> parent residue (otherwise they score 0)
_ALIASES = {"MSE": "MET", "SEC": "CYS", "SEP": "SER", "TPO": "THR", "PTR": "TYR",
            "MLY": "LYS", "KCX": "LYS", "HYP": "PRO", "CSO": "CYS", "CME": "CYS"}


def _kd_scores(sequence) -> np.ndarray:
    return np.array([KYTE_DOOLITTLE.get(res, KYTE_DOOLITTLE.get(_ALIASES.get(res, ""), 0.0))
                     for res in sequence], dtype=float)


def _windowed_mean(scores: np.ndarray, window: int, fragments) -> np.ndarray:
    """Centred sliding mean that never reaches across a fragment boundary (the window is
    truncated at fragment ends)."""
    n = len(scores)
    out = np.empty(n)
    half = window // 2
    for f0, f1 in fragments:
        seg = scores[f0:f1 + 1]
        c = np.concatenate([[0.0], np.cumsum(seg)])
        i = np.arange(len(seg))
        lo = np.maximum(0, i - half)
        hi = np.minimum(len(seg), i + half + 1)
        out[f0:f1 + 1] = (c[hi] - c[lo]) / (hi - lo)
    return out


def _smooth_hydrophobicity(sequence, window: int = SMOOTH_WINDOW, fragments=None) -> np.ndarray:
    """Windowed Kyte-Doolittle mean. With ``fragments`` (inclusive position ranges of
    covalently continuous stretches) the window never reaches across a chain break,
    so residues on both sides of an unresolved loop are not averaged together."""
    scores = _kd_scores(sequence)
    n = len(scores)
    if n == 0:
        return scores
    return _windowed_mean(scores, window, fragments or [(0, n - 1)])


def _split_long(s: int, e: int, edge_profile: np.ndarray, min_length: int,
                max_length: int) -> list[tuple[int, int]]:
    """A hydrophobic run longer than one helix can hold is (almost always) two helices
    joined by a short loop that the edge profile did not resolve: cut it at its least
    hydrophobic point, keeping both halves >= min_length, recursively."""
    if e - s + 1 <= max_length:
        return [(s, e)]
    lo, hi = s + min_length, e - min_length
    if hi < lo:
        return [(s, e)]
    k = lo + int(np.argmin(edge_profile[lo:hi + 1]))
    return (_split_long(s, k - 1, edge_profile, min_length, max_length)
            + _split_long(k + 1, e, edge_profile, min_length, max_length))


def _find_tm_segments(sequence, fragments=None,
                      threshold: float = TM_HYDRO_THRESHOLD, window: int = SMOOTH_WINDOW,
                      edge_threshold: float = TM_EDGE_THRESHOLD,
                      edge_window: int = TM_EDGE_WINDOW,
                      min_length: int = TM_SEGMENT_MIN_LENGTH,
                      max_length: int = TM_MAX_LENGTH) -> list[tuple[int, int]]:
    """TM segments (inclusive positions) by hysteresis on two Kyte-Doolittle profiles.

    * DETECTION - the classic Kyte-Doolittle criterion: a 19-residue window whose mean
      exceeds 1.6 means a TM helix is present somewhere under that window.
    * EXTENT - taken from a short-window profile (``edge_window``): the segment is the
      maximal run where the short profile stays >= ``edge_threshold`` and that contains
      at least one detection position. The 19-residue profile cannot place the ends
      (it is still above 1.6 while its window already covers several polar flank
      residues) and it does not dip between two helices joined by a short loop.
    * Runs longer than ``max_length`` are split at their least hydrophobic point
      (`_split_long`); pieces shorter than ``min_length`` are dropped.

    Nothing is merged or extended across a fragment boundary.
    """
    scores = _kd_scores(sequence)
    n = len(scores)
    if n == 0:
        return []
    fragments = fragments or [(0, n - 1)]
    detect = _windowed_mean(scores, window, fragments)
    edge = _windowed_mean(scores, edge_window, fragments)
    segments = []
    for f0, f1 in fragments:
        i = f0
        while i <= f1:
            if edge[i] < edge_threshold:
                i += 1
                continue
            j = i
            while j + 1 <= f1 and edge[j + 1] >= edge_threshold:
                j += 1
            if np.any(detect[i:j + 1] > threshold):
                for s, e in _split_long(i, j, edge, min_length, max_length):
                    if e - s + 1 >= min_length:
                        segments.append((s, e))
            i = j + 1
    return segments


def _find_hydrophobic_segments(hydro, threshold, min_length, merge_gap=TM_MERGE_GAP,
                               window: int = SMOOTH_WINDOW):
    """LEGACY - no longer used by SequenceTMProvider (kept for other callers).

    Positions (inclusive) of windows whose mean exceeds ``threshold``, widened by half a
    window. Known problems: the widening overshoots each end by 2-5 residues (the mean
    passes the threshold while the window still covers polar flanks), helices joined by
    loops shorter than ~12 residues are merged, and ``min_length`` never filters anything
    (every widened run is >= ``window`` long). Use `_find_tm_segments`."""
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
    the orchestrator (alternation + positive-inside rule).

    Known limitation: hydropathy cannot see beta-barrel strands (they alternate
    hydrophobic/polar residues), so outer-membrane proteins get no TM segments."""

    LABELER = "Kyte-Doolittle_Sequence"

    def predict_tm(self, file_path: Path, params: Optional[TMParams] = None,
                   frame: Optional[ResidueFrame] = None, **kwargs) -> TMPrediction:
        frame = frame if frame is not None else load_residue_frame(file_path, kwargs.get("chain_id"))
        n = len(frame)
        if n == 0:
            return TMPrediction(labeler=self.LABELER, warnings=["no amino-acid residues in chain"])

        segments = _find_tm_segments(frame.names, fragments=frame.fragments())

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