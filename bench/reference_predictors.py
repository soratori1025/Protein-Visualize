"""
Reference predictors implemented for the comparison, not for the application.

The topology benchmark of the manuscript compares the providers of ProteinLab with
each other.  A reader is entitled to ask how they compare with *published* methods
that the authors did not write.  This module provides two such references, both
re-implementations of published rules and both labelled as such wherever they are
reported:

* **TopPred-style** - the original sliding-window hydropathy rule of von Heijne (1992):
  a transmembrane segment is predicted when two or more *consecutive* 19-residue
  windows of the Kyte-Doolittle scale are both favourable; the segment then runs from
  the first to the last of those windows.  The consecutive-window requirement is what
  distinguishes it from a single-peak threshold and is the reason it does not produce
  the isolated hydrophobic stretches that a naive threshold does.

* **TOPCONS-style consensus** - the consensus rule of Fukasaka and Tomii (2016):
  the predictions of several methods are combined by majority vote, after a
  dissimilarity filter removes a method whose prediction agrees with the others on
  less than half of its own segments.  Here the constituents are the TopPred-style
  predictor, TMHMM 2.0 and the Kyte-Doolittle scan, which are the three predictors
  that can be run without an external service.  **This is not TOPCONS**; it is the
  same consensus rule applied to a different set of constituents, and it is reported
  as "a TOPCONS-style consensus" throughout.

Neither predictor is part of the application.  They exist so that the manuscript can
compare the method with something outside it.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from app.core.constants import KYTE_DOOLITTLE  # noqa: E402
from app.services.topology.providers.tm.tmhmm import TMHMMProvider  # noqa: E402

# von Heijne (1992) TopPred: window and threshold of the original method.
TOPPRED_WINDOW = 19
TOPPRED_THRESHOLD = 1.6
# A consensus segment must be at least this long to be called (the same floor the
# application uses for a single provider).
CONSENSUS_MIN_LENGTH = 13
# The dissimilarity filter of TOPCONS: a method is dropped when fewer than this
# fraction of its segments survive in the majority vote.
CONSENSUS_SIMILARITY = 0.5

_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def hydropathy_profile(sequence: str, window: int = TOPPRED_WINDOW) -> list[float]:
    """Mean Kyte-Doolittle hydropathy over a sliding window (0 outside the chain).

    The application's scale is keyed by the three-letter residue name, so the
    one-letter sequence is converted first.
    """
    values = [KYTE_DOOLITTLE.get(_ONE_TO_THREE.get(c, ""), 0.0) for c in sequence]
    n = len(values)
    profile = []
    for i in range(n):
        lo = max(0, i - window // 2)
        hi = min(n, lo + window)
        lo = max(0, hi - window)
        profile.append(sum(values[lo:hi]) / max(1, hi - lo))
    return profile


def toppred_segments(sequence: str) -> list[tuple[int, int]]:
    """The TopPred rule: two consecutive favourable windows make a segment.

    Runs of favourable windows that are separated by less than one window are
    merged, which is what the original method does: a single dip in the profile
    does not split a transmembrane helix into two.
    """
    profile = hydropathy_profile(sequence)
    favourable = [p >= TOPPRED_THRESHOLD for p in profile]
    n = len(favourable)

    # 1. maximal runs of favourable windows
    runs = []
    i = 0
    while i < n:
        if not favourable[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and favourable[j + 1]:
            j += 1
        runs.append((i, j))
        i = j + 1

    # 2. a run of a single window is not a transmembrane segment
    runs = [(a, b) for a, b in runs if b - a + 1 >= 2]
    if not runs:
        return []

    # 3. merge runs separated by less than one window
    merged = [list(runs[0])]
    for start, end in runs[1:]:
        if start - merged[-1][1] <= TOPPRED_WINDOW:
            merged[-1][1] = end
        else:
            merged.append([start, end])

    # 4. the segment runs from the first to the last of its windows
    return [(max(0, a - TOPPRED_WINDOW // 2), min(n - 1, b + TOPPRED_WINDOW // 2))
            for a, b in merged]


def tmhmm_segments(file_path: Path, sequence: str,
                    chain_id: str | None = None) -> list[tuple[int, int]]:
    """TMHMM 2.0 helices, in the frame of the analysed chain."""
    try:
        prediction = TMHMMProvider().predict_tm(file_path, chain_id=chain_id)
    except Exception:
        return []
    return [(int(a), int(b)) for a, b in (prediction.segments or [])]


def _segments_from_mask(mask: Sequence[bool], min_length: int) -> list[tuple[int, int]]:
    out = []
    start = None
    for i, flag in enumerate(mask):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            if i - start >= min_length:
                out.append((start, i - 1))
            start = None
    if start is not None and len(mask) - start >= min_length:
        out.append((start, len(mask) - 1))
    return out


def _segments_agree(a: tuple[int, int], b: tuple[int, int], fraction: float = 0.5) -> bool:
    """Do two predicted segments overlap by at least `fraction` of the smaller one?"""
    lo = max(a[0], b[0])
    hi = min(a[1], b[1])
    if hi < lo:
        return False
    overlap = hi - lo + 1
    smaller = min(a[1] - a[0] + 1, b[1] - b[0] + 1)
    return overlap >= fraction * smaller


def consensus_segments(sequence: str, providers: dict[str, list[tuple[int, int]]],
                        n: int, min_length: int = CONSENSUS_MIN_LENGTH,
                        similarity: float = CONSENSUS_SIMILARITY) -> list[tuple[int, int]]:
    """The TOPCONS consensus rule: dissimilarity filter, then majority vote.

    `providers` maps a method name to its predicted segments in frame positions.
    """
    if not providers:
        return []
    clean = {name: [(int(s), int(e)) for s, e in segs
                    if int(s) <= int(e) and int(s) < n] for name, segs in providers.items()}

    # --- dissimilarity filter: drop a method that agrees with the others on less
    #     than `similarity` of its own segments
    keep = set(clean)
    for name in list(clean):
        others = [segs for other, segs in clean.items() if other != name and other in keep]
        if not others:
            continue
        mine = clean[name]
        if not mine:
            keep.discard(name)
            continue
        agreed = sum(1 for s in mine
                      if any(_segments_agree(s, other_segment)
                             for other in others for other_segment in other))
        if agreed < similarity * len(mine):
            keep.discard(name)

    # --- majority vote over the surviving methods
    votes = [0] * n
    for name in keep:
        for start, end in clean[name]:
            for i in range(max(0, start), min(n, end + 1)):
                votes[i] += 1
    needed = len(keep) // 2 + 1
    mask = [votes[i] >= needed for i in range(n)]
    return _segments_from_mask(mask, min_length)


def reference_predictors(sequence: str, file_path: Path, n: int,
                         hydropathy_segments: Sequence[tuple[int, int]] | None = None,
                         chain_id: str | None = None) -> dict[str, list]:
    """All reference predictors for one protein, in frame positions.

    `hydropathy_segments` is the application's own sequence scan; it is passed in
    because it lives on the residue frame, which only the caller has.  `chain_id` is
    the chain that is being analysed - the TMHMM reference has to be computed on the
    same chain, or its numbering is a different one.
    """
    toppred = toppred_segments(sequence)
    tmhmm = tmhmm_segments(file_path, sequence, chain_id)
    hydropathy = list(hydropathy_segments or [])
    providers = {"toppred": toppred, "tmhmm": tmhmm, "hydropathy": hydropathy}
    return {
        "toppred": toppred,
        "tmhmm": tmhmm,
        "topcons_style_consensus": consensus_segments(sequence, providers, n),
    }
