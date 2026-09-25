"""
3D slab geometry TM block + the numeric core it needs.

Membrane placement (the ONLY thing this module decides):
  A fixed-thickness lipid slab is fitted to the C-alpha cloud by searching
  membrane-normal orientations on a Fibonacci sphere and, for each orientation,
  sliding a slab to maximise the two-sided hydrophobic contrast (see `_best_slab`).
  When the file already carries the bilayer (OPM / PPM / memembed DUM atoms) that
  placement is used instead. A residue is "Transmembrane" iff its C-alpha lies in
  the slab; the others are Cytoplasmic / Extracellular by face + positive-inside.

`GeometryTMProvider` exposes that as a TM block for the orchestrator, which merges it
with DSSP/STRIDE in the one consensus flow. The old standalone pipelines
(`predict_topology_structure`, `predict_topology_ss_first`) had their own SS rules
(snapping TM runs to helix ends, fusing partial elements); they now run the same
consensus flow and are kept only so older callers keep working.

Residue indexing: every per-residue list here is indexed by ResidueFrame POSITION.

Known limitations (state them in a paper):
  * One chain of the first model is analysed.
  * The hydrophobic thickness is a fixed parameter (default 30 A).
  * `membrane_score` is a heuristic (mean hydrophobic weight inside the slab); below
    MIN_MEMBRANE_SCORE the protein is reported as soluble.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

from app.core.constants import (
    KYTE_DOOLITTLE, MEMBRANE_THICKNESS, N_AXIS_SAMPLES, N_CENTER_SAMPLES, SMOOTH_WINDOW,
    MAX_JITTER_LEN, JITTER_MARGIN, MIN_TM_CORE, MIN_FACE_RESIDUES, MIN_FACE_FRACTION,
)
from app.schemas.topology import TMParams
from app.services.topology import labels as L
from app.services.topology.providers.tm.base import TMPrediction, TMProvider
from app.services.topology.membrane import membrane_from_file
from app.services.topology.residues import ResidueFrame, load_residue_frame

# =============================================================================
# STAGE 1 helpers - pure geometry / hydrophobicity (numpy only)
# =============================================================================
def _fibonacci_sphere(n: int) -> np.ndarray:
    """n near-uniform points on the unit sphere, used as candidate membrane
    normals."""
    phi = np.pi * (3.0 - np.sqrt(5.0))
    i = np.arange(n)
    y = 1.0 - (i / float(n - 1)) * 2.0 if n > 1 else np.array([0.0])
    r = np.sqrt(np.clip(1.0 - y * y, 0.0, None))
    theta = phi * i
    return np.stack([np.cos(theta) * r, y, np.sin(theta) * r], axis=1)


def _membrane_weights(residue_names, window: int = SMOOTH_WINDOW) -> np.ndarray:
    """Smoothed, SIGNED hydrophobic weight per residue (Kyte-Doolittle).

    The sign is essential and must NOT be clipped: hydrophilic residues carry
    negative weight, so a slab that swallows charged loops is penalised. This is
    what stops the fit from choosing a degenerate orientation whose slab simply
    contains the whole protein. See `_best_slab` for why the two-pointer sweep is
    still correct with signed weights.
    """
    raw = np.array([KYTE_DOOLITTLE.get(name, 0.0) for name in residue_names])
    half = window // 2
    return np.array([
        raw[max(0, i - half): min(len(raw), i + half + 1)].mean()
        for i in range(len(raw))
    ])


def _best_slab(coords: np.ndarray, weights: np.ndarray,
               thickness: float = MEMBRANE_THICKNESS,
               n_axes: int = N_AXIS_SAMPLES,
               n_centers: int = N_CENTER_SAMPLES):
    """Fit the membrane slab.

    Returns (axis, center, mean_in) where `axis` is the membrane normal, `center`
    is the slab centre projected on that axis (both in the centroid frame), and
    `mean_in` is the mean signed hydrophobicity inside the winning slab (used as a
    membrane-confidence score downstream).

    Objective: over all candidate normals and slab positions, maximise the
    TWO-SIDED hydrophobic contrast

        score = min( mean_in - mean_below ,  mean_in - mean_above )

    i.e. the slab must be more hydrophobic than the protein BOTH above and below it,
    scored on a slab of exactly `thickness` with >= MIN_FACE_RESIDUES on each face.
    The min-of-both-flanks is what centres the slab on the hydrophobic core: a slab
    shoved off-centre has one flank made of hydrophobic core residues, so that flank's
    contrast (and hence the min) collapses. See the inline comment for why the simpler
    objectives (max sum inside; max mean_in - mean_out; even a one-sided t-statistic)
    each fail on proteins narrower than the bilayer. Weights are signed (see
    `_membrane_weights`). A short-peptide fallback (nothing protrudes on both faces)
    keeps a plain max-sum slab so the function always returns a slab.
    """
    axes = _fibonacci_sphere(n_axes)
    centroid = coords.mean(axis=0)
    centered = coords - centroid
    n = len(coords)
    half = thickness / 2.0
    face_min = max(MIN_FACE_RESIDUES, int(np.ceil(MIN_FACE_FRACTION * n)))

    # Primary objective: TWO-SIDED hydrophobic contrast, min(mean_in - mean_below,
    # mean_in - mean_above), scored on a slab of EXACTLY `thickness` (identical to
    # the slab used later for classification) with >= MIN_FACE_RESIDUES on each face.
    # Why not a simpler objective:
    #   * max SUM inside -> engulfs the whole protein when it is narrower than the
    #     bilayer is thick (picks an in-plane axis).
    #   * max (mean_in - mean_out) -> maximised by excluding a few maximally
    #     hydrophilic residues; a near-engulfing diagonal wins.
    #   * one-sided separation (mean_in vs ALL outside, even a t-statistic) -> lets
    #     the slab slide off-centre, keeping just MIN_FACE residues on one face and
    #     dumping every hydrophilic residue on the other; the core ends up at the
    #     slab edge, merging adjacent crossings.
    # Requiring the core to beat BOTH flanks pins the slab onto a hydrophobic band
    # that is flanked by hydrophilic protein on each side - the defining profile of
    # a membrane-spanning region.
    best_score = -np.inf
    best_axis, best_center, best_mean_in = axes[0], 0.0, 0.0

    # Fallback for a structure too small to protrude on both faces (short peptide):
    # keep the plain max-hydrophobic-sum slab so the function always returns a slab.
    fb_sum = -np.inf
    fb_axis, fb_center, fb_mean_in = axes[0], 0.0, 0.0

    for axis in axes:
        z = centered @ axis                                   # (n,)
        centers = np.linspace(z.min(), z.max(), n_centers)    # (C,)
        d = z[None, :] - centers[:, None]                     # (C, n)
        inside = np.abs(d) <= half
        below_mask = d < -half
        above_mask = d > half
        n_in = inside.sum(axis=1).astype(float)
        n_below = below_mask.sum(axis=1).astype(float)
        n_above = above_mask.sum(axis=1).astype(float)
        sum_in = (inside * weights[None, :]).sum(axis=1)
        sum_below = (below_mask * weights[None, :]).sum(axis=1)
        sum_above = (above_mask * weights[None, :]).sum(axis=1)

        # unconstrained fallback: densest hydrophobic slab on this axis
        fi = int(np.argmax(sum_in))
        if n_in[fi] > 0 and sum_in[fi] > fb_sum:
            fb_sum = float(sum_in[fi])
            fb_axis, fb_center, fb_mean_in = axis, float(centers[fi]), float(sum_in[fi] / n_in[fi])

        valid = (n_below >= face_min) & (n_above >= face_min) & (n_in >= MIN_TM_CORE)
        if not valid.any():
            continue

        mean_in = sum_in / np.maximum(n_in, 1.0)
        mean_below = sum_below / np.maximum(n_below, 1.0)
        mean_above = sum_above / np.maximum(n_above, 1.0)
        score = np.where(valid,
                         np.minimum(mean_in - mean_below, mean_in - mean_above),
                         -np.inf)

        ci = int(np.argmax(score))
        if score[ci] > best_score:
            best_score = float(score[ci])
            best_axis, best_center, best_mean_in = axis, float(centers[ci]), float(mean_in[ci])

    if best_score == -np.inf:
        return fb_axis, fb_center, fb_mean_in
    return best_axis, best_center, best_mean_in


def _classify_by_slab(coords: np.ndarray, weights: np.ndarray,
                      thickness: float = MEMBRANE_THICKNESS):
    """Stage-1 decision. Every residue whose C-alpha lies inside the fitted slab is
    Transmembrane; residues above/below become Side_A/Side_B. Returns
    (classifications, d, membrane_score, axis, center) where `d` is each residue's
    signed distance from the slab centre along the membrane normal."""
    centroid = coords.mean(axis=0)
    axis, center, _ = _best_slab(coords, weights, thickness)
    z = (coords - centroid) @ axis
    half = thickness / 2.0
    d = z - center

    classifications = []
    inside_weights = []
    for wi, di in zip(weights, d):
        if abs(di) <= half:
            classifications.append("Transmembrane")
            inside_weights.append(wi)
        elif di > half:
            classifications.append("Side_A")
        else:
            classifications.append("Side_B")

    membrane_score = float(np.mean(inside_weights)) if inside_weights else 0.0
    return classifications, d, membrane_score, axis, float(center)


def _smooth_flickers(classifications, d, half,
                     max_jitter_len: int = MAX_JITTER_LEN,
                     jitter_margin: float = JITTER_MARGIN):
    """Merge ONLY boundary-jitter dips back into the membrane.

    A short Side run sitting between two Transmembrane runs is reclassified as
    Transmembrane iff every residue in it stays within `jitter_margin` of the slab
    face (`|d| <= half + jitter_margin`) - i.e. the chain merely grazed the hard
    slab edge. A Side run that reaches farther out is a genuine connecting loop and
    is preserved, so two distinct membrane crossings joined by a short (2-4 residue)
    turn - common in transporters and beta-barrels - are NEVER merged into one. This
    is deliberately conservative: it never merges two TM runs across a real loop, and
    it leaves spurious short TM specks to `_drop_short_tm`.
    """
    n = len(classifications)
    if n == 0:
        return classifications
    out = list(classifications)

    runs = []
    start = 0
    for i in range(1, n + 1):
        if i == n or out[i] != out[start]:
            runs.append((start, i - 1, out[start]))
            start = i

    for idx in range(1, len(runs) - 1):
        s, e, t = runs[idx]
        if t == "Transmembrane":
            continue
        if runs[idx - 1][2] != "Transmembrane" or runs[idx + 1][2] != "Transmembrane":
            continue
        if e - s + 1 > max_jitter_len:
            continue
        if all(abs(d[k]) <= half + jitter_margin for k in range(s, e + 1)):
            for k in range(s, e + 1):
                out[k] = "Transmembrane"
    return out


def _drop_short_tm(classifications, min_core: int = MIN_TM_CORE, d=None):
    """Demote TM runs whose in-slab core is shorter than `min_core` back to a side
    (spurious slab dips by a terminal tail or a sharp turn). The side is the one both
    neighbours share; if they differ (or at a terminus) and `d` is given, the face the
    run's C-alphas are closer to is used."""
    n = len(classifications)
    out = list(classifications)
    i = 0
    while i < n:
        if out[i] != "Transmembrane":
            i += 1
            continue
        j = i
        while j < n and out[j] == "Transmembrane":
            j += 1
        if j - i < min_core:
            left = out[i - 1] if i > 0 and out[i - 1] != "Transmembrane" else None
            right = out[j] if j < n and out[j] != "Transmembrane" else None
            if left and left == right:
                side = left
            elif d is not None:
                side = "Side_A" if float(np.mean(d[i:j])) >= 0 else "Side_B"
            else:
                side = left or right or "Side_A"
            for k in range(i, j):
                out[k] = side
        i = j
    return out


def _tm_runs(classifications) -> list[tuple[int, int]]:
    runs, i, n = [], 0, len(classifications)
    while i < n:
        if not L.is_tm(classifications[i]):
            i += 1
            continue
        j = i
        while j < n and L.is_tm(classifications[j]):
            j += 1
        runs.append((i, j - 1))
        i = j
    return runs


# =============================================================================
# Small helpers
# =============================================================================
def _apply_positive_inside_rule(descriptions, residues_data):
    """Rename Side_A/Side_B to Cytoplasmic/Extracellular by the positive-inside
    rule (the cytoplasmic face is enriched in Lys/Arg). TM labels pass through."""
    return L.apply_positive_inside_rule(descriptions, [r["name"] for r in residues_data])


def _groups_for_runs(n, runs):
    groups = [-1] * n
    for g, (s, e) in enumerate(runs):
        for k in range(s, e + 1):
            groups[k] = g
    return groups


def _sides_by_sign(d) -> list[str]:
    """Nearest slab face for every residue. (The old rule `Side_A if d > half else
    Side_B` sent every in-slab residue on the A half to side B.)"""
    return ["Side_A" if di >= 0 else "Side_B" for di in d]



def _load_structure_residues(file_path: Path, chain_id: str | None = None):
    """Kept for older callers: (residues_data, ca_coords) of the analysed chain."""
    frame = load_residue_frame(file_path, chain_id)
    return frame.residues_data(), frame.coords



def _helix_axis_normal(coords: np.ndarray, k: int = 4, breaks=None) -> np.ndarray:
    """Estimate the membrane normal from geometry alone: the dominant direction of
    the chain, taken as the top eigenvector of the scatter of local CA(i+k)-CA(i)
    vectors. Transmembrane helices are the longest straight runs so they dominate,
    and the outer-product scatter is sign-invariant, so up- and down-helices
    reinforce the same axis. Vectors that span a chain break (unresolved residues)
    are ignored - they point anywhere."""
    if len(coords) <= k:
        return np.array([0.0, 0.0, 1.0])
    dirs = coords[k:] - coords[:-k]
    if breaks is not None:
        br = np.asarray(breaks, dtype=int)
        cum = np.cumsum(br)
        spans_break = (cum[k:] - cum[:-k]) > 0         # any break in (i, i+k]
        dirs = dirs[~spans_break]
        if len(dirs) == 0:
            return np.array([0.0, 0.0, 1.0])
    norm = np.linalg.norm(dirs, axis=1, keepdims=True)
    norm[norm == 0] = 1.0
    dirs = dirs / norm
    _, vec = np.linalg.eigh(dirs.T @ dirs)
    axis = vec[:, -1]
    return axis / np.linalg.norm(axis)



# =============================================================================
# Older entry points - now the one consensus flow
# =============================================================================
def _run_consensus(file_path: Path, labeler: str, params: Optional[TMParams],
                   chain_id: Optional[str]) -> dict:
    """3D slab TM block + DSSP/STRIDE through the orchestrator (the only flow)."""
    from app.services.topology.orchestrator import TopologyOrchestrator   # lazy
    from app.services.topology.providers.ss.dssp import DSSPProvider
    from app.services.topology.providers.ss.stride import STRIDEProvider

    name = (labeler or "DSSP").strip().upper()
    ss = None if name in ("__NONE__", "NONE") else (
        STRIDEProvider() if name == "STRIDE" else DSSPProvider())
    response = TopologyOrchestrator(GeometryTMProvider(), ss).execute(
        Path(file_path), params=params, chain_id=chain_id)
    return response.model_dump()


def predict_topology_structure(file_path: Path, labeler: str = "DSSP",
                               params: TMParams | None = None,
                               chain_id: str | None = None) -> dict:
    """Kept for older callers: 3D slab + DSSP/STRIDE, consensus flow. Returns the
    TopologyResponse as a dict."""
    return _run_consensus(file_path, labeler, params, chain_id)


def predict_topology_ss_first(file_path: Path, labeler: str = "DSSP",
                              params: TMParams | None = None,
                              chain_id: str | None = None) -> dict:
    """Kept for older callers: identical to `predict_topology_structure` (there is
    one flow now)."""
    return _run_consensus(file_path, labeler, params, chain_id)


# =============================================================================
# Orchestrator provider: slab geometry as a TM block
# =============================================================================
class GeometryTMProvider(TMProvider):
    """Stage 1 of the slab pipeline as a TMProvider: per-position labels
    (Transmembrane / Cytoplasmic / Extracellular) on the shared ResidueFrame."""

    LABELER = "3D_Geometry"

    def predict_tm(self, file_path: Path, params: Optional[TMParams] = None,
                   frame: Optional[ResidueFrame] = None, **kwargs) -> TMPrediction:
        params = params if params is not None else TMParams()
        frame = frame if frame is not None else load_residue_frame(file_path, kwargs.get("chain_id"))
        n = len(frame)
        if n < 10:
            # (the old code returned TMPrediction(boundaries=[]) -> pydantic ValidationError)
            return TMPrediction(labeler=self.LABELER,
                                warnings=[f"only {n} residues - too few to fit a membrane slab"])

        residues_data = frame.residues_data()
        weights = _membrane_weights(frame.names)
        labeler = self.LABELER
        placed = membrane_from_file(frame)
        if placed is not None:
            # The file already carries the bilayer (OPM / PPM / memembed DUM atoms):
            # use that placement instead of re-fitting a hydrophobic slab.
            half = placed.half_thickness
            d = placed.depth
            inside = np.abs(d) <= half
            classifications = ["Transmembrane" if ins else ("Side_A" if di > 0 else "Side_B")
                               for ins, di in zip(inside, d)]
            membrane_score = float(weights[inside].mean()) if inside.any() else 0.0
            axis = placed.normal
            labeler = f"{self.LABELER} (membrane from file)"
        else:
            half = params.membrane_thickness / 2.0
            classifications, d, membrane_score, axis, _ = _classify_by_slab(
                frame.coords, weights, thickness=params.membrane_thickness)
        normal = [round(float(a), 4) for a in axis]

        if membrane_score < params.min_membrane_score:
            labels = _apply_positive_inside_rule(_sides_by_sign(d), residues_data)
            return TMPrediction(
                labels=labels, segments=[], membrane_score=membrane_score,
                membrane_normal=normal, labeler=labeler,
                warnings=[f"membrane_score {membrane_score:.2f} < {params.min_membrane_score} "
                          "- treated as soluble (no TM segments)"])

        classifications = _smooth_flickers(classifications, d, half)
        classifications = _drop_short_tm(classifications, d=d)
        labels = _apply_positive_inside_rule(classifications, residues_data)
        segments = _tm_runs(labels)
        return TMPrediction(
            regions=L.build_regions(frame, labels, _groups_for_runs(n, segments)),
            labels=labels,
            segments=segments,
            membrane_normal=normal,
            membrane_score=membrane_score,
            labeler=labeler,
        )
