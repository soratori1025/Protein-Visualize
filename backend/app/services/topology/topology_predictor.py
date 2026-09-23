"""
Unified transmembrane (TM) topology predictor.

The pipeline is deliberately split into two stages with a strict separation of
responsibilities, so that the reason a residue is called "in the membrane" is
always defensible:

  STAGE 1 - MEMBRANE PLACEMENT  (3D geometry + hydrophobicity, the ONLY decider)
    A fixed-thickness lipid slab is fitted to the C-alpha cloud by searching
    membrane-normal orientations on a Fibonacci sphere and, for each orientation,
    sliding a slab to maximise the *non-negative* hydrophobic weight it contains.
    A residue is transmembrane iff its C-alpha falls inside that slab. This is the
    only place where "is this in the membrane?" is answered, and it is answered
    from real PDB coordinates (this is the fix for false-positive TM calls made by
    pure sequence hydrophobicity scanning).

  STAGE 2 - SECONDARY-STRUCTURE LABELLING  (DSSP or STRIDE, DESCRIPTIVE ONLY)
    For each membrane-spanning run found in stage 1, the majority DSSP/STRIDE code
    decides Alpha helix / Beta strand / Loop, and the run boundaries are snapped
    outward to the true ends of the overlapping SS element. Snapping is capped so
    that at least one loop residue always remains between two TM runs, so two
    distinct crossings can never be merged. DSSP/STRIDE never adds or removes a TM
    call.

Residue indexing: the structure is parsed ONCE into a ResidueFrame (see
``app.services.topology.residues``) and every per-residue list in this module is
indexed by frame POSITION. DSSP/STRIDE output is mapped onto those positions by
residue key verified against residue identity (sequence-alignment fallback), never
by looking up a bare author residue number.

Sides are named Cytoplasmic / Extracellular geometrically (which face of the slab
a loop sits on) plus the positive-inside rule - NOT by topological alternation.
This keeps re-entrant / half-membrane loops correct, which matters for the
transporter/GPCR families this tool targets.

Design notes / known limitations (state these in the paper):
  * One chain of the first model is analysed (the requested chain, else the first
    chain that contains amino acids). For an oligomeric bundle, fit per chain or
    merge chains upstream before calling this.
  * The hydrophobic-core thickness is a fixed parameter (default 30 Angstrom); it
    is not optimised per structure.
  * `membrane_score` is a heuristic confidence (mean hydrophobic weight inside the
    slab). Non-membrane proteins are reported with no TM regions when it falls
    below MIN_MEMBRANE_SCORE; expose the score so a caller can apply its own threshold.

Heavy imports (Bio.PDB, the DSSP/STRIDE adapters) are loaded lazily inside the
functions that need them, so the numeric core can be imported and unit-tested with
numpy alone.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from app.core.constants import (
    KYTE_DOOLITTLE, POSITIVE_RESIDUES, HELIX_CODES, STRAND_CODES,
    MEMBRANE_THICKNESS, N_AXIS_SAMPLES, N_CENTER_SAMPLES, SMOOTH_WINDOW,
    MAX_JITTER_LEN, JITTER_MARGIN, MIN_TM_CORE, MIN_FACE_RESIDUES,
    MIN_FACE_FRACTION, MAX_SNAP, MIN_MEMBRANE_SCORE, MIN_EXTRA_SS_LEN,
    MIN_TM_ELEMENT_IN_SLAB, FULL_CROSS_FRAC, BROKEN_GAP_MAX, MIN_CROSS_SPAN_FRAC
)
from app.schemas.topology import TMParams
from app.services.topology import labels as L
from app.services.topology.providers.tm.base import TMPrediction, TMProvider
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
# STAGE 2 helpers - DSSP / STRIDE labelling (descriptive only)
# =============================================================================
def _run_labeler(file_path: Path, labeler: str, frame: ResidueFrame):
    """Per-POSITION DSSP/STRIDE codes on `frame` ('C' where the tool said nothing).

    Returns (codes, warnings). The old version returned {residue_number: code} for
    ALL chains keyed by whatever the tool printed (str or int), which was then looked
    up with int author numbers - chains overwrote each other and string keys never
    matched."""
    from app.services.topology.providers.ss.dssp import DSSPProvider
    from app.services.topology.providers.ss.stride import STRIDEProvider

    provider = STRIDEProvider() if (labeler or "DSSP").upper() == "STRIDE" else DSSPProvider()
    result = provider.assign(file_path, frame)
    if result.matched == 0:
        return [], result.warnings
    return [c if c is not None else "C" for c in result.codes], result.warnings


def _classify_code(code: str):
    if code in HELIX_CODES:
        return "Alpha"
    if code in STRAND_CODES:
        return "Beta"
    return None


def _ss_label(code: str) -> str:
    """DSSP/STRIDE code -> coarse SS class for drawing (Helix / Strand / Coil)."""
    return L.ss_word(L.coarse_ss(code))


def _label_extramembrane(descriptions, ss_codes: Sequence[str],
                         min_ss_len: int = MIN_EXTRA_SS_LEN, breaks=None):
    """Append secondary-structure detail to the extramembrane (Cytoplasmic /
    Extracellular) residues, so the loops of a topology snake-plot can show their
    own helices and strands (e.g. EL2, EL3a/3b in a transporter diagram) instead of
    a bare line. `ss_codes` is indexed by frame position."""
    return L.label_extramembrane_ss(descriptions, [L.coarse_ss(c) for c in ss_codes],
                                    min_len=min_ss_len, breaks=breaks)


def _region_side_ss(description: str):
    """(side, ss) tags derived from a region description (kept for older callers)."""
    _, side, ss = L.region_fields(description)
    return side, ss


def _label_and_snap(classifications, ss_codes: Sequence[str], max_snap: int = MAX_SNAP,
                    breaks=None):
    """Label each TM run Alpha/Beta/Loop by majority SS code and snap boundaries
    outward along the same SS element. Each side may grow by at most
    min(max_snap, (gap - 1) // 2) residues, so at least one loop residue always
    separates two crossings, and never across a chain break."""
    n = len(classifications)
    out = list(classifications)
    runs = _tm_runs(classifications)

    for idx, (start, end) in enumerate(runs):
        codes = ss_codes[start:end + 1]
        alpha = sum(1 for c in codes if c in HELIX_CODES)
        beta = sum(1 for c in codes if c in STRAND_CODES)
        if alpha == 0 and beta == 0:
            label, target = "Transmembrane Loop", None
        elif alpha >= beta:
            label, target = "Transmembrane Alpha Helix", HELIX_CODES
        else:
            label, target = "Transmembrane Beta Strand", STRAND_CODES

        s, e = start, end
        if target is not None:
            prev_end = runs[idx - 1][1] if idx > 0 else -1
            next_start = runs[idx + 1][0] if idx + 1 < len(runs) else n
            left_gap = start - prev_end - 1
            right_gap = next_start - end - 1
            # towards a terminus: half the tail (as before); between two runs: keep >= 1
            # loop residue (the old `gap // 2` on both sides could consume an even gap)
            left_budget = min(max_snap, left_gap // 2 if idx == 0 else (left_gap - 1) // 2)
            right_budget = min(max_snap, right_gap // 2 if idx + 1 == len(runs)
                               else (right_gap - 1) // 2)

            moved = 0
            while (s > 0 and moved < left_budget and ss_codes[s - 1] in target
                   and not (breaks is not None and breaks[s])):
                s -= 1
                moved += 1
            moved = 0
            while (e < n - 1 and moved < right_budget and ss_codes[e + 1] in target
                   and not (breaks is not None and breaks[e + 1])):
                e += 1
                moved += 1

        for k in range(s, e + 1):
            out[k] = label
    return out


def _apply_positive_inside_rule(descriptions, residues_data):
    """Rename Side_A/Side_B to Cytoplasmic/Extracellular by the positive-inside
    rule (the cytoplasmic face is enriched in Lys/Arg). TM labels pass through."""
    return L.apply_positive_inside_rule(descriptions, [r["name"] for r in residues_data])


def _build_regions(frame: ResidueFrame, descriptions, group_ids=None):
    """Per-position descriptions -> region dicts (author numbering + insertion codes)."""
    return [r.model_dump() for r in L.build_regions(frame, descriptions, group_ids)]


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


# =============================================================================
# Structure loading (lazy Bio.PDB)
# =============================================================================
def _load_structure_residues(file_path: Path, chain_id: str | None = None):
    """Kept for older callers: (residues_data, ca_coords) of the analysed chain."""
    frame = load_residue_frame(file_path, chain_id)
    return frame.residues_data(), frame.coords


def _empty_result(labeler, params, name, frame=None, warnings=None, score=0.0):
    return {"uniprot_id": "CALCULATED", "protein_name": name,
            "gene_name": "", "organism": "Computed",
            "membrane_score": round(float(score), 3), "labeler": labeler,
            "parameters_used": params.to_response_dict(), "regions": [],
            "chain_id": frame.chain_id if frame is not None else None,
            "warnings": list(warnings or [])}


# =============================================================================
# Public entry points
# =============================================================================
def predict_topology_structure(file_path: Path, labeler: str = "DSSP",
                               params: TMParams | None = None,
                               chain_id: str | None = None) -> dict:
    """Unified structure-based prediction: slab geometry gates TM (stage 1),
    DSSP/STRIDE labels alpha/beta/loop (stage 2)."""
    if params is None:
        params = TMParams()
    thickness = params.membrane_thickness

    frame = load_residue_frame(file_path, chain_id)
    residues_data = frame.residues_data()
    if len(frame) < 10:
        return _empty_result(labeler, params, "Too few residues", frame)

    weights = _membrane_weights(frame.names)
    breaks = frame.chain_breaks()
    warns: list[str] = []

    # STAGE 1 - geometry decides TM.
    classifications, d, membrane_score, axis, center = _classify_by_slab(
        frame.coords, weights, thickness)
    classifications = _smooth_flickers(classifications, d, thickness / 2.0)
    classifications = _drop_short_tm(classifications, d=d)

    # Guard: if the best slab is not convincingly hydrophobic, treat as soluble.
    if membrane_score < params.min_membrane_score:
        descriptions = _apply_positive_inside_rule(_sides_by_sign(d), residues_data)
        return {
            "uniprot_id": "CALCULATED",
            "protein_name": "No membrane slab detected (likely soluble)",
            "gene_name": "", "organism": "Computed",
            "membrane_score": round(membrane_score, 3), "labeler": labeler,
            "parameters_used": params.to_response_dict(),
            "regions": _build_regions(frame, descriptions),
            "chain_id": frame.chain_id, "warnings": warns,
        }

    # STAGE 2 - DSSP/STRIDE labels the already-fixed TM runs.
    ss_codes, labeler_used = [], labeler
    if labeler == "__none__":
        labeler_used = "geometry-only (no SS tool)"
    else:
        try:
            ss_codes, ss_warns = _run_labeler(file_path, labeler, frame)
            warns.extend(ss_warns)
        except Exception as error:  # noqa: BLE001 - tool may be missing at runtime
            print(f"[topology_predictor] {labeler} failed, geometry-only labels. {error}")
            warns.append(f"{labeler} unavailable ({error}); geometry-only labels")
        if not ss_codes:
            labeler_used = f"{labeler} (unavailable)"

    if ss_codes:
        descriptions = _label_and_snap(classifications, ss_codes, breaks=breaks)
        descriptions = _apply_positive_inside_rule(descriptions, residues_data)
        # Give extramembrane loops their own helix/strand/coil elements to draw.
        descriptions = _label_extramembrane(descriptions, ss_codes, breaks=breaks)
    else:
        # No SS tool: keep TM calls, but we cannot tell alpha from beta.
        descriptions = ["Transmembrane (helical, unverified)"
                        if c == "Transmembrane" else c for c in classifications]
        descriptions = _apply_positive_inside_rule(descriptions, residues_data)

    groups = _groups_for_runs(len(frame), _tm_runs(descriptions))
    return {
        "uniprot_id": "CALCULATED",
        "protein_name": f"Structure-based TM prediction (slab-fit + {labeler_used})",
        "gene_name": "", "organism": "Computed",
        "membrane_score": round(membrane_score, 3),
        "membrane_normal": [round(float(a), 4) for a in axis],
        "labeler": labeler_used,
        "parameters_used": params.to_response_dict(),
        "regions": _build_regions(frame, descriptions, groups),
        "chain_id": frame.chain_id, "warnings": warns,
    }


# =============================================================================
# SS-element-first pipeline  (more accurate boundaries on real polytopic proteins)
# =============================================================================
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


def _ss_runs(ss_codes: Sequence[str], breaks=None):
    """Maximal contiguous runs of one coarse SS class over the frame positions,
    split at chain breaks. Returns [(start_idx, end_idx, 'Helix'|'Strand'|'Coil'), ...]."""
    classes = [_ss_label(c) for c in ss_codes]
    runs = []
    i, n = 0, len(classes)
    while i < n:
        j = i + 1
        while j < n and classes[j] == classes[i] and not (breaks is not None and breaks[j]):
            j += 1
        runs.append((i, j - 1, classes[i]))
        i = j
    return runs


def _membrane_center_geometric(z, runs, half, loose: int = 6):
    """Locate the membrane centre along a fixed normal from GEOMETRY, not
    hydrophobicity (which is ~uniform along a helix and cannot localise a centre).
    Every TM helix is centred on the membrane, so the mean position of the residues
    that belong to membrane-embedded SS elements is the membrane centre. Two or three
    fixed-point passes converge: start from the median of all SS-element residues,
    classify which elements sit in the slab, recentre on their residues, repeat."""
    ss_idx = [i for (a, b, c) in runs if c != "Coil" for i in range(a, b + 1)]
    if not ss_idx:
        return float(np.median(z))
    center = float(np.median(z[ss_idx]))
    for _ in range(4):
        in_slab = np.abs(z - center) <= half
        tm_idx = []
        for a, b, c in runs:
            if c != "Coil" and int(in_slab[a:b + 1].sum()) >= loose:
                tm_idx.extend(range(a, b + 1))
        if not tm_idx:
            break
        new_center = float(np.mean(z[tm_idx]))
        if abs(new_center - center) < 0.1:
            center = new_center
            break
        center = new_center
    return center


def _membrane_crossings(runs, z, center, half, thickness, params: TMParams):
    """Turn DSSP/STRIDE secondary-structure elements into transmembrane crossings.

    An element (helix or strand) is transmembrane when >= params.min_tm_element_in_slab
    of its residues fall inside the slab. A crossing is one full-height element, OR two
    (rarely more) short partial elements fused across a *short, intramembrane* break
    whose two parts continue in the same direction and together span the bilayer -
    that is a discontinuous helix (TM1a/1b, TM6a/6b), which stays one crossing.
    Two antiparallel helices joined by a loop that leaves the membrane never fuse,
    because the break is not intramembrane. Returns [(start_idx, end_idx, cls), ...].
    """
    d = z - center
    in_slab = np.abs(d) <= half
    full = params.full_cross_frac * thickness

    elems = []
    for a, b, cls in runs:
        if cls == "Coil":
            continue
        if int(in_slab[a:b + 1].sum()) < params.min_tm_element_in_slab:
            continue  # extramembrane helix/strand (a loop feature), not a crossing
        zr = z[a:b + 1]
        elems.append({"a": a, "b": b, "cls": cls, "span": float(zr.max() - zr.min())})

    crossings = []
    i = 0
    while i < len(elems):
        e = elems[i]
        a, b = e["a"], e["b"]
        if e["span"] >= full:
            crossings.append((a, b, e["cls"]))
            i += 1
            continue
        # partial element: try to fuse with the next partial(s) into one crossing
        j = i
        while j + 1 < len(elems):
            nxt = elems[j + 1]
            gap_lo, gap_hi = elems[j]["b"] + 1, nxt["a"]
            gap = gap_hi - gap_lo
            gap_intramembrane = gap == 0 or bool(np.all(in_slab[gap_lo:gap_hi]))
            comb_span = float(z[a:nxt["b"] + 1].max() - z[a:nxt["b"] + 1].min())
            if (nxt["span"] < full and gap <= params.broken_gap_max
                    and gap_intramembrane and comb_span <= thickness * 1.5):
                b = nxt["b"]
                j += 1
                if comb_span >= full:
                    break
            else:
                break
        crossings.append((a, b, e["cls"]))
        i = j + 1

    # Keep only elements that actually span the bilayer. A helix confined to one
    # leaflet (interfacial or re-entrant) clears min_tm_element_in_slab but does not
    # cross, so it must be dropped here rather than counted as a TM segment.
    min_span = params.min_cross_span_frac * thickness
    return [(a, b, c) for (a, b, c) in crossings
            if float(z[a:b + 1].max() - z[a:b + 1].min()) >= min_span]


def predict_topology_ss_first(file_path: Path, labeler: str = "DSSP",
                              params: TMParams | None = None,
                              chain_id: str | None = None) -> dict:
    """Recommended structure-based predictor for real polytopic proteins.

    Order is inverted vs `predict_topology_structure`: DSSP/STRIDE define the actual
    secondary-structure ELEMENTS (accurate boundaries, and two adjacent helices are
    never one blob), and the membrane slab only decides which elements cross it. The
    normal is estimated from helix geometry, so a weakly-hydrophobic transporter
    core no longer tilts the fit. Falls back to the slab-first pipeline if no SS
    tool is available."""
    if params is None:
        params = TMParams()
    thickness = params.membrane_thickness

    frame = load_residue_frame(file_path, chain_id)
    residues_data = frame.residues_data()
    if len(frame) < 10:
        return _empty_result(labeler, params, "Too few residues", frame)

    ss_codes, labeler_used, warns = [], labeler, []
    try:
        ss_codes, warns = _run_labeler(file_path, labeler, frame)
    except Exception as error:  # noqa: BLE001
        print(f"[topology_predictor] {labeler} failed; falling back to slab-first. {error}")
    if not ss_codes:
        return predict_topology_structure(file_path, labeler="__none__", params=params,
                                          chain_id=frame.chain_id)

    weights = _membrane_weights(frame.names)
    breaks = frame.chain_breaks()
    coords = frame.coords

    normal = _helix_axis_normal(coords, breaks=breaks)
    centroid = coords.mean(axis=0)
    z = (coords - centroid) @ normal
    half = thickness / 2.0

    runs = _ss_runs(ss_codes, breaks)
    center = _membrane_center_geometric(z, runs, half)
    d = z - center
    inside = np.abs(d) <= half
    mean_in = float(weights[inside].mean()) if inside.any() else 0.0
    crossings = _membrane_crossings(runs, z, center, half, thickness, params)

    if not crossings or mean_in < params.min_membrane_score:
        descriptions = _apply_positive_inside_rule(_sides_by_sign(d), residues_data)
        descriptions = _label_extramembrane(descriptions, ss_codes, breaks=breaks)
        return {
            "uniprot_id": "CALCULATED",
            "protein_name": "No transmembrane crossings detected (likely soluble)",
            "gene_name": "", "organism": "Computed",
            "membrane_score": round(mean_in, 3), "labeler": labeler_used,
            "parameters_used": params.to_response_dict(),
            "regions": _build_regions(frame, descriptions),
            "chain_id": frame.chain_id, "warnings": warns,
        }

    # Per-residue labels: TM crossings first, then sides, then extramembrane SS.
    descriptions = _sides_by_sign(d)
    groups = [-1] * len(frame)
    for g, (a, b, cls) in enumerate(crossings):
        label = "Transmembrane Alpha Helix" if cls == "Helix" else "Transmembrane Beta Strand"
        for k in range(a, b + 1):
            descriptions[k] = label
            groups[k] = g
    descriptions = _apply_positive_inside_rule(descriptions, residues_data)
    descriptions = _label_extramembrane(descriptions, ss_codes, breaks=breaks)

    return {
        "uniprot_id": "CALCULATED",
        "protein_name": f"Structure-based TM prediction (SS-element-first + {labeler_used})",
        "gene_name": "", "organism": "Computed",
        "membrane_score": round(mean_in, 3),
        "membrane_normal": [round(float(a), 4) for a in normal],
        "labeler": labeler_used,
        "parameters_used": params.to_response_dict(),
        "regions": _build_regions(frame, descriptions, groups),
        "chain_id": frame.chain_id, "warnings": warns,
    }


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
        half = params.membrane_thickness / 2.0
        classifications, d, membrane_score, axis, _ = _classify_by_slab(
            frame.coords, weights, thickness=params.membrane_thickness)
        normal = [round(float(a), 4) for a in axis]

        if membrane_score < params.min_membrane_score:
            labels = _apply_positive_inside_rule(_sides_by_sign(d), residues_data)
            return TMPrediction(
                labels=labels, segments=[], membrane_score=membrane_score,
                membrane_normal=normal, labeler=self.LABELER,
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
            labeler=self.LABELER,
        )
