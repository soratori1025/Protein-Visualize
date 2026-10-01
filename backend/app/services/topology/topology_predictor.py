"""
3D slab geometry TM block + the numeric core it needs.

Membrane placement (the ONLY thing this module decides):
  A fixed-thickness lipid slab is fitted to the C-alpha cloud by searching
  membrane-normal orientations on a Fibonacci sphere and, for each orientation,
  sliding a slab to maximise the two-sided hydrophobic contrast (see `_best_slab`).
  When the file already carries the bilayer (OPM / PPM / memembed DUM atoms) that
  placement is used instead. A residue starts as "Transmembrane" iff its C-alpha lies
  in the slab; the others are Cytoplasmic / Extracellular by face + positive-inside.

Post-processing of the in-slab runs (each step is chain-break aware - runs are never
joined across unresolved residues):
  1. `_smooth_flickers`   - heal boundary jitter, but never across a real turn.
  2. `_split_at_turns`    - a run that climbs to one face and comes back down is two
                            crossings (or a re-entrant loop), not one segment.
  3. `_drop_short_tm`     - demote in-slab specks shorter than MIN_TM_CORE.
  4. `_drop_non_crossing` - demote runs that do not traverse the slab (interfacial /
                            re-entrant) or are strongly hydrophilic (e.g. a disordered
                            AlphaFold tail that happens to pass through the slab plane).

`GeometryTMProvider` exposes that as a TM block for the orchestrator, which merges it
with DSSP/STRIDE in the one consensus flow. The old standalone pipelines
(`predict_topology_structure`, `predict_topology_ss_first`) had their own SS rules
(snapping TM runs to helix ends, fusing partial elements); they now run the same
consensus flow and are kept only so older callers keep working.

Residue indexing: every per-residue list here is indexed by ResidueFrame POSITION.

Known limitations (state them in a paper):
  * One chain of the first model is classified. The slab is fitted on that chain only
    unless the caller passes the other chains as `context_frames`. The orientation term
    in `_best_slab` keeps narrow chains from being engulfed by a tilted slab; a warning
    is still emitted if most of the chain ends up inside the slab.
  * Beta-barrel strands alternate hydrophobic/polar residues, so their mean
    Kyte-Doolittle weight is near 0 and `membrane_score` can fall below
    MIN_MEMBRANE_SCORE: outer-membrane proteins may be reported as soluble.
  * The hydrophobic thickness is a fixed parameter (default 30 A).
  * `membrane_score` is a heuristic (mean hydrophobic weight inside the slab); below
    MIN_MEMBRANE_SCORE the protein is reported as soluble.
  * In-plane interfacial helices that continue straight out of a TM helix inside the
    slab edge are not split off here (DSSP/STRIDE element boundaries do that).
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

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

TM_LABEL = "Transmembrane"

# --- local tuning constants (move to app.core.constants if you want them in the UI) ---
MAX_CA_CA = 4.2               # A; consecutive C-alphas farther apart are NOT bonded (gap)
TURN_TRAVEL = 3.0             # A; min travel along the normal on each side of a Side dip
                              #    for it to count as a real turn (not jitter)
TURN_LOOKAHEAD = 6            # residues used to measure that travel
SPLIT_TRAVEL_FRAC = 0.4       # a run is split at an apex when both arms travel
                              #    >= SPLIT_TRAVEL_FRAC * thickness along the normal
APEX_BAND = 1.5               # A; residues this close to the apex depth are demoted
MIN_RUN_SPAN_FRAC = 0.5       # a TM run must span >= this fraction of the thickness
MIN_RUN_HYDROPATHY = -1.0     # mean raw Kyte-Doolittle of a TM run; below = hydrophilic
MIN_OPEN_RUN_SPAN_FRAC = 0.25 # same for a run cut by a chain break / terminus
ENGULF_FRACTION = 0.80        # warn when this much of the chain sits inside the slab
DIRECTION_SPAN = 4            # residues; local chain direction = CA(i+2) - CA(i-2)
SLAB_FIT_WINDOW = 7           # residues; hydropathy smoothing used to FIT the slab. The
                              #    19-residue SMOOTH_WINDOW blurs the hydrophobic band
                              #    (worse normal/centre); it is still used for
                              #    membrane_score, where it keeps soluble cores out.

# =============================================================================
# Residue-name handling for the hydrophobicity table
# =============================================================================
_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
_ONE_TO_THREE = {v: k for k, v in _THREE_TO_ONE.items()}
# modified / force-field-specific residue names seen in PDB, OPM and MD files
_ALIASES = {
    "MSE": "MET", "SEC": "CYS", "CYX": "CYS", "CYM": "CYS", "HSD": "HIS", "HSE": "HIS",
    "HSP": "HIS", "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HISD": "HIS", "HISE": "HIS",
    "HISH": "HIS", "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "ARN": "ARG", "SEP": "SER",
    "TPO": "THR", "PTR": "TYR", "MLY": "LYS", "KCX": "LYS", "PYL": "LYS", "HYP": "PRO",
}
_KD_TABLE = {str(k).strip().upper(): float(v) for k, v in KYTE_DOOLITTLE.items()}


def _kd_value(name) -> Optional[float]:
    """Kyte-Doolittle value for a residue name whatever the table is keyed by
    (three-letter or one-letter, any case), with common modified residues mapped to
    their parent. None if unknown."""
    key = str(name).strip().upper()
    key = _ALIASES.get(key, key)
    for cand in (key, _THREE_TO_ONE.get(key), _ONE_TO_THREE.get(key)):
        if cand is not None and cand in _KD_TABLE:
            return _KD_TABLE[cand]
    return None


def _kd_raw(residue_names) -> tuple[np.ndarray, list[str]]:
    """Raw per-residue Kyte-Doolittle values (unknown names -> 0.0) + the unknown names."""
    vals, unknown = [], []
    for name in residue_names:
        v = _kd_value(name)
        if v is None:
            unknown.append(str(name))
            v = 0.0
        vals.append(v)
    return np.asarray(vals, dtype=float), unknown


# =============================================================================
# Chain continuity
# =============================================================================
def _chain_breaks(coords: np.ndarray) -> np.ndarray:
    """breaks[i] is True iff residue i is NOT bonded to residue i-1 (unresolved residues in
    between). breaks[0] is False. Decided from geometry (C-alpha distance), which is what
    matters here and is robust to numbering gaps/insertion codes."""
    n = len(coords)
    out = np.zeros(n, dtype=bool)
    if n > 1:
        out[1:] = np.linalg.norm(np.diff(np.asarray(coords, float), axis=0), axis=1) > MAX_CA_CA
    return out


def _fragments(n: int, breaks=None) -> list[tuple[int, int]]:
    """Continuous pieces (start, end inclusive) of the chain."""
    if n == 0:
        return []
    if breaks is None:
        return [(0, n - 1)]
    starts = [0] + [i for i in range(1, n) if breaks[i]]
    ends = [s - 1 for s in starts[1:]] + [n - 1]
    return list(zip(starts, ends))


def _runs_where(mask: Sequence[bool], breaks=None) -> list[tuple[int, int]]:
    """Maximal runs (start, end inclusive) of True in `mask`, split at chain breaks."""
    n = len(mask)
    runs = []
    i = 0
    while i < n:
        if not mask[i]:
            i += 1
            continue
        j = i + 1
        while j < n and mask[j] and not (breaks is not None and breaks[j]):
            j += 1
        runs.append((i, j - 1))
        i = j
    return runs


def _side_for(d_values) -> str:
    return "Side_A" if float(np.mean(d_values)) >= 0 else "Side_B"


def _local_directions(coords: np.ndarray, breaks=None, span: int = DIRECTION_SPAN) -> np.ndarray:
    """Unit vector of the local chain direction at every residue, CA(i+span/2) -
    CA(i-span/2), never taken across a chain break (the window is shifted inside the
    fragment; one-residue fragments get a zero vector). For an alpha-helix the i -> i+4
    vector is almost parallel to the helix axis; for a strand it follows the strand."""
    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    out = np.zeros((n, 3))
    if breaks is None:
        breaks = _chain_breaks(coords)
    for s, e in _fragments(n, breaks):
        if e == s:
            continue
        idx = np.arange(s, e + 1)
        a = np.clip(idx - span // 2, s, e)
        b = np.clip(a + span, s, e)
        a = np.clip(b - span, s, e)
        v = coords[b] - coords[a]
        norm = np.linalg.norm(v, axis=1, keepdims=True)
        norm[norm == 0] = 1.0
        out[s:e + 1] = v / norm
    return out


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


def _membrane_weights(residue_names, window: int = SMOOTH_WINDOW, breaks=None) -> np.ndarray:
    """Smoothed, SIGNED hydrophobic weight per residue (Kyte-Doolittle).

    The sign is essential and must NOT be clipped: hydrophilic residues carry
    negative weight, so a slab that swallows charged loops is penalised. This is
    what stops the fit from choosing a degenerate orientation whose slab simply
    contains the whole protein.

    The sliding window never crosses a chain break: residues on either side of an
    unresolved stretch are not sequence neighbours.
    """
    raw, _ = _kd_raw(residue_names)
    n = len(raw)
    out = np.empty(n, dtype=float)
    half = window // 2
    for s, e in _fragments(n, breaks):
        seg = raw[s:e + 1]
        c = np.concatenate([[0.0], np.cumsum(seg)])
        idx = np.arange(len(seg))
        lo = np.maximum(0, idx - half)
        hi = np.minimum(len(seg), idx + half + 1)
        out[s:e + 1] = (c[hi] - c[lo]) / (hi - lo)
    return out


def _best_slab(coords: np.ndarray, weights: np.ndarray,
               thickness: float = MEMBRANE_THICKNESS,
               n_axes: int = N_AXIS_SAMPLES,
               n_centers: int = N_CENTER_SAMPLES,
               directions: Optional[np.ndarray] = None):
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
    contrast (and hence the min) collapses. Weights are signed (see
    `_membrane_weights`). A short-peptide fallback (nothing protrudes on both faces)
    keeps a plain max-sum slab so the function always returns a slab.

    Orientation term (`directions`, unit local chain directions, see
    `_local_directions`): a positive contrast is multiplied by the mean |cos| between
    the in-slab chain directions and the candidate normal. Without it, a chain that is
    narrower than the slab is thick AND has its termini on opposite faces (odd number
    of TM segments: 3, 5, compact 7-TM, one subunit of an oligomer) is fitted with the
    slab tilted ~70 deg: the slab swallows the whole TM domain lying on its side and
    uses the two polar terminal tails as its flanks, which beats the true contrast.
    Membrane-spanning helices/strands run roughly along the normal (tilt < ~40 deg), so
    that degenerate slab has a low alignment and loses. With directions=None the old
    contrast-only objective is used.
    """
    axes = _fibonacci_sphere(n_axes)
    centroid = coords.mean(axis=0)
    centered = coords - centroid
    n = len(coords)
    half = thickness / 2.0
    face_min = max(MIN_FACE_RESIDUES, int(np.ceil(MIN_FACE_FRACTION * n)))
    align = None if directions is None else np.asarray(directions, dtype=float)

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
        contrast = np.minimum(mean_in - mean_below, mean_in - mean_above)
        if align is not None:
            cosang = np.abs(align @ axis)                                  # (n,)
            mean_align = (inside * cosang[None, :]).sum(axis=1) / np.maximum(n_in, 1.0)
            contrast = np.where(contrast > 0, contrast * mean_align, contrast)
        score = np.where(valid, contrast, -np.inf)

        ci = int(np.argmax(score))
        if score[ci] > best_score:
            best_score = float(score[ci])
            best_axis, best_center, best_mean_in = axis, float(centers[ci]), float(mean_in[ci])

    if best_score == -np.inf:
        return fb_axis, fb_center, fb_mean_in
    return best_axis, best_center, best_mean_in


def _classify_by_slab(coords: np.ndarray, weights: np.ndarray,
                      thickness: float = MEMBRANE_THICKNESS,
                      fit_coords: Optional[np.ndarray] = None,
                      fit_weights: Optional[np.ndarray] = None,
                      fit_directions: Optional[np.ndarray] = None):
    """Stage-1 decision. Every residue whose C-alpha lies inside the fitted slab is
    Transmembrane; residues above/below become Side_A/Side_B. Returns
    (classifications, d, membrane_score, axis, center) where `d` is each residue's
    signed distance from the slab centre along the membrane normal.

    The slab is fitted on (`fit_coords`, `fit_weights`) when given - e.g. every chain
    of the assembly - and the residues of `coords` are classified against it. Chains
    can simply be concatenated: the jump between two chains is > MAX_CA_CA, so it is
    treated as a break everywhere."""
    if fit_coords is None or fit_weights is None:
        fit_coords, fit_weights = coords, weights
    fit_coords = np.asarray(fit_coords, dtype=float)
    if fit_directions is None:
        fit_directions = _local_directions(fit_coords, _chain_breaks(fit_coords))
    centroid = fit_coords.mean(axis=0)
    axis, center, _ = _best_slab(fit_coords, fit_weights, thickness,
                                 directions=fit_directions)
    z = (coords - centroid) @ axis
    half = thickness / 2.0
    d = z - center

    classifications = []
    inside_weights = []
    for wi, di in zip(weights, d):
        if abs(di) <= half:
            classifications.append(TM_LABEL)
            inside_weights.append(wi)
        elif di > half:
            classifications.append("Side_A")
        else:
            classifications.append("Side_B")

    membrane_score = float(np.mean(inside_weights)) if inside_weights else 0.0
    return classifications, d, membrane_score, axis, float(center)


def _smooth_flickers(classifications, d, half,
                     max_jitter_len: int = MAX_JITTER_LEN,
                     jitter_margin: float = JITTER_MARGIN,
                     breaks=None,
                     turn_travel: float = TURN_TRAVEL,
                     lookahead: int = TURN_LOOKAHEAD):
    """Merge ONLY boundary-jitter dips back into the membrane.

    A short Side run sitting between two Transmembrane runs is reclassified as
    Transmembrane iff every residue in it stays within `jitter_margin` of the slab
    face (`|d| <= half + jitter_margin`) - i.e. the chain merely grazed the slab edge -
    AND it is not a turn. A dip is a turn when the chain travels towards the face before
    it and away from the face after it (>= `turn_travel` A each way over `lookahead`
    residues): that is two crossings joined by a short loop (helical hairpins in
    transporters, short ECLs in GPCRs, beta hairpins) and must stay two segments. The
    old rule merged those whenever the loop hugged the face. A dip that spans a chain
    break is never healed.
    """
    n = len(classifications)
    if n == 0:
        return classifications
    out = list(classifications)
    d = np.asarray(d, dtype=float)

    runs = []
    start = 0
    for i in range(1, n + 1):
        if i == n or out[i] != out[start] or (breaks is not None and breaks[i]):
            runs.append((start, i - 1, out[start]))
            start = i

    for idx in range(1, len(runs) - 1):
        s, e, t = runs[idx]
        if t == TM_LABEL:
            continue
        ps, pe, pt = runs[idx - 1]
        ns, ne, nt = runs[idx + 1]
        if pt != TM_LABEL or nt != TM_LABEL:
            continue
        if pe + 1 != s or e + 1 != ns:
            continue
        if breaks is not None and (breaks[s] or breaks[ns] or np.any(breaks[s:e + 1])):
            continue
        if e - s + 1 > max_jitter_len:
            continue
        if not all(abs(d[k]) <= half + jitter_margin for k in range(s, e + 1)):
            continue
        before = d[pe] - d[max(ps, pe - lookahead)]
        after = d[min(ne, ns + lookahead)] - d[ns]
        if abs(before) >= turn_travel and abs(after) >= turn_travel \
                and np.sign(before) != np.sign(after):
            continue                                    # genuine turn: keep the loop
        for k in range(s, e + 1):
            out[k] = TM_LABEL
    return out


def _split_at_turns(classifications, d, half, breaks=None,
                    travel_frac: float = SPLIT_TRAVEL_FRAC, apex_band: float = APEX_BAND):
    """Split an in-slab run where the chain climbs towards one face and comes back.

    When the connecting loop of a hairpin stays inside the slab (short loop, slab a few
    A off-centre) there is no Side dip at all and the two crossings arrive here as one
    long run. A run is cut at its apex when BOTH arms travel >= travel_frac * thickness
    along the normal; the residues within `apex_band` of the apex depth are demoted to
    the nearer face. A re-entrant loop (dips half-way in and returns to the same face)
    is cut the same way and its two halves are then rejected by `_drop_non_crossing`.
    """
    out = list(classifications)
    d = np.asarray(d, dtype=float)
    min_travel = travel_frac * 2.0 * half
    stack = _runs_where([c == TM_LABEL for c in out], breaks)
    while stack:
        s, e = stack.pop()
        if e - s < 2:
            continue
        seg = d[s:e + 1]
        pre_min, pre_max = np.minimum.accumulate(seg), np.maximum.accumulate(seg)
        suf_min = np.minimum.accumulate(seg[::-1])[::-1]
        suf_max = np.maximum.accumulate(seg[::-1])[::-1]
        top = np.minimum(seg - pre_min, seg - suf_min)      # up to the apex and down again
        bottom = np.minimum(pre_max - seg, suf_max - seg)   # down to the apex and up again
        depth = np.maximum(top, bottom)
        m = int(np.argmax(depth))
        if depth[m] < min_travel:
            continue
        apex = seg[m]
        lo = m
        while lo > 0 and abs(seg[lo - 1] - apex) <= apex_band:
            lo -= 1
        hi = m
        while hi < len(seg) - 1 and abs(seg[hi + 1] - apex) <= apex_band:
            hi += 1
        side = "Side_A" if apex >= 0 else "Side_B"
        for k in range(s + lo, s + hi + 1):
            out[k] = side
        if lo > 0:
            stack.append((s, s + lo - 1))
        if hi < len(seg) - 1:
            stack.append((s + hi + 1, e))
    return out


def _drop_short_tm(classifications, min_core: int = MIN_TM_CORE, d=None, breaks=None):
    """Demote TM runs whose in-slab core is shorter than `min_core` back to a side
    (spurious slab dips by a terminal tail or a sharp turn). The side is the one both
    neighbours share; if they differ (or at a terminus) and `d` is given, the face the
    run's C-alphas are closer to is used. Runs are split at chain breaks."""
    n = len(classifications)
    out = list(classifications)
    for i, j_last in _runs_where([c == TM_LABEL for c in out], breaks):
        j = j_last + 1
        if j - i >= min_core:
            continue
        left = out[i - 1] if i > 0 and out[i - 1] != TM_LABEL else None
        right = out[j] if j < n and out[j] != TM_LABEL else None
        if left and left == right:
            side = left
        elif d is not None:
            side = _side_for(np.asarray(d)[i:j])
        else:
            side = left or right or "Side_A"
        for k in range(i, j):
            out[k] = side
    return out


def _drop_non_crossing(classifications, d, half, raw_kd, breaks=None,
                       min_span_frac: float = MIN_RUN_SPAN_FRAC,
                       min_hydropathy: float = MIN_RUN_HYDROPATHY,
                       open_span_frac: float = MIN_OPEN_RUN_SPAN_FRAC):
    """Demote in-slab runs that are not membrane crossings. Returns (labels, n_dropped).

    * span: a crossing must travel >= min_span_frac * thickness along the normal. Runs
      that lie along the interface or dip in and back out do not. Runs cut short by a
      chain break or a chain terminus only need open_span_frac * thickness (the rest of
      the crossing is missing) - enough to drop the tip of an in-plane helix that starts
      right after an unresolved stretch and grazes the core (Piezo1 1513-1518).
    * hydropathy: a run whose mean raw Kyte-Doolittle value is below min_hydropathy is a
      polar strand that merely passes through the slab plane (typical for low-pLDDT
      AlphaFold tails). The threshold is lenient enough for beta-barrel strands.
    """
    n = len(classifications)
    out = list(classifications)
    d = np.asarray(d, dtype=float)
    dropped = 0
    for s, e in _runs_where([c == TM_LABEL for c in out], breaks):
        open_start = s > 0 and not (breaks is not None and breaks[s])
        open_end = e < n - 1 and not (breaks is not None and breaks[e + 1])
        span = float(d[s:e + 1].max() - d[s:e + 1].min())
        need = min_span_frac if (open_start and open_end) else open_span_frac
        too_flat = span < need * 2.0 * half
        too_polar = float(np.mean(raw_kd[s:e + 1])) < min_hydropathy
        if too_flat or too_polar:
            side = _side_for(d[s:e + 1])
            for k in range(s, e + 1):
                out[k] = side
            dropped += 1
    return out, dropped


def _tm_runs(classifications, breaks=None) -> list[tuple[int, int]]:
    """TM segments (start, end inclusive); never joined across a chain break."""
    return _runs_where([L.is_tm(c) for c in classifications], breaks)


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
    (Transmembrane / Cytoplasmic / Extracellular) on the shared ResidueFrame.

    Optional kwarg `context_frames`: ResidueFrames of the OTHER chains of the same model.
    When given, the slab is fitted on the whole assembly and only this chain is
    classified - the fix for oligomers whose single subunit is narrower than the bilayer.
    """

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

        warnings: list[str] = []
        residues_data = frame.residues_data()
        coords = np.asarray(frame.coords, dtype=float)
        breaks = _chain_breaks(coords)
        raw_kd, unknown = _kd_raw(frame.names)
        weights = _membrane_weights(frame.names, breaks=breaks)          # membrane_score
        fit_w = _membrane_weights(frame.names, window=SLAB_FIT_WINDOW, breaks=breaks)
        if unknown:
            uniq = sorted(set(unknown))
            warnings.append(f"{len(unknown)} residue(s) not in the Kyte-Doolittle table "
                            f"({', '.join(uniq[:6])}{'...' if len(uniq) > 6 else ''}) - "
                            "weighted 0")
        if breaks.any():
            warnings.append(f"{int(breaks.sum())} chain break(s) (unresolved residues); "
                            "TM segments are not joined across them")

        labeler = self.LABELER
        placed = membrane_from_file(frame)
        if placed is not None:
            # The file already carries the bilayer (OPM / PPM / memembed DUM atoms):
            # use that placement instead of re-fitting a hydrophobic slab.
            half = placed.half_thickness
            d = np.asarray(placed.depth, dtype=float)
            inside = np.abs(d) <= half
            classifications = [TM_LABEL if ins else ("Side_A" if di > 0 else "Side_B")
                               for ins, di in zip(inside, d)]
            membrane_score = float(weights[inside].mean()) if inside.any() else 0.0
            axis = placed.normal
            labeler = f"{self.LABELER} (membrane from file)"
        else:
            half = params.membrane_thickness / 2.0
            fit_coords, fit_weights = coords, fit_w
            context = kwargs.get("context_frames") or []
            if context:
                fit_coords = np.vstack([coords] + [np.asarray(o.coords, dtype=float)
                                                   for o in context])
                fit_names = list(frame.names) + [nm for o in context for nm in o.names]
                fit_weights = _membrane_weights(fit_names, window=SLAB_FIT_WINDOW,
                                                breaks=_chain_breaks(fit_coords))
                labeler = f"{self.LABELER} (fitted on {1 + len(context)} chains)"
            classifications, d, membrane_score, axis, _ = _classify_by_slab(
                coords, weights, thickness=params.membrane_thickness,
                fit_coords=fit_coords, fit_weights=fit_weights)
            frac_in = float(np.mean(np.abs(d) <= half))
            if not context and frac_in >= ENGULF_FRACTION:
                warnings.append(f"{frac_in:.0%} of the chain lies inside the fitted slab - the "
                                "chain may be narrower than the bilayer and the membrane "
                                "orientation unreliable; fit on the whole assembly")
        normal = [round(float(a), 4) for a in axis]

        if membrane_score < params.min_membrane_score:
            labels = _apply_positive_inside_rule(_sides_by_sign(d), residues_data)
            return TMPrediction(
                labels=labels, segments=[], membrane_score=membrane_score,
                membrane_normal=normal, labeler=labeler,
                warnings=warnings + [
                    f"membrane_score {membrane_score:.2f} < {params.min_membrane_score} "
                    "- treated as soluble (no TM segments)"])

        classifications = _smooth_flickers(classifications, d, half, breaks=breaks)
        classifications = _split_at_turns(classifications, d, half, breaks=breaks)
        classifications = _drop_short_tm(classifications, d=d, breaks=breaks)
        classifications, n_dropped = _drop_non_crossing(classifications, d, half, raw_kd,
                                                        breaks=breaks)
        if n_dropped:
            warnings.append(f"{n_dropped} in-slab run(s) dropped: they do not cross the "
                            "membrane or are strongly hydrophilic")
        labels = _apply_positive_inside_rule(classifications, residues_data)
        segments = _tm_runs(labels, breaks=breaks)
        return TMPrediction(
            regions=L.build_regions(frame, labels, _groups_for_runs(n, segments)),
            labels=labels,
            segments=segments,
            membrane_normal=normal,
            membrane_score=membrane_score,
            labeler=labeler,
            warnings=warnings,
        )