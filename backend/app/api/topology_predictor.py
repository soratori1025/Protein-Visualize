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
    outward to the true ends of the overlapping SS element. Snapping is capped at
    half the gap to the neighbouring TM run, so two distinct crossings can never be
    merged. DSSP/STRIDE never adds or removes a TM call.

Sides are named Cytoplasmic / Extracellular geometrically (which face of the slab
a loop sits on) plus the positive-inside rule - NOT by topological alternation.
This keeps re-entrant / half-membrane loops correct, which matters for the
transporter/GPCR families this tool targets.

Design notes / known limitations (state these in the paper):
  * Only the first chain of the first model is analysed. For an oligomeric bundle,
    fit per chain or merge chains upstream before calling this.
  * The hydrophobic-core thickness is a fixed parameter (default 30 Angstrom); it
    is not optimised per structure.
  * `membrane_score` is a heuristic confidence (mean hydrophobic weight inside the
    slab). Non-membrane proteins are reported with regions=[] when it falls below
    MIN_MEMBRANE_SCORE; expose the score so a caller can apply its own threshold.

Heavy imports (Bio.PDB, the DSSP/STRIDE adapters) are loaded lazily inside the
functions that need them, so the numeric core can be imported and unit-tested with
numpy alone.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np

# --- Kyte-Doolittle hydropathy (higher = more hydrophobic) --------------------
KYTE_DOOLITTLE = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
    "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
    "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
    "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}
POSITIVE_RESIDUES = {"LYS", "ARG"}  # positive-inside rule

# --- Secondary-structure code sets (shared by DSSP and STRIDE) ----------------
HELIX_CODES = {"H", "G", "I"}          # alpha, 3-10, pi
STRAND_CODES = {"E", "B", "b"}         # extended strand / beta bridge

# --- Membrane geometry parameters ---------------------------------------------
MEMBRANE_THICKNESS = 30.0   # Angstrom, typical hydrophobic-core thickness
N_AXIS_SAMPLES = 500        # candidate membrane-normal directions
N_CENTER_SAMPLES = 160      # slab positions scanned along each candidate normal
SMOOTH_WINDOW = 19          # residues, hydropathy smoothing window (~1 TM helix)
MAX_JITTER_LEN = 3          # max length of a slab-edge jitter dip that may be merged
JITTER_MARGIN = 3.0         # Angstrom a jitter residue may sit beyond the slab face
MIN_TM_CORE = 5             # drop TM runs whose in-slab core is shorter than this
MIN_FACE_RESIDUES = 5       # min residues required OUTSIDE the slab on each face
MIN_FACE_FRACTION = 0.08    # ...and each face must hold >= this fraction of residues
                            # (balance: rejects a slab shoved to one extreme, which
                            #  is how a single hydrophilic sliver gamed the contrast)
MAX_SNAP = 4                # max residues a boundary may snap past the slab edge
MIN_MEMBRANE_SCORE = 0.5    # below this the structure is treated as non-membrane
MIN_EXTRA_SS_LEN = 3        # extramembrane helix/strand shorter than this -> Coil

# --- SS-element-first pipeline (predict_topology_ss_first) ---------------------
MIN_TM_ELEMENT_IN_SLAB = 4  # a DSSP element needs this many residues in the slab
FULL_CROSS_FRAC = 0.66      # element whose z-span >= this*thickness crosses fully
BROKEN_GAP_MAX = 9          # max intramembrane break between the two halves of one
                            # discontinuous helix (TM1a/1b, TM6a/6b in a LeuT fold)
MIN_CROSS_SPAN_FRAC = 0.45  # a crossing must reach across at least this*thickness of
                            # the bilayer; short helices that sit in one leaflet
                            # (interfacial / re-entrant) are not TM crossings


# =============================================================================
# User-configurable parameter bundle
# =============================================================================
@dataclass
class TMParams:
    """All biology-dependent thresholds that affect TM detection.

    Each field carries a sensible default derived from literature.  The values
    can be overridden per-call via the API query string so that users can tune
    the predictor to their specific membrane type or protein family without
    editing source code.

    References for default values:
      * membrane_thickness: Mitra et al., Biochemistry 2004; OPM database
        (Lomize et al., Proteins 2006). Typical range 25–32 Å.
      * min_membrane_score: heuristic; calibrated on OPM soluble-vs-membrane
        classification.  Expose the score so callers can apply their own cut.
      * Kyte-Doolittle scale: Kyte & Doolittle, J. Mol. Biol. 157:105–132, 1982.
      * Positive-inside rule: Von Heijne, J. Mol. Biol. 225:487–494, 1992.
    """
    # --- Biology-dependent thresholds (user-tunable) --------------------------
    membrane_thickness: float = MEMBRANE_THICKNESS
    min_tm_element_in_slab: int = MIN_TM_ELEMENT_IN_SLAB
    min_cross_span_frac: float = MIN_CROSS_SPAN_FRAC
    full_cross_frac: float = FULL_CROSS_FRAC
    broken_gap_max: int = BROKEN_GAP_MAX
    min_membrane_score: float = MIN_MEMBRANE_SCORE

    def to_response_dict(self) -> dict:
        """Serialise into a dict suitable for inclusion in the API response
        under the key ``parameters_used``, so users / papers know exactly
        which values produced the result."""
        return {
            "membrane_thickness_angstrom": self.membrane_thickness,
            "min_tm_element_in_slab": self.min_tm_element_in_slab,
            "min_cross_span_frac": self.min_cross_span_frac,
            "full_cross_frac": self.full_cross_frac,
            "broken_gap_max": self.broken_gap_max,
            "min_membrane_score": self.min_membrane_score,
        }


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


def _drop_short_tm(classifications, min_core: int = MIN_TM_CORE):
    """Demote TM runs whose in-slab core is shorter than `min_core` back to the
    neighbouring side (spurious slab dips by a terminal tail or a sharp turn)."""
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
            side = out[i - 1] if i > 0 and out[i - 1] != "Transmembrane" else (
                out[j] if j < n and out[j] != "Transmembrane" else "Side_A")
            for k in range(i, j):
                out[k] = side
        i = j
    return out


# =============================================================================
# STAGE 2 helpers - DSSP / STRIDE labelling (descriptive only)
# =============================================================================
def _run_labeler(file_path: Path, labeler: str) -> dict:
    """Return {residue_number: ss_code} from DSSP or STRIDE. Imported lazily."""
    name = (labeler or "DSSP").upper()
    if name == "STRIDE":
        from protein_engine.secondary_structure.stride import STRIDEMethod
        result = STRIDEMethod().assign(file_path).to_dict()
    else:
        from protein_engine.secondary_structure.dssp import DSSPMethod
        result = DSSPMethod().assign(file_path).to_dict()
    # First chain only, matching _load_structure_residues.
    ss_map = {}
    for r in result.get("residues", []):
        ss_map[r["residue_number"]] = r["code"]
    return ss_map


def _classify_code(code: str):
    if code in HELIX_CODES:
        return "Alpha"
    if code in STRAND_CODES:
        return "Beta"
    return None


def _ss_label(code: str) -> str:
    """DSSP/STRIDE code -> coarse SS class for drawing (Helix / Strand / Coil)."""
    if code in HELIX_CODES:
        return "Helix"
    if code in STRAND_CODES:
        return "Strand"
    return "Coil"


def _label_extramembrane(descriptions, residues_data, ss_map,
                         min_ss_len: int = MIN_EXTRA_SS_LEN):
    """Append secondary-structure detail to the extramembrane (Cytoplasmic /
    Extracellular) residues, so the loops of a topology snake-plot can show their
    own helices and strands (e.g. EL2, EL3a/3b in a transporter diagram) instead of
    a bare line. TM labels and any un-renamed sides pass through untouched.

    Each side residue's description becomes "<side> <Helix|Strand|Coil>", which
    `_build_regions` then splits into contiguous drawable elements. Short helix/
    strand specks (< min_ss_len) are demoted to Coil so a 1-2 residue DSSP flicker
    in a loop does not become a cylinder.
    """
    n = len(descriptions)
    out = list(descriptions)

    raw: list[str | None] = [None] * n
    for i in range(n):
        if out[i] in ("Cytoplasmic", "Extracellular"):
            raw[i] = _ss_label(ss_map.get(residues_data[i]["id"], "C"))

    # demote short helix/strand runs (within one side) to Coil
    i = 0
    while i < n:
        if raw[i] in ("Helix", "Strand"):
            j = i
            while j < n and raw[j] == raw[i]:
                j += 1
            if j - i < min_ss_len:
                for k in range(i, j):
                    raw[k] = "Coil"
            i = j
        else:
            i += 1

    for i in range(n):
        if raw[i] is not None:
            out[i] = f"{out[i]} {raw[i]}"
    return out


def _region_side_ss(description: str):
    """Derive explicit (side, ss) tags from a region description so the frontend
    needs no string parsing. side in {membrane, Cytoplasmic, Extracellular};
    ss in {Helix, Strand, Loop, Coil, None}."""
    if "Transmembrane" in description:
        side = "membrane"
    elif "Cytoplasmic" in description:
        side = "Cytoplasmic"
    elif "Extracellular" in description:
        side = "Extracellular"
    else:
        side = None
    if "Helix" in description:
        ss = "Helix"
    elif "Strand" in description:
        ss = "Strand"
    elif "Loop" in description:
        ss = "Loop"
    elif "Coil" in description:
        ss = "Coil"
    else:
        ss = None
    return side, ss


def _label_and_snap(classifications, residues_data, ss_map, max_snap: int = MAX_SNAP):
    """Label each TM run Alpha/Beta/Loop by majority SS code and snap boundaries
    outward along the same SS element, capped at half the gap to the neighbouring
    TM run so distinct crossings never merge."""
    n = len(classifications)
    out = list(classifications)

    # Collect TM run intervals [start, end] inclusive.
    runs = []
    i = 0
    while i < n:
        if classifications[i] != "Transmembrane":
            i += 1
            continue
        j = i
        while j < n and classifications[j] == "Transmembrane":
            j += 1
        runs.append((i, j - 1))
        i = j

    for idx, (start, end) in enumerate(runs):
        codes = [ss_map.get(residues_data[k]["id"], "C") for k in range(start, end + 1)]
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
            left_budget = min(max_snap, (start - prev_end - 1) // 2)
            right_budget = min(max_snap, (next_start - end - 1) // 2)

            moved = 0
            while (s > 0 and moved < left_budget
                   and ss_map.get(residues_data[s - 1]["id"], "C") in target):
                s -= 1
                moved += 1
            moved = 0
            while (e < n - 1 and moved < right_budget
                   and ss_map.get(residues_data[e + 1]["id"], "C") in target):
                e += 1
                moved += 1

        for k in range(s, e + 1):
            out[k] = label
    return out


def _apply_positive_inside_rule(descriptions, residues_data):
    """Rename Side_A/Side_B to Cytoplasmic/Extracellular by the positive-inside
    rule (the cytoplasmic face is enriched in Lys/Arg). TM labels pass through."""
    def _pos_density(side):
        idx = [i for i, c in enumerate(descriptions) if c == side]
        if not idx:
            return 0.0
        pos = sum(1 for i in idx if residues_data[i]["name"] in POSITIVE_RESIDUES)
        return pos / len(idx)

    if "Side_A" not in descriptions and "Side_B" not in descriptions:
        return list(descriptions)

    cyto = "Side_A" if _pos_density("Side_A") >= _pos_density("Side_B") else "Side_B"
    rename = {cyto: "Cytoplasmic",
              ("Side_B" if cyto == "Side_A" else "Side_A"): "Extracellular"}
    return [rename.get(c, c) for c in descriptions]


def _build_regions(residues_data, descriptions):
    def _emit(desc, start_id, end_id):
        base = "Transmembrane" if "Transmembrane" in desc else "Topological domain"
        side, ss = _region_side_ss(desc)
        return {"type": base, "start": start_id, "end": end_id,
                "description": desc, "side": side, "ss": ss}

    regions = []
    current, start_id = descriptions[0], residues_data[0]["id"]
    for i in range(1, len(descriptions)):
        if descriptions[i] != current:
            regions.append(_emit(current, start_id, residues_data[i - 1]["id"]))
            current, start_id = descriptions[i], residues_data[i]["id"]
    regions.append(_emit(current, start_id, residues_data[-1]["id"]))
    return regions


# =============================================================================
# Structure loading (lazy Bio.PDB)
# =============================================================================
def _load_structure_residues(file_path: Path, chain_id: str | None = None):
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        is_cif = file_path.suffix.lower() in (".cif", ".mmcif")
        parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(file_path))

    model = next(iter(structure))
    chain = model[chain_id] if chain_id else next(iter(model))

    residues_data, ca_coords = [], []
    for residue in chain:
        if "CA" in residue and residue.get_resname() in KYTE_DOOLITTLE:
            ca_coords.append(residue["CA"].get_coord())
            residues_data.append({"id": residue.get_id()[1], "name": residue.get_resname()})

    coords = np.array(ca_coords, dtype=float) if ca_coords else np.empty((0, 3))
    return residues_data, coords


# =============================================================================
# Public entry points
# =============================================================================
def predict_topology_structure(file_path: Path, labeler: str = "DSSP",
                               params: TMParams | None = None) -> dict:
    """Unified structure-based prediction: slab geometry gates TM (stage 1),
    DSSP/STRIDE labels alpha/beta/loop (stage 2)."""
    if params is None:
        params = TMParams()
    thickness = params.membrane_thickness

    residues_data, coords = _load_structure_residues(file_path)
    if len(residues_data) < 10:
        return {"uniprot_id": "PREDICTED", "protein_name": "Too few residues",
                "gene_name": "", "organism": "Computed",
                "membrane_score": 0.0, "labeler": labeler,
                "parameters_used": params.to_response_dict(), "regions": []}

    names = [r["name"] for r in residues_data]
    weights = _membrane_weights(names)

    # STAGE 1 - geometry decides TM.
    classifications, d, membrane_score, axis, center = _classify_by_slab(
        coords, weights, thickness)
    classifications = _smooth_flickers(classifications, d, thickness / 2.0)
    classifications = _drop_short_tm(classifications)

    # Guard: if the best slab is not convincingly hydrophobic, treat as soluble.
    if membrane_score < params.min_membrane_score:
        descriptions = _apply_positive_inside_rule(
            ["Side_A" if c == "Transmembrane" else c for c in classifications],
            residues_data)
        return {
            "uniprot_id": "CALCULATED",
            "protein_name": "No membrane slab detected (likely soluble)",
            "gene_name": "", "organism": "Computed",
            "membrane_score": round(membrane_score, 3), "labeler": labeler,
            "parameters_used": params.to_response_dict(),
            "regions": _build_regions(residues_data, descriptions),
        }

    # STAGE 2 - DSSP/STRIDE labels the already-fixed TM runs.
    ss_map, labeler_used = {}, labeler
    if labeler == "__none__":
        labeler_used = "geometry-only (no SS tool)"
    else:
        try:
            ss_map = _run_labeler(file_path, labeler)
        except Exception as error:  # noqa: BLE001 - tool may be missing at runtime
            print(f"[topology_predictor] {labeler} failed, geometry-only labels. {error}")
            labeler_used = f"{labeler} (unavailable)"

    if ss_map:
        descriptions = _label_and_snap(classifications, residues_data, ss_map)
        descriptions = _apply_positive_inside_rule(descriptions, residues_data)
        # Give extramembrane loops their own helix/strand/coil elements to draw.
        descriptions = _label_extramembrane(descriptions, residues_data, ss_map)
    else:
        # No SS tool: keep TM calls, but we cannot tell alpha from beta.
        descriptions = ["Transmembrane (helical, unverified)"
                        if c == "Transmembrane" else c for c in classifications]
        descriptions = _apply_positive_inside_rule(descriptions, residues_data)

    regions = _build_regions(residues_data, descriptions)

    return {
        "uniprot_id": "CALCULATED",
        "protein_name": f"Structure-based TM prediction (slab-fit + {labeler_used})",
        "gene_name": "", "organism": "Computed",
        "membrane_score": round(membrane_score, 3),
        "membrane_normal": [round(float(a), 4) for a in axis],
        "labeler": labeler_used,
        "parameters_used": params.to_response_dict(),
        "regions": regions,
    }


# =============================================================================
# SS-element-first pipeline  (more accurate boundaries on real polytopic proteins)
# =============================================================================
def _helix_axis_normal(coords: np.ndarray, k: int = 4) -> np.ndarray:
    """Estimate the membrane normal from geometry alone: the dominant direction of
    the chain, taken as the top eigenvector of the scatter of local CA(i+k)-CA(i)
    vectors. Transmembrane helices are the longest straight runs so they dominate,
    and the outer-product scatter is sign-invariant, so up- and down-helices
    reinforce the same axis. Unlike a hydrophobicity fit this does not degrade when
    the TM core is weakly hydrophobic (a substrate pore lined with polar residues),
    which is exactly where the slab-contrast fit fails."""
    if len(coords) <= k:
        return np.array([0.0, 0.0, 1.0])
    dirs = coords[k:] - coords[:-k]
    norm = np.linalg.norm(dirs, axis=1, keepdims=True)
    norm[norm == 0] = 1.0
    dirs = dirs / norm
    _, vec = np.linalg.eigh(dirs.T @ dirs)
    axis = vec[:, -1]
    return axis / np.linalg.norm(axis)


def _ss_runs(ss_map: dict, residues_data):
    """Maximal contiguous runs of one coarse SS class over the sequence.
    Returns [(start_idx, end_idx, 'Helix'|'Strand'|'Coil'), ...]."""
    classes = [_ss_label(ss_map.get(r["id"], "C")) for r in residues_data]
    runs = []
    i = 0
    n = len(classes)
    while i < n:
        j = i
        while j < n and classes[j] == classes[i]:
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
                              params: TMParams | None = None) -> dict:
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

    residues_data, coords = _load_structure_residues(file_path)
    if len(residues_data) < 10:
        return {"uniprot_id": "PREDICTED", "protein_name": "Too few residues",
                "gene_name": "", "organism": "Computed",
                "membrane_score": 0.0, "labeler": labeler,
                "parameters_used": params.to_response_dict(), "regions": []}

    ss_map, labeler_used = {}, labeler
    try:
        ss_map = _run_labeler(file_path, labeler)
    except Exception as error:  # noqa: BLE001
        print(f"[topology_predictor] {labeler} failed; falling back to slab-first. {error}")
    if not ss_map:
        return predict_topology_structure(file_path, labeler="__none__", params=params)

    names = [r["name"] for r in residues_data]
    weights = _membrane_weights(names)
    n = len(residues_data)

    normal = _helix_axis_normal(coords)
    centroid = coords.mean(axis=0)
    z = (coords - centroid) @ normal
    half = thickness / 2.0

    runs = _ss_runs(ss_map, residues_data)
    center = _membrane_center_geometric(z, runs, half)
    d = z - center
    inside = np.abs(d) <= half
    mean_in = float(weights[inside].mean()) if inside.any() else 0.0
    crossings = _membrane_crossings(runs, z, center, half, thickness, params)

    if not crossings or mean_in < params.min_membrane_score:
        descriptions = _apply_positive_inside_rule(
            ["Side_A" if di > half else "Side_B" for di in d], residues_data)
        descriptions = _label_extramembrane(descriptions, residues_data, ss_map)
        return {
            "uniprot_id": "CALCULATED",
            "protein_name": "No transmembrane crossings detected (likely soluble)",
            "gene_name": "", "organism": "Computed",
            "membrane_score": round(mean_in, 3), "labeler": labeler_used,
            "parameters_used": params.to_response_dict(),
            "regions": _build_regions(residues_data, descriptions),
        }

    # Per-residue labels: TM crossings first, then sides, then extramembrane SS.
    descriptions = ["Side_A" if di > half else "Side_B" for di in d]
    for a, b, cls in crossings:
        label = "Transmembrane Alpha Helix" if cls == "Helix" else "Transmembrane Beta Strand"
        for k in range(a, b + 1):
            descriptions[k] = label
    descriptions = _apply_positive_inside_rule(descriptions, residues_data)
    descriptions = _label_extramembrane(descriptions, residues_data, ss_map)
    regions = _build_regions(residues_data, descriptions)

    return {
        "uniprot_id": "CALCULATED",
        "protein_name": f"Structure-based TM prediction (SS-element-first + {labeler_used})",
        "gene_name": "", "organism": "Computed",
        "membrane_score": round(mean_in, 3),
        "membrane_normal": [round(float(a), 4) for a in normal],
        "labeler": labeler_used,
        "parameters_used": params.to_response_dict(),
        "regions": regions,
    }


def predict_topology(file_path: Path, algorithm: str = "dssp_slab",
                     params: TMParams | None = None) -> dict:
    """Router-facing dispatcher.

      dssp_ss    -> SS-element-first + DSSP   (recommended, most accurate)
      stride_ss  -> SS-element-first + STRIDE
      dssp_slab  -> slab geometry + DSSP labelling
      stride_slab-> slab geometry + STRIDE labelling
      kd_slab    -> slab geometry only, no SS tool (cannot distinguish beta)
    """
    if params is None:
        params = TMParams()
    if algorithm == "kd_slab":
        return predict_topology_structure(file_path, labeler="__none__", params=params)
    if algorithm == "stride_slab":
        return predict_topology_structure(file_path, labeler="STRIDE", params=params)
    if algorithm == "dssp_slab":
        return predict_topology_structure(file_path, labeler="DSSP", params=params)
    if algorithm == "stride_ss":
        return predict_topology_ss_first(file_path, labeler="STRIDE", params=params)
    # default + "dssp_ss"
    return predict_topology_ss_first(file_path, labeler="DSSP", params=params)