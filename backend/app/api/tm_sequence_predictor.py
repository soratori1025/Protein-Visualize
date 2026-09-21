import numpy as np
from pathlib import Path
import warnings
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from protein_engine.secondary_structure.dssp import DSSPMethod

KYTE_DOOLITTLE = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
    'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
    'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
    'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
}
POSITIVE_RESIDUES = {"LYS", "ARG"}

TM_WINDOW = 19

# --- Confident pass: pure hydrophobicity, no DSSP involved. This is the primary,
# high-precision signal. ---
TM_HYDRO_THRESHOLD = 1.6
TM_MIN_LENGTH = 15
TM_MERGE_GAP = 3

# --- Rescue pass: catches real TM helices/strands whose average hydrophobicity is
# too weak to clear the confident threshold (common for transporters/GPCRs with
# polar residues lining a pore), but ONLY when DSSP shows a long, uninterrupted
# helix or strand there. Hydrophobicity still gates it (RESCUE_HYDRO_THRESHOLD);
# DSSP is corroborating evidence, not the primary decision. ---
RESCUE_HYDRO_THRESHOLD = 0.6
RESCUE_MIN_ALPHA_LENGTH = 12
RESCUE_MIN_BETA_LENGTH = 5


def _smooth_hydrophobicity(sequence, window: int = TM_WINDOW) -> np.ndarray:
    scores = [KYTE_DOOLITTLE.get(res, 0.0) for res in sequence]
    half = window // 2
    out = []
    for i in range(len(scores)):
        lo, hi = max(0, i - half), min(len(scores), i + half + 1)
        out.append(float(np.mean(scores[lo:hi])))
    return np.array(out)


def _find_hydrophobic_segments(hydro, threshold, min_length, merge_gap=TM_MERGE_GAP):
    """Generic sliding-threshold peak finder over a hydropathy trace."""
    raw_segments = []
    in_seg = False
    start = 0
    for i, h in enumerate(hydro):
        if h > threshold:
            if not in_seg:
                in_seg = True
                start = i
        else:
            if in_seg:
                in_seg = False
                raw_segments.append((start, i - 1))
    if in_seg:
        raw_segments.append((start, len(hydro) - 1))

    merged = []
    for seg in raw_segments:
        if merged and seg[0] - merged[-1][1] - 1 <= merge_gap:
            merged[-1] = (merged[-1][0], seg[1])
        else:
            merged.append(seg)

    return [(s, e) for s, e in merged if (e - s + 1) >= min_length]


def _find_rescued_segments(hydro, ss_map, residues_data):
    """
    Second, lower-confidence hydrophobicity pass: only promoted to a TM call when
    DSSP confirms a long, uninterrupted helix or strand overlapping it. This is what
    catches real TM segments the confident pass misses because part of the helix is
    lined with polar/charged residues (pores, substrate-binding transporters, GPCRs).
    """
    n = len(residues_data)
    runs = []  # contiguous DSSP alpha/beta runs: (start, end, type)
    current_type = None
    run_start = None

    for i in range(n):
        code = ss_map.get(residues_data[i]['id'], 'C')
        this_type = 'Alpha' if code in ('H', 'G', 'I') else ('Beta' if code in ('E', 'B') else None)
        if this_type == current_type and this_type is not None:
            continue
        if current_type is not None:
            runs.append((run_start, i - 1, current_type))
        current_type = this_type
        run_start = i
    if current_type is not None:
        runs.append((run_start, n - 1, current_type))

    rescued = []
    for start, end, stype in runs:
        length = end - start + 1
        min_len = RESCUE_MIN_ALPHA_LENGTH if stype == 'Alpha' else RESCUE_MIN_BETA_LENGTH
        if length < min_len:
            continue
        avg_hydro = float(np.mean(hydro[start:end + 1]))
        if avg_hydro > RESCUE_HYDRO_THRESHOLD:
            rescued.append((start, end))

    return rescued


def _union_segments(*segment_lists, merge_gap=0):
    """Union multiple (start, end) segment lists, merging anything that overlaps
    or touches within merge_gap."""
    all_segs = sorted(seg for segs in segment_lists for seg in segs)
    merged = []
    for seg in all_segs:
        if merged and seg[0] - merged[-1][1] - 1 <= merge_gap:
            merged[-1] = (merged[-1][0], max(merged[-1][1], seg[1]))
        else:
            merged.append(seg)
    return merged


def _refine_with_dssp(tm_segments, ss_map, residues_data):
    """
    Descriptive/refining only, never decides TM vs non-TM:
      (a) classify each segment as Alpha or Beta (majority DSSP code within it)
      (b) snap its boundaries out to the edges of the overlapping SS element.
    """
    refined = []
    n = len(residues_data)

    for start, end in tm_segments:
        codes_in_segment = [ss_map.get(residues_data[i]['id'], 'C') for i in range(start, end + 1)]
        alpha_count = sum(1 for c in codes_in_segment if c in ('H', 'G', 'I'))
        beta_count = sum(1 for c in codes_in_segment if c in ('E', 'B'))

        if alpha_count == 0 and beta_count == 0:
            refined.append((start, end, None))
            continue

        stype = 'Alpha' if alpha_count >= beta_count else 'Beta'
        target_codes = ('H', 'G', 'I') if stype == 'Alpha' else ('E', 'B')

        new_start = start
        while new_start > 0 and ss_map.get(residues_data[new_start - 1]['id'], 'C') in target_codes:
            new_start -= 1

        new_end = end
        while new_end < n - 1 and ss_map.get(residues_data[new_end + 1]['id'], 'C') in target_codes:
            new_end += 1

        refined.append((new_start, new_end, stype))

    return refined


def _merge_typed_overlaps(typed_segments):
    typed_segments = sorted(typed_segments, key=lambda t: t[0])
    merged = []
    for seg in typed_segments:
        if merged and seg[0] <= merged[-1][1] + 1:
            prev = merged[-1]
            merged[-1] = (prev[0], max(prev[1], seg[1]), prev[2] or seg[2])
        else:
            merged.append(seg)
    return merged


def predict_topology_from_sequence(file_path: Path):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)
        parser = MMCIFParser() if file_path.suffix.lower() in ('.cif', '.mmcif') else PDBParser()
        structure = parser.get_structure('protein', str(file_path))

    model = next(iter(structure))
    chain = next(iter(model))

    residues_data = []
    for residue in chain:
        if 'CA' in residue and residue.get_resname() in KYTE_DOOLITTLE:
            residues_data.append({'id': residue.get_id()[1], 'name': residue.get_resname()})

    if not residues_data:
        return {"uniprot_id": "PREDICTED", "regions": []}

    sequence = [r['name'] for r in residues_data]

    # === STEP 1: hydrophobicity (Kyte-Doolittle / TMHMM-style) decides WHERE TM regions are ===
    hydro = _smooth_hydrophobicity(sequence)
    confident_segments = _find_hydrophobic_segments(hydro, TM_HYDRO_THRESHOLD, TM_MIN_LENGTH)

    ss_map = {}
    dssp_success = False
    try:
        ss_res = DSSPMethod().assign(file_path).to_dict()
        for r in ss_res.get('residues', []):
            ss_map[r['residue_number']] = r['code']
        dssp_success = True
    except Exception as e:
        print(f"DSSP failed for sequence predictor. Falling back to pure hydrophobicity TM calls. Error: {e}")

    # === STEP 1b: rescue segments that are real TM but too weakly hydrophobic to
    # clear the confident threshold on their own — only counted with DSSP backing ===
    rescued_segments = []
    if dssp_success:
        rescued_segments = _find_rescued_segments(hydro, ss_map, residues_data)

    tm_segments = _union_segments(confident_segments, rescued_segments)

    print(
        f"[tm_sequence_predictor] confident={len(confident_segments)} "
        f"rescued={len(rescued_segments)} union={len(tm_segments)} "
        f"dssp_success={dssp_success} residues={len(residues_data)}"
    )

    # === STEP 2: DSSP only labels/refines segments already found above ===
    if dssp_success and tm_segments:
        tm_segments_typed = _refine_with_dssp(tm_segments, ss_map, residues_data)
        tm_segments_typed = _merge_typed_overlaps(tm_segments_typed)
    else:
        tm_segments_typed = [(s, e, 'Alpha') for s, e in tm_segments]

    # --- Build per-residue classification from the final TM segments ---
    classifications = ["Side_A"] * len(sequence)
    region_details = ["Side_A"] * len(sequence)

    for start, end, stype in tm_segments_typed:
        label = ("Transmembrane Alpha Helix" if stype == 'Alpha'
                 else "Transmembrane Beta Strand" if stype == 'Beta'
                 else "Transmembrane")
        for i in range(start, end + 1):
            classifications[i] = "Transmembrane"
            region_details[i] = label

    # --- Assign loops to Side_A / Side_B alternatingly across each TM crossing ---
    current_side = "Side_A"
    for i in range(len(sequence)):
        if classifications[i] == "Transmembrane":
            if i > 0 and classifications[i - 1] != "Transmembrane":
                current_side = "Side_B" if current_side == "Side_A" else "Side_A"
        else:
            classifications[i] = current_side
            region_details[i] = current_side

    # --- Positive-inside rule ---
    def _pos_density(side):
        idx = [i for i, c in enumerate(classifications) if c == side]
        if not idx:
            return 0.0
        pos = sum(1 for i in idx if sequence[i] in POSITIVE_RESIDUES)
        return pos / len(idx)

    cyto_side = "Side_A" if _pos_density("Side_A") >= _pos_density("Side_B") else "Side_B"
    side_label = {cyto_side: "Cytoplasmic",
                  ("Side_B" if cyto_side == "Side_A" else "Side_A"): "Extracellular"}

    for i in range(len(region_details)):
        if region_details[i] in ["Side_A", "Side_B"]:
            region_details[i] = side_label.get(region_details[i], region_details[i])

    # --- Collapse into contiguous regions for output ---
    regions, current_type, start_id = [], region_details[0], residues_data[0]['id']
    for i in range(1, len(region_details)):
        if region_details[i] != current_type:
            end_id = residues_data[i - 1]['id']
            base_type = "Transmembrane" if "Transmembrane" in current_type else "Topological domain"
            regions.append({"type": base_type,
                             "start": start_id, "end": end_id, "description": current_type})
            current_type, start_id = region_details[i], residues_data[i]['id']

    end_id = residues_data[-1]['id']
    base_type = "Transmembrane" if "Transmembrane" in current_type else "Topological domain"
    regions.append({"type": base_type,
                     "start": start_id, "end": end_id, "description": current_type})

    return {
        "uniprot_id": "CALCULATED",
        "protein_name": "Prediction (Kyte-Doolittle TM scan, DSSP-refined + rescue)",
        "gene_name": "",
        "organism": "Computed",
        "regions": regions,
    }