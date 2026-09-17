import numpy as np

from protein_engine.secondary_structure.dssp import DSSPMethod

from pathlib import Path
import warnings
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

KYTE_DOOLITTLE = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
    'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
    'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
    'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
}
POSITIVE_RESIDUES = {"LYS", "ARG"}
MEMBRANE_THICKNESS = 30.0   # Å, bề dày lõi hydrophobic điển hình của bilayer
N_AXIS_SAMPLES = 300        # số hướng ứng viên

def _fibonacci_sphere(n: int) -> np.ndarray:
    """Sinh n điểm gần-đều trên mặt cầu đơn vị (dùng làm trục ứng viên)."""
    phi = np.pi * (3.0 - np.sqrt(5.0))
    pts = []
    for i in range(n):
        y = 1 - (i / float(n - 1)) * 2
        r = np.sqrt(max(0.0, 1 - y * y))
        theta = phi * i
        pts.append((np.cos(theta) * r, y, np.sin(theta) * r))
    return np.array(pts)

def _smooth_hydrophobicity(residues, window: int = 19) -> np.ndarray:
    scores = [KYTE_DOOLITTLE.get(r['name'], 0.0) for r in residues]
    half = window // 2
    out = []
    for i in range(len(scores)):
        lo, hi = max(0, i - half), min(len(scores), i + half + 1)
        out.append(float(np.mean(scores[lo:hi])))
    return np.array(out)

def _best_slab(coords: np.ndarray, weights: np.ndarray,
               thickness: float = MEMBRANE_THICKNESS,
               n_axes: int = N_AXIS_SAMPLES) -> tuple[np.ndarray, float, float]:
    axes = _fibonacci_sphere(n_axes)
    centroid = coords.mean(axis=0)
    centered = coords - centroid

    best_score = -np.inf
    best: tuple[np.ndarray, float] | None = None
    for axis in axes:
        z = centered @ axis
        order = np.argsort(z)
        z_sorted, w_sorted = z[order], weights[order]
        cum = np.concatenate([[0.0], np.cumsum(w_sorted)])

        lo = 0
        for hi in range(len(z_sorted)):
            while z_sorted[hi] - z_sorted[lo] > thickness:
                lo += 1
            inside = cum[hi + 1] - cum[lo]
            if inside > best_score:
                best_score = inside
                best = (axis, 0.5 * (z_sorted[lo] + z_sorted[hi]))

    assert best is not None, "no axis candidate produced a slab — coords/weights empty?"
    return best[0], best[1], best_score

def predict_topology_kd(file_path: Path):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)
        parser = MMCIFParser() if file_path.suffix.lower() in ('.cif', '.mmcif') else PDBParser()
        structure = parser.get_structure('protein', str(file_path))

    model = next(iter(structure))
    chain = next(iter(model))

    residues_data, ca_coords = [], []
    for residue in chain:
        if 'CA' in residue and residue.get_resname() in KYTE_DOOLITTLE:
            ca_coords.append(residue['CA'].get_coord())
            residues_data.append({'id': residue.get_id()[1], 'name': residue.get_resname()})

    if len(ca_coords) < 10:
        return {"uniprot_id": "PREDICTED", "regions": []}

    coords = np.array(ca_coords)
    hydro = _smooth_hydrophobicity(residues_data, window=19)

    axis, center, _ = _best_slab(coords, hydro)
    centroid = coords.mean(axis=0)
    z = (coords - centroid) @ axis
    half = MEMBRANE_THICKNESS / 2.0

    classifications = []
    for i in range(len(residues_data)):
        if abs(z[i] - center) <= half and hydro[i] > np.median(hydro):
            classifications.append("Transmembrane")
        elif z[i] - center > half:
            classifications.append("Side_A")
        else:
            classifications.append("Side_B")

    # positive-inside rule: bên nào giàu Lys/Arg hơn -> Cytoplasmic
    def _pos_density(side):
        idx = [i for i, c in enumerate(classifications) if c == side]
        if not idx:
            return 0.0
        pos = sum(1 for i in idx if residues_data[i]['name'] in POSITIVE_RESIDUES)
        return pos / len(idx)

    cyto_side = "Side_A" if _pos_density("Side_A") >= _pos_density("Side_B") else "Side_B"
    side_label = {cyto_side: "Cytoplasmic",
                  ("Side_B" if cyto_side == "Side_A" else "Side_A"): "Extracellular"}
    classifications = [side_label.get(c, c) for c in classifications]

    # --- phần merge/filter region: giữ nguyên logic bạn đã có ---
    regions, current_type, start_id = [], classifications[0], residues_data[0]['id']
    for i in range(1, len(classifications)):
        if classifications[i] != current_type:
            end_id = residues_data[i - 1]['id']
            desc = "Calculated Helical TM region" if current_type == "Transmembrane" else current_type
            regions.append({"type": "Transmembrane" if current_type == "Transmembrane" else "Topological domain",
                             "start": start_id, "end": end_id, "description": desc})
            current_type, start_id = classifications[i], residues_data[i]['id']
    end_id = residues_data[-1]['id']
    desc = "Calculated Helical TM region" if current_type == "Transmembrane" else current_type
    regions.append({"type": "Transmembrane" if current_type == "Transmembrane" else "Topological domain",
                     "start": start_id, "end": end_id, "description": desc})

    filtered = []
    for r in regions:
        if r['type'] == 'Transmembrane' and (r['end'] - r['start']) < 10:
            r['type'], r['description'] = 'Topological domain', 'Loop'
        filtered.append(r)

    merged = []
    if filtered:
        cur = filtered[0].copy()
        for r in filtered[1:]:
            if cur['type'] == r['type'] and cur['description'] == r['description']:
                cur['end'] = r['end']
            else:
                merged.append(cur)
                cur = r.copy()
        merged.append(cur)

    return {
        "uniprot_id": "CALCULATED",
        "protein_name": "Structure-based Prediction (Slab-fit + Kyte-Doolittle + positive-inside rule)",
        "gene_name": "",
        "organism": "Computed",
        "regions": merged,
    }


def predict_topology_dssp(file_path: Path):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)
        parser = MMCIFParser() if file_path.suffix.lower() in ('.cif', '.mmcif') else PDBParser()
        structure = parser.get_structure('protein', str(file_path))

    model = next(iter(structure))
    chain = next(iter(model))

    residues_data, ca_coords = [], []
    for residue in chain:
        if 'CA' in residue and residue.get_resname() in KYTE_DOOLITTLE:
            ca_coords.append(residue['CA'].get_coord())
            residues_data.append({'id': residue.get_id()[1], 'name': residue.get_resname()})

    if len(ca_coords) < 10:
        return {"uniprot_id": "PREDICTED", "regions": []}

    coords = np.array(ca_coords)
    hydro = _smooth_hydrophobicity(residues_data, window=19)

    axis, center, _ = _best_slab(coords, hydro)
    centroid = coords.mean(axis=0)
    z = (coords - centroid) @ axis

    # Run DSSP
    try:
        ss_res = DSSPMethod().assign(file_path).to_dict()
    except Exception as e:
        print(f"DSSP failed, falling back to KD. Error: {e}")
        return predict_topology_kd(file_path)
    
    ss_map = {}
    for r in ss_res.get('residues', []):
        ss_map[r['residue_number']] = r['code']

    classifications = ["Side_A"] * len(residues_data)
    
    # Extract helices
    helices = []
    current_helix = []
    for i, r in enumerate(residues_data):
        code = ss_map.get(r['id'], 'C')
        if code in ['H', 'G', 'I']:
            current_helix.append(i)
        else:
            if current_helix:
                helices.append(current_helix)
                current_helix = []
    if current_helix:
        helices.append(current_helix)

    # Classify helices
    half = MEMBRANE_THICKNESS / 2.0
    for h in helices:
        z_vals = [z[idx] - center for idx in h]
        z_min = min(z_vals)
        z_max = max(z_vals)
        # Check if helix spans a significant portion of the hydrophobic core
        if z_min < 5 and z_max > -5 and (z_max - z_min) > 10:
            for idx in h:
                classifications[idx] = "Transmembrane"

    # Classify remaining residues as Side_A or Side_B based on their Z coordinate relative to center
    for i in range(len(residues_data)):
        if classifications[i] != "Transmembrane":
            if z[i] - center > 0:
                classifications[i] = "Side_A"
            else:
                classifications[i] = "Side_B"

    # Positive-inside rule
    def _pos_density(side):
        idx = [i for i, c in enumerate(classifications) if c == side]
        if not idx:
            return 0.0
        pos = sum(1 for i in idx if residues_data[i]['name'] in POSITIVE_RESIDUES)
        return pos / len(idx)

    cyto_side = "Side_A" if _pos_density("Side_A") >= _pos_density("Side_B") else "Side_B"
    side_label = {cyto_side: "Cytoplasmic",
                  ("Side_B" if cyto_side == "Side_A" else "Side_A"): "Extracellular"}
    classifications = [side_label.get(c, c) for c in classifications]

    regions, current_type, start_id = [], classifications[0], residues_data[0]['id']
    for i in range(1, len(classifications)):
        if classifications[i] != current_type:
            end_id = residues_data[i - 1]['id']
            desc = "Calculated Helical TM region" if current_type == "Transmembrane" else current_type
            regions.append({"type": "Transmembrane" if current_type == "Transmembrane" else "Topological domain",
                             "start": start_id, "end": end_id, "description": desc})
            current_type, start_id = classifications[i], residues_data[i]['id']
    end_id = residues_data[-1]['id']
    desc = "Calculated Helical TM region" if current_type == "Transmembrane" else current_type
    regions.append({"type": "Transmembrane" if current_type == "Transmembrane" else "Topological domain",
                     "start": start_id, "end": end_id, "description": desc})

    return {
        "uniprot_id": "CALCULATED",
        "protein_name": "Structure-based Prediction (DSSP + Slab-fit)",
        "gene_name": "",
        "organism": "Computed",
        "regions": regions,
    }

def predict_topology(file_path: Path, algorithm: str = "dssp_slab"):
    if algorithm == "kd_slab":
        return predict_topology_kd(file_path)
    return predict_topology_dssp(file_path)
