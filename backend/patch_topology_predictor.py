import sys
import re

with open("app/api/topology_predictor.py", "r", encoding="utf-8") as f:
    content = f.read()

content = content.replace("def predict_topology(file_path: Path):", "def predict_topology_kd(file_path: Path):")

import_statements = """
from protein_engine.secondary_structure.dssp import DSSPMethod
"""

# Insert imports at the top
content = content.replace("import numpy as np", "import numpy as np\n" + import_statements)

dssp_algorithm = """
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
        z_start = z[h[0]] - center
        z_end = z[h[-1]] - center
        # Check if helix spans across the membrane center
        if (z_start < -3 and z_end > 3) or (z_start > 3 and z_end < -3):
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
"""

content = content + "\n\n" + dssp_algorithm

with open("app/api/topology_predictor.py", "w", encoding="utf-8") as f:
    f.write(content)
print("Updated topology_predictor.py successfully")
