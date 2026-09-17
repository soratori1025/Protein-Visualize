import numpy as np
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

def _smooth_hydrophobicity(sequence, window: int = 19) -> np.ndarray:
    scores = [KYTE_DOOLITTLE.get(res, 0.0) for res in sequence]
    half = window // 2
    out = []
    for i in range(len(scores)):
        lo, hi = max(0, i - half), min(len(scores), i + half + 1)
        out.append(float(np.mean(scores[lo:hi])))
    return np.array(out)

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
    hydro = _smooth_hydrophobicity(sequence, window=19)
    
    classifications = ["Side_A"] * len(sequence)
    
    # Find TM regions (simple peak finding)
    tm_regions = []
    in_tm = False
    start_tm = 0
    threshold = 1.0 # KD threshold for TM
    
    for i, h in enumerate(hydro):
        if h > threshold:
            if not in_tm:
                in_tm = True
                start_tm = i
        else:
            if in_tm:
                in_tm = False
                end_tm = i - 1
                if end_tm - start_tm >= 10: # Minimum TM length
                    tm_regions.append((start_tm, end_tm))

    for start, end in tm_regions:
        for i in range(start, end + 1):
            classifications[i] = "Transmembrane"
            
    # Assign loops to Side_A and Side_B alternatingly
    current_side = "Side_A"
    for i in range(len(sequence)):
        if classifications[i] == "Transmembrane":
            if i > 0 and classifications[i-1] != "Transmembrane":
                current_side = "Side_B" if current_side == "Side_A" else "Side_A"
        else:
            classifications[i] = current_side

    # Positive-inside rule
    def _pos_density(side):
        idx = [i for i, c in enumerate(classifications) if c == side]
        if not idx:
            return 0.0
        pos = sum(1 for i in idx if sequence[i] in POSITIVE_RESIDUES)
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
        "protein_name": "Sequence-based Prediction (Kyte-Doolittle HMM)",
        "gene_name": "",
        "organism": "Computed",
        "regions": regions,
    }
