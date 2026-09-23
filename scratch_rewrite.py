import sys

with open("backend/app/services/topology/providers/tm/geometry.py", "r", encoding="utf-8") as f:
    content = f.read()

# I only need up to _drop_short_tm (inclusive)
import re
match = re.search(r'def _assign_labels\(.*?\n\n\n', content, flags=re.DOTALL)
if match:
    core_logic = content[:match.start()]
else:
    core_logic = content

new_class = """
from app.services.topology.providers.tm.base import TMProvider, TMPrediction, TMBoundary
import warnings
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

class GeometryTMProvider(TMProvider):
    def predict_tm(self, file_path: Path, params: TMParams = None, **kwargs) -> TMPrediction:
        if params is None:
            params = TMParams()
            
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            parser = MMCIFParser() if file_path.suffix.lower() in ('.cif', '.mmcif') else PDBParser()
            structure = parser.get_structure('protein', str(file_path))

        model = next(iter(structure))
        chain = next(iter(model))

        residues_data = []
        ca_coords = []
        for residue in chain:
            if 'CA' in residue and residue.get_resname() in KYTE_DOOLITTLE:
                ca_coords.append(residue['CA'].get_coord())
                residues_data.append({'id': residue.get_id()[1], 'name': residue.get_resname()})

        if len(ca_coords) < 10:
            return TMPrediction(boundaries=[])

        coords = np.array(ca_coords)
        weights = _membrane_weights([r['name'] for r in residues_data])

        classifications, d, membrane_score, axis, center = _classify_by_slab(
            coords, weights, thickness=params.membrane_thickness
        )
        
        if membrane_score < params.min_membrane_score:
            return TMPrediction(boundaries=[], membrane_score=membrane_score, membrane_normal=list(axis))

        half = params.membrane_thickness / 2.0
        classifications = _smooth_flickers(classifications, d, half)
        classifications = _drop_short_tm(classifications)

        # Convert classifications ("Transmembrane") to boundaries
        boundaries = []
        in_tm = False
        start_idx = -1
        for i, c in enumerate(classifications):
            if c == "Transmembrane":
                if not in_tm:
                    in_tm = True
                    start_idx = i
            else:
                if in_tm:
                    boundaries.append(TMBoundary(start=residues_data[start_idx]['id'], end=residues_data[i-1]['id']))
                    in_tm = False
        
        if in_tm:
            boundaries.append(TMBoundary(start=residues_data[start_idx]['id'], end=residues_data[-1]['id']))

        return TMPrediction(
            boundaries=boundaries, 
            membrane_normal=list(axis), 
            membrane_score=membrane_score, 
            labeler="3D_Geometry"
        )
"""

final_content = core_logic + "\n" + new_class

with open("backend/app/services/topology/providers/tm/geometry.py", "w", encoding="utf-8") as f:
    f.write(final_content)
