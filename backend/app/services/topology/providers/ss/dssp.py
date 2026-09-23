from pathlib import Path
from typing import Dict
from app.services.topology.providers.ss.base import SSProvider
from protein_engine.secondary_structure.dssp import DSSPMethod

class DSSPProvider(SSProvider):
    def get_secondary_structure(self, file_path: Path) -> Dict[int, str]:
        ss_res = DSSPMethod().assign(file_path).to_dict()
        ss_map = {}
        for r in ss_res.get('residues', []):
            try:
                res_id = int(r['residue_number'])
                ss_map[res_id] = r['code']
            except ValueError:
                pass # Ignore non-integer IDs or unparseable
        return ss_map
