from pathlib import Path
from typing import Dict
from app.services.topology.providers.ss.base import SSProvider
from protein_engine.secondary_structure.stride import STRIDEMethod

class STRIDEProvider(SSProvider):
    def get_secondary_structure(self, file_path: Path) -> Dict[int, str]:
        ss_res = STRIDEMethod().assign(file_path).to_dict()
        ss_map = {}
        for r in ss_res.get('residues', []):
            try:
                res_id = int(r['residue_number'])
                ss_map[res_id] = r['code']
            except ValueError:
                pass
        return ss_map
