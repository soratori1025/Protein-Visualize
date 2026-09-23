from pathlib import Path
from typing import List, Mapping

from app.services.topology.providers.ss.base import SSProvider


class STRIDEProvider(SSProvider):
    labeler_name = "STRIDE"

    def raw_residues(self, file_path: Path) -> List[Mapping]:
        from protein_engine.secondary_structure.stride import STRIDEMethod   # lazy
        return STRIDEMethod().assign(file_path).to_dict().get("residues", []) or []
