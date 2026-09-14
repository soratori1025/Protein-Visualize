from abc import ABC, abstractmethod
from pathlib import Path
from protein_engine.secondary_structure.result import SecondaryStructureResult

class SecondaryStructureMethod(ABC):
    name: str
    @abstractmethod
    def assign(self, structure_path: Path, model_id: int = 0) -> SecondaryStructureResult:
        pass
