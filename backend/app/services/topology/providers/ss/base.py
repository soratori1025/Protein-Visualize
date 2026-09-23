from abc import ABC, abstractmethod
from typing import Dict
from pathlib import Path

class SSProvider(ABC):
    @abstractmethod
    def get_secondary_structure(self, file_path: Path) -> Dict[int, str]:
        """
        Returns a mapping from residue integer ID to its secondary structure code.
        Typically, codes are 'H' (Helix), 'E' (Strand), 'C' (Coil/Loop).
        """
        pass
