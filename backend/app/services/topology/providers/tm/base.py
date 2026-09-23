from abc import ABC, abstractmethod
from typing import List, Tuple, Dict, Any, Optional
from pathlib import Path
from pydantic import BaseModel

class TMBoundary(BaseModel):
    start: int
    end: int

class TMPrediction(BaseModel):
    boundaries: List[TMBoundary]
    membrane_normal: Optional[List[float]] = None
    membrane_score: float = 0.0
    labeler: str = ""

class TMProvider(ABC):
    @abstractmethod
    def predict_tm(self, file_path: Path, **kwargs) -> TMPrediction:
        """
        Returns TM boundaries (start and end residue IDs) along with any metadata 
        (e.g., membrane_normal for geometry-based methods).
        """
        pass
