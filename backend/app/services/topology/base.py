from abc import ABC, abstractmethod
from typing import Optional, Any
from pathlib import Path
from app.schemas.topology import TMParams, TopologyResponse

class BaseTopologyPredictor(ABC):
    """
    Abstract Base Class for Transmembrane Topology Predictors.
    Follows Open-Closed Principle allowing new predictors (e.g. sequence-based, AlphaFold-based)
    to be added without modifying existing code.
    """
    
    @abstractmethod
    def predict(self, file_path: Path, labeler: str = "DSSP", params: Optional[TMParams] = None) -> TopologyResponse:
        """
        Predict the transmembrane topology of a protein.
        
        Args:
            file_path: Path to the protein structure or sequence file.
            labeler: The secondary structure labeler to use (DSSP or STRIDE).
            params: User-configurable biological thresholds.
            
        Returns:
            TopologyResponse containing the predicted regions.
        """
        pass
