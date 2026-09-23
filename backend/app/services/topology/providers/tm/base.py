from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Optional, Tuple

from pydantic import BaseModel, Field

from app.schemas.topology import TMParams, TopologyRegion
from app.services.topology.residues import ResidueFrame


class TMPrediction(BaseModel):
    """Output of a TM provider.

    ``labels`` / ``segments`` are the authoritative output: they are indexed by
    ResidueFrame POSITION (0..n-1) on the frame the orchestrator passed in, so no
    residue-number conversion (and no offset bug) can happen between blocks.
    ``regions`` (author numbering) is kept for display and for older callers; the
    orchestrator only falls back to it for third-party providers without labels."""
    regions: List[TopologyRegion] = Field(default_factory=list)
    labels: Optional[List[str]] = None                    # one per frame position
    segments: Optional[List[Tuple[int, int]]] = None      # TM segments, inclusive positions
    membrane_normal: Optional[List[float]] = None
    membrane_score: float = 0.0
    labeler: str = ""
    warnings: List[str] = Field(default_factory=list)


class TMProvider(ABC):
    @abstractmethod
    def predict_tm(self, file_path: Path, params: Optional[TMParams] = None,
                   frame: Optional[ResidueFrame] = None, **kwargs) -> TMPrediction:
        """Predict TM segments for the chain described by ``frame``.

        ``frame`` is the shared residue index built once by the orchestrator. A
        provider must NOT re-parse the file with its own residue filter when a frame
        is given; if it is None the provider builds one itself (standalone use)."""
