from dataclasses import dataclass
from typing import Optional
from pydantic import BaseModel, Field
from app.core.constants import (
    MEMBRANE_THICKNESS,
    MIN_MEMBRANE_SCORE,
    MIN_TM_ELEMENT_IN_SLAB,
    MIN_CROSS_SPAN_FRAC,
    FULL_CROSS_FRAC,
    BROKEN_GAP_MAX,
)

@dataclass
class TMParams:
    """All biology-dependent thresholds that affect TM detection."""
    membrane_thickness: float = MEMBRANE_THICKNESS
    min_membrane_score: float = MIN_MEMBRANE_SCORE
    min_tm_element_in_slab: int = MIN_TM_ELEMENT_IN_SLAB
    min_cross_span_frac: float = MIN_CROSS_SPAN_FRAC
    full_cross_frac: float = FULL_CROSS_FRAC
    broken_gap_max: int = BROKEN_GAP_MAX
    def to_response_dict(self) -> dict:
        return {
            "membrane_thickness": self.membrane_thickness,
            "min_tm_element_in_slab": self.min_tm_element_in_slab,
            "min_cross_span_frac": self.min_cross_span_frac,
            "full_cross_frac": self.full_cross_frac,
            "broken_gap_max": self.broken_gap_max,
            "min_membrane_score": self.min_membrane_score,
        }

class TopologyRegion(BaseModel):
    type: str
    start: int                       # author residue number (first observed residue)
    end: int                         # author residue number (last observed residue)
    description: str
    side: Optional[str] = None
    ss: Optional[str] = None
    # Insertion codes of the boundary residues ("100A"); None when blank. Additive,
    # so existing clients that only read start/end keep working.
    start_icode: Optional[str] = None
    end_icode: Optional[str] = None

class TopologyResponse(BaseModel):
    uniprot_id: str
    protein_name: str
    gene_name: str
    organism: str
    membrane_score: float = 0.0
    labeler: str
    membrane_normal: Optional[list[float]] = None
    parameters_used: Optional[dict] = None
    regions: list[TopologyRegion] = Field(default_factory=list)
    # Additive fields: which chain was analysed, and anything the caller should know
    # (SS tool missing, residue numbering re-mapped, features outside the model ...).
    chain_id: Optional[str] = None
    warnings: list[str] = Field(default_factory=list)
