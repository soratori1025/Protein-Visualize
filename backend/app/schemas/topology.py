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
    """All biology-dependent thresholds that affect TM detection (consensus flow)."""
    membrane_thickness: float = MEMBRANE_THICKNESS      # 3D slab fit + estimated bilayer
    min_membrane_score: float = MIN_MEMBRANE_SCORE      # 3D slab: below -> soluble protein
    min_cross_span_frac: float = MIN_CROSS_SPAN_FRAC    # hairpin guard / bilayer centring
    full_cross_frac: float = FULL_CROSS_FRAC            # hairpin guard: a run crosses alone
    treat_turn_as_helix: bool = False                   # T/S next to helices (orchestrator)
    # Read by the removed tm_then_ss / ss_then_tm flows only. Kept so old requests still
    # validate; the service reports them as ignored.
    min_tm_element_in_slab: int = MIN_TM_ELEMENT_IN_SLAB
    broken_gap_max: int = BROKEN_GAP_MAX

    def to_response_dict(self) -> dict:
        return {
            "membrane_thickness": self.membrane_thickness,
            "min_tm_element_in_slab": self.min_tm_element_in_slab,
            "min_cross_span_frac": self.min_cross_span_frac,
            "full_cross_frac": self.full_cross_frac,
            "broken_gap_max": self.broken_gap_max,
            "min_membrane_score": self.min_membrane_score,
            "treat_turn_as_helix": self.treat_turn_as_helix,
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
    # 'high' | 'medium' | 'low': agreement of the TM block, the membrane geometry and
    # DSSP/STRIDE over the region's residues (heuristic, not a calibrated probability).
    confidence: Optional[str] = None
    # Transmembrane regions only: 1-based crossing number. A discontinuous crossing
    # (TM1a/1b ...) is reported as several regions with the SAME crossing number:
    # part "a", the unwound stretch (part None, description "Transmembrane Unwound"),
    # part "b".
    crossing: Optional[int] = None
    part: Optional[str] = None
    # Membrane role (TM_CROSSING, BROKEN_TM, UNWOUND, REENTRANT, INTERFACIAL,
    # EXTRAMEMBRANE, SIGNAL) is separate from secondary structure (`ss`).
    membrane_role: Optional[str] = None
    # TM1, TM1a, TM1 unwound, EL2, EL3a, IL1, N-term, C-term, RE1 ...
    topology_label: Optional[str] = None
    parent_tm: Optional[str] = None
    # Evidence for every fragment junction inside this crossing (why 1a/1b or not).
    transitions: Optional[list[dict]] = None


class ResidueAnnotation(BaseModel):
    """One row of the residue-level evidence matrix; regions are derived from these."""
    index: int                          # position in the analysed chain (0-based)
    residue_number: int
    insertion_code: Optional[str] = None
    aa: str
    ss_raw: Optional[str] = None        # DSSP/STRIDE 8-state code (H G I E B T S C)
    ss: Optional[str] = None            # coarse: H / E / C
    tm_evidence: str                    # what the TM block said (Transmembrane, Cytoplasmic, ...)
    depth: Optional[float] = None       # signed distance from the bilayer mid-plane (A)
    zone: Optional[str] = None          # CORE / EDGE / OUT
    plddt: Optional[float] = None       # predicted models only (B-factor column)
    label: str                          # final label (same vocabulary as region descriptions)
    confidence: Optional[str] = None


class ConsensusResidue(BaseModel):
    """One entry of the consensus residue map (the only flow).
    Only residues whose SS flag is on (helix; strand in a beta TM segment; T/S next to
    a helix with treat_turn_as_helix) appear."""
    index: int                          # position in the analysed chain (0-based)
    residue_number: int
    insertion_code: Optional[str] = None
    aa: str
    label: str                          # TM_in | TM_C | TM_E | Turn_in | Turn_C | Turn_E
    ss_raw: str                         # raw DSSP/STRIDE code (H G I E; T/S for promoted turns)
    tm_segment: Optional[int] = None    # TM_in: 1-based TM segment of the TM block (UniProt TM feature)
    crossing: Optional[int] = None      # TM_in: 1-based crossing it is drawn in (None = not drawn)
    part: Optional[str] = None          # "a" / "b" when that crossing is drawn broken (TM1a / TM1b)


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
    # 'alpha_helical' | 'beta_barrel' | 'beta' | 'mixed' | 'irregular' | None
    domain_type: Optional[str] = None
    # {"source", "normal", "half_thickness"} of the membrane placement used
    membrane: Optional[dict] = None
    residues: Optional[list[ResidueAnnotation]] = None
    # residue map TM_in / TM_C / TM_E (see ConsensusResidue); None without an SS result
    consensus_map: Optional[list[ConsensusResidue]] = None
