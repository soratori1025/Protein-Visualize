from typing import Optional
from fastapi import Query
from app.schemas.topology import TMParams
from app.services.topology.base import BaseTopologyPredictor

def get_tm_params(
    thickness: Optional[float] = Query(
        None, ge=20.0, le=40.0,
        description="Hydrophobic core thickness in Angstrom. "
                    "Bacterial IM ~27, eukaryotic PM ~30, ER ~25. (Mitra 2004; OPM)"),
    min_tm_element: Optional[int] = Query(
        None, ge=2, le=15,
        description="Min residues of an SS element that must fall inside the slab "
                    "to be considered transmembrane. Lower = more sensitive."),
    min_cross_span: Optional[float] = Query(
        None, ge=0.2, le=0.9,
        description="Min fraction of membrane thickness a crossing must span. "
                    "Lower admits shallower crossings (e.g. TM12 of hSERT)."),
    full_cross_frac: Optional[float] = Query(
        None, ge=0.3, le=1.0,
        description="Fraction of thickness above which a single element counts "
                    "as a full crossing (no fusion needed)."),
    broken_gap_max: Optional[int] = Query(
        None, ge=3, le=20,
        description="Max gap (residues) between two partial helices that can be "
                    "fused into one crossing (broken/discontinuous helix)."),
    min_membrane_score: Optional[float] = Query(
        None, ge=0.0, le=3.0,
        description="Min mean hydrophobicity inside the slab to call membrane. "
                    "Below this the structure is treated as soluble."),
) -> TMParams:
    """Dependency to extract TMParams from query string."""
    overrides = {}
    if thickness is not None: overrides["membrane_thickness"] = thickness
    if min_tm_element is not None: overrides["min_tm_element_in_slab"] = min_tm_element
    if min_cross_span is not None: overrides["min_cross_span_frac"] = min_cross_span
    if full_cross_frac is not None: overrides["full_cross_frac"] = full_cross_frac
    if broken_gap_max is not None: overrides["broken_gap_max"] = broken_gap_max
    if min_membrane_score is not None: overrides["min_membrane_score"] = min_membrane_score
    return TMParams(**overrides)

def get_topology_predictor(algorithm: str = "dssp_slab") -> BaseTopologyPredictor:
    """Dependency to inject the correct topology predictor based on algorithm."""
    if algorithm == "tmhmm_seq":
        from app.services.topology.sequence_predictor import SequenceTopologyPredictor
        return SequenceTopologyPredictor()
    else:
        from app.services.topology.structure_predictor import StructureTopologyPredictor
        return StructureTopologyPredictor()
