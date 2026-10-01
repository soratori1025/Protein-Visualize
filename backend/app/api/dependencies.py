from typing import Optional
from fastapi import Query
from app.schemas.topology import TMParams

def get_tm_params(
    thickness: Optional[float] = Query(None, ge=20.0, le=40.0),
    min_tm_element: Optional[int] = Query(None, ge=2, le=15),
    min_cross_span: Optional[float] = Query(None, ge=0.2, le=0.9),
    full_cross_frac: Optional[float] = Query(None, ge=0.3, le=1.0),
    broken_gap_max: Optional[int] = Query(None, ge=3, le=20),
    min_membrane_score: Optional[float] = Query(None, ge=0.0, le=3.0),
    treat_turn_as_helix: Optional[bool] = Query(False)
) -> TMParams:
    overrides = {}
    if thickness is not None: overrides["membrane_thickness"] = thickness
    if min_tm_element is not None: overrides["min_tm_element_in_slab"] = min_tm_element
    if min_cross_span is not None: overrides["min_cross_span_frac"] = min_cross_span
    if full_cross_frac is not None: overrides["full_cross_frac"] = full_cross_frac
    if broken_gap_max is not None: overrides["broken_gap_max"] = broken_gap_max
    if min_membrane_score is not None: overrides["min_membrane_score"] = min_membrane_score
    if treat_turn_as_helix is not None: overrides["treat_turn_as_helix"] = treat_turn_as_helix
    return TMParams(**overrides)

def get_topology_orchestrator(
    tm_algo: str = Query("3d_slab_geom", description="TM Algorithm: kyte_doolittle_seq, 3d_slab_geom, uniprot_api, 3d_energy, tmhmm"),
    ss_algo: str = Query("dssp", description="SS Algorithm: dssp, stride, none"),
    flow_type: str = Query("ss_then_tm", description="Flow Type: ss_then_tm, tm_then_ss, parallel_merge")
):
    from app.services.topology.orchestrator import TopologyOrchestrator
    from app.services.topology.service import TM_ALGORITHMS, SS_ALGORITHMS
    from fastapi import HTTPException
    
    tm_key = tm_algo.strip().lower()
    ss_key = ss_algo.strip().lower()
    
    if tm_key not in TM_ALGORITHMS:
        raise HTTPException(status_code=400, detail=f"Unknown tm_algo '{tm_algo}'")
    if ss_key not in SS_ALGORITHMS:
        raise HTTPException(status_code=400, detail=f"Unknown ss_algo '{ss_algo}'")
        
    # Select TM Provider
    tm_provider = TM_ALGORITHMS[tm_key]()
        
    # Select SS Provider
    ss_factory = SS_ALGORITHMS[ss_key]
    ss_provider = ss_factory() if ss_factory else None

    return TopologyOrchestrator(tm_provider, ss_provider, flow_type)
