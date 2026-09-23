from typing import Optional
from fastapi import Query
from app.schemas.topology import TMParams

def get_tm_params(
    thickness: Optional[float] = Query(None, ge=20.0, le=40.0),
    min_tm_element: Optional[int] = Query(None, ge=2, le=15),
    min_cross_span: Optional[float] = Query(None, ge=0.2, le=0.9),
    full_cross_frac: Optional[float] = Query(None, ge=0.3, le=1.0),
    broken_gap_max: Optional[int] = Query(None, ge=3, le=20),
    min_membrane_score: Optional[float] = Query(None, ge=0.0, le=3.0)
) -> TMParams:
    overrides = {}
    if thickness is not None: overrides["membrane_thickness"] = thickness
    if min_tm_element is not None: overrides["min_tm_element_in_slab"] = min_tm_element
    if min_cross_span is not None: overrides["min_cross_span_frac"] = min_cross_span
    if full_cross_frac is not None: overrides["full_cross_frac"] = full_cross_frac
    if broken_gap_max is not None: overrides["broken_gap_max"] = broken_gap_max
    if min_membrane_score is not None: overrides["min_membrane_score"] = min_membrane_score
    return TMParams(**overrides)

def get_topology_orchestrator(
    tm_algo: str = Query("3d_slab_geom", description="TM Algorithm: kyte_doolittle_seq, 3d_slab_geom, uniprot_api"),
    ss_algo: str = Query("dssp", description="SS Algorithm: dssp, stride, none"),
    flow_type: str = Query("ss_then_tm", description="Flow Type: ss_then_tm, tm_then_ss, parallel_merge")
):
    from app.services.topology.orchestrator import TopologyOrchestrator
    
    # Select TM Provider
    if tm_algo == "kyte_doolittle_seq":
        from app.services.topology.providers.tm.sequence import SequenceTMProvider
        tm_provider = SequenceTMProvider()
    elif tm_algo == "uniprot_api":
        from app.services.topology.providers.tm.uniprot import UniprotTMProvider
        tm_provider = UniprotTMProvider()
    else:
        from app.services.topology.providers.tm.geometry import GeometryTMProvider
        tm_provider = GeometryTMProvider()
        
    # Select SS Provider
    ss_provider = None
    if ss_algo == "dssp":
        from app.services.topology.providers.ss.dssp import DSSPProvider
        ss_provider = DSSPProvider()
    elif ss_algo == "stride":
        from app.services.topology.providers.ss.stride import STRIDEProvider
        ss_provider = STRIDEProvider()

    return TopologyOrchestrator(tm_provider, ss_provider, flow_type)
