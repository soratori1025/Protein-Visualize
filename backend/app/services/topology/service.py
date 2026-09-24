"""
Request-level entry point for "Calculated topology": turns the query parameters the
frontend sends (tm_algo, ss_algo, flow_type, uniprot_id, chain_id and the advanced
TM parameters) into providers + TMParams, and runs the orchestrator.

Keeping this mapping here (not in the FastAPI route) means the frontend/backend
contract is in one testable place. The route only resolves the uploaded file and
turns ValueError into HTTP 400:

    @router.get("/predict-topology/{filename}", response_model=TopologyResponse)
    def predict_topology_endpoint(filename: str, tm_algo: str = "3d_slab_geom",   # flow_type also
                                  # accepts the aliases tm_first / structure_guided / consensus
                                  ss_algo: str = "dssp", flow_type: str = "ss_then_tm",
                                  uniprot_id: str | None = None, chain_id: str | None = None,
                                  thickness: float | None = None, broken_gap_max: int | None = None,
                                  min_membrane_score: float | None = None,
                                  min_tm_element: int | None = None, min_cross_span: float | None = None,
                                  full_cross_frac: float | None = None,
                                  include_residues: bool = True):
        path = <resolve uploaded filename as the other endpoints do>
        try:
            return predict_topology(path, tm_algo=tm_algo, ss_algo=ss_algo, flow_type=flow_type,
                                    uniprot_id=uniprot_id, chain_id=chain_id, thickness=thickness,
                                    broken_gap_max=broken_gap_max, min_membrane_score=min_membrane_score,
                                    min_tm_element=min_tm_element, min_cross_span=min_cross_span,
                                    full_cross_frac=full_cross_frac,
                                    include_residues=include_residues)
        except ValueError as error:
            raise HTTPException(status_code=400, detail=str(error))
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Callable, Optional

from app.schemas.topology import TMParams, TopologyResponse
from app.services.topology.orchestrator import FLOW_ALIASES, FlowType, TopologyOrchestrator
from app.services.topology.providers.ss.base import SSProvider
from app.services.topology.providers.ss.dssp import DSSPProvider
from app.services.topology.providers.ss.stride import STRIDEProvider
from app.services.topology.providers.tm.base import TMProvider
from app.services.topology.providers.tm.sequence import SequenceTMProvider
from app.services.topology.providers.tm.uniprot import UniprotTMProvider


def _geometry_provider() -> TMProvider:
    # imported lazily: the predictor module is large and pulls in the numeric core
    from app.services.topology.topology_predictor import GeometryTMProvider
    return GeometryTMProvider()


TM_ALGORITHMS: dict[str, Callable[[], TMProvider]] = {
    "3d_slab_geom": _geometry_provider,
    "kyte_doolittle_seq": SequenceTMProvider,
    "uniprot_api": UniprotTMProvider,
}
SS_ALGORITHMS: dict[str, Optional[Callable[[], SSProvider]]] = {
    "dssp": DSSPProvider,
    "stride": STRIDEProvider,
    "none": None,
}

# Accepted ranges for user-tunable parameters (same bounds as the frontend inputs).
PARAM_BOUNDS = {
    "membrane_thickness": (10.0, 60.0),
    "min_membrane_score": (-5.0, 5.0),
    "broken_gap_max": (0, 30),
    "min_tm_element_in_slab": (1, 40),
    "min_cross_span_frac": (0.05, 1.0),
    "full_cross_frac": (0.1, 1.5),
}
# Crossing-geometry parameters are read only by the structure-guided flow
# (ss_then_tm) and by the SS-element-first predictor.
STRUCTURE_GUIDED_ONLY = {"min_tm_element_in_slab", "min_cross_span_frac", "full_cross_frac"}


def build_params(thickness: Optional[float] = None, min_membrane_score: Optional[float] = None,
                 broken_gap_max: Optional[int] = None, min_tm_element: Optional[int] = None,
                 min_cross_span: Optional[float] = None,
                 full_cross_frac: Optional[float] = None) -> TMParams:
    """TMParams from optional overrides; None keeps the default from app.core.constants."""
    overrides = {
        "membrane_thickness": thickness,
        "min_membrane_score": min_membrane_score,
        "broken_gap_max": broken_gap_max,
        "min_tm_element_in_slab": min_tm_element,
        "min_cross_span_frac": min_cross_span,
        "full_cross_frac": full_cross_frac,
    }
    clean = {}
    for field, value in overrides.items():
        if value is None:
            continue
        lo, hi = PARAM_BOUNDS[field]
        if not (lo <= value <= hi):
            raise ValueError(f"{field}={value} is outside the accepted range [{lo}, {hi}]")
        clean[field] = value
    return replace(TMParams(), **clean)


def predict_topology(file_path: Path, tm_algo: str = "3d_slab_geom", ss_algo: str = "dssp",
                     flow_type: str = "ss_then_tm", uniprot_id: Optional[str] = None,
                     chain_id: Optional[str] = None, include_residues: bool = True,
                     **param_overrides) -> TopologyResponse:
    tm_key = (tm_algo or "").strip().lower()
    ss_key = (ss_algo or "none").strip().lower()
    if tm_key not in TM_ALGORITHMS:
        raise ValueError(f"Unknown tm_algo '{tm_algo}'; expected one of {sorted(TM_ALGORITHMS)}")
    if ss_key not in SS_ALGORITHMS:
        raise ValueError(f"Unknown ss_algo '{ss_algo}'; expected one of {sorted(SS_ALGORITHMS)}")
    flow_type = FLOW_ALIASES.get(flow_type, flow_type)
    if flow_type not in {f.value for f in FlowType}:
        raise ValueError(f"Unknown flow_type '{flow_type}'; expected one of "
                         f"{[f.value for f in FlowType] + list(FLOW_ALIASES)}")

    params = build_params(**param_overrides)
    ss_factory = SS_ALGORITHMS[ss_key]
    orchestrator = TopologyOrchestrator(TM_ALGORITHMS[tm_key](),
                                        ss_factory() if ss_factory else None, flow_type,
                                        include_residues=include_residues)
    kwargs = {"uniprot_id": uniprot_id.strip()} if uniprot_id and uniprot_id.strip() else {}
    response = orchestrator.execute(Path(file_path), params=params,
                                    chain_id=(chain_id or None), **kwargs)
    ignored = [name for name, value in param_overrides.items()
               if value is not None and _field_for(name) in STRUCTURE_GUIDED_ONLY
               and flow_type != FlowType.SS_THEN_TM.value]
    if ignored:
        response.warnings.append(f"parameters {ignored} only affect the structure-guided flow "
                                 "(ss_then_tm); ignored here")
    return response


def _field_for(query_name: str) -> str:
    return {"thickness": "membrane_thickness", "min_tm_element": "min_tm_element_in_slab",
            "min_cross_span": "min_cross_span_frac"}.get(query_name, query_name)
