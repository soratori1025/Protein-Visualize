"""
Request-level entry point for "Calculated topology": turns the query parameters the
frontend sends (tm_algo, ss_algo, uniprot_id, chain_id and the advanced TM parameters)
into providers + TMParams, and runs the orchestrator. There is ONE flow (consensus);
``flow_type`` is still accepted so old clients keep working.

Keeping this mapping here (not in the FastAPI route) means the frontend/backend
contract is in one testable place. The route only resolves the uploaded file and
turns ValueError into HTTP 400:

    @router.get("/predict-topology/{filename}", response_model=TopologyResponse)
    def predict_topology_endpoint(filename: str, tm_algo: str = "3d_slab_geom",
                                  ss_algo: str = "dssp", flow_type: str | None = None,
                                  uniprot_id: str | None = None, chain_id: str | None = None,
                                  thickness: float | None = None,
                                  min_membrane_score: float | None = None,
                                  min_cross_span: float | None = None,
                                  full_cross_frac: float | None = None,
                                  treat_turn_as_helix: bool | None = None,
                                  broken_gap_max: int | None = None,      # ignored (old flows)
                                  min_tm_element: int | None = None,      # ignored (old flows)
                                  include_residues: bool = True):
        path = <resolve uploaded filename as the other endpoints do>
        try:
            return predict_topology(path, tm_algo=tm_algo, ss_algo=ss_algo, flow_type=flow_type,
                                    uniprot_id=uniprot_id, chain_id=chain_id, thickness=thickness,
                                    min_membrane_score=min_membrane_score,
                                    min_cross_span=min_cross_span, full_cross_frac=full_cross_frac,
                                    treat_turn_as_helix=treat_turn_as_helix,
                                    broken_gap_max=broken_gap_max, min_tm_element=min_tm_element,
                                    include_residues=include_residues)
        except ValueError as error:
            raise HTTPException(status_code=400, detail=str(error))
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Callable, Optional

from app.schemas.topology import TMParams, TopologyResponse
from app.services.topology.orchestrator import (
    CONSENSUS_NAMES, FLOW, REMOVED_FLOWS, TopologyOrchestrator,
)
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
# Read only by the removed tm_then_ss / ss_then_tm flows: accepted, reported as ignored.
IGNORED_PARAMS = {"broken_gap_max", "min_tm_element_in_slab"}


def build_params(thickness: Optional[float] = None, min_membrane_score: Optional[float] = None,
                 broken_gap_max: Optional[int] = None, min_tm_element: Optional[int] = None,
                 min_cross_span: Optional[float] = None,
                 full_cross_frac: Optional[float] = None,
                 treat_turn_as_helix: Optional[bool] = None) -> TMParams:
    """TMParams from optional overrides; None keeps the default from app.core.constants."""
    overrides = {
        "membrane_thickness": thickness,
        "min_membrane_score": min_membrane_score,
        "broken_gap_max": broken_gap_max,
        "min_tm_element_in_slab": min_tm_element,
        "min_cross_span_frac": min_cross_span,
        "full_cross_frac": full_cross_frac,
    }
    clean: dict = {}
    for field, value in overrides.items():
        if value is None:
            continue
        lo, hi = PARAM_BOUNDS[field]
        if not (lo <= value <= hi):
            raise ValueError(f"{field}={value} is outside the accepted range [{lo}, {hi}]")
        clean[field] = value
    if treat_turn_as_helix is not None:
        clean["treat_turn_as_helix"] = bool(treat_turn_as_helix)
    return replace(TMParams(), **clean)


def predict_topology(file_path: Path, tm_algo: str = "3d_slab_geom", ss_algo: str = "dssp",
                     flow_type: Optional[str] = None, uniprot_id: Optional[str] = None,
                     chain_id: Optional[str] = None, include_residues: bool = True,
                     **param_overrides) -> TopologyResponse:
    tm_key = (tm_algo or "").strip().lower()
    ss_key = (ss_algo or "none").strip().lower()
    if tm_key not in TM_ALGORITHMS:
        raise ValueError(f"Unknown tm_algo '{tm_algo}'; expected one of {sorted(TM_ALGORITHMS)}")
    if ss_key not in SS_ALGORITHMS:
        raise ValueError(f"Unknown ss_algo '{ss_algo}'; expected one of {sorted(SS_ALGORITHMS)}")
    flow = (flow_type or FLOW).strip()
    if flow not in CONSENSUS_NAMES | REMOVED_FLOWS:
        raise ValueError(f"Unknown flow_type '{flow_type}'; the only flow is '{FLOW}'")

    params = build_params(**param_overrides)
    ss_factory = SS_ALGORITHMS[ss_key]
    orchestrator = TopologyOrchestrator(TM_ALGORITHMS[tm_key](),
                                        ss_factory() if ss_factory else None, flow,
                                        include_residues=include_residues)
    kwargs = {"uniprot_id": uniprot_id.strip()} if uniprot_id and uniprot_id.strip() else {}
    response = orchestrator.execute(Path(file_path), params=params,
                                    chain_id=(chain_id or None), **kwargs)
    ignored = [name for name, value in param_overrides.items()
               if value is not None and _field_for(name) in IGNORED_PARAMS]
    if ignored:
        response.warnings.append(f"parameters {ignored} belonged to the removed flows; "
                                 "ignored by the consensus flow")
    return response


def _field_for(query_name: str) -> str:
    return {"thickness": "membrane_thickness", "min_tm_element": "min_tm_element_in_slab",
            "min_cross_span": "min_cross_span_frac"}.get(query_name, query_name)
