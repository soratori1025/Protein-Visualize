import subprocess
from pathlib import Path
from shutil import which
from typing import Optional

import numpy as np
from fastapi import APIRouter, Depends, HTTPException, Query

from protein_engine.secondary_structure.dssp import DSSPMethod
from protein_engine.secondary_structure.stride import STRIDEMethod
from app.schemas.topology import TMParams
from app.services.topology.orchestrator import TopologyOrchestrator
from app.services.topology.providers.tm.uniprot import ACCESSION_RE, fetch_uniprot_entry
from app.services.topology.residues import load_residue_frame, map_to_reference
from app.api.dependencies import get_tm_params, get_topology_orchestrator

router = APIRouter(prefix="/api/secondary-structure", tags=["secondary-structure"])
ROOT = Path(__file__).resolve().parents[3]
UPLOAD_DIR = ROOT / "backend" / "data" / "uploads"

TOPOLOGY_FEATURES = ("Transmembrane", "Topological domain", "Intramembrane")


def _uploaded(filename: str) -> Path:
    path = UPLOAD_DIR / Path(filename).name
    if not path.is_file():
        raise HTTPException(status_code=404, detail="Uploaded structure not found")
    return path


def _as_dict(model) -> dict:
    # pydantic v2 (model_dump) and v1 (dict)
    return model.model_dump() if hasattr(model, "model_dump") else model.dict()


@router.get("/capabilities")
def secondary_structure_capabilities() -> dict:
    dssp_exec = DSSPMethod().executable
    stride_exec = STRIDEMethod().executable
    return {
        "DSSP": {"available": which(dssp_exec) is not None or Path(dssp_exec).is_file(), "executable": dssp_exec},
        "STRIDE": {"available": which(stride_exec) is not None or Path(stride_exec).is_file(), "executable": stride_exec},
    }


@router.post("/{filename}")
def assign_secondary_structure(filename: str, method: str = "DSSP") -> dict:
    path = _uploaded(filename)

    normalized_method = method.upper()
    if normalized_method == "DSSP":
        adapter = DSSPMethod()
    elif normalized_method == "STRIDE":
        adapter = STRIDEMethod()
    else:
        raise HTTPException(status_code=400, detail="method must be DSSP or STRIDE")

    try:
        return adapter.assign(path).to_dict()
    except FileNotFoundError as error:
        raise HTTPException(status_code=503, detail=f"{normalized_method} executable is not installed") from error
    except OSError as error:
        raise HTTPException(status_code=503, detail=f"Unable to execute {normalized_method}: {error}") from error
    except subprocess.CalledProcessError as error:
        detail = error.stderr.strip() if error.stderr else f"{normalized_method} failed with exit code {error.returncode}"
        raise HTTPException(status_code=422, detail=detail) from error
    except RuntimeError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


def _map_regions_to_structure(path: Path, chain_id: Optional[str], sequence: str,
                              regions: list[dict]) -> dict:
    """Add the structure's own residue numbers to every UniProt region.

    UniProt positions are NOT author residue numbers: deposited files skip unresolved
    residues and sometimes renumber stretches (8ZU3 / 9VMX: author 767-857 = UniProt
    789-879). Painting UniProt start/end directly onto author numbers puts TM17-TM19 of
    PIEZO1 on the wrong residues and makes TM19 look unresolved. Positions are matched by
    sequence alignment (the same mapping the UniProt TM provider uses)."""
    frame = load_residue_frame(path, chain_id)
    mapping = map_to_reference(frame, [(u, "") for u in range(1, len(sequence) + 1)],
                               list(sequence), offsets=(0,), source="UniProt",
                               allow_order_fallback=False, note_remap=False)
    if mapping.matched == 0:
        raise HTTPException(
            status_code=422,
            detail=f"UniProt sequence does not match chain {frame.chain_id} of {path.name} "
                   f"(identity {mapping.identity or 0:.0%}). Wrong UniProt ID or chain?")
    pos_of_u = np.full(len(sequence) + 2, -1, dtype=int)
    for p, j in enumerate(mapping.ref_index):
        if j is not None:
            pos_of_u[j + 1] = p

    shifted = 0
    for region in regions:
        u0, u1 = max(1, int(region["start"])), min(len(sequence), int(region["end"]))
        hits = pos_of_u[u0:u1 + 1] if u0 <= u1 else np.array([], dtype=int)
        hits = hits[hits >= 0]
        region["resolved_residues"] = int(hits.size)
        region["total_residues"] = max(0, u1 - u0 + 1)
        region["resolved"] = bool(hits.size)
        if hits.size:
            first, last = frame.residues[int(hits.min())], frame.residues[int(hits.max())]
            region["structure_start"] = first.resseq
            region["structure_end"] = last.resseq
            region["structure_start_icode"] = first.icode or None
            region["structure_end_icode"] = last.icode or None
            shifted += (first.resseq, last.resseq) != (u0, u1) and hits.size == u1 - u0 + 1
        else:
            region["structure_start"] = region["structure_end"] = None
            region["structure_start_icode"] = region["structure_end_icode"] = None

    tms = [r for r in regions if r["type"] == "Transmembrane"]
    warnings = list(mapping.warnings)
    if shifted:
        warnings.append(f"{shifted} fully resolved region(s) have author numbers different from "
                        "UniProt positions: draw them with structure_start/structure_end")
    return {
        "chain_id": frame.chain_id,
        "mapping_method": mapping.method,
        "mapping_identity": mapping.identity,
        "resolved_tm_count": sum(r["resolved"] for r in tms),
        "tm_count": len(tms),
        "mapping_warnings": warnings,
    }


@router.get("/uniprot/{uniprot_id}")
def get_uniprot_topology(
    uniprot_id: str,
    filename: Optional[str] = Query(None, description="Uploaded structure: map the regions onto its residue numbers"),
    chain_id: Optional[str] = Query(None, description="Chain of that structure (default: first protein chain)"),
) -> dict:
    """UniProt topology in UniProt numbering ("start"/"end"). With ``filename`` every region
    also gets "structure_start"/"structure_end" (author numbers of the uploaded file, null
    when not resolved) - the frontend must use those to paint the 3D structure."""
    clean_id = uniprot_id.strip().upper()
    if not ACCESSION_RE.match(clean_id):
        raise HTTPException(status_code=400, detail=f"'{clean_id}' is not a valid UniProt accession")
    data = fetch_uniprot_entry(clean_id)        # cached; HTTP errors already mapped to 404/502/504

    protein_name = (
        data.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value", clean_id)
    )
    genes = data.get("genes", [])
    gene_name = genes[0].get("geneName", {}).get("value", "") if genes else ""
    organism = data.get("organism", {}).get("scientificName", "")

    topology_regions = []
    for f in data.get("features", []):
        ftype = f.get("type")
        if ftype not in TOPOLOGY_FEATURES:
            continue
        loc = f.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        desc = f.get("description", "")

        helix_name = ""
        if "Name=" in desc:
            helix_name = desc.split("Name=", 1)[1].split(";")[0].split(".")[0].strip()

        if start is not None and end is not None:
            topology_regions.append({
                "type": ftype,
                "start": start,
                "end": end,
                "description": desc,
                "name": helix_name,
            })

    result = {
        "uniprot_id": clean_id,
        "protein_name": protein_name,
        "gene_name": gene_name,
        "organism": organism,
        "regions": topology_regions,
    }
    if filename:
        sequence = (data.get("sequence") or {}).get("value", "")
        if not sequence:
            raise HTTPException(status_code=502, detail=f"UniProt entry {clean_id} has no sequence")
        try:
            result["structure"] = _map_regions_to_structure(_uploaded(filename), chain_id,
                                                            sequence, topology_regions)
        except ValueError as error:                       # e.g. unknown chain
            raise HTTPException(status_code=400, detail=str(error)) from error
    return result


@router.get("/predict-topology/{filename}")
def predict_topology_from_structure(
    filename: str,
    uniprot_id: Optional[str] = Query(None, description="Required only if tm_algo=uniprot_api"),
    chain_id: Optional[str] = Query(None, description="Chain to analyse (default: first protein chain)"),
    params: TMParams = Depends(get_tm_params),
    orchestrator: TopologyOrchestrator = Depends(get_topology_orchestrator)
) -> dict:
    path = _uploaded(filename)

    try:
        response = orchestrator.execute(path, params=params, uniprot_id=uniprot_id,
                                        chain_id=(chain_id or None))
        return _as_dict(response)
    except HTTPException:
        raise                                   # keep the provider's 400/404/422/502/504
    except ValueError as error:                 # bad chain / parameters
        raise HTTPException(status_code=400, detail=str(error)) from error
    except Exception as error:
        raise HTTPException(status_code=500, detail=str(error)) from error


@router.get("/predict-topology/{filename}/membrane-normal")
def get_membrane_normal(
    filename: str,
    chain_id: Optional[str] = Query(None, description="Chain to analyse (default: first protein chain)"),
) -> dict:
    path = _uploaded(filename)

    try:
        # We always use the geometry provider to find the membrane normal, regardless of the chosen algorithm
        from app.services.topology.providers.tm.geometry import GeometryTMProvider
        geom_predictor = GeometryTMProvider()

        prediction = geom_predictor.predict_tm(path, chain_id=chain_id)
        if prediction.membrane_normal is None:
            raise HTTPException(status_code=400, detail="No membrane normal detected (likely soluble)")

        return {
            "membrane_normal": prediction.membrane_normal,
            "membrane_score": prediction.membrane_score
        }
    except HTTPException:
        raise
    except ValueError as error:
        raise HTTPException(status_code=400, detail=str(error)) from error
    except Exception as error:
        raise HTTPException(status_code=500, detail=str(error)) from error