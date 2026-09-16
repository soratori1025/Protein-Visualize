import json
from pathlib import Path
from shutil import which
import subprocess
import urllib.request
from urllib.error import HTTPError, URLError
from fastapi import APIRouter, HTTPException
from protein_engine.secondary_structure.dssp import DSSPMethod
from protein_engine.secondary_structure.stride import STRIDEMethod
from .topology_predictor import predict_topology

router = APIRouter(prefix="/api/secondary-structure", tags=["secondary-structure"])
ROOT = Path(__file__).resolve().parents[3]
UPLOAD_DIR = ROOT / "backend" / "data" / "uploads"


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
    path = UPLOAD_DIR / Path(filename).name
    if not path.is_file():
        raise HTTPException(status_code=404, detail="Uploaded structure not found")

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


@router.get("/uniprot/{uniprot_id}")
def get_uniprot_topology(uniprot_id: str) -> dict:
    clean_id = uniprot_id.strip().upper()
    url = f"https://rest.uniprot.org/uniprotkb/{clean_id}.json"
    req = urllib.request.Request(url, headers={"User-Agent": "ProteinVisualizeApp/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=10) as response:
            if response.status != 200:
                raise HTTPException(status_code=response.status, detail="Failed to fetch data from UniProt API")
            data = json.loads(response.read().decode("utf-8"))
    except HTTPError as err:
        if err.code == 404:
            raise HTTPException(status_code=404, detail=f"UniProt ID {clean_id} not found") from err
        raise HTTPException(status_code=502, detail=f"UniProt API error: {err.reason}") from err
    except URLError as err:
        raise HTTPException(status_code=504, detail=f"Network error connecting to UniProt API: {err.reason}") from err

    protein_name = (
        data.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value", clean_id)
    )
    genes = data.get("genes", [])
    gene_name = genes[0].get("geneName", {}).get("value", "") if genes else ""
    organism = data.get("organism", {}).get("scientificName", "")

    features = data.get("features", [])
    topology_regions = []

    for f in features:
        ftype = f.get("type")
        if ftype in ("Transmembrane", "Topological domain", "Intramembrane"):
            loc = f.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            desc = f.get("description", "")

            helix_name = ""
            if "Name=" in desc:
                parts = desc.split("Name=")
                if len(parts) > 1:
                    helix_name = parts[1].split(";")[0].split(".")[0].strip()

            if start is not None and end is not None:
                topology_regions.append({
                    "type": ftype,
                    "start": start,
                    "end": end,
                    "description": desc,
                    "name": helix_name,
                })

    return {
        "uniprot_id": clean_id,
        "protein_name": protein_name,
        "gene_name": gene_name,
        "organism": organism,
        "regions": topology_regions,
    }

@router.get("/predict-topology/{filename}")
def predict_topology_from_structure(filename: str) -> dict:
    path = UPLOAD_DIR / Path(filename).name
    if not path.is_file():
        raise HTTPException(status_code=404, detail="Uploaded structure not found")
    
    try:
        return predict_topology(path)
    except Exception as error:
        raise HTTPException(status_code=500, detail=str(error)) from error
