from pathlib import Path
from shutil import which
import subprocess
from fastapi import APIRouter, HTTPException
from protein_engine.secondary_structure.dssp import DSSPMethod
from protein_engine.secondary_structure.stride import STRIDEMethod

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