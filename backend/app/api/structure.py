import sys
from pathlib import Path

from fastapi import APIRouter, File, UploadFile

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from app.services.protein_service import ProteinService
from protein_engine.parser.pdb import parse_pdb

router = APIRouter(prefix="/api/structure", tags=["structure"])


@router.post("/upload")
async def upload_structure(file: UploadFile = File(...)):
    saved_path = ProteinService.save_uploaded_file(await file.read(), file.filename)
    structure = parse_pdb(str(saved_path))

    summary = []
    for model in structure:
        chains = []
        for chain in model:
            residues = list(chain.get_residues())
            chains.append({
                "id": chain.id,
                "residue_count": len(residues),
            })
        summary.append({
            "id": model.id,
            "chains": chains,
        })

    return {
        "filename": file.filename,
        "saved_path": str(saved_path),
        "models": summary,
    }
