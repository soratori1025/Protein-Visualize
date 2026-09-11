import sys
from pathlib import Path

from fastapi import APIRouter, File, HTTPException, UploadFile
from fastapi.responses import FileResponse
from fastapi import HTTPException

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from app.services.protein_service import ProteinService
from protein_engine.parser.pdb import parse_pdb

router = APIRouter(prefix="/api/structure", tags=["structure"])
UPLOAD_DIR = Path("data/uploads")


@router.get("/file/{filename}")
def get_structure_file(filename: str):
    path = UPLOAD_DIR / Path(filename).name
    if not path.is_file():
        raise HTTPException(status_code=404, detail="Uploaded structure not found")
    return FileResponse(path)


@router.post("/upload")
async def upload_structure(file: UploadFile = File(...)):
    filename = Path(file.filename or "structure.pdb").name
    if Path(filename).suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
        raise HTTPException(status_code=400, detail="Only PDB and mmCIF files are supported")
    saved_path = ProteinService.save_uploaded_file(await file.read(), filename)
    return {
        "filename": filename,
        "saved_path": str(saved_path),
        "models": build_structure_summary(saved_path),
    }


def build_structure_summary(path: Path) -> list[dict]:
    structure = parse_pdb(str(path))
    summary = []
    for model in structure:
        chains = []
        for chain in model:
            residues = []
            for residue in chain.get_residues():
                if residue.id[0].strip():
                    continue
                atoms = [
                    {
                        "name": atom.name,
                        "element": atom.element or atom.name[0],
                        "x": round(float(atom.coord[0]), 3),
                        "y": round(float(atom.coord[1]), 3),
                        "z": round(float(atom.coord[2]), 3),
                    }
                    for atom in residue.get_atoms()
                ]
                residues.append({
                    "id": residue.id[1],
                    "name": residue.resname.strip(),
                    "atoms": atoms,
                })
            chains.append({
                "id": chain.id,
                "residue_count": len(residues),
                "sequence": "".join(THREE_TO_ONE.get(residue["name"], "X") for residue in residues),
                "residues": residues,
            })
        summary.append({
            "id": model.id,
            "chains": chains,
        })

    return summary


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
