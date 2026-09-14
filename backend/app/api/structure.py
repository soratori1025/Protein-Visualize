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

import re
import urllib.request
import json

UNIPROT_REGEX = re.compile(r'\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])\b', re.IGNORECASE)

def extract_uniprot_id(path: Path) -> str | None:
    # 1. Match from filename (e.g. P31645.pdb, AF-P23975-F1.pdb, P00533_model.cif)
    stem = path.stem.upper()
    match = UNIPROT_REGEX.search(stem)
    if match:
        return match.group(1).upper()

    # 2. Match from PDB DBREF / mmCIF headers inside file content
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("DBREF"):
                    parts = line.split()
                    for idx, part in enumerate(parts):
                        if part in ("UNP", "SWS") and idx + 1 < len(parts):
                            return parts[idx + 1].upper()
                if "_struct_ref.pdbx_db_accession" in line or "_struct_ref.db_code" in line:
                    parts = line.split()
                    if len(parts) >= 2:
                        val = parts[-1].strip("'\" ")
                        if len(val) >= 6:
                            return val.upper()
                # If DBREF is missing, check if line contains UNP accession code
                if "UNP" in line or "UNIPROT" in line.upper():
                    unp_match = UNIPROT_REGEX.search(line)
                    if unp_match:
                        return unp_match.group(1).upper()
    except Exception:
        pass

    # 3. Fallback: check 4-letter PDB ID in filename (e.g. 7r2v.pdb, 6dzz.pdb)
    if len(path.stem) == 4 and path.stem.isalnum():
        pdb_id = path.stem.lower()
        try:
            url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
            req = urllib.request.Request(url, headers={"User-Agent": "ProteinVisualizeApp/1.0"})
            with urllib.request.urlopen(req, timeout=3) as resp:
                if resp.status == 200:
                    data = json.loads(resp.read().decode("utf-8"))
                    # Search structure_keywords or primary citation
                    text_str = json.dumps(data)
                    m = UNIPROT_REGEX.search(text_str)
                    if m:
                        return m.group(1).upper()
        except Exception:
            pass

    return None


@router.post("/upload")
async def upload_structure(file: UploadFile = File(...)):
    filename = Path(file.filename or "structure.pdb").name
    if Path(filename).suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
        raise HTTPException(status_code=400, detail="Only PDB and mmCIF files are supported")
    saved_path = ProteinService.save_uploaded_file(await file.read(), filename)
    uniprot_id = extract_uniprot_id(saved_path)
    return {
        "filename": filename,
        "saved_path": str(saved_path),
        "uniprot_id": uniprot_id,
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
