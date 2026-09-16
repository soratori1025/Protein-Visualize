import re
import sys
import json
import urllib.request
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


# A UniProt accession is exactly 6 characters: [OPQ]+digit+3 alnum+digit, or
# a letter from A-N/R-Z + digit + 3 alnum + digit. Anchored so it must match a
# WHOLE token, not just appear somewhere inside a longer alnum run.
UNIPROT_REGEX = re.compile(r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])$', re.IGNORECASE)

# A PDB ID is exactly 4 characters: digit followed by 3 alphanumeric characters.
PDB_ID_REGEX = re.compile(r'^[0-9][A-Z0-9]{3}$', re.IGNORECASE)


def _tokenize(name: str) -> list[str]:
    """Split a filename stem into alnum tokens on any non-alphanumeric separator
    (underscore, dash, space, dot, ...). e.g. "hNET_8WTW_ions" -> ["hNET", "8WTW", "ions"]."""
    return [t for t in re.split(r'[^A-Za-z0-9]+', name) if t]


def _lookup_uniprot_from_pdb_id(pdb_id: str) -> str | None:
    """Resolve a 4-character PDB ID to its primary UniProt accession via the
    RCSB Data API. Returns None on any failure (network, missing mapping, etc.)
    so the caller can fall back gracefully instead of erroring the upload."""
    try:
        entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
        req = urllib.request.Request(entry_url, headers={"User-Agent": "ProteinVisualizeApp/1.0"})
        with urllib.request.urlopen(req, timeout=5) as resp:
            if resp.status != 200:
                return None
            entry_data = json.loads(resp.read().decode("utf-8"))

        entity_ids = entry_data.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
        if not entity_ids:
            return None

        # Check each polymer entity (not just the first) in case the primary
        # chain of interest isn't entity 1.
        for entity_id in entity_ids:
            entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id.lower()}/{entity_id}"
            req = urllib.request.Request(entity_url, headers={"User-Agent": "ProteinVisualizeApp/1.0"})
            try:
                with urllib.request.urlopen(req, timeout=5) as resp:
                    if resp.status != 200:
                        continue
                    entity_data = json.loads(resp.read().decode("utf-8"))
            except Exception:
                continue

            # Preferred: direct UniProt accession list on the entity container identifiers.
            uniprot_ids = (
                entity_data.get("rcsb_polymer_entity_container_identifiers", {}).get("uniprot_ids") or []
            )
            if uniprot_ids:
                return str(uniprot_ids[0]).upper()

            # Fallback: SIFTS sequence-alignment records.
            for align in entity_data.get("rcsb_polymer_entity_align", []) or []:
                if align.get("reference_database_name") == "UniProt":
                    acc = align.get("reference_database_accession")
                    if acc:
                        return str(acc).upper()
    except Exception:
        pass
    return None


def extract_uniprot_id(path: Path) -> str | None:
    stem = path.stem
    tokens = _tokenize(stem)

    # 1. A token that is itself a valid UniProt accession
    #    (e.g. "P31645.pdb", "hSERT_P31645_model.cif")
    for token in tokens:
        if UNIPROT_REGEX.match(token):
            return token.upper()

    # 2. DBREF / mmCIF headers inside the file content, if the uploader kept
    #    them (works even when the filename itself has no useful hints).
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
                if "UNP" in line or "UNIPROT" in line.upper():
                    for token in re.split(r'[^A-Za-z0-9]+', line):
                        if UNIPROT_REGEX.match(token):
                            return token.upper()
    except Exception:
        pass

    # 3. Any token that looks like a 4-character PDB ID (e.g. "8WTW" inside
    #    "hNET_8WTW_ions.pdb") -> resolve it to its UniProt accession via RCSB.
    for token in tokens:
        if PDB_ID_REGEX.match(token):
            resolved = _lookup_uniprot_from_pdb_id(token)
            if resolved:
                return resolved

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