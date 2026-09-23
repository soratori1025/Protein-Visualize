from pathlib import Path

from fastapi import APIRouter, HTTPException, Query

from app.api.structure import build_structure_summary
from protein_engine.analysis.metrics import analyze_chain

router = APIRouter(prefix="/api/analysis", tags=["analysis"])
ROOT = Path(__file__).resolve().parents[2]


@router.get("/{filename}/chain/{chain_id}")
def analyze_structure_chain(
    filename: str,
    chain_id: str,
    contact_cutoff: float = Query(default=8.0, ge=3.0, le=20.0),
) -> dict:
    path = ROOT / "data" / "uploads" / Path(filename).name
    if not path.is_file():
        raise HTTPException(status_code=404, detail="Uploaded structure not found")
    summary = build_structure_summary(path)
    chain = next((item for item in summary[0]["chains"] if item["id"] == chain_id), None)
    if chain is None:
        raise HTTPException(status_code=404, detail=f"Chain {chain_id} not found")
    return {"input": {"filename": path.name, "chain_id": chain_id}, "result": analyze_chain(chain, contact_cutoff)}