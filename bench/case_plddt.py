"""Reproduce the pLDDT-filter claim on a specific structure (e.g. 9VMX)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

import app.services.topology.membrane_energy as ME
from app.schemas.topology import TMParams
from app.services.topology.providers.tm.energy import EnergyTMProvider
from app.services.topology.residues import load_residue_frame

pdb = sys.argv[1] if len(sys.argv) > 1 else "9VMX"
chain = sys.argv[2] if len(sys.argv) > 2 else "A"
path = Path("bench/results/alphafold") / f"{pdb}.cif"
if not path.exists():
    import urllib.request
    url = f"https://files.rcsb.org/download/{pdb}.cif"
    req = urllib.request.Request(url, headers={"User-Agent": "ProteinLab-benchmark/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(r.read())

frame = load_residue_frame(path, chain)
print(f"{pdb} chain {frame.chain_id}: {len(frame)} residues, "
      f"predicted_model={frame.predicted_model}")

for label, value in (("with pLDDT filter", ME.MIN_PLDDT), ("without pLDDT filter", -1.0)):
    original = ME.MIN_PLDDT
    ME.MIN_PLDDT = value
    try:
        r = EnergyTMProvider().predict_tm(path, params=TMParams(), frame=frame)
    finally:
        ME.MIN_PLDDT = original
    print(f"  {label:22s} nTM={len(r.segments or []):3d}  {r.labeler}")
    for w in r.warnings:
        print("      -", w)
