"""Case-study numbers for the manuscript (Piezo1 Q92508 and the KcsA control)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))

from app.services.topology.service import predict_topology

DATA = Path(__file__).resolve().parents[1] / "backend" / "data" / "uploads"

for fname, chain in (("8YEZ_Q92508.cif", "A"), ("8ZU3.cif", "A"),
                     ("hSERT_5I6X_ions.pdb", "A")):
    p = DATA / fname
    for algo in ("3d_energy", "3d_slab_geom"):
        r = predict_topology(p, tm_algo=algo, ss_algo="none", flow_type="consensus",
                             chain_id=chain)
        tms = [reg for reg in r.regions if reg.type == "Transmembrane"]
        print(f"{fname:22s} {algo:16s} chain={r.chain_id} nTM={len(tms)}")
        print(f"   {r.labeler}")
        for reg in tms:
            print(f"     TM {reg.start}-{reg.end} conf={reg.confidence}")
        for w in r.warnings:
            if "curved" in w or "dropped" in w:
                print("   warn:", w)
    print("-" * 90)
