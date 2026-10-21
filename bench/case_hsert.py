"""Case study: human serotonin transporter (5I6X, UniProtKB P31645).

Prints the computed transmembrane segments next to the curated UniProtKB features,
both in author residue numbering, so that the manuscript can quote the comparison.
"""
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

warnings.filterwarnings("ignore")
import app.main  # noqa: E402,F401  - puts tools/bin (mkdssp) on PATH

from app.services.topology.providers.tm.uniprot import fetch_uniprot_entry  # noqa: E402
from app.services.topology.residues import load_residue_frame, map_to_reference  # noqa: E402
from app.services.topology.service import predict_topology  # noqa: E402

path = ROOT / "backend" / "data" / "uploads" / "hSERT_5I6X_ions.pdb"
acc = "P31645"

frame = load_residue_frame(path, "A")
entry = fetch_uniprot_entry(acc)
useq = entry["sequence"]["value"]
mapping = map_to_reference(frame, [(u, "") for u in range(1, len(useq) + 1)], list(useq),
                            offsets=[0], source="UniProt", allow_order_fallback=False)
pos_of_u = {j + 1: p for p, j in enumerate(mapping.ref_index) if j is not None}

curated = []
for f in entry.get("features", []):
    if f.get("type") != "Transmembrane":
        continue
    u0 = int(f["location"]["start"]["value"])
    u1 = int(f["location"]["end"]["value"])
    hits = [pos_of_u[u] for u in range(u0, u1 + 1) if u in pos_of_u]
    if hits:
        lo, hi = min(hits), max(hits)
        # report in author residue numbering, like every region of the result
        curated.append((frame.residues[lo].resseq, frame.residues[hi].resseq,
                        f.get("description", "")))

r = predict_topology(path, tm_algo="3d_energy", ss_algo="dssp", flow_type="consensus",
                     chain_id="A")
computed = []
for reg in r.regions:
    if reg.type == "Transmembrane":
        computed.append((reg.start, reg.end, reg.confidence))

print(r.labeler)
print(f"computed {len(computed)} TM segments; curated {len(curated)} TM helices\n")
print(f"{'computed (author)':>22s}   {'curated (author)':>22s}   dstart  dend")
for i in range(max(len(computed), len(curated))):
    c = f"{computed[i][0]}-{computed[i][1]}" if i < len(computed) else "-"
    u = f"{curated[i][0]}-{curated[i][1]}" if i < len(curated) else "-"
    ds = dd = ""
    if i < len(computed) and i < len(curated):
        ds = str(computed[i][0] - curated[i][0])
        dd = str(computed[i][1] - curated[i][1])
    print(f"{c:>22s}   {u:>22s}   {ds:>6s} {dd:>5s}")
print()
print("warnings:")
for w in r.warnings:
    print("  -", w)
