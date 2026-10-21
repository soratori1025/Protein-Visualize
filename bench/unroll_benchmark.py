"""
Validation of the 2-D unrolling.

    python bench/unroll_benchmark.py --results bench/results

The unrolling is one of the three claimed contributions of the paper and, until now,
had no quantitative validation.  This harness runs it over the evaluation cohort and
reports the geometric distortion it introduces:

  * the closure error of every inter-transmembrane loop (RMSD of the closed end point
    from the native one, in angstrom);
  * the number of loops that could not be closed and were drawn schematically;
  * the path deviation of the rebuilt loops from their guide curve;
  * the heavy-atom clashes the layout introduced;
  * the peptide bonds it stretched and the disulfides it broke;
  * the longest C-N bond in the output.

The unrolling itself is TypeScript and runs in the browser; it is executed here by the
*production module* under Node's native TypeScript type stripping, so what is measured
is what ships.  The rigid 'aligned' mode is run on the same inputs as a control: it
deforms nothing, so every quantity it reports is the input's own.

The consensus map that the unrolling needs comes from the same backend call the
application makes, so the coordinate frame, the chain and the crossings are the ones a
user would see.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

import numpy as np  # noqa: E402

import app.main  # noqa: F401,E402  - puts tools/bin (mkdssp, stride) on PATH
from app.schemas.topology import TMParams  # noqa: E402
from app.services.topology.service import predict_topology  # noqa: E402

UA = {"User-Agent": "ProteinLab-benchmark/1.0 (academic validation)"}
NODE = "node"


def download(url: str, dest: Path) -> bool:
    if dest.exists() and dest.stat().st_size > 1000:
        return True
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=180) as resp:
            data = resp.read()
    except Exception:
        return False
    if len(data) < 1000:
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    return True


def consensus_map(response):
    """The per-residue consensus map in the shape the frontend expects."""
    out = []
    for m in response.consensus_map:
        out.append({
            "residue_number": m.residue_number,
            "insertion_code": m.insertion_code or "",
            "aa": m.aa,
            "label": m.label,
            "crossing": m.crossing,
            "part": m.part,
        })
    return out


def run_one(job_path: Path, out_path: Path) -> bool:
    # --experimental-transform-types (not --experimental-strip-types): the production
    # module uses TypeScript parameter properties, which the strip-only mode rejects.
    cmd = [NODE, "--experimental-transform-types",
           str(ROOT / "bench" / "unroll_validate.ts"), str(job_path), str(out_path)]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    except Exception:
        return False
    if res.returncode != 0:
        tail = [l for l in (res.stderr or "").splitlines() if l.strip()]
        print("   node failed:", tail[-1] if tail else "?")
        return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="bench/results")
    ap.add_argument("--limit", type=int, default=40)
    args = ap.parse_args()
    results = Path(args.results)
    records = json.loads((results / "benchmark.json").read_text())[: args.limit]
    structs = results / "structures"
    structs.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="unroll_"))

    rows = []
    for r in records:
        pdb = r["pdb_id"]
        path = structs / f"{pdb}.cif"
        if not download(f"https://files.rcsb.org/download/{pdb}.cif", path):
            print(f"[unroll] {pdb}: download failed")
            continue
        try:
            resp = predict_topology(path, tm_algo="3d_energy", ss_algo="dssp",
                                    flow_type="consensus", chain_id=r["chain"])
        except Exception as exc:
            print(f"[unroll] {pdb}: topology failed ({exc})")
            continue
        cmap = consensus_map(resp)
        crossings = {m["crossing"] for m in cmap if m.get("crossing") is not None}
        if len(crossings) < 2:
            continue
        mem = resp.membrane
        membrane = None
        if isinstance(mem, dict):
            normal = mem.get("normal")
            half = mem.get("half_thickness")
            if normal and half:
                membrane = {"normal": list(normal), "half_thickness": float(half)}
        job = {
            "structure": str(path),
            "chainId": r["chain"],
            "consensusMap": cmap,
            "membrane": membrane,
        }
        job_path = tmp / f"{pdb}.job.json"
        out_path = tmp / f"{pdb}.out.json"
        job_path.write_text(json.dumps(job), encoding="utf-8")
        if not run_one(job_path, out_path):
            continue
        rep = json.loads(out_path.read_text())
        u = rep["unrolled"]
        a = rep["aligned"]
        row = {
            "pdb_id": pdb, "accession": r["accession"], "chain": r["chain"],
            "n_crossings": len(crossings),
            "n_loops": len(u["loops"]),
            "n_loops_closed": sum(1 for l in u["loops"] if l.get("closed")),
            "n_loops_schematic": sum(1 for l in u["loops"] if l.get("schematic")),
            "closure_rmsd": [l.get("closureRmsd") for l in u["loops"]
                             if l.get("closureRmsd") is not None],
            "path_deviation": [l.get("pathDeviation") for l in u["loops"]
                               if l.get("pathDeviation") is not None],
            "clashes": u["clashes"],
            "broken_peptide_bonds": u["brokenPeptideBonds"],
            "broken_disulfides": u["brokenDisulfides"],
            "max_peptide_bond": u["maxPeptideBond"],
            "n_warnings": len(u["warnings"]),
            "aligned_max_peptide_bond": a["maxPeptideBond"],
            "aligned_clashes": a["clashes"],
            "aligned_broken_peptide_bonds": a["brokenPeptideBonds"],
        }
        rows.append(row)
        print(f"[unroll] {pdb} {r['chain']}: {row['n_crossings']} crossings, "
              f"{row['n_loops']} loops ({row['n_loops_schematic']} schematic), "
              f"clashes {row['clashes']}, max C-N {row['max_peptide_bond']}")
        # written after every structure so that a long run is not lost
        (results / "unroll_benchmark.json").write_text(json.dumps(rows, indent=2))

    (results / "unroll_benchmark.json").write_text(json.dumps(rows, indent=2))

    all_rms = [v for row in rows for v in row["closure_rmsd"]]
    all_dev = [v for row in rows for v in row["path_deviation"]]
    n_loops = sum(r["n_loops"] for r in rows)
    n_schematic = sum(r["n_loops_schematic"] for r in rows)
    n_closed = sum(r["n_loops_closed"] for r in rows)
    print(f"\nstructures analysed      : {len(rows)}")
    print(f"loops                    : {n_loops}")
    print(f"  closed                 : {n_closed}")
    print(f"  schematic (not closed) : {n_schematic}")
    if all_rms:
        print(f"closure RMSD (A)         : median {np.median(all_rms):.2f}, "
              f"mean {np.mean(all_rms):.2f}, p90 {np.percentile(all_rms, 90):.2f}, "
              f"max {max(all_rms):.2f}")
        print(f"  within 0.3 A           : {sum(1 for v in all_rms if v <= 0.3)}/{len(all_rms)}")
        print(f"  within 1.0 A           : {sum(1 for v in all_rms if v <= 1.0)}/{len(all_rms)}")
    if all_dev:
        print(f"path deviation (A)       : median {np.median(all_dev):.2f}, "
              f"p90 {np.percentile(all_dev, 90):.2f}")
    if rows:
        print(f"new clashes per structure: median "
              f"{np.median([r['clashes'] for r in rows]):.1f}")
        print(f"stretched bonds          : {sum(r['broken_peptide_bonds'] for r in rows)}")
        print(f"broken disulfides        : {sum(r['broken_disulfides'] for r in rows)}")


if __name__ == "__main__":
    main()
