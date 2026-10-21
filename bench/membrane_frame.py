"""
Membrane-frame statistics for the evaluation cohort.

    python bench/membrane_frame.py --results bench/results

Re-runs the 3D Energy provider over the coordinate files already downloaded by the
benchmark and records, for each protein: the transfer free energy, the fitted
hydrophobic thickness, whether a curved (spherical) mid-surface was accepted and its
radius, and the number of chains that took part in the fit.  Writes
`membrane_frame.json` next to the benchmark results.
"""
from __future__ import annotations

import argparse
import itertools
import json
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

import numpy as np  # noqa: E402

from app.schemas.topology import TMParams  # noqa: E402
from app.services.topology.providers.tm.energy import EnergyTMProvider  # noqa: E402
from app.services.topology.residues import load_residue_frame  # noqa: E402

CURVED_RE = re.compile(r"curved R=([0-9.]+) A")
THICK_RE = re.compile(r"hydrophobic thickness ([0-9.]+) A")
CHAINS_RE = re.compile(r"([0-9]+) chains")


def one(job):
    acc, pdb_id, chain, path = job
    try:
        frame = load_residue_frame(Path(path), chain or None)
        pred = EnergyTMProvider().predict_tm(Path(path), params=TMParams(), frame=frame)
    except Exception as exc:
        return {"accession": acc, "pdb_id": pdb_id, "error": f"{type(exc).__name__}: {exc}"}
    labeler = pred.labeler or ""
    m = CURVED_RE.search(labeler)
    t = THICK_RE.search(labeler)
    c = CHAINS_RE.search(labeler)
    return {
        "accession": acc, "pdb_id": pdb_id, "chain": chain,
        "dG_transfer": round(-float(pred.membrane_score), 2),
        "hydrophobic_thickness_A": float(t.group(1)) if t else None,
        "curved": bool(m),
        "curvature_radius_A": float(m.group(1)) if m else None,
        "n_chains_in_fit": int(c.group(1)) if c else 1,
        "n_tm_segments": len(pred.segments or []),
        "labeler": labeler,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="bench/results")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    results = Path(args.results)
    records = json.loads((results / "benchmark.json").read_text())
    jobs = []
    for r in records:
        path = results / "structures" / f"{r['pdb_id']}.cif"
        if path.exists():
            jobs.append((r["accession"], r["pdb_id"], r["chain"], str(path)))
    print(f"[frame] {len(jobs)} structures")

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(one, jobs))
    rows = [r for r in rows if "error" not in r]
    (results / "membrane_frame.json").write_text(json.dumps(rows, indent=2))

    dg = np.array([r["dG_transfer"] for r in rows], float)
    th = np.array([r["hydrophobic_thickness_A"] for r in rows if r["hydrophobic_thickness_A"]], float)
    curved = [r for r in rows if r["curved"]]
    print(f"dG_transfer      n={dg.size}  min={dg.min():.2f}  max={dg.max():.2f}  "
          f"median={np.median(dg):.2f}  mean={dg.mean():.2f}")
    print(f"thickness (A)    n={th.size}  min={th.min():.1f}  max={th.max():.1f}  "
          f"median={np.median(th):.1f}  mean={th.mean():.1f}")
    print(f"curved           {len(curved)} / {len(rows)}")
    for r in curved:
        print(f"    {r['accession']} {r['pdb_id']} chain {r['chain']} "
              f"R={r['curvature_radius_A']} A  dG={r['dG_transfer']}")
    multi = [r for r in rows if r["n_chains_in_fit"] > 1]
    print(f"multi-chain fits {len(multi)} / {len(rows)}")


if __name__ == "__main__":
    main()
