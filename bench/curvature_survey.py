"""How often the spherical mid-surface is accepted, over the topology cohort.

    python bench/curvature_survey.py --out bench/results --workers 8

The implicit-membrane model fits a plane and, when the plane is clearly worse
than a sphere, a spherical mid-surface (PPM 3.0 style).  The topology benchmark
does not record which of the two was accepted, so the curvature claim of the
manuscript needs its own pass.  This script re-runs the energy provider on the
150 proteins of ``benchmark.json`` and writes, for every one of them, whether
the sphere was accepted and at which radius.

The output is ``curvature_survey.json``.  Nothing is downloaded: the structures
are the ones the topology benchmark already committed under
``bench/results/structures/``, and the seed and the cohort are the ones of that
benchmark.  Re-running is deterministic.
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

import numpy as np  # noqa: E402

from app.schemas.topology import TMParams  # noqa: E402
from app.services.topology.providers.tm.energy import EnergyTMProvider  # noqa: E402
from app.services.topology.residues import load_residue_frame  # noqa: E402

# The cohort is the one of the topology benchmark; the image of it is committed.
COHORT = "benchmark.json"
STRUCTURES = "structures"


def survey_one(job) -> dict | None:
    """Run the energy provider on one cohort member.  Module level for pickling."""
    out_dir, rec = job
    out = Path(out_dir)
    pdb_id, chain = rec["pdb_id"], rec["chain"]
    path = out / STRUCTURES / f"{pdb_id}.cif"
    if not path.exists():
        print(f"[survey] {rec['accession']} / {pdb_id}: structure missing")
        return None
    try:
        frame = load_residue_frame(path, chain or None)
    except Exception as exc:  # noqa: BLE001 - the benchmark logs and moves on
        print(f"[survey] {rec['accession']} / {pdb_id}: parse failed ({exc})")
        return None
    if len(frame) < 50:
        return None
    try:
        pred = EnergyTMProvider().predict_tm(path, params=TMParams(), frame=frame)
    except Exception as exc:  # noqa: BLE001
        print(f"[survey] {rec['accession']} / {pdb_id}: provider failed ({exc})")
        return None
    radius = getattr(pred, "membrane_radius", None)
    half = getattr(pred, "half_thickness", None)
    return {
        "accession": rec["accession"],
        "pdb_id": pdb_id,
        "chain": chain,
        "protein_name": rec.get("protein_name", ""),
        "n_residues": int(len(frame)),
        "n_tm_computed": int(len(pred.segments or [])),
        "curved": radius is not None,
        "radius_A": float(radius) if radius is not None else None,
        "dG_transfer_kcal_per_mol": round(-float(pred.membrane_score), 2),
        "hydrophobic_thickness_A": round(2 * float(half), 2) if half is not None else None,
        "placement_quality": getattr(pred, "placement_quality", None),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "bench" / "results")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    cohort = json.loads((args.out / COHORT).read_text())
    rows: list[dict] = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i, row in enumerate(pool.map(survey_one,
                                         [(str(args.out), r) for r in cohort], chunksize=4)):
            if row is not None:
                rows.append(row)
            if (i + 1) % 25 == 0:
                print(f"[survey] {i + 1}/{len(cohort)}")

    curved = [r for r in rows if r["curved"]]
    radii = [r["radius_A"] for r in curved]
    summary = {
        "cohort": COHORT,
        "n_cohort": len(rows),
        "n_curved": len(curved),
        "radius_min_A": min(radii) if radii else None,
        "radius_max_A": max(radii) if radii else None,
        "radius_median_A": float(np.median(radii)) if radii else None,
        "curved_entries": [f"{r['pdb_id']}({r['chain']}) {r['radius_A']:.0f} A "
                           f"{r['protein_name']}" for r in curved],
    }
    payload = {"summary": summary, "rows": rows}
    (args.out / "curvature_survey.json").write_text(json.dumps(payload, indent=2))

    print(f"\n[survey] {summary['n_curved']} of {summary['n_cohort']} structures accepted "
          f"the spherical mid-surface")
    for entry in summary["curved_entries"]:
        print("   ", entry)
    if radii:
        print(f"[survey] radii {summary['radius_min_A']:.0f}"
              f"-{summary['radius_max_A']:.0f} A, median "
              f"{summary['radius_median_A']:.0f} A")
    print(f"[survey] written to {args.out / 'curvature_survey.json'}")


if __name__ == "__main__":
    main()
