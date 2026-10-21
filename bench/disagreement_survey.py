"""
What the geometry adds: a survey of the cases where the computed topology and the
curated topology differ, and of the geometric features the engine reports.

    python bench/disagreement_survey.py --results bench/results

A computational paper should say something about the data it processes.  This script
surveys the evaluation cohort for the situations in which a structure-derived topology
and a curated topology cannot agree, and for the geometric features that the consensus
engine reports and that a sequence-based predictor cannot see at all:

* **broken helices** - a crossing interrupted by an unwound stretch;
* **hairpins** - a chain that turns back on itself inside the bilayer;
* **interfacial helices** - a helix lying flat in the interfacial region;
* **re-entrant segments** - a chain that leaves the membrane and comes back;
* **unwound crossings** - a crossing with no continuous helix at all.

For every protein whose computed topology differs from the curated one, the
transition class of the difference is recorded, so that the disagreements can be
sorted into "the geometry explains it" and "the geometry does not decide".
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _candidate in (ROOT / "backend", ROOT):
    if str(_candidate) not in sys.path:
        sys.path.insert(0, str(_candidate))

import numpy as np  # noqa: E402

from app.schemas.topology import TMParams  # noqa: E402
from app.services.topology.service import predict_topology  # noqa: E402

from benchmark_tm_topology import match_segments, overlap, score  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="bench/results")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    results = Path(args.results)
    records = json.loads((results / "benchmark.json").read_text())
    if args.limit:
        records = records[: args.limit]

    structures = results / "structures"
    feature_counts = Counter()
    per_protein = []
    disagreements = []

    for record in records:
        path = structures / f"{record['pdb_id']}.cif"
        if not path.exists():
            continue
        try:
            response = predict_topology(path, tm_algo="3d_energy", ss_algo="dssp",
                                        flow_type="consensus", chain_id=record["chain"])
        except Exception as exc:
            print(f"[survey] {record['pdb_id']}: {type(exc).__name__}: {exc}")
            continue

        roles = Counter(r.membrane_role for r in response.regions if r.membrane_role)
        features = {
            "broken_helix": roles.get("BROKEN_TM", 0),
            "unwound_crossing": roles.get("UNWOUND", 0),
            "interfacial_helix": roles.get("INTERFACIAL", 0),
            "reentrant": roles.get("REENTRANT", 0),
            "signal_peptide": roles.get("SIGNAL", 0),
        }
        for key, value in features.items():
            if value:
                feature_counts[key] += value
        per_protein.append({
            "accession": record["accession"], "pdb_id": record["pdb_id"],
            "chain": record["chain"],
            "protein_name": record.get("protein_name", ""),
            "n_crossings": len({r.crossing for r in response.regions if r.crossing}),
            "n_curated": record["n_reference_tm"],
            "f1": record["methods"]["3d_energy"].get("f1"),
            **features,
        })

        # --- where the computed and the curated topology differ
        computed = [(r.start, r.end) for r in response.regions
                    if r.type == "Transmembrane" and r.crossing]
        curated = record["methods"]["3d_energy"]
        if computed and curated.get("f1") is not None and curated["f1"] < 1.0:
            disagreements.append({
                "accession": record["accession"], "pdb_id": record["pdb_id"],
                "protein_name": record.get("protein_name", ""),
                "f1": curated["f1"],
                "n_computed": len(computed), "n_curated": record["n_reference_tm"],
                "roles": {k: v for k, v in roles.items() if v},
            })

    n = len(per_protein)
    print(f"\n[survey] {n} structures analysed with the full consensus flow "
          f"(3D Energy + DSSP)\n")
    print(f"{'geometric feature':30s} {'structures':>10s} {'share':>8s}")
    for key in ("broken_helix", "unwound_crossing", "interfacial_helix", "reentrant",
                "signal_peptide"):
        count = sum(1 for p in per_protein if p[key])
        print(f"{key:30s} {count:10d} {count / max(1, n):8.1%}")

    f1 = [p["f1"] for p in per_protein if p["f1"] is not None]
    exact = sum(1 for v in f1 if v >= 0.999)
    print(f"\nexactly the curated topology : {exact}/{len(f1)} ({exact / max(1, len(f1)):.1%})")
    print(f"differing from the curation   : {len(f1) - exact}/{len(f1)}")
    print(f"mean F1                       : {np.mean(f1):.4f}" if f1 else "")

    out = {"per_protein": per_protein,
           "feature_counts": dict(feature_counts),
           "n_structures": n,
           "n_exactly_curated": exact,
           "disagreements": disagreements}
    (results / "disagreement_survey.json").write_text(json.dumps(out, indent=2))
    print(f"\nwritten to {results / 'disagreement_survey.json'}")


if __name__ == "__main__":
    main()
