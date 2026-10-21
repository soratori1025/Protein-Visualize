"""
Paired comparison of the topology providers on the evaluation cohort.

    python bench/compare_methods.py --results bench/results

For every pair of providers the per-protein F1 values are compared with a paired
Wilcoxon signed-rank test (two-sided) and a bootstrap confidence interval for the
mean difference (10 000 resamples of the proteins, seeded).  The per-protein pairing
matters: all providers see the same structures and the same reference topology, so
the protein is the natural experimental unit and the comparison is paired.

Because the comparison produces one test per pair of providers, the raw *P* values
are also corrected with the Benjamini-Hochberg procedure, and the correction is
reported next to them.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

ORDER = ["3d_energy", "3d_slab_geom", "tmhmm", "kyte_doolittle_seq", "toppred",
         "topcons_style_consensus", "uniprot_api"]
BOOTSTRAP = 10000


def per_protein_f1(records, method):
    out = {}
    for r in records:
        m = r.get("methods", {}).get(method) or {}
        if "f1" in m:
            out[r["accession"]] = m["f1"]
    return out


def compare(a: dict, b: dict, rng):
    keys = sorted(set(a) & set(b))
    if len(keys) < 3:
        return None
    da = np.array([a[k] for k in keys], float)
    db = np.array([b[k] for k in keys], float)
    d = da - db
    nz = d[d != 0]
    if len(nz) == 0:
        p = 1.0
    else:
        try:
            p = float(wilcoxon(da, db, zero_method="wilcox", alternative="two-sided").pvalue)
        except ValueError:
            p = float("nan")
    boot = np.array([rng.choice(d, size=len(d), replace=True).mean() for _ in range(BOOTSTRAP)])
    return {
        "n_paired": len(keys),
        "mean_difference": round(float(d.mean()), 4),
        "median_difference": round(float(np.median(d)), 4),
        "ci95_low": round(float(np.percentile(boot, 2.5)), 4),
        "ci95_high": round(float(np.percentile(boot, 97.5)), 4),
        "n_better": int((d > 0).sum()),
        "n_worse": int((d < 0).sum()),
        "n_tied": int((d == 0).sum()),
        "wilcoxon_p": None if np.isnan(p) else round(p, 4),
    }


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg step-up adjusted *P* values (monotone, clipped at 1)."""
    n = len(p_values)
    if n == 0:
        return []
    order = np.argsort(p_values)
    adjusted = np.empty(n, float)
    running = 1.0
    for rank in range(n - 1, -1, -1):
        index = order[rank]
        value = p_values[index] * n / (rank + 1)
        running = min(running, value)
        adjusted[index] = min(1.0, running)
    return adjusted.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="bench/results")
    ap.add_argument("--seed", type=int, default=20261021)
    args = ap.parse_args()
    results = Path(args.results)
    records = json.loads((results / "benchmark.json").read_text())
    rng = np.random.default_rng(args.seed)

    present = [m for m in ORDER
               if any(m in (r.get("methods") or {}) for r in records)]
    f1 = {m: per_protein_f1(records, m) for m in present}
    print(f"{'method':26s} {'n proteins':>10s} {'mean F1':>9s} {'median F1':>10s}")
    for m in present:
        vals = np.array(list(f1[m].values()), float)
        if vals.size:
            print(f"{m:26s} {vals.size:10d} {vals.mean():9.4f} {np.median(vals):10.4f}")

    print("\npaired comparisons (F1 difference, method A minus method B)")
    out = []
    for i, a in enumerate(present):
        for b in present[i + 1:]:
            res = compare(f1[a], f1[b], rng)
            if res is None:
                continue
            res["method_a"] = a
            res["method_b"] = b
            out.append(res)

    # multiplicity: one test per pair -> Benjamini-Hochberg over the whole family
    raw = [r["wilcoxon_p"] for r in out]
    adjusted = benjamini_hochberg([p if p is not None else 1.0 for p in raw])
    for r, q in zip(out, adjusted):
        r["wilcoxon_p_bh"] = round(float(q), 4)

    for r in out:
        print(f"  {r['method_a']:24s} - {r['method_b']:24s} "
              f"dF1={r['mean_difference']:+.4f} "
              f"[{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}]  "
              f"W={r['n_better']}/{r['n_worse']}/{r['n_tied']}  "
              f"P={r['wilcoxon_p']}  P(BH)={r['wilcoxon_p_bh']}")

    (results / "paired_comparisons.json").write_text(json.dumps(out, indent=2))
    print(f"\nwritten to {results / 'paired_comparisons.json'}")


if __name__ == "__main__":
    main()
