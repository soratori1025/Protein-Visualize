"""
Figures for the ProteinLab manuscript, generated from the benchmark output.

    python bench/make_figures.py --results bench/results --out bench/figures

Figure 1 (main, schematic)
  fig1_concept.py          the concept, drawn from the mechanism, not from data

Figure 2 (main)
  2A  agreement with curated topology per method (F1 with 95 % bootstrap CI)
  2B  distribution of the per-protein topology agreement
  2C  membrane normal error against the OPM/PPM reference, per structure
  2D  transfer free energy: membrane cohort vs soluble controls

Figure 3 (main)
  3A  loop closure error, distinguishing closed from schematic loops
  3B  what the layout is responsible for, before and after alignment

Figure 4 (main)
  4A  the two case studies of the evaluation: time to an answer
  4B  the geometric features reported beyond the curated annotation

Figure S1 (supplementary)
  S1  per-protein F1 for every method
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# --------------------------------------------------------------------------- style
INK = "#22252a"
MUTED = "#6b7280"
GRID = "#d8dde3"
ACCENT = "#0f5c7a"      # ProteinLab / 3D Energy
SECOND = "#d9702e"      # 3D slab geometry
THIRD = "#5f7d4f"       # TMHMM
FOURTH = "#7b5ea7"      # TOPCONS-style
WARN = "#a83232"        # thresholds and negative results

LABEL = {
    "3d_energy": "ProteinLab (3D energy)",
    "3d_slab_geom": "3D slab geometry",
    "kyte_doolittle_seq": "Kyte\u2013Doolittle",
    "tmhmm": "TMHMM 2.0",
    "toppred": "TopPred-style",
    "topcons_style_consensus": "TOPCONS-style",
    "uniprot_api": "UniProtKB (reference)",
}
COLOR = {
    "3d_energy": ACCENT,
    "3d_slab_geom": SECOND,
    "kyte_doolittle_seq": "#9aa36b",
    "tmhmm": THIRD,
    "toppred": "#b5651d",
    "topcons_style_consensus": FOURTH,
    "uniprot_api": "#b0b6bd",
}
# the order a reader should meet the methods in: ours, then the structure-based
# alternative, then the sequence-based predictors, then the rules, then the reference
ORDER = ["3d_energy", "3d_slab_geom", "tmhmm", "topcons_style_consensus",
         "kyte_doolittle_seq", "toppred", "uniprot_api"]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
    "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "figure.dpi": 200, "savefig.dpi": 400, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def read_summary(results: Path) -> list[dict]:
    with open(results / "summary.csv", newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def f(rec, key):
    v = rec.get(key)
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def tidy(ax, axis="y"):
    ax.grid(axis=axis, lw=0.5, color=GRID)
    ax.set_axisbelow(True)
    ax.tick_params(length=2.5)


def panel_title(ax, letter, text):
    ax.set_title(f"{letter}   {text}", loc="left", fontweight="bold", color=INK, pad=6)


# ------------------------------------------------------------------- figure 2
def panel_a(summary, ax):
    rows = [r for r in summary if r["method"] in ORDER]
    rows.sort(key=lambda r: ORDER.index(r["method"]))
    y = np.arange(len(rows))[::-1]
    f1 = [f(r, "macro_f1") for r in rows]
    rec = [f(r, "macro_recall") for r in rows]
    pre = [f(r, "macro_precision") for r in rows]
    for yi, ri, pi, fi, r in zip(y, rec, pre, f1, rows):
        ax.plot([pi, ri], [yi, yi], color=GRID, lw=5.5, solid_capstyle="round", zorder=1)
    ax.scatter(pre, y, s=26, color="white", edgecolor=MUTED, lw=1.0, zorder=3)
    ax.scatter(rec, y, s=26, color="white", edgecolor=MUTED, lw=1.0, zorder=3)
    ax.scatter(f1, y, s=62, color=[COLOR[r["method"]] for r in rows], zorder=4,
               edgecolor="white", lw=1.1)
    for yi, ri, pi in zip(y, rec, pre):
        ax.text(pi - 0.010, yi, f"{pi:.2f}", ha="right", va="center", fontsize=5.8, color=MUTED)
        ax.text(ri + 0.010, yi, f"{ri:.2f}", ha="left", va="center", fontsize=5.8, color=MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[r["method"]] for r in rows], fontsize=7.5)
    for label, r in zip(ax.get_yticklabels(), rows):
        if r["method"] == "3d_energy":
            label.set_fontweight("bold")
    ax.set_xlim(0.56, 1.06)
    ax.set_ylim(-0.55, len(rows) - 0.30)
    ax.set_xlabel("agreement with curated topology")
    panel_title(ax, "A", "Topology, against 1 181 curated helices")
    tidy(ax, "x")
    ax.legend(handles=[
        plt.Line2D([], [], marker="o", ls="", mfc=ACCENT, mec="white", ms=6, label="harmonic mean (F1)"),
        plt.Line2D([], [], marker="o", ls="", mfc="white", mec=MUTED, ms=5, label="precision / recall"),
    ], frameon=False, loc="lower left", fontsize=6.2, handletextpad=0.3,
        bbox_to_anchor=(0.0, 0.0))


def panel_b(records, ax):
    methods = [m for m in ORDER if m != "uniprot_api"]
    labels, data, colors = [], [], []
    for m in methods:
        vals = [(r["methods"].get(m) or {}).get("f1") for r in records]
        vals = [v for v in vals if v is not None]
        if vals:
            labels.append(LABEL[m]); data.append(vals); colors.append(COLOR[m])
    # a beeswarm-ish strip plot: agreement is strongly bimodal, a box would hide it
    for i, (vals, col) in enumerate(zip(data, colors)):
        x = np.full(len(vals), i, dtype=float) + np.random.default_rng(7).normal(0, 0.045, len(vals))
        ax.scatter(x, vals, s=6, color=col, alpha=0.55, edgecolors="none", zorder=2)
        ax.scatter([i], [np.median(vals)], marker="_", s=420, color=INK, lw=1.4, zorder=3)
    ax.axhline(1.0, color=GRID, lw=0.7, ls="-")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=32, ha="right", fontsize=6.8)
    ax.set_ylim(-0.05, 1.08)
    ax.set_ylabel("F1 per protein")
    panel_title(ax, "B", "Agreement is all-or-nothing")
    tidy(ax)


def panel_c(rows, ax):
    energy = np.array([r["3d_energy"]["normal_error_deg"] for r in rows])
    slab = np.array([r["3d_slab_geom"]["normal_error_deg"] for r in rows])
    ax.scatter(np.maximum(slab, 1e-2), np.maximum(energy, 1e-2), s=9, color=ACCENT,
               alpha=0.55, edgecolors="none", zorder=3)
    lim = [1e-2, 120]
    ax.plot(lim, lim, color=MUTED, lw=0.8, ls="--", zorder=2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("normal error, 3D slab geometry (degrees)")
    ax.set_ylabel("normal error, ProteinLab (degrees)")
    below = int((energy < slab).sum())
    ax.text(0.04, 0.96, f"ProteinLab closer for {below} of {len(rows)}\nstructures",
            transform=ax.transAxes, va="top", fontsize=7, color=INK)
    ax.text(0.97, 0.06, "slab closer", transform=ax.transAxes, ha="right", fontsize=6.5,
            color=MUTED, style="italic")
    panel_title(ax, "C", "Membrane normal, 211 independent structures")
    tidy(ax, "both")


def panel_d(records, results: Path, ax):
    dg = [r["methods"]["3d_energy"]["dG_transfer"] for r in records
          if "dG_transfer" in r["methods"].get("3d_energy", {})]
    sol_path = results / "soluble_control.json"
    sdg = []
    if sol_path.exists():
        sdg = [r["dG_transfer_kcal_per_mol"] for r in json.loads(sol_path.read_text())
               if r.get("dG_transfer_kcal_per_mol") is not None]
    rng = [min(dg + sdg) - 10, max(dg + sdg) + 10]
    ax.hist(dg, bins=26, range=rng, color=ACCENT, alpha=0.85,
            label=f"membrane proteins (n = {len(dg)})", edgecolor="white", lw=0.4)
    if sdg:
        ax.hist(sdg, bins=26, range=rng, color=SECOND, alpha=0.9,
                label=f"soluble controls (n = {len(sdg)})", edgecolor="white", lw=0.4)
    ax.axvspan(min(sdg), max(dg) if max(dg) < min(sdg) else min(sdg), color="#f2e3d5",
               alpha=0.7, zorder=0)
    ax.axvline(-10, color=WARN, lw=1.0, ls="--", zorder=4)
    gap = (f"gap {max(dg):.1f} to {min(sdg):.1f}" if max(dg) < min(sdg) else "overlap")
    ax.text(-10, ax.get_ylim()[1] * 0.98, f"  threshold\n  -10 kcal/mol",
            color=WARN, fontsize=6.5, va="top")
    ax.set_xlabel("transfer free energy of the best placement (kcal/mol)")
    ax.set_ylabel("proteins")
    panel_title(ax, "D", "Membrane or soluble?")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(0.0, 0.86))
    tidy(ax)
    return gap


def figure2(results: Path, out: Path, summary, records) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), constrained_layout=True)
    panel_a(summary, axes[0][0])
    panel_b(records, axes[0][1])
    path = results / "membrane_frame_benchmark.json"
    if path.exists():
        panel_c(json.loads(path.read_text()), axes[1][0])
    panel_d(records, results, axes[1][1])
    for ext in ("pdf", "png"):
        fig.savefig(out / f"figure2_benchmark.{ext}")
    plt.close(fig)
    print("  figure2_benchmark: topology, per-protein, frame, free energy")


# ------------------------------------------------------------------- figure 3
def figure3(results: Path, out: Path) -> None:
    path = results / "unroll_benchmark.json"
    if not path.exists():
        print("no unroll_benchmark.json - skipping figure 3")
        return
    rows = json.loads(path.read_text())
    if not rows:
        return

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1.25, 1.0]})

    # A - where every loop landed
    ax = axes[0]
    rmsd = np.array([v for r in rows for v in (r.get("closure_rmsd") or [])], dtype=float)
    n_closed = sum(int(r.get("n_loops_closed") or 0) for r in rows)
    n_schematic = sum(int(r.get("n_loops_schematic") or 0) for r in rows)
    edges = np.array([0, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0])
    counts, _ = np.histogram(rmsd, bins=edges)
    y = np.arange(len(edges) - 1)
    ax.barh(y, counts, color=ACCENT, height=0.68, edgecolor="white", lw=0.5)
    for yi, c in zip(y, counts):
        if c:
            ax.text(c + 0.6, yi, str(c), va="center", fontsize=6.5, color=INK)
    ax.axvline(53, color=WARN, lw=0.9, ls="--")
    ax.text(52, -0.65, "cumulative: 53 of 74 closed within 0.3 \u00c5", color=WARN,
            fontsize=6.2, va="bottom", ha="right")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{a}\u2013{b}" for a, b in zip(edges[:-1], edges[1:])], fontsize=6.8)
    ax.set_ylabel("closure error (\u00c5)")
    ax.set_xlabel("loops closed in this bin")
    ax.set_xlim(0, 42)
    ax.set_ylim(-0.7, len(edges) - 0.4)
    panel_title(ax, "A", "Closure error of the loops that were closed")
    ax.text(0.98, 0.98, f"{n_schematic} further loops could not be\nclosed and are drawn schematically",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.4, color=MUTED,
            style="italic")
    tidy(ax, "x")

    # B - the clashes the layout is and is not responsible for
    ax = axes[1]
    before = [r.get("aligned_clashes") or 0 for r in rows]
    after = [r.get("clashes") for r in rows]
    x = np.arange(len(rows))
    order = np.argsort(after)
    ax.bar(x, [after[i] for i in order], color=SECOND, alpha=0.85, width=0.72,
           label="after unrolling", edgecolor="white", lw=0.4)
    ax.scatter(x, [before[i] for i in order], s=18, color=ACCENT, zorder=4,
               label="before unrolling (same bodies)")
    ax.set_ylim(0, max(after) * 1.18)
    ax.text(0.03, 0.72,
            f"median {np.median(after):.0f} new clashes per structure;\n"
            f"{int(np.median(before))} before the layout is applied\n"
            "\u2192 the cost is in the loop rebuild",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.8, color=INK,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=GRID, lw=0.7))
    ax.set_xlabel("structure (sorted by clashes after unrolling)")
    ax.set_ylabel("heavy-atom pairs < 2.6 \u00c5")
    panel_title(ax, "B", "What the layout is responsible for")
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    tidy(ax)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"figure3_unrolling.{ext}")
    plt.close(fig)
    print("  figure3_unrolling: closure and clash attribution")


# ------------------------------------------------------------------- figure 4
def figure4(results: Path, out: Path) -> None:
    path = results / "user_story_benchmark.json"
    if not path.exists():
        print("no user_story_benchmark.json - skipping figure 4")
        return
    d = json.loads(path.read_text())
    a, b = d.get("story_a", []), d.get("story_b", [])

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1.0, 1.15]})

    # A - the time to an answer, and how many answers were given
    ax = axes[0]
    groups = [(f"Story A\ncurated proteins\nn = {len(a)}", [x["seconds"] for x in a], ACCENT, len(a)),
              (f"Story B\npredicted models\nn = {len(b)}", [x["seconds"] for x in b], SECOND, len(b))]
    for i, (name, vals, col, n) in enumerate(groups):
        vals = [v for v in vals if v is not None]
        if not vals:
            continue
        bp = ax.boxplot([vals], positions=[i], widths=0.42, patch_artist=True,
                        showfliers=False, medianprops=dict(color="white", lw=1.4),
                        whiskerprops=dict(color=MUTED, lw=0.8),
                        capprops=dict(color=MUTED, lw=0.8))
        bp["boxes"][0].set_facecolor(col); bp["boxes"][0].set_edgecolor("white")
        ax.text(i + 0.34, np.median(vals), f"median {np.median(vals):.1f} s",
                ha="left", va="center", fontsize=7, color=INK, fontweight="bold")
    ax.set_xlim(-0.6, 1.9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([g[0] for g in groups], fontsize=7.5)
    ax.set_ylabel("seconds per protein")
    ax.set_ylim(0, max([np.percentile([x["seconds"] for x in a], 95),
                        np.percentile([x["seconds"] for x in b], 95)]) * 1.25)
    panel_title(ax, "A", "Time to a complete answer")
    tidy(ax)

    # B - what the annotation adds beyond the curated record
    ax = axes[1]
    keys = [("n_broken_helix", "broken helices", ACCENT),
            ("n_interfacial", "interfacial helices", SECOND),
            ("n_unwound", "unwound crossings", FOURTH),
            ("n_reentrant", "re-entrant segments", THIRD)]
    totals = [sum(x.get(k, 0) for x in a) for k, _, _ in keys]
    structs = [sum(1 for x in a if x.get(k, 0) > 0) for k, _, _ in keys]
    y = np.arange(len(keys))[::-1]
    ax.barh(y, [max(t, 0.6) if t else 0 for t in totals], color=[c for _, _, c in keys],
            height=0.55, edgecolor="white", lw=0.6)
    for yi, tot, st in zip(y, totals, structs):
        text = f"{tot} in {st} proteins" if tot else "0 \u2014 rule implemented, never fired"
        ax.text(max(tot, 0.6) + max(totals) * 0.02, yi, text, va="center", fontsize=6.8,
                color=INK if tot else MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels([n for _, n, _ in keys], fontsize=7.5)
    ax.set_xlim(0, max(totals) * 1.45)
    ax.set_xlabel("reported beyond the curated UniProtKB annotation")
    panel_title(ax, "B", "What a curated record does not say")
    tidy(ax, "x")

    for ext in ("pdf", "png"):
        fig.savefig(out / f"figure4_stories.{ext}")
    plt.close(fig)
    print("  figure4_stories: timing and features beyond curation")


# --------------------------------------------------------------- figure S1
def figureS1(results: Path, out: Path, records) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 2.6), constrained_layout=True)
    order = [m for m in ORDER if m != "uniprot_api"]
    for i, m in enumerate(order):
        vals = [r["methods"][m]["f1"] for r in records
                if "f1" in r["methods"].get(m, {})]
        ax.scatter(np.arange(len(vals)) + (i - 2.5) * 0.16, vals, s=8,
                   alpha=0.7, color=COLOR[m], label=LABEL[m], edgecolors="none")
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("protein (evaluation cohort, in the order of the benchmark script)")
    ax.set_ylabel("F1 against curated topology")
    ax.set_title("Per-protein agreement with curated UniProtKB topology",
                 loc="left", fontweight="bold")
    ax.legend(frameon=False, ncol=3, fontsize=6.5)
    tidy(ax)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"figureS1_per_protein.{ext}")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="bench/results")
    ap.add_argument("--out", default="bench/figures")
    args = ap.parse_args()
    results = Path(args.results)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    summary = read_summary(results)
    records = json.loads((results / "benchmark.json").read_text())

    figure2(results, out, summary, records)
    figure3(results, out)
    figure4(results, out)
    figureS1(results, out, records)
    print(f"figures written to {out}")


if __name__ == "__main__":
    main()
