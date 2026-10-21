"""
Figures for the ProteinLab manuscript, generated from the benchmark output.

    python bench/make_figures.py --results bench/results --out bench/figures

Figure 3 (main text)
  3A  agreement with curated topology per method (F1, recall, precision)
  3B  distribution of the worse boundary error per method
  3C  per-protein comparison of the two structure-based providers
  3D  transfer free energy: membrane cohort vs soluble controls

Figure S1 (supplementary)
  S1A per-protein F1 for every method
  S1B scatter of predicted vs curated helix counts
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

LABEL = {
    "3d_energy": "3D Energy",
    "3d_slab_geom": "3D slab geometry",
    "kyte_doolittle_seq": "Kyte-Doolittle",
    "tmhmm": "TMHMM 2.0",
    "toppred": "TopPred-style",
    "topcons_style_consensus": "TOPCONS-style",
    "uniprot_api": "UniProtKB (reference)",
}
COLOR = {
    "3d_energy": "#1f6f8b",
    "3d_slab_geom": "#e07a3f",
    "kyte_doolittle_seq": "#8a9a5b",
    "tmhmm": "#6b4c9a",
    "toppred": "#b5651d",
    "topcons_style_consensus": "#3f7d8c",
    "uniprot_api": "#bbbbbb",
}
ORDER = ["3d_energy", "3d_slab_geom", "tmhmm", "topcons_style_consensus",
         "kyte_doolittle_seq", "toppred", "uniprot_api"]
plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
    "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "figure.dpi": 200, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
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


def panel_a(summary, ax):
    rows = [r for r in summary if r["method"] in ORDER]
    rows.sort(key=lambda r: ORDER.index(r["method"]))
    x = np.arange(len(rows))
    w = 0.26
    ax.bar(x - w, [f(r, "macro_recall") for r in rows], w, label="recall",
           color="#4c9f70", edgecolor="white", linewidth=0.6)
    ax.bar(x, [f(r, "macro_f1") for r in rows], w, label="F1",
           color="#1f6f8b", edgecolor="white", linewidth=0.6)
    ax.bar(x + w, [f(r, "macro_precision") for r in rows], w, label="precision",
           color="#e07a3f", edgecolor="white", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([LABEL[r["method"]] for r in rows], rotation=20, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("agreement with curated topology")
    ax.set_title("A  Agreement with curated topology", loc="left", fontweight="bold")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.30))
    ax.grid(axis="y", lw=0.4, alpha=0.35)
    ax.set_axisbelow(True)


def panel_b(records, ax):
    present = {m for r in records for m in r.get("methods", {})}
    methods = [m for m in ORDER if m in present and m != "uniprot_api"]
    data = []
    for m in methods:
        vals = []
        for r in records:
            mm = r["methods"].get(m) or {}
            v = mm.get("median_max_boundary")
            if v is not None:
                vals.append(v)
        data.append(vals)
    bp = ax.boxplot(data, tick_labels=[LABEL[m] for m in methods], patch_artist=True,
                    widths=0.55, showfliers=False, medianprops=dict(color="black", lw=1.0))
    for patch, m in zip(bp["boxes"], methods):
        patch.set_facecolor(COLOR[m])
        patch.set_edgecolor("white")
        patch.set_linewidth(0.8)
    for key in ("whiskers", "caps"):
        for artist in bp[key]:
            artist.set_color("#555555")
            artist.set_linewidth(0.7)
    ax.axhline(5, color="#b03a3a", lw=0.9, ls="--")
    ax.text(0.02, 5.3, "5 residues", color="#b03a3a", fontsize=6.5, ha="left")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=28, ha="right")
    ax.set_ylabel("worse boundary error (residues)")
    ax.set_title("B  Boundary error", loc="left", fontweight="bold")
    ax.grid(axis="y", lw=0.4, alpha=0.35)
    ax.set_axisbelow(True)


def panel_c(records, ax):
    a, b = [], []
    for r in records:
        ma = r["methods"].get("3d_energy") or {}
        mb = r["methods"].get("3d_slab_geom") or {}
        if "f1" in ma and "f1" in mb:
            a.append(ma["f1"])
            b.append(mb["f1"])
    ax.scatter(a, b, s=14, c="#1f6f8b", alpha=0.75, edgecolors="none")
    lim = [0, 1.02]
    ax.plot(lim, lim, color="#888888", lw=0.8, ls="--")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("F1, 3D Energy")
    ax.set_ylabel("F1, 3D slab geometry")
    ax.set_title("C  Structure-based providers, per protein", loc="left", fontweight="bold")
    ax.grid(lw=0.4, alpha=0.35)
    ax.set_axisbelow(True)


def panel_d(records, results: Path, ax):
    dg = [r["methods"]["3d_energy"]["dG_transfer"] for r in records
          if "dG_transfer" in r["methods"].get("3d_energy", {})]
    sol_path = results / "soluble_control.json"
    sdg = []
    if sol_path.exists():
        sdg = [r["dG_transfer_kcal_per_mol"] for r in json.loads(sol_path.read_text())
               if r.get("dG_transfer_kcal_per_mol") is not None]
    rng = [min(dg + sdg) - 10, max(dg + sdg) + 10]
    ax.hist(dg, bins=40, range=rng, color="#1f6f8b", alpha=0.85,
            label=f"membrane cohort (n = {len(dg)})")
    if sdg:
        ax.hist(sdg, bins=40, range=rng, color="#e07a3f", alpha=0.85,
                label=f"soluble controls (n = {len(sdg)})")
    ax.axvline(-10, color="#b03a3a", lw=1.0, ls="--")
    ax.text(-10, ax.get_ylim()[1] * 0.92, " decision threshold\n -10 kcal/mol",
            color="#b03a3a", fontsize=6.5, va="top")
    ax.set_xlabel("transfer free energy of the best placement (kcal/mol)")
    ax.set_ylabel("proteins")
    ax.set_title("D  Membrane versus soluble discrimination", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.5)
    ax.grid(axis="y", lw=0.4, alpha=0.35)
    ax.set_axisbelow(True)


def figure4_frame(results: Path, out: Path) -> None:
    """Figure 4A / Supplementary Figure S2: the fitted frame against OPM/PPM."""
    path = results / "membrane_frame_benchmark.json"
    if not path.exists():
        print("no membrane_frame_benchmark.json - skipping figure 4A")
        return
    rows = json.loads(path.read_text())
    if not rows:
        return
    energy = [r["3d_energy"]["normal_error_deg"] for r in rows
              if "3d_energy" in r and "normal_error_deg" in r["3d_energy"]]
    slab = [r["3d_slab_geom"]["normal_error_deg"] for r in rows
            if "3d_slab_geom" in r and "normal_error_deg" in r["3d_slab_geom"]]
    offset = [r["3d_energy"]["midplane_offset_A"] for r in rows
              if "3d_energy" in r and "midplane_offset_A" in r["3d_energy"]]
    thickness = [r["3d_energy"]["thickness_A"] for r in rows
                 if "3d_energy" in r and "thickness_A" in r["3d_energy"]]
    reference = [r["reference_thickness_A"] for r in rows if "reference_thickness_A" in r]

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7), constrained_layout=True)

    ax = axes[0]
    ax.hist(energy, bins=24, color=COLOR["3d_energy"], alpha=0.9, label="3D Energy")
    if slab:
        ax.hist(slab, bins=24, color=COLOR["3d_slab_geom"], alpha=0.65, label="3D slab geometry")
    ax.axvline(5, color="#b03a3a", lw=0.9, ls="--")
    ax.text(5.4, ax.get_ylim()[1] * 0.55, "5 degrees", color="#b03a3a", fontsize=6.5)
    ax.set_xlabel("normal error (degrees)")
    ax.set_ylabel("structures")
    ax.set_title("A  Membrane normal", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.5)
    ax.grid(axis="y", lw=0.4, alpha=0.35)
    ax.set_axisbelow(True)

    ax = axes[1]
    ax.hist(offset, bins=24, color=COLOR["3d_energy"], alpha=0.9)
    ax.axvline(1.0, color="#b03a3a", lw=0.9, ls="--")
    ax.text(1.6, ax.get_ylim()[1] * 0.55, "1 A", color="#b03a3a", fontsize=6.5)
    ax.set_xlabel("mid-plane offset (A)")
    ax.set_ylabel("structures")
    ax.set_title("B  Mid-plane", loc="left", fontweight="bold")
    ax.grid(axis="y", lw=0.4, alpha=0.35)
    ax.set_axisbelow(True)

    ax = axes[2]
    ax.hist(thickness, bins=24, color=COLOR["3d_energy"], alpha=0.9, label="fitted")
    if reference:
        ax.hist(reference, bins=24, color="#bbbbbb", alpha=0.8, label="OPM/PPM reference")
    ax.set_xlabel("hydrophobic thickness (A)")
    ax.set_ylabel("structures")
    ax.set_title("C  Bilayer thickness", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.5)
    ax.grid(axis="y", lw=0.4, alpha=0.35)
    ax.set_axisbelow(True)

    fig.savefig(out / "figure4_frame.pdf")
    fig.savefig(out / "figure4_frame.png")
    plt.close(fig)
    print(f"  figure4_frame: n = {len(energy)}")


def figure4_unrolling(results: Path, out: Path) -> None:
    """Figure 4B / Supplementary Figure S3: the geometric audit of the layout."""
    path = results / "unroll_benchmark.json"
    if not path.exists():
        print("no unroll_benchmark.json - skipping figure 4B")
        return
    rows = json.loads(path.read_text())
    if not rows:
        return
    rmsd = [v for r in rows for v in r.get("closure_rmsd", [])]
    schematic = [r["n_loops_schematic"] for r in rows]
    loops = [r["n_loops"] for r in rows]
    clashes = [r["clashes"] for r in rows]

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7), constrained_layout=True)

    ax = axes[0]
    if rmsd:
        ax.hist(rmsd, bins=30, color="#1f6f8b", alpha=0.9)
        ax.axvline(0.3, color="#b03a3a", lw=0.9, ls="--")
        ax.text(0.4, ax.get_ylim()[1] * 0.92, "0.3 A", color="#b03a3a", fontsize=6.5)
        ax.set_yscale("log")
    ax.set_xlabel("closure error (A)")
    ax.set_ylabel("loops")
    ax.set_title("A  Loop closure error", loc="left", fontweight="bold")
    ax.grid(axis="y", lw=0.4, alpha=0.35)
    ax.set_axisbelow(True)

    ax = axes[1]
    fraction = [s / l if l else 0.0 for s, l in zip(schematic, loops)]
    ax.hist(fraction, bins=20, color="#e07a3f", alpha=0.9)
    ax.set_xlabel("fraction of loops drawn schematically")
    ax.set_ylabel("structures")
    ax.set_title("B  Loops not closed", loc="left", fontweight="bold")
    ax.grid(axis="y", lw=0.4, alpha=0.35)
    ax.set_axisbelow(True)

    ax = axes[2]
    ax.hist(clashes, bins=20, color="#6b4c9a", alpha=0.9)
    ax.set_xlabel("new heavy-atom clashes introduced")
    ax.set_ylabel("structures")
    ax.set_title("C  New steric clashes", loc="left", fontweight="bold")
    ax.grid(axis="y", lw=0.4, alpha=0.35)
    ax.set_axisbelow(True)

    fig.savefig(out / "figure4_unrolling.pdf")
    fig.savefig(out / "figure4_unrolling.png")
    plt.close(fig)
    print(f"  figure4_unrolling: {len(rows)} structures")


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

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8), constrained_layout=True)
    panel_a(summary, axes[0][0])
    panel_b(records, axes[0][1])
    panel_c(records, axes[1][0])
    panel_d(records, results, axes[1][1])
    fig.savefig(out / "figure3_benchmark.png")
    fig.savefig(out / "figure3_benchmark.pdf")
    plt.close(fig)

    # supplementary per-protein panel
    fig, ax = plt.subplots(figsize=(7.2, 2.6), constrained_layout=True)
    order = [m for m in ORDER if m != "uniprot_api"]
    for i, m in enumerate(order):
        vals = [r["methods"][m]["f1"] for r in records
                if "f1" in r["methods"].get(m, {})]
        ax.scatter(np.arange(len(vals)) + (i - 1.5) * 0.18, vals, s=9,
                   alpha=0.65, color=COLOR[m], label=LABEL[m], edgecolors="none")
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("protein (evaluation cohort, in the order of the benchmark script)")
    ax.set_ylabel("F1 against curated topology")
    ax.set_title("Per-protein agreement with curated UniProtKB topology",
                 loc="left", fontweight="bold")
    ax.legend(frameon=False, ncol=4, fontsize=6.5)
    ax.grid(axis="y", lw=0.4, alpha=0.35)
    ax.set_axisbelow(True)
    fig.savefig(out / "figureS1_per_protein.png")
    fig.savefig(out / "figureS1_per_protein.pdf")
    plt.close(fig)

    figure4_frame(results, out)
    figure4_unrolling(results, out)
    print(f"figures written to {out}")


if __name__ == "__main__":
    main()
