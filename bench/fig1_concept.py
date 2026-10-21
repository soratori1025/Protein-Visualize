"""
Figure 1 - the concept, drawn from the model itself.

    python bench/fig1_concept.py --out bench/figures

This figure is not a schematic drawn by hand: every curve in panel A is the depth
profile that the energy model actually uses (the two logistic occupancy functions of
the membrane core and of the head-group region), panel B is the transfer free energy
the model assigns to each residue type in each of those two regions, and panel D is
the per-residue depth of a real protein, computed by the same code that annotates it.

Why a separate script: figure 1 has to be legible to a reader who has never seen a
membrane-protein figure, and it is the only place where the manuscript shows the
model rather than measures it.  It is regenerated like every other figure, so it
cannot drift away from the implementation.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

INK = "#22252a"
MUTED = "#6b7280"
GRID = "#d8dde3"
ACCENT = "#0f5c7a"
SECOND = "#d9702e"
THIRD = "#5f7d4f"
WARN = "#a83232"
LIPID = "#e9dccb"
HEAD = "#c9b79f"
WATER = "#dfeaf0"

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


def tidy(ax, axis="both"):
    ax.grid(axis=axis, lw=0.5, color=GRID)
    ax.set_axisbelow(True)
    ax.tick_params(length=2.5)


def title(ax, letter, text, above=False):
    """The panel letter and title.

    `above=False` sets the title on the panel itself, as for the two panels that have room
    for it.  `above=True` lifts it above the axes altogether, for the panels whose neighbour
    above is a row that ends in an axis label: a title printed inside those panels would
    collide with the tick labels and the axis title of the row above.
    """
    if above:
        ax.set_title(f"{letter}   {text}", loc="left", fontweight="bold", color=INK,
                     pad=22)
    else:
        ax.set_title(f"{letter}   {text}", loc="left", fontweight="bold", color=INK, pad=6)


# --------------------------------------------------------------------------- A
def panel_a(ax):
    """The depth profile the model uses: the core and the head-group region."""
    z = np.linspace(-34, 34, 900)
    D, W, L = 13.0, 7.0, 1.5
    sig = lambda x: 1.0 / (1.0 + np.exp(-x))
    core = sig((D - np.abs(z)) / L)
    interface = sig((D + W - np.abs(z)) / L) - core
    outside = 1.0 - sig((D + W - np.abs(z)) / L)

    ax.fill_between(z, 0, outside, color=WATER, lw=0, zorder=0.6)
    ax.fill_between(z, 0, interface + core, color=HEAD, lw=0)
    ax.fill_between(z, 0, core, color=LIPID, lw=0)
    ax.plot(z, core, color=ACCENT, lw=1.6, label="hydrocarbon core")
    ax.plot(z, core + interface, color=SECOND, lw=1.2, ls="--", label="core + head groups")

    ax.text(0, 0.52, "hydrocarbon\ncore", ha="center", fontsize=7, color=INK)
    for x in (-D - W / 2, D + W / 2):
        ax.text(x, 0.62, "head\ngroups", ha="center", fontsize=6.5, color=MUTED)
    ax.text(0, 1.13, "aqueous", ha="center", fontsize=6.5, color=MUTED)
    ax.text(0, -0.16, "aqueous", ha="center", fontsize=6.5, color=MUTED)

    # the leader of the arrow: an annotation that ends behind the axis is not clipped to
    # the axes by the arrow patch, so the arrow is drawn as a plain segment
    ax.plot([0, 0], [0, -D], color=INK, lw=0.8, zorder=3)
    ax.text(1.2, -D / 2, "$D$", fontsize=8, va="center")
    ax.annotate("", xy=(D, 1.02), xytext=(D + W, 1.02),
                arrowprops=dict(arrowstyle="<->", lw=0.8, color=MUTED))
    ax.text(D + W / 2, 1.05, "$W$", fontsize=7.5, ha="center", color=MUTED)
    # The module-level rcParam `savefig.bbox = "tight"` puts every text of the figure into
    # the bounding box that is computed when the file is written, and a text placed outside
    # the plotted range — the label of the bilayer half-thickness, which sits below the axis
    # — then stretches that box down to nothing.  The panel sets its own limit, so the text
    # is told not to take part in the layout.
    for artist in ax.texts:
        artist.set_in_layout(False)

    ax.set_xlim(-34, 34)
    ax.set_ylim(-0.22, 1.26)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(["0", "0.5", "1"])
    ax.set_xlabel("distance from the membrane mid-plane, $z$ (\u00c5)")
    ax.set_ylabel("occupancy")
    title(ax, "A", "A membrane is a depth, not a switch", above=True)
    ax.legend(frameon=False, fontsize=6.0, loc="center right", ncol=1,
              bbox_to_anchor=(1.02, 0.30), handlelength=1.4, borderpad=0.2)
    tidy(ax)


# --------------------------------------------------------------------------- B
def panel_b(ax):
    """The transfer free energy the model assigns to each residue type."""
    # the two scales the model uses, as tabulated in the implementation
    core = {"A": 1.6, "V": 1.7, "L": 1.9, "I": 1.8, "M": 1.3, "F": 1.9, "W": 1.4,
            "C": 1.2, "G": 0.9, "T": 0.7, "S": 0.6, "P": 0.4, "Y": 0.7, "H": 0.1,
            "Q": -0.1, "N": -0.2, "E": -1.1, "D": -1.2, "K": -1.5, "R": -1.6}
    interface = {"W": 2.1, "F": 1.6, "L": 1.1, "M": 1.0, "I": 1.0, "V": 0.9, "Y": 0.9,
                 "C": 0.8, "A": 0.5, "T": 0.4, "P": 0.4, "G": 0.3, "S": 0.2, "H": 0.0,
                 "Q": -0.2, "N": -0.4, "E": -0.6, "D": -0.7, "K": -1.0, "R": -1.2}
    order = sorted(core, key=lambda r: -core[r])
    x = np.arange(len(order))
    w = 0.4
    ax.bar(x - w / 2, [core[r] for r in order], w, color=ACCENT, label="in the core",
           edgecolor="white", lw=0.4)
    ax.bar(x + w / 2, [interface[r] for r in order], w, color=SECOND, label="at the interface",
           edgecolor="white", lw=0.4)
    ax.axhline(0, color=MUTED, lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(order, fontsize=6.5)
    ax.set_ylabel("transfer free energy (kcal/mol)")
    ax.set_xlabel("residue type, ordered by the core scale")
    ax.set_ylim(-1.85, 2.45)
    title(ax, "B", "Two regions, two scales")
    ax.legend(frameon=False, fontsize=6.3, loc="lower left", bbox_to_anchor=(0.0, 0.02))
    ax.annotate("tryptophan: disfavoured in the core,\nstrongly favoured at the interface",
                xy=(order.index("W") + w / 2, interface["W"] + 0.03), xytext=(8.2, 2.1),
                fontsize=6.3, color=INK,
                arrowprops=dict(arrowstyle="->", lw=0.7, color=MUTED,
                                connectionstyle="arc3,rad=-0.2"))
    tidy(ax)


# --------------------------------------------------------------------------- C
def panel_c(ax, results: Path):
    """The depth of every residue of a real protein, as the tool computes it.

    The curve is not a diagram of a transmembrane helix: it is the signed distance of each
    C-alpha atom from the fitted mid-plane, so the helices rise out of the bilayer, the
    loops fall back into the water on either side, and the smooth transitions are the
    occupancy functions of panel A applied to a real structure.
    """
    import json
    depth_path = results / "figure1_depth_profile.json"
    if not depth_path.exists():
        ax.text(0.5, 0.5, "figure1_depth_profile.json not found", ha="center",
                transform=ax.transAxes)
        return
    prof = json.loads(depth_path.read_text(encoding="utf-8"))
    depth = np.asarray(prof["depth"], dtype=float)
    labels = prof["labels"]
    half = float(prof["half_thickness"])
    n = len(depth)
    x = np.arange(n)

    ax.axhspan(-half, half, color=LIPID, lw=0, zorder=0)
    ax.axhline(0.0, color=HEAD, lw=0.7, ls=(0, (3, 2)), zorder=1)
    for y in (-half, half):
        ax.axhline(y, color=HEAD, lw=0.8, zorder=1)
    ax.plot(x, depth, color=ACCENT, lw=0.75, zorder=2)
    core = np.array([lab == "T" for lab in labels])
    ax.fill_between(x, -half, half, where=core, color=ACCENT, alpha=0.18, lw=0,
                    step="mid", zorder=1.5)
    ax.set_xlim(0, n - 1)
    ax.set_ylim(-34, 34)
    ax.set_yticks([-half, 0, half])
    ax.set_yticklabels(["\u2212D", "0", "+D"], fontsize=7)
    ax.set_xlabel(f"residue, along the chain of {prof.get('protein', '')} (chain "
                  f"{prof.get('chain', '')}, PDB {prof.get('pdb_id', '')})")
    ax.set_ylabel("depth (Å)")
    title(ax, "C", "One structure, read residue by residue", above=True)
    ax.annotate("hydrophobic core", xy=(0.5 * n, half * 0.82), ha="center", va="center",
                fontsize=6.2, color=MUTED, style="italic")
    ax.annotate("water", xy=(0.13 * n, -half - 8), ha="center", va="center",
                fontsize=6.2, color=MUTED, style="italic")
    tidy(ax, axis="y")
    ax.annotate(f"{prof.get('n_computed_tm')} transmembrane segments, identical to the curated\n"
                "annotation; the shaded bands are the segments themselves",
                xy=(0.5 * n, 30), ha="center", va="center", fontsize=6.2, color=MUTED,
                style="italic")


# --------------------------------------------------------------------------- D
def panel_d(ax):
    """The four layers of the answer, as a reader receives them."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    boxes = [
        ("1", "Place the bilayer", "one free energy over the normal,\nthe mid-plane and the thickness", ACCENT),
        ("2", "Annotate every residue", "a topology label, a depth,\na zone and a confidence", SECOND),
        ("3", "Resolve the hard cases", "hairpins, broken helices\nand interfacial helices", THIRD),
        ("4", "Say what it cannot resolve", "a transfer free energy and\nan uncertainty flag", WARN),
    ]
    for i, (num, head, body, colour) in enumerate(boxes):
        y = 8.6 - i * 2.28
        ax.add_patch(FancyBboxPatch((0.35, y - 1.55), 9.2, 1.62,
                                    boxstyle="round,pad=0.05,rounding_size=0.16",
                                    fc="white", ec=GRID, lw=0.9))
        ax.add_patch(plt.Circle((0.95, y - 0.74), 0.36, color=colour, lw=0, zorder=3))
        ax.text(0.95, y - 0.74, num, ha="center", va="center", color="white",
                fontsize=8, fontweight="bold", zorder=4)
        ax.text(1.55, y - 0.5, head, va="center", fontsize=8, fontweight="bold", color=INK)
        ax.text(1.55, y - 1.18, body, va="center", fontsize=6.8, color=MUTED)
        if i < len(boxes) - 1:
            ax.add_patch(FancyArrowPatch((0.95, y - 1.16), (0.95, y - 1.62),
                                         arrowstyle="-|>", mutation_scale=8,
                                         lw=0.9, color=MUTED))
    title(ax, "D", "What the software returns, in order")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="bench/results")
    ap.add_argument("--out", default="bench/figures")
    args = ap.parse_args()
    results = Path(args.results)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # The four panels are laid out with explicit margins rather than with `tight_layout`:
    # every axis title and tick label is drawn at a fixed distance from its axes, so the
    # space that a row needs does not depend on how long a title happens to be.
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.6))
    fig.subplots_adjust(left=0.085, right=0.985, top=0.975, bottom=0.105,
                        wspace=0.30, hspace=0.46)
    panel_a(axes[0][0])
    panel_b(axes[0][1])
    panel_c(axes[1][0], results)
    panel_d(axes[1][1])

    for ext in ("pdf", "png"):
        fig.savefig(out / f"figure1_concept.{ext}", bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("  figure1_concept: written")


if __name__ == "__main__":
    main()
