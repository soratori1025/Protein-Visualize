"""
The depth profile of one structure, for panel C of Figure 1.

    python bench/fig1_depth_profile.py --out bench/results/figure1_depth_profile.json

Panel C of Figure 1 is the one place in the article where a real structure is shown
residue by residue, and until now it showed a strip of labels rather than the model: the
per-residue results file stores how many residues carry each label, not which residue
carries which.  This script runs the production pipeline - the same `predict_topology`
that the web application calls - on the structure that stories A and B use, and writes the
signed depth of every residue that the model computed, together with its label.

The depth is the distance of each C-alpha atom from the fitted mid-plane, measured along
the fitted normal: it is what `MembraneFrame.depth` holds inside the service, and the
script checks the reconstruction against the depths that the service reports before
writing anything.  The check is a hard failure, not a warning: the figure must not claim a
profile that the tool does not compute.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from app.services.topology.residues import load_residue_frame  # noqa: E402
from app.services.topology.service import predict_topology  # noqa: E402

PDB = "7YK7"
CHAIN = "R"
STRUCTS = ROOT / "bench" / "results" / "structures"

LABEL = {"Transmembrane": "T", "Extracellular": "O", "Cytoplasmic": "I",
         "Interfacial": "F", "Extramembrane": "E"}


def structure() -> Path:
    path = STRUCTS / f"{PDB}.cif"
    if path.exists():
        return path
    STRUCTS.mkdir(parents=True, exist_ok=True)
    import urllib.request
    urllib.request.urlretrieve(f"https://files.rcsb.org/download/{PDB}.cif", path)
    return path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="bench/results/figure1_depth_profile.json")
    args = ap.parse_args()

    path = structure()
    response = predict_topology(path, tm_algo="3d_energy", ss_algo="none", chain_id=CHAIN,
                                include_residues=True)
    reported = np.array([r.depth for r in response.residues], dtype=float)
    normal = np.array(response.membrane["normal"], dtype=float)
    half = float(response.membrane["half_thickness"])

    frame = load_residue_frame(path, chain_id=CHAIN)
    if len(frame.coords) != len(reported):
        raise SystemExit(f"{len(frame.coords)} residues in the frame, {len(reported)} reported")

    # the mid-plane is the translation along the normal; it is recovered from the data
    # because the response reports the normal and the thickness, not the offset
    projection = frame.coords @ normal
    offsets = np.arange(-60.0, 60.0, 0.01)
    errors = [np.median(np.abs((projection - (projection.mean() + o)) - reported))
              for o in offsets]
    mid = projection.mean() + offsets[int(np.argmin(errors))]
    depth = projection - mid
    if np.max(np.abs(depth - reported)) > 0.05:
        raise SystemExit("the reconstructed depth does not match the service: "
                         f"max |error| = {np.max(np.abs(depth - reported)):.3f} A")

    labels = "".join(LABEL.get(r.label, "?") for r in response.residues)
    runs, start = [], 0
    for i in range(1, len(labels) + 1):
        if i == len(labels) or labels[i] != labels[start]:
            runs.append((labels[start], i - start))
            start = i
    # a transmembrane segment is a run that is long enough to cross the bilayer; the short
    # runs are the flanking residues that the boundary placement puts inside the membrane
    n_computed = sum(1 for lab, length in runs if lab == "T" and length >= 15)
    curated = next((r for r in json.loads(
        (ROOT / "bench" / "results" / "user_story_benchmark.json").read_text(
            encoding="utf-8"))["story_a"] if r.get("pdb_id") == PDB), {})
    out = {
        "pdb_id": PDB,
        "chain": CHAIN,
        "protein": curated.get("protein", "Relaxin-3 receptor 2"),
        "accession": curated.get("accession"),
        "normal": [round(float(x), 4) for x in normal],
        "midplane_point": [round(float(x), 3) for x in normal * mid],
        "half_thickness": round(half, 2),
        "n_residues": len(depth),
        "n_computed_tm": n_computed,
        "n_curated_tm": curated.get("n_curated_tm"),
        "f1_vs_curation": curated.get("f1_vs_curation"),
        "depth": [round(float(x), 2) for x in depth],
        "resseq": [r.resseq for r in frame.residues],
        "labels": labels,
        "provenance": ("bench/fig1_depth_profile.py: predict_topology(tm_algo='3d_energy') "
                       "on the downloaded structure, depth = CA @ normal - mid-plane"),
    }
    out_path = ROOT / args.out
    out_path.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"  {out_path}: {len(depth)} residues, {n_computed} transmembrane segments, "
          f"max |reconstruction error| = {np.max(np.abs(depth - reported)):.3f} A")


if __name__ == "__main__":
    main()
