"""Case-study numbers for the manuscript: Piezo1 (Q92508), the KcsA control, and hSERT.

    python bench/case_piezo1.py

Writes ``bench/results/piezo1_case.json``, which the curvature subsection of the
manuscript quotes: for each structure, the topology count, the energy of the best
planar placement, the energy of the best spherical placement and the fitted radius.
The structures are read from ``backend/data/uploads/``; they are downloaded from
the RCSB PDB and are not redistributed with the source archive.

The plane and sphere are computed here exactly as ``membrane_energy.fit_structure``
computes them, so the numbers do not depend on which of the two the acceptance rule
happens to choose.  The acceptance rule itself (``CURVED_MIN_GAIN`` and
``CURVED_MIN_GAIN_REL``) is reported with every structure so that a reader can see
whether the sphere would have been accepted.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

import numpy as np  # noqa: E402

from app.services.topology import membrane_energy as ME  # noqa: E402
from app.services.topology.service import predict_topology  # noqa: E402

DATA = ROOT / "backend" / "data" / "uploads"
OUT = ROOT / "bench" / "results" / "piezo1_case.json"

CASES = (("8YEZ_Q92508.cif", "A", "PIEZO1, 8YEZ"),
         ("8ZU3.cif", "A", "PIEZO1, 8ZU3"),
         ("hSERT_5I6X_ions.pdb", "A", "hSERT, 5I6X"))


def plane_and_sphere(path: Path) -> dict:
    """Both placements on the same weights, without the acceptance rule."""
    residues, predicted = ME.read_structure(path)
    warnings: list[str] = []
    residues = ME.confident_residues(residues, predicted, warnings)
    rsa, centroids, _, _ = ME.residue_exposure(residues)
    atoms = np.array([p for r in residues for p in r.xyz])
    pore = ME.pore_weight(ME.escape_fraction(centroids, atoms))
    names = [r.resname for r in residues]
    weights = rsa * pore
    planar = ME.place_membrane(centroids, weights, names)
    curved = ME.place_membrane_curved(centroids, weights, names, planar)
    gain = None
    if curved is not None:
        gain = planar.energy - curved.energy
    accepted = bool(curved is not None and gain is not None
                    and gain >= max(ME.CURVED_MIN_GAIN,
                                    ME.CURVED_MIN_GAIN_REL * abs(planar.energy)))
    return {
        "n_residues": len(residues),
        "n_chains": len({r.chain for r in residues}),
        "plane_energy_kcal_per_mol": round(float(planar.energy), 2),
        "sphere_energy_kcal_per_mol": round(float(curved.energy), 2) if curved else None,
        "sphere_radius_A": round(float(curved.radius), 1) if curved else None,
        "gain_kcal_per_mol": round(float(gain), 2) if gain is not None else None,
        "accepted": accepted,
    }


def main() -> None:
    rows = []
    for fname, chain, label in CASES:
        path = DATA / fname
        if not path.exists():
            print(f"[case] {fname}: missing (download it from the RCSB PDB)")
            continue
        result = predict_topology(path, tm_algo="3d_energy", ss_algo="none",
                                  flow_type="consensus", chain_id=chain)
        tms = [reg for reg in result.regions if reg.type == "Transmembrane"]
        row = {"file": fname, "label": label, "chain": result.chain_id,
               "n_tm": len(tms), **plane_and_sphere(path)}
        rows.append(row)
        print(f"{fname:22s} chain={result.chain_id} nTM={len(tms)} "
              f"plane={row['plane_energy_kcal_per_mol']:.1f} "
              f"sphere={row['sphere_energy_kcal_per_mol']} "
              f"R={row['sphere_radius_A']} accepted={row['accepted']}")
        print(f"   {result.labeler}")

    OUT.write_text(json.dumps({"structures": rows,
                               "acceptance_rule": {
                                   "min_gain_kcal_per_mol": ME.CURVED_MIN_GAIN,
                                   "min_gain_relative": ME.CURVED_MIN_GAIN_REL}},
                              indent=2))
    print(f"written to {OUT}")


if __name__ == "__main__":
    main()
