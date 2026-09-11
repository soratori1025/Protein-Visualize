from __future__ import annotations

from collections import Counter
from math import sqrt
from typing import Any

from Bio.SeqUtils.ProtParam import ProteinAnalysis

HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def analyze_chain(chain: dict[str, Any], contact_cutoff: float = 8.0) -> dict[str, Any]:
    sequence = chain["sequence"]
    analysis = ProteinAnalysis(sequence) if sequence else None
    ca_atoms = []
    for residue in chain["residues"]:
        atom = next((item for item in residue["atoms"] if item["name"] == "CA"), None)
        if atom:
            ca_atoms.append((residue["id"], atom))

    contacts = []
    for index, (first_id, first_atom) in enumerate(ca_atoms):
        for second_id, second_atom in ca_atoms[index + 1:]:
            distance = _distance(first_atom, second_atom)
            if distance <= contact_cutoff:
                contacts.append({"source": first_id, "target": second_id, "distance": round(distance, 3)})

    composition = Counter(sequence)
    return {
        "chain_id": chain["id"],
        "algorithm": "protein_engine.analysis.metrics",
        "version": "0.1.0",
        "parameters": {"contact_cutoff": contact_cutoff},
        "sequence": {"length": len(sequence), "composition": dict(sorted(composition.items())), "unknown_residues": composition.get("X", 0)},
        "physicochemical": _physicochemical(analysis, sequence),
        "geometry": {"ca_atoms": len(ca_atoms), "contacts": len(contacts), "contact_density": round(len(contacts) / max(len(ca_atoms), 1), 4)},
        "contacts": contacts[:5000],
    }


def _physicochemical(analysis: ProteinAnalysis | None, sequence: str) -> dict[str, float | None]:
    if not analysis:
        return {"molecular_weight": None, "isoelectric_point": None, "aromaticity": None, "instability_index": None, "mean_hydropathy": None}
    mean_hydropathy = sum(HYDROPATHY.get(letter, 0.0) for letter in sequence) / max(len(sequence), 1)
    return {"molecular_weight": round(analysis.molecular_weight(), 3), "isoelectric_point": round(analysis.isoelectric_point(), 3), "aromaticity": round(analysis.aromaticity(), 4), "instability_index": round(analysis.instability_index(), 3), "mean_hydropathy": round(mean_hydropathy, 3)}


def _distance(first: dict[str, float], second: dict[str, float]) -> float:
    return sqrt(sum((first[axis] - second[axis]) ** 2 for axis in ("x", "y", "z")))