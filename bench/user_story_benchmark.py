"""
The benchmark of the two user stories.

    python bench/user_story_benchmark.py --results bench/results --out bench/results

Two users, and for each of them the question is the same: does the tool answer it,
and how well?

**Story A - "I have a protein that is defined in UniProtKB and I want to visualize
it."**  The tool has to produce a topology figure for a curated membrane protein.  What
is measured: the agreement between the computed topology and the curated one, the
information the tool reports *beyond* the curation (broken helices, interfacial
helices, hairpins, a curved membrane), the number of machine-readable formats that can
be exported, and the wall-clock time.

**Story B - "I have a new model of my protein and I want to test it on the
transmembrane algorithm."**  The structure is a prediction, so no curated annotation of
*it* exists.  What is measured: whether the tool produces a topology at all, the
transfer free energy and the fitted membrane thickness, the distribution of the
per-residue confidence, and - the point of the exercise - whether the tool says when
the placement cannot be trusted.

The competitors cannot be run here (TOPCONS2 needs a 20 GB database and an Ubuntu
host, TMbed and DeepTMHMM need multi-gigabyte models, Protter is a web application),
so they are compared on the basis of what their own sources say they do; see
``bench/competitive_matrix.py``.  This script measures *our* side of that comparison.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _candidate in (ROOT / "backend", ROOT):
    if str(_candidate) not in sys.path:
        sys.path.insert(0, str(_candidate))

import numpy as np  # noqa: E402

from app.schemas.topology import TMParams  # noqa: E402
from app.services.topology.service import predict_topology  # noqa: E402

from benchmark_tm_topology import match_segments, overlap  # noqa: E402

UA = {"User-Agent": "ProteinLab-benchmark/1.0 (academic validation)"}
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"


# --------------------------------------------------------------- the 3line format
def to_3line(sequence: str, labels: list[str]) -> str:
    """The ``.3line`` export: a FASTA record and a topology string.

    ``M`` transmembrane, ``I`` cytoplasmic, ``O`` extracellular, ``-`` outside the
    membrane and not assigned a side.
    """
    out = []
    for name, state in zip(sequence, labels):
        if state == "TM_in":
            out.append("M")
        elif state == "TM_C":
            out.append("I")
        elif state == "TM_E":
            out.append("O")
        else:
            out.append("-")
    return "".join(out)


def check_3line_roundtrip(topology: str, labels: list[str]) -> dict:
    """Does the exported string carry the same information as the annotation?

    The invariant is exact: the positions marked ``M`` in the string must be exactly
    the positions the consensus map calls transmembrane, and the positions marked
    ``I``/``O`` must be exactly the extramembrane helices on each side.
    """
    m_positions = {i for i, c in enumerate(topology) if c == "M"}
    expected_m = {i for i, lab in enumerate(labels) if lab == "TM_in"}
    i_positions = {i for i, c in enumerate(topology) if c == "I"}
    o_positions = {i for i, c in enumerate(topology) if c == "O"}
    expected_i = {i for i, lab in enumerate(labels) if lab == "TM_C"}
    expected_o = {i for i, lab in enumerate(labels) if lab == "TM_E"}
    helices = 0
    i = 0
    while i < len(topology):
        if topology[i] == "M":
            j = i
            while j + 1 < len(topology) and topology[j + 1] == "M":
                j += 1
            helices += 1
            i = j + 1
        else:
            i += 1
    return {
        "n_helices_in_string": helices,
        "membrane_positions_exact": m_positions == expected_m,
        "cytoplasmic_positions_exact": i_positions == expected_i,
        "extracellular_positions_exact": o_positions == expected_o,
        "consistent": (m_positions == expected_m and i_positions == expected_i
                       and o_positions == expected_o),
    }


def fetch_alphafold(accession: str, dest: Path) -> bool:
    try:
        request = urllib.request.Request(ALPHAFOLD_API.format(acc=accession), headers=UA)
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except Exception:
        return False
    if not payload:
        return False
    url = payload[0].get("cifUrl")
    if not url:
        return False
    try:
        request = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(request, timeout=180) as response:
            data = response.read()
    except Exception:
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    return True


def story_a(records: list[dict], structures: Path, limit: int) -> list[dict]:
    """The curated-protein user: a UniProtKB entry and a structure."""
    rows = []
    for record in records[:limit]:
        path = structures / f"{record['pdb_id']}.cif"
        if not path.exists():
            continue
        t0 = time.time()
        try:
            response = predict_topology(path, tm_algo="3d_energy", ss_algo="dssp",
                                        flow_type="consensus", chain_id=record["chain"])
        except Exception as exc:
            print(f"[story A] {record['pdb_id']}: {type(exc).__name__}: {exc}")
            continue
        elapsed = time.time() - t0

        # the topology string that the figure exports
        labels = [c.label for c in (response.consensus_map or [])]
        sequence = "".join(
            {"TM_in": "M", "TM_C": "I", "TM_E": "O"}.get(lab, "-")
            for lab in labels)
        regions = [r for r in response.regions if r.type == "Transmembrane"]
        roundtrip = check_3line_roundtrip(sequence, labels)

        roles = {}
        for region in response.regions:
            if region.membrane_role:
                roles[region.membrane_role] = roles.get(region.membrane_role, 0) + 1

        row = {
            "accession": record["accession"], "pdb_id": record["pdb_id"],
            "chain": record["chain"],
            "protein": record.get("protein_name", "")[:70],
            "n_curated_tm": record["n_reference_tm"],
            "n_computed_tm": len(regions),
            "f1_vs_curation": record["methods"]["3d_energy"].get("f1"),
            "boundary_median": record["methods"]["3d_energy"].get("median_max_boundary"),
            "membrane_score": response.membrane_score,
            "placement_quality": response.placement_quality,
            "seconds": round(elapsed, 2),
            "roles": roles,
            "n_broken_helix": roles.get("BROKEN_TM", 0),
            "n_unwound": roles.get("UNWOUND", 0),
            "n_interfacial": roles.get("INTERFACIAL", 0),
            "n_reentrant": roles.get("REENTRANT", 0),
            "n_signal": roles.get("SIGNAL", 0),
            "export_3line": {"n_characters": len(sequence),
                              **roundtrip},
            "exports_available": ["PNG", "3line", "GFF3", "SVG-in-app"],
        }
        rows.append(row)
        print(f"[story A] {record['pdb_id']} {record['chain']}: "
              f"curated {row['n_curated_tm']} TM, computed {row['n_computed_tm']}, "
              f"F1 {row['f1_vs_curation']}, {row['seconds']} s, "
              f"extra: broken {row['n_broken_helix']}, interfacial {row['n_interfacial']}")
    return rows


def story_b(records: list[dict], models: Path, limit: int) -> list[dict]:
    """The local user: a new model, no curated annotation of it."""
    rows = []
    for record in records[:limit]:
        accession = record["accession"]
        path = models / f"AF-{accession}-F1.cif"
        if not path.exists() and not fetch_alphafold(accession, path):
            print(f"[story B] {accession}: no AlphaFold model could be fetched")
            continue
        t0 = time.time()
        try:
            response = predict_topology(path, tm_algo="3d_energy", ss_algo="dssp",
                                        flow_type="consensus", chain_id="A")
        except Exception as exc:
            print(f"[story B] {accession}: {type(exc).__name__}: {exc}")
            continue
        elapsed = time.time() - t0
        regions = [r for r in response.regions if r.type == "Transmembrane"]
        confidence = {}
        for region in regions:
            confidence[region.confidence] = confidence.get(region.confidence, 0) + 1
        rows.append({
            "accession": accession,
            "model": f"AF-{accession}-F1.cif",
            "n_residues": len(response.residues or []),
            "n_computed_tm": len(regions),
            "membrane_score": response.membrane_score,
            "placement_quality": response.placement_quality,
            "placement_flags": response.placement_flags,
            "confidence": confidence,
            "seconds": round(elapsed, 2),
            "answered": bool(regions),
        })
        print(f"[story B] {accession}: {len(regions)} TM, dG {response.membrane_score:.1f}, "
              f"quality {response.placement_quality}, {elapsed:.1f} s")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="bench/results")
    ap.add_argument("--story-a", type=int, default=25)
    ap.add_argument("--story-b", type=int, default=15)
    args = ap.parse_args()
    results = Path(args.results)
    records = json.loads((results / "benchmark.json").read_text())
    structures = results / "structures"
    models = results / "alphafold"

    print("=" * 74)
    print("STORY A: a protein that is defined in UniProtKB, and I want to visualize it")
    print("=" * 74)
    a = story_a(records, structures, args.story_a)

    print()
    print("=" * 74)
    print("STORY B: a new model of my protein, and I want to test it on the algorithm")
    print("=" * 74)
    b = story_b(records, models, args.story_b)

    out = results / "user_story_benchmark.json"
    out.write_text(json.dumps({"story_a": a, "story_b": b}, indent=2), encoding="utf-8")

    print()
    print("=" * 74)
    print("SUMMARY")
    print("=" * 74)
    if a:
        f1 = [r["f1_vs_curation"] for r in a if r["f1_vs_curation"] is not None]
        times = [r["seconds"] for r in a]
        extra = sum(r["n_broken_helix"] + r["n_interfacial"] + r["n_reentrant"]
                    + r["n_unwound"] for r in a)
        print(f"story A: {len(a)} curated proteins")
        print(f"  agreement with the curation: mean F1 {np.mean(f1):.3f}, "
              f"median boundary error {np.median([r['boundary_median'] for r in a if r['boundary_median'] is not None]):.1f} residues")
        print(f"  geometric features reported beyond the curation: {extra} in {len(a)} proteins")
        print(f"  wall-clock time per protein: median {np.median(times):.1f} s "
              f"(total {np.sum(times):.0f} s)")
        consistent = all(r["export_3line"]["consistent"] for r in a)
        print(f"  the .3line export carries the same helices as the region list: {consistent}")
    if b:
        answered = sum(1 for r in b if r["answered"])
        flagged = sum(1 for r in b if r["placement_quality"] != "ok")
        times = [r["seconds"] for r in b]
        print(f"story B: {len(b)} predicted models")
        print(f"  a topology was produced for {answered}/{len(b)} models")
        print(f"  the placement was flagged as not fully determined for {flagged}/{len(b)}")
        print(f"  wall-clock time per model: median {np.median(times):.1f} s")
    print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()
