"""
Case study: an AlphaFold model with no experimental structure.

    python bench/case_alphafold.py --accession Q92508

Downloads the AlphaFold DB model of the given accession, runs the 3D Energy provider
on it, and repeats the run with the pLDDT filter disabled, so that the effect of the
filter on a predicted model can be measured rather than asserted.  The number of
transmembrane segments is compared with the curated UniProtKB annotation.
"""
from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT / "backend"))

import app.services.topology.membrane_energy as ME  # noqa: E402
from app.schemas.topology import TMParams  # noqa: E402
from app.services.topology.providers.tm.energy import EnergyTMProvider  # noqa: E402
from app.services.topology.providers.tm.uniprot import fetch_uniprot_entry  # noqa: E402
from app.services.topology.residues import load_residue_frame  # noqa: E402

UA = {"User-Agent": "ProteinLab-benchmark/1.0 (academic validation)"}


def download(url: str, dest: Path) -> bool:
    if dest.exists() and dest.stat().st_size > 1000:
        return True
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=180) as resp:
        data = resp.read()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", default="Q92508")
    ap.add_argument("--out", default="bench/results/alphafold")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    api = f"https://alphafold.ebi.ac.uk/api/prediction/{args.accession}"
    req = urllib.request.Request(api, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as resp:
        pred = json.loads(resp.read().decode("utf-8"))[0]
    path = out / f"AF-{args.accession}-F1.cif"
    print("downloading", pred["cifUrl"])
    download(pred["cifUrl"], path)

    frame = load_residue_frame(path, None)
    print(f"chain {frame.chain_id}, {len(frame)} residues, "
          f"predicted_model={frame.predicted_model}, "
          f"B-factor range {frame.b_factors.min():.1f}-{frame.b_factors.max():.1f}")

    entry = fetch_uniprot_entry(args.accession)
    curated = [f for f in entry.get("features", []) if f.get("type") == "Transmembrane"]
    print(f"curated UniProtKB transmembrane helices: {len(curated)}")

    provider = EnergyTMProvider()
    r_with = provider.predict_tm(path, params=TMParams(), frame=frame)
    print(f"with the pLDDT filter:    {len(r_with.segments or [])} TM segments, "
          f"{r_with.labeler}")
    for w in r_with.warnings:
        print("   -", w)

    original = ME.MIN_PLDDT
    ME.MIN_PLDDT = -1.0        # disable the filter
    try:
        r_without = EnergyTMProvider().predict_tm(path, params=TMParams(), frame=frame)
    finally:
        ME.MIN_PLDDT = original
    print(f"without the pLDDT filter: {len(r_without.segments or [])} TM segments, "
          f"{r_without.labeler}")
    for w in r_without.warnings:
        print("   -", w)

    (out / f"alphafold_{args.accession}.json").write_text(json.dumps({
        "accession": args.accession,
        "model": pred["cifUrl"],
        "n_residues": len(frame),
        "n_curated_tm": len(curated),
        "n_tm_with_plddt_filter": len(r_with.segments or []),
        "n_tm_without_plddt_filter": len(r_without.segments or []),
        "labeler_with": r_with.labeler,
        "labeler_without": r_without.labeler,
        "warnings_with": r_with.warnings,
        "warnings_without": r_without.warnings,
    }, indent=2))


if __name__ == "__main__":
    main()
