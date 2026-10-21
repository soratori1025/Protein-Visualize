"""
Quick scan: how does the 3D Energy provider behave on AlphaFold models?

    python bench/case_alphafold_scan.py

Downloads the AlphaFold DB model of each accession in ``TARGETS`` and reports the
transfer free energy, the fitted thickness and the number of transmembrane segments
against the curated UniProtKB annotation.  It is a diagnostic, not a published result;
the two accessions quoted in the manuscript are analysed in detail by
``case_alphafold.py``.
"""
from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))
sys.path.insert(0, str(ROOT / "bench"))

from app.schemas.topology import TMParams  # noqa: E402
from app.services.topology.providers.tm.energy import EnergyTMProvider  # noqa: E402
from app.services.topology.providers.tm.uniprot import fetch_uniprot_entry  # noqa: E402
from app.services.topology.residues import load_residue_frame  # noqa: E402

import case_alphafold as ca  # noqa: E402

TARGETS = ["P31645", "P35354", "P05067", "P48751", "Q92508", "P00533", "P08100"]


def main() -> None:
    out = ROOT / "bench" / "results" / "alphafold"
    out.mkdir(parents=True, exist_ok=True)
    for accession in TARGETS:
        try:
            api = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
            request = urllib.request.Request(api, headers=ca.UA)
            with urllib.request.urlopen(request, timeout=60) as response:
                prediction = json.loads(response.read().decode())[0]
            path = out / f"AF-{accession}-F1.cif"
            ca.download(prediction["cifUrl"], path)
            frame = load_residue_frame(path, None)
            result = EnergyTMProvider().predict_tm(path, params=TMParams(), frame=frame)
            curated = len([f for f in fetch_uniprot_entry(accession).get("features", [])
                           if f.get("type") == "Transmembrane"])
            print(f"{accession} res={len(frame):5d} curated={curated:3d} "
                  f"computed={len(result.segments or []):3d}  {result.labeler}")
        except Exception as exc:
            print(f"{accession} ERR {type(exc).__name__}: {exc}")


if __name__ == "__main__":
    main()
