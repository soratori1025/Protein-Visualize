"""
An independent cohort for the membrane-frame benchmark.

    python bench/frame_cohort.py --out bench/frame --limit 250

The frame comparison in `membrane_frame_benchmark.py` is limited by the OPM coverage of
the topology cohort (58 of 150 structures).  This script builds a *separate* cohort for
the frame benchmark alone, drawn from the RCSB search API rather than from UniProtKB:
a random sample of PDB entries whose own keywords classify them as membrane proteins,
for which an OPM-oriented file with a reference bilayer exists.

The two cohorts are independent, which matters: the frame benchmark then measures the
placement on structures that were not used to choose or tune anything.

Output: `bench/frame/frame_cohort.json`, one record per structure with the deposited
and the OPM file paths, and `bench/frame/structures/` with the downloaded files.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bench"))

from benchmark_tm_topology import download  # noqa: E402

OPM_URL = "https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb"
RCSB_URL = "https://files.rcsb.org/download/{pdb}.cif"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
UA = {"User-Agent": "ProteinLab-benchmark/1.0 (academic validation)"}
# entries the RCSB keyword search returns for membrane proteins
MEMBRANE_QUERY = {
    "query": {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "struct_keywords.pdbx_keywords",
            "operator": "contains_phrase",
            "value": "membrane protein",
        },
    },
    "return_type": "entry",
    "request_options": {"paginate": {"start": 0, "rows": 10000}},
}


def rcsb_membrane_entries() -> list[str]:
    """PDB identifiers whose keyword classification includes 'membrane protein'."""
    body = json.dumps(MEMBRANE_QUERY).encode("utf-8")
    request = urllib.request.Request(SEARCH_URL, data=body,
                                     headers={**UA, "Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = json.loads(response.read().decode("utf-8"))
    return [hit["identifier"] for hit in payload.get("result_set", [])]


def build_one(job):
    pdb_id, out_dir = job
    structs = Path(out_dir)
    opm_path = structs / f"{pdb_id}.pdb"
    if not download(OPM_URL.format(pdb=pdb_id.lower()), opm_path):
        return None
    deposited = structs / f"{pdb_id}.cif"
    if not download(RCSB_URL.format(pdb=pdb_id.upper()), deposited):
        return None
    return {"pdb_id": pdb_id, "accession": "", "chain": "",
            "opm_file": str(opm_path), "structure_file": str(deposited)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="bench/frame")
    ap.add_argument("--limit", type=int, default=250)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--seed", type=int, default=20261022)
    args = ap.parse_args()

    out = Path(args.out)
    structs = out / "structures"
    structs.mkdir(parents=True, exist_ok=True)

    identifiers = rcsb_membrane_entries()
    random.Random(args.seed).shuffle(identifiers)
    print(f"[frame] {len(identifiers)} membrane-protein entries in the RCSB keyword search")

    # --- download the coordinate files (network bound, many files)
    records = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(build_one, [(pdb, str(structs)) for pdb in identifiers]):
            if result is None:
                continue
            records.append(result)
            if len(records) % 25 == 0:
                print(f"[frame] {len(records)} structures with an OPM reference")
            if len(records) >= args.limit:
                break
    records = records[: args.limit]
    (out / "frame_cohort.json").write_text(json.dumps(records, indent=2))
    print(f"[frame] {len(records)} structures written to {out / 'frame_cohort.json'}")


if __name__ == "__main__":
    main()
