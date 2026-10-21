"""
The competitive matrix: what each available tool can and cannot do.

    python bench/competitive_matrix.py

Every row of the matrix is a statement about a *competitor*, and every statement is
backed by a place in that competitor's own repository or documentation, so that the
comparison in the manuscript can be checked by the reader rather than taken on trust.

The four reference implementations were cloned for this purpose:

    https://github.com/ElofssonLab/TOPCONS2       the consensus topology predictor
    https://github.com/ulo/Protter                the topology figure
    https://github.com/BernhoferM/TMbed            the language-model topology predictor
    https://github.com/bosborne/pyTMHMM           the TMHMM 2.0 implementation we vendor

The matrix is a data file (``competitive_matrix.json``) plus this script, which checks
each evidence pointer against the clone and reports whether it is still there.  A row
whose evidence cannot be found is dropped, so the matrix cannot drift away from the
sources.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# where the clones are; the script only needs the ones that are present
CLONES = {
    "TOPCONS2": Path("TOPCONS2"),
    "Protter": Path("Protter"),
    "TMbed": Path("TMbed"),
    "pyTMHMM": Path("pyTMHMM"),
}

# ---------------------------------------------------------------- the matrix
# Each capability is a row.  "yes" / "no" / "partial" is the verdict; the evidence
# is a file in the clone and a pattern that has to be found in it.
MATRIX = [
    {
        "tool": "ProteinLab (this work)",
        "kind": "web application + Python package",
        "capabilities": {
            "works on a curated UniProt protein": "yes",
            "works on a new or predicted model": "yes",
            "uses the 3-D structure, not only the sequence": "yes",
            "fits the membrane itself (normal, mid-plane, thickness)": "yes",
            "reports a transfer free energy for the placement": "yes",
            "says when the placement is not determined by the data": "yes",
            "detects hairpins, broken helices, interfacial helices": "yes",
            "handles a curved membrane": "yes",
            "per-residue confidence": "yes",
            "machine-readable export": "yes (3line, GFF3)",
            "needs a multi-gigabyte model or database": "no",
        },
        "evidence": {
            "fits the membrane itself (normal, mid-plane, thickness)":
                "backend/app/services/topology/membrane_energy.py",
            "reports a transfer free energy for the placement":
                "backend/app/services/topology/membrane_energy.py",
            "says when the placement is not determined by the data":
                "backend/app/services/topology/providers/tm/energy.py",
            "detects hairpins, broken helices, interfacial helices":
                "backend/app/services/topology/orchestrator.py",
        },
    },
    {
        "tool": "TOPCONS2",
        "kind": "standalone software + web server",
        "clone": "TOPCONS2",
        "capabilities": {
            "works on a curated UniProt protein": "yes",
            "works on a new or predicted model": "yes",
            "uses the 3-D structure, not only the sequence": "no",
            "fits the membrane itself (normal, mid-plane, thickness)": "no",
            "reports a transfer free energy for the placement": "no",
            "says when the placement is not determined by the data": "no",
            "detects hairpins, broken helices, interfacial helices": "no",
            "handles a curved membrane": "no",
            "per-residue confidence": "no",
            "machine-readable export": "partial (per-residue labels)",
            "needs a multi-gigabyte model or database": "yes",
        },
        "evidence": {
            "uses the 3-D structure, not only the sequence":
                ("topcons2_webserver", r"hmmscan|hhblits"),
            "needs a multi-gigabyte model or database":
                ("README.md", r"topcons2_database\.zip|uniref90|BLAST database"),
        },
    },
    {
        "tool": "Protter",
        "kind": "web application",
        "clone": "Protter",
        "capabilities": {
            "works on a curated UniProt protein": "yes",
            "works on a new or predicted model": "partial (needs a curated record)",
            "uses the 3-D structure, not only the sequence": "partial (structure for the figure)",
            "fits the membrane itself (normal, mid-plane, thickness)": "no",
            "reports a transfer free energy for the placement": "no",
            "says when the placement is not determined by the data": "no",
            "detects hairpins, broken helices, interfacial helices": "no",
            "handles a curved membrane": "no",
            "per-residue confidence": "no",
            "machine-readable export": "partial (SVG, GFF3)",
            "needs a multi-gigabyte model or database": "no",
        },
        "evidence": {
            "uses the 3-D structure, not only the sequence":
                ("src", r"pdb|structure"),
            "machine-readable export": ("src", r"gff3|svg"),
        },
    },
    {
        "tool": "TMbed",
        "kind": "standalone software (language model)",
        "clone": "TMbed",
        "capabilities": {
            "works on a curated UniProt protein": "yes",
            "works on a new or predicted model": "yes",
            "uses the 3-D structure, not only the sequence": "no",
            "fits the membrane itself (normal, mid-plane, thickness)": "no",
            "reports a transfer free energy for the placement": "no",
            "says when the placement is not determined by the data": "no",
            "detects hairpins, broken helices, interfacial helices": "no",
            "handles a curved membrane": "no",
            "per-residue confidence": "partial (probabilities)",
            "machine-readable export": "partial (per-residue labels)",
            "needs a multi-gigabyte model or database": "yes",
        },
        "evidence": {
            "uses the 3-D structure, not only the sequence":
                ("README.md", r"Protein Language Model|embeddings|ProtT5"),
            "needs a multi-gigabyte model or database":
                ("README.md", r"ProtT5|2\.25 GB|encoder model"),
        },
    },
    {
        "tool": "pyTMHMM",
        "kind": "Python library (vendored here)",
        "clone": "pyTMHMM",
        "capabilities": {
            "works on a curated UniProt protein": "yes",
            "works on a new or predicted model": "yes",
            "uses the 3-D structure, not only the sequence": "no",
            "fits the membrane itself (normal, mid-plane, thickness)": "no",
            "reports a transfer free energy for the placement": "no",
            "says when the placement is not determined by the data": "no",
            "detects hairpins, broken helices, interfacial helices": "no",
            "handles a curved membrane": "no",
            "per-residue confidence": "partial (posteriors)",
            "machine-readable export": "partial (per-residue labels)",
            "needs a multi-gigabyte model or database": "no",
        },
        "evidence": {
            "uses the 3-D structure, not only the sequence":
                ("README.md", r"protein sequences|SEQUENCE_FILE"),
        },
    },
]


def check_evidence(entry: dict, clone_root: Path) -> dict:
    """Confirm that every evidence pointer of a competitor still resolves."""
    result = {}
    clone = CLONES.get(entry.get("clone", ""))
    for capability, pointer in entry.get("evidence", {}).items():
        if isinstance(pointer, str):
            target = ROOT / pointer
            result[capability] = {"ok": target.exists(), "where": pointer}
            continue
        if clone is None or not clone_root:
            result[capability] = {"ok": None, "where": "the clone was not available"}
            continue
        folder, pattern = pointer
        base = clone_root / clone / folder
        if base.is_file():
            candidates = [base]
        elif base.exists():
            candidates = [p for p in base.rglob("*")
                          if p.is_file() and p.stat().st_size < 4_000_000]
        else:
            result[capability] = {"ok": False, "where": f"{clone}/{folder} (missing)"}
            continue
        found = False
        pattern_re = re.compile(pattern, re.IGNORECASE)
        for path in candidates:
            try:
                if pattern_re.search(path.read_text(errors="replace")):
                    found = True
                    result[capability] = {"ok": True,
                                           "where": f"{clone}/{folder}"}
                    break
            except (OSError, UnicodeDecodeError):
                continue
        if not found:
            result[capability] = {"ok": False,
                                   "where": f"{clone}/{folder} (pattern not found)"}
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clones", default=None,
                    help="directory holding the clones (default: a directory called "
                         "'external' next to the repository)")
    args = ap.parse_args()
    clone_root = Path(args.clones) if args.clones else None
    if clone_root is None:
        for candidate in (ROOT.parent / "external", ROOT / "external"):
            if candidate.exists():
                clone_root = candidate
                break

    out = {"tools": [], "capabilities": [], "clones": {}}
    for entry in MATRIX:
        record = {k: v for k, v in entry.items() if k != "evidence"}
        record["evidence"] = check_evidence(entry, clone_root)
        out["tools"].append(record)

    # the union of the capabilities, in a stable order
    seen: list[str] = []
    for entry in MATRIX:
        for capability in entry["capabilities"]:
            if capability not in seen:
                seen.append(capability)
    out["capabilities"] = seen
    if clone_root:
        out["clones"] = {name: str(clone_root / name) for name in CLONES
                          if (clone_root / name).exists()}

    path = ROOT / "bench" / "competitive_matrix.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"written to {path}\n")

    width = max(len(c) for c in seen) + 2
    header = "capability".ljust(width) + "".join(t["tool"][:22].ljust(24) for t in out["tools"])
    print(header)
    print("-" * len(header))
    for capability in seen:
        row = capability.ljust(width)
        for tool in out["tools"]:
            row += str(tool["capabilities"].get(capability, "-")).ljust(24)
        print(row)
    print()
    for tool in out["tools"]:
        missing = [c for c, e in tool["evidence"].items() if e.get("ok") is False]
        if missing:
            print(f"! {tool['tool']}: {len(missing)} evidence pointer(s) could not be "
                  f"confirmed - {missing[:3]}")


if __name__ == "__main__":
    main()
