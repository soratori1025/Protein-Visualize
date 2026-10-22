"""
Build the source archive that accompanies the manuscript as a supplementary file.

    python bench/make_source_archive.py --out manuscript/ProteinLab-source-archive.zip

BMC Bioinformatics asks that the software be included with the submitted manuscript as
additional files, so that it remains available even if the repository moves.  The archive
contains everything needed to run the package and to reproduce every number in the
paper, and nothing that belongs to the paper itself: no manuscript, no supplementary
information, no working document.

What is in
  backend/            the FastAPI service, the topology engine, the unit tests
  protein_engine/     the shared structure-parsing and secondary-structure library
  frontend/src/       the React + TypeScript application
  bench/              the validation harness, its results and the figures
  tools/              the prebuilt DSSP and STRIDE executables for Windows
  LICENSE, CITATION.cff, README.md, .github/  the licence, the citation metadata, CI

What is not in
  .git, node_modules, .venv, the downloaded coordinate files, the manuscript, the
  supplementary information, and the TMHMM 2.0 model file, which DTU distributes under
  an academic licence that does not permit redistribution (the user installs it with
  backend/scripts/install_tmhmm_model.py).
"""
from __future__ import annotations

import argparse
import hashlib
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

INCLUDE_DIRS = ["backend", "protein_engine", "frontend/src", "bench", ".github"]
INCLUDE_FILES = ["LICENSE", "CITATION.cff", "README.md", ".gitignore"]
# The prebuilt Windows executables and the sample coordinate files are not part of the
# archive: the executables are third-party binaries (DSSP and STRIDE have their own
# distributions, and the README says how to install them), and the sample structures are
# downloaded by the benchmark.  The archive therefore stays small enough to be attached
# to a submission.
EXCLUDE_DIRS = ["tools", "backend/data/uploads"]
EXCLUDE_PARTS = {
    "__pycache__", "node_modules", ".venv", "dist", "structures", "opm", "alphafold",
    ".pytest_cache", "results_n64", "figures_test", ".mypy_cache", "data",
}
EXCLUDE_SUFFIXES = {".pyc", ".pyo", ".log", ".tsbuildinfo"}
MAX_FILE_MB = 5.0
# Third-party files that may not be redistributed.  The guard below refuses to build an
# archive that contains one, so a future copy cannot slip in unnoticed.
EXCLUDE_FILES = {
    "backend/app/services/topology/providers/tm/_tmhmm/TMHMM2.0.model",
}
NEVER_REDISTRIBUTE = {"TMHMM2.0.model"}


def wanted(absolute: Path) -> bool:
    relative = absolute.relative_to(ROOT)
    if relative.as_posix() in EXCLUDE_FILES:
        return False
    if any(part in EXCLUDE_PARTS for part in relative.parts):
        return False
    if any(str(relative).startswith(d) for d in EXCLUDE_DIRS):
        return False
    if absolute.suffix in EXCLUDE_SUFFIXES:
        return False
    if absolute.stat().st_size > MAX_FILE_MB * 1e6:
        return False
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="manuscript/ProteinLab-source-archive.zip")
    args = ap.parse_args()
    out = Path(args.out)
    if not out.is_absolute():
        out = ROOT / out
    out.parent.mkdir(parents=True, exist_ok=True)

    files: list[Path] = []
    for name in INCLUDE_DIRS:
        base = ROOT / name
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if path.is_file() and wanted(path):
                files.append(path)
    for name in INCLUDE_FILES:
        path = ROOT / name
        if path.exists():
            files.append(path)

    leaked = sorted(p.relative_to(ROOT).as_posix() for p in files
                    if p.name in NEVER_REDISTRIBUTE)
    if leaked:
        raise SystemExit("refusing to build: not redistributable: " + ", ".join(leaked))

    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for path in files:
            arcname = Path("ProteinLab") / path.relative_to(ROOT)
            zf.write(path, arcname)

    digest = hashlib.sha256(out.read_bytes()).hexdigest()
    print(f"archive: {out}")
    print(f"  files : {len(files)}")
    print(f"  size  : {out.stat().st_size / 1e6:.1f} MB")
    print(f"  sha256: {digest}")
    print()
    print("Upload this file to Zenodo with the 'anonymous' option (or an equivalent")
    print("anonymous repository) to obtain the persistent identifier for the")
    print("Availability and requirements section, and attach the same file to the")
    print("submission as an additional file.")


if __name__ == "__main__":
    main()
