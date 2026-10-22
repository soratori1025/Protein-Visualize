"""Install the TMHMM 2.0 model file that the vendored port reads.

DTU Health Tech distributes TMHMM 2.0 under an academic licence that does not permit
redistribution, so ProteinLab does **not** ship the model.  Download TMHMM 2.0 from

    https://services.healthtech.dtu.dk/services/TMHMM-2.0/

unpack it, and run

    python backend/scripts/install_tmhmm_model.py path/to/TMHMM2.0.model

The file is parsed once to make sure it is the model the code expects, then copied to
``<repository>/.tmhmm/TMHMM2.0.model``.  That folder is ignored by git and is found
automatically by the provider, so no environment variable is needed.  Alternatively,
point ``TMHMM_MODEL`` at the file wherever it lives.

    python backend/scripts/install_tmhmm_model.py --check

reports whether a usable model is already installed.  Without the model the TMHMM
provider answers HTTP 503 and every other provider keeps working; the benchmark tables
in the paper can still be reproduced from the cached predictions in ``bench/results``.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[1]
REPOSITORY = BACKEND.parent
TARGET = REPOSITORY / ".tmhmm" / "TMHMM2.0.model"

sys.path.insert(0, str(BACKEND))

from app.services.topology.providers.tm._tmhmm import MODEL_ENV, model_path  # noqa: E402
from app.services.topology.providers.tm._tmhmm.model import parse  # noqa: E402


def _read(path: Path):
    """Parse the model, raising a readable error if this is not a TMHMM 2.0 model."""
    try:
        _, model = parse(str(path))
    except Exception as err:                      # noqa: BLE001 - report, do not traceback
        raise SystemExit(f"{path} is not a readable TMHMM 2.0 model file: {err}") from err
    states = len(model[0]) if isinstance(model, tuple) else None
    return states


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("source", nargs="?",
                    help="TMHMM2.0.model from the DTU TMHMM 2.0 package, or the "
                         "package directory that contains it")
    ap.add_argument("--check", action="store_true",
                    help="report whether a model is installed, and exit")
    args = ap.parse_args()

    current = model_path()
    if args.check:
        if current.is_file():
            print(f"model: {current}")
            print(f"states: {_read(current)}")
            return 0
        print(f"model: not installed (looked for {current})")
        print("download TMHMM 2.0 from DTU and pass the model file to this script")
        return 1

    if not args.source:
        ap.error("give the path to TMHMM2.0.model (or use --check)")

    source = Path(args.source).expanduser()
    if source.is_dir():
        source = source / "TMHMM2.0.model"
    if not source.is_file():
        raise SystemExit(f"no such file: {source}")

    states = _read(source)
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, TARGET)
    print(f"installed: {TARGET}")
    print(f"states   : {states}")
    print()
    print("The file is git-ignored.  The TMHMM provider now runs; without it the")
    print("provider reports HTTP 503 and every other provider is unaffected.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
