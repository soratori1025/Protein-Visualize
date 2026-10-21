"""Command-line interface to the ProteinLab topology engine.

    proteinlab topology 1bl8.cif --algorithm 3d_energy --chain A
    proteinlab topology 1bl8.cif --algorithm 3d_energy --json out.json
    proteinlab secondary-structure 1bl8.cif --method dssp
    proteinlab capabilities

The command line exposes exactly the pipeline the web application runs, so that every
number in the paper can be reproduced without a browser and without the API.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# The DSSP and STRIDE executables are not Python packages; make the bundled ones
# discoverable exactly as the API does.
import app.main  # noqa: F401  (side effect: PATH and LIBCIFPP_DATA_DIR)

from app.schemas.topology import TMParams
from app.services.topology.service import predict_topology

ALGORITHMS = ["3d_energy", "3d_slab_geom", "kyte_doolittle_seq", "tmhmm", "uniprot_api"]
SS_ALGORITHMS = ["dssp", "stride", "none"]


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("structure", help="PDB or mmCIF file")
    parser.add_argument("--chain", default=None, help="chain identifier to analyse")
    parser.add_argument("--uniprot", default=None, help="UniProtKB accession")
    parser.add_argument("--thickness", type=float, default=None,
                        help="membrane thickness, 10-60 A")
    parser.add_argument("--min-membrane-score", type=float, default=None)
    parser.add_argument("--min-cross-span", type=float, default=None)
    parser.add_argument("--full-cross-frac", type=float, default=None)
    parser.add_argument("--treat-turn-as-helix", action="store_true")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="proteinlab",
        description="Structure-derived transmembrane topology annotation.")
    sub = parser.add_subparsers(dest="command", required=True)

    topology = sub.add_parser("topology", help="annotate the transmembrane topology")
    _add_common(topology)
    topology.add_argument("--algorithm", default="3d_energy", choices=ALGORITHMS)
    topology.add_argument("--secondary-structure", default="dssp", choices=SS_ALGORITHMS)
    topology.add_argument("--json", default=None, help="write the full result as JSON")

    secondary = sub.add_parser("secondary-structure", help="assign secondary structure")
    _add_common(secondary)
    secondary.add_argument("--method", default="dssp", choices=["dssp", "stride"])
    secondary.add_argument("--json", default=None)

    sub.add_parser("capabilities", help="report which external tools are installed")
    return parser


def _overrides(args: argparse.Namespace) -> dict:
    """Map the command-line options onto the names `build_params` expects."""
    mapping = {
        "thickness": "thickness",
        "min_membrane_score": "min_membrane_score",
        "min_cross_span": "min_cross_span",
        "full_cross_frac": "full_cross_frac",
    }
    out = {}
    for attribute, name in mapping.items():
        value = getattr(args, attribute, None)
        if value is not None:
            out[name] = value
    if getattr(args, "treat_turn_as_helix", False):
        out["treat_turn_as_helix"] = True
    return out


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "capabilities":
        from protein_engine.secondary_structure.dssp import DSSPMethod
        from protein_engine.secondary_structure.stride import STRIDEMethod
        from shutil import which

        for name, method in (("DSSP", DSSPMethod), ("STRIDE", STRIDEMethod)):
            resolved = method().executable
            found = bool(resolved) and (Path(resolved).is_file() or which(resolved))
            print(f"{name:8s} {'available' if found else 'not installed':16s} {resolved}")
        return 0

    structure = Path(args.structure)
    if not structure.exists():
        print(f"no such file: {structure}", file=sys.stderr)
        return 2

    if args.command == "secondary-structure":
        from app.services.topology.providers.ss.dssp import DSSPProvider
        from app.services.topology.providers.ss.stride import STRIDEProvider

        provider = DSSPProvider() if args.method == "dssp" else STRIDEProvider()
        assignment = provider.assign(structure)
        payload = {"method": provider.name, "coverage": assignment.coverage,
                   "residues": assignment.codes}
        if args.json:
            Path(args.json).write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"{provider.name}: {assignment.coverage:.0%} of the residues assigned")
        return 0

    response = predict_topology(structure, tm_algo=args.algorithm,
                                ss_algo=args.secondary_structure, flow_type="consensus",
                                chain_id=args.chain, uniprot_id=args.uniprot,
                                **_overrides(args))
    if args.json:
        Path(args.json).write_text(json.dumps(response.model_dump(mode="json"), indent=2),
                                   encoding="utf-8")
    print(response.labeler)
    for region in response.regions:
        if region.type == "Transmembrane":
            print(f"  TM {region.start}-{region.end}  side={region.side}  "
                  f"confidence={region.confidence}")
    for warning in response.warnings:
        print(f"  ! {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
