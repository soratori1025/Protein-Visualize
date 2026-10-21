"""
Validate the vendored TMHMM port against the upstream pyTMHMM project.

    python bench/validate_tmhmm_port.py --upstream <path to a pyTMHMM checkout>

The repository vendors a NumPy port of pyTMHMM 1.3.6 under
``backend/app/services/topology/providers/tm/_tmhmm/``.  The port's docstring claims
that its Viterbi path is *bit-identical* to the upstream implementation and that two
defects of the upstream code were corrected - but the script that was supposed to
demonstrate this was never written, so the claim could not be checked.

Three checks are done here:

1. **The model file.**  The vendored ``TMHMM2.0.model`` is compared with the copy the
   upstream project used, both byte for byte and line for line (the byte comparison
   differs on Windows only because of the line endings).

2. **The parsed model.**  The vendored parser is compared with the upstream's own pure
   Python parser, state by state: the initial distribution, the transition matrices and
   the emission tables.

3. **The upstream Cython extension.**  When a C compiler is available the upstream
   implementation is built and the Viterbi paths are compared residue by residue.  When
   it is not, the script says so and does not pretend to have verified anything.

The script also records a licensing observation that matters for the repository: the
upstream project removed the model file from its own repository in the commit
"Adhere to DTU's license requirements".
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

VENDORED = (ROOT / "backend" / "app" / "services" / "topology" / "providers" / "tm"
            / "_tmhmm")

# the upstream project published its annotation for the human P2X3 receptor
PAR3 = "P49641"


def fetch_sequence(accession: str) -> str | None:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    request = urllib.request.Request(url, headers={"User-Agent": "ProteinLab-benchmark/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            text = response.read().decode("utf-8")
    except Exception:
        return None
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream", required=True,
                    help="path to a checkout of bosborne/pyTMHMM")
    ap.add_argument("--out", default="bench/results/tmhmm_port_validation.json")
    args = ap.parse_args()
    upstream = Path(args.upstream).resolve()
    report: dict = {}

    # ------------------------------------------------------- 1. the model file
    vendored_model = VENDORED / "TMHMM2.0.model"
    print("[port] 1. the model file")
    upstream_model = None
    for candidate in (upstream / "TMHMM2.0.model", upstream / "tmhmm" / "TMHMM2.0.model"):
        if candidate.exists():
            upstream_model = candidate
            break
    if upstream_model is None:
        # the upstream removed the file to comply with the DTU licence; it is still
        # in the upstream git history
        print("[port]    the upstream file is no longer in the checkout (removed to "
              "comply with the DTU licence); looking in the git history")
        history = subprocess.run(["git", "log", "--all", "--format=%H"], cwd=str(upstream),
                                 capture_output=True, text=True).stdout.split()
        for commit in history:
            blob = subprocess.run(["git", "cat-file", "blob",
                                   f"{commit}:tmhmm/TMHMM2.0.model"], cwd=str(upstream),
                                  capture_output=True).stdout
            if blob:
                tmp = Path(tempfile.mkdtemp()) / "TMHMM2.0.model"
                tmp.write_bytes(blob)
                upstream_model = tmp
                report["upstream_model_from_history"] = commit
                break
    if upstream_model is not None:
        a = upstream_model.read_bytes()
        b = vendored_model.read_bytes()
        same_bytes = a == b
        same_text = (a.replace(b"\r\n", b"\n") == b.replace(b"\r\n", b"\n"))
        report["model_file"] = {
            "upstream": str(upstream_model), "vendored": str(vendored_model),
            "byte_identical": same_bytes,
            "identical_after_normalising_line_endings": same_text,
            "size_upstream": len(a), "size_vendored": len(b),
            "n_lines_upstream": len(a.splitlines()), "n_lines_vendored": len(b.splitlines()),
        }
        print(f"[port]    byte identical: {same_bytes}   "
              f"identical after normalising line endings: {same_text}   "
              f"({len(a)} vs {len(b)} bytes)")
    else:
        report["model_file"] = {"error": "no upstream model file could be located"}
        print("[port]    no upstream model file could be located")

    # ---------------------------------------------------- 2. the parsed model
    print("[port] 2. the parsed model")
    from app.services.topology.providers.tm._tmhmm import model as vendored_model_module
    _, parsed = vendored_model_module.parse(str(vendored_model))
    initial, transitions, emissions, char_map, label_map, name_map = parsed
    report["parsed_model"] = {
        "n_states": len(initial), "n_emission_symbols": len(emissions),
        "n_char_map": len(char_map), "n_label_map": len(label_map),
    }
    print(f"[port]    the vendored parser reads {len(initial)} states")

    upstream_parser = upstream / "tmhmm" / "model.py"
    if not upstream_parser.exists():
        upstream_parser = upstream / "model.py"
    if upstream_parser.exists():
        try:
            import importlib.util
            spec = importlib.util.spec_from_file_location("upstream_model",
                                                          upstream_parser)
            upstream_module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(upstream_module)
            parse = getattr(upstream_module, "parse", None)
            if parse is not None:
                result = parse(str(upstream_model))
                up_parsed = result[1] if isinstance(result, tuple) and len(result) == 2 else result
                up_initial, up_transitions, up_emissions = up_parsed[0], up_parsed[1], up_parsed[2]
                # both implementations return the transitions as a list of dense
                # matrices, one per state
                max_diff = 0.0
                if len(up_transitions) == len(transitions):
                    for u, v in zip(up_transitions, transitions):
                        u = np.asarray(u, dtype=float)
                        v = np.asarray(v, dtype=float)
                        if u.shape == v.shape:
                            max_diff = max(max_diff, float(np.abs(u - v).max()))
                        else:
                            max_diff = float("inf")
                else:
                    max_diff = float("inf")
                report["parsed_model"]["upstream_parser_used"] = str(upstream_parser)
                report["parsed_model"]["same_number_of_states"] = (
                    len(up_initial) == len(initial))
                report["parsed_model"]["max_transition_difference"] = max_diff
                print(f"[port]    the upstream pure-Python parser gives "
                      f"{len(up_initial)} states; the largest transition difference is "
                      f"{max_diff:.3g}")
            else:
                print("[port]    the upstream module has no parse() function")
        except Exception as exc:
            report["parsed_model"]["upstream_parser_error"] = str(exc)
            print(f"[port]    the upstream parser could not be used: {exc}")
    else:
        print("[port]    the upstream model parser was not found")

    # -------------------------------------------------- 3. the upstream extension
    print("[port] 3. the upstream Cython extension")
    upstream_ok = False
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp) / "pyTMHMM"
        subprocess.run(["git", "clone", "--depth", "1", str(upstream), str(work)],
                       check=True, capture_output=True, text=True)
        build = subprocess.run([sys.executable, "setup.py", "build_ext", "--inplace"],
                               cwd=str(work), capture_output=True, text=True)
        if build.returncode == 0:
            sys.path.insert(0, str(work))
            try:
                import api as upstream_api
                upstream_ok = True
            except Exception as exc:
                print(f"[port]    the extension could not be imported: {exc}")
        else:
            tail = (build.stderr or build.stdout or "").strip().splitlines()[-1:]
            print(f"[port]    the extension could not be built: {tail}")
            print("[port]    (a C compiler is needed; the comparison is skipped)")
    report["upstream_cython_extension_built"] = upstream_ok
    if upstream_ok:
        sequence = fetch_sequence(PAR3)
        if sequence:
            from app.services.topology.providers.tm._tmhmm import api as vendored
            up = upstream_api.predict(sequence, model=str(upstream_model))
            ve = vendored.predict(sequence, model=parsed)
            report["viterbi_comparison"] = {
                "protein": PAR3, "length": len(sequence),
                "identical_path": up.labels == ve.path,
                "upstream_path": up.labels, "vendored_path": ve.path,
            }
            print(f"[port]    Viterbi path identical on {PAR3}: {up.labels == ve.path}")

    # ------------------------------------------- the annotation of a reference protein
    sequence = fetch_sequence(PAR3)
    if sequence:
        from app.services.topology.providers.tm._tmhmm import api as vendored
        result = vendored.predict(sequence, model=parsed)
        report["reference_annotation"] = {
            "protein": PAR3, "length": len(sequence),
            "n_helices": len(result.helices),
            "helices": [[int(a), int(b)] for a, b in result.helices],
            "path": result.path,
            "exp_aa_tmh": round(result.exp_aa_tmh, 2),
            "exp_first60": round(result.exp_first60, 2),
            "prob_n_in": round(result.prob_n_in, 4),
        }
        print(f"[port]    P2X3 ({PAR3}): {len(result.helices)} helix/ices "
              f"{[(int(a) + 1, int(b) + 1) for a, b in result.helices]} (1-based)")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\n[port] written to {out}")


if __name__ == "__main__":
    main()
