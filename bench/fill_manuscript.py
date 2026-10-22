"""
Fill the manuscript placeholders from the committed benchmark results.

    python bench/fill_manuscript.py --results bench/results \
        --template manuscript/ProteinLab_manuscript_submission.tex \
        --out manuscript/ProteinLab_manuscript_submission.tex

Every number in the manuscript is read from a file produced by a committed script, so
the paper cannot silently disagree with the data.  The script reports every placeholder
it could not fill, so that nothing is left to memory.

The placeholders are the ``\\todo{...}`` macros of the LaTeX source: the author
metadata and the citation checks keep their marker, and every numeric placeholder is
replaced by the value from the result files.

The competitive matrix (Table 5) is not a placeholder: it is a generated file of its own,
``manuscript/parts/table_5_matrix.tex``, which the manuscript inputs and which
``bench/make_matrix_table.py`` regenerates from ``bench/competitive_matrix.json``.

Formatting rules
----------------
* a count is an integer;
* a fraction is written with the number of decimals that the source carries;
* a *P* value is written as "< 0.001" below that threshold, otherwise with three
  decimals;
* a length keeps the unit of its source.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np

PLACEHOLDER = re.compile(r"\\todo\{(?:\\detokenize\{([^{}]*)\}|([^{}]*))\}")


def fmt(value, decimals=None, unit="", percent=False, pvalue=False):
    if value is None:
        return "n/a"
    if pvalue:
        if value < 0.001:
            return "$< 0.001$"
        return f"${value:.3f}$"
    if percent:
        return f"{value:.1f} %"
    if decimals is None:
        if isinstance(value, float):
            decimals = 3 if abs(value) < 1 else 1
            text = f"{value:.{decimals}f}"
        else:
            text = str(value)
    else:
        text = f"{value:.{decimals}f}"
    return f"{text}{unit}"


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def collect(results: Path) -> dict:
    """Every value the manuscript needs, read from the result files."""
    out: dict = {}
    summary = {r["method"]: r for r in read_csv(results / "summary.csv")}
    out["_summary"] = summary

    def g(method, key, cast=float):
        try:
            return cast(summary[method][key])
        except (KeyError, TypeError, ValueError):
            return None

    # ---- table 1, one row per method
    for method in summary:
        for key, name in (("n_ref", "nref"), ("macro_recall", "recall"),
                           ("macro_precision", "precision"), ("macro_f1", "f1"),
                           ("residue_recall", "residue"),
                           ("median_max_boundary_error", "bnd"),
                           ("frac_within_5_residues", "within5"),
                           ("mean_runtime_s", "time")):
            out[f"T|{method}|{name}"] = g(method, key)

    for method, key in (("3d_energy", "F1_ENERGY"), ("3d_slab_geom", "F1_SLAB"),
                        ("tmhmm", "F1_TMHMM"), ("toppred", "F1_TOPPRED"),
                        ("topcons_style_consensus", "F1_CONSENSUS"),
                        ("kyte_doolittle_seq", "F1_KD")):
        out[key] = g(method, "macro_f1")

    # ---- paired comparisons
    paired_path = results / "paired_comparisons.json"
    paired = json.loads(paired_path.read_text()) if paired_path.exists() else []
    wanted = [("3d_energy", "tmhmm", "energy_vs_tmhmm"),
              ("3d_energy", "kyte_doolittle_seq", "energy_vs_kd"),
              ("tmhmm", "kyte_doolittle_seq", "tmhmm_vs_kd"),
              ("3d_energy", "3d_slab_geom", "energy_vs_slab")]
    for a, b, prefix in wanted:
        for rec in paired:
            if rec["method_a"] == a and rec["method_b"] == b:
                out[f"{prefix}_DIFF"] = rec["mean_difference"]
                out[f"{prefix}_LO"] = rec["ci95_low"]
                out[f"{prefix}_HI"] = rec["ci95_high"]
                out[f"{prefix}_P"] = rec["wilcoxon_p"]
                out[f"{prefix}_PBH"] = rec.get("wilcoxon_p_bh")
                out[f"{prefix}_TIES"] = rec["n_tied"]
                break

    # ---- transfer free energy
    records = json.loads((results / "benchmark.json").read_text())
    out["_records"] = records
    dg = [r["methods"]["3d_energy"]["dG_transfer"] for r in records
          if "dG_transfer" in r["methods"].get("3d_energy", {})]
    out["DG_MIN"] = min(dg) if dg else None
    out["DG_MAX"] = max(dg) if dg else None
    out["DG_MED"] = float(np.median(dg)) if dg else None
    out["DG_CALLED"] = f"{sum(1 for v in dg if v <= -10)} / {len(dg)}" if dg else None

    soluble_path = results / "soluble_control_summary.json"
    if soluble_path.exists():
        soluble = json.loads(soluble_path.read_text())
        out["SDG_MIN"] = soluble.get("dG_min")
        out["SDG_MAX"] = soluble.get("dG_max")
        out["SDG_MED"] = soluble.get("dG_median")
        out["SDG_CALLED_S"] = f"{soluble.get('n_called_membrane')} / {soluble.get('n')}"
        out["SDG_N"] = soluble.get("n")

    # ---- the membrane frame
    frame_path = results / "membrane_frame_benchmark.json"
    if frame_path.exists():
        frame = json.loads(frame_path.read_text())
        out["_frame"] = frame
        energy = [r["3d_energy"]["normal_error_deg"] for r in frame
                  if "normal_error_deg" in r.get("3d_energy", {})]
        slab = [r["3d_slab_geom"]["normal_error_deg"] for r in frame
                if "normal_error_deg" in r.get("3d_slab_geom", {})]
        offset = [r["3d_energy"]["midplane_offset_A"] for r in frame
                  if "midplane_offset_A" in r.get("3d_energy", {})]
        thickness = [r["3d_energy"]["thickness_A"] for r in frame
                     if "thickness_A" in r.get("3d_energy", {})]
        reference = [r["reference_thickness_A"] for r in frame
                     if "reference_thickness_A" in r]
        rmsd = [r["superposition_rmsd_A"] for r in frame
                if "superposition_rmsd_A" in r]
        out["FRAME_N"] = len(frame)
        out["FRAME_N_OPM"] = len(energy)
        out["FRAME_NORMAL_DEG"] = float(np.median(energy)) if energy else None
        out["FRAME_NORMAL_MEAN"] = float(np.mean(energy)) if energy else None
        out["FRAME_NORMAL_P90"] = float(np.percentile(energy, 90)) if energy else None
        out["FRAME_NORMAL_MAX"] = float(np.max(energy)) if energy else None
        out["FRAME_OFFSET_A"] = float(np.median(offset)) if offset else None
        out["FRAME_OFFSET_P90"] = float(np.percentile(offset, 90)) if offset else None
        out["SLAB_NORMAL_DEG"] = float(np.median(slab)) if slab else None
        out["SLAB_NORMAL_MEAN"] = float(np.mean(slab)) if slab else None
        out["SLAB_NORMAL_P90"] = float(np.percentile(slab, 90)) if slab else None
        out["THICKNESS_MEDIAN"] = float(np.median(thickness)) if thickness else None
        out["THICKNESS_OPM"] = float(np.median(reference)) if reference else None
        out["THICKNESS_GAP"] = (round(out["THICKNESS_OPM"] - out["THICKNESS_MEDIAN"], 1)
                                if out.get("THICKNESS_OPM") and out.get("THICKNESS_MEDIAN")
                                else None)
        out["FRAME_RMSD"] = float(np.median(rmsd)) if rmsd else None
        out["FRAME_RMSD_MAX"] = float(np.max(rmsd)) if rmsd else None
        out["FR1"] = out.get("FRAME_NORMAL_DEG")
        out["FR2"] = out.get("FRAME_NORMAL_MEAN")
        out["FR3"] = out.get("FRAME_NORMAL_P90")
        out["FR4"] = out.get("FRAME_NORMAL_MAX")
        out["FR5"] = out.get("FRAME_OFFSET_A")
        out["FR6"] = out.get("FRAME_OFFSET_P90")
        out["FR7"] = out.get("THICKNESS_MEDIAN")
        out["FR8"] = out.get("THICKNESS_OPM")
        out["FR9"] = out.get("SLAB_NORMAL_DEG")
        out["FR10"] = out.get("SLAB_NORMAL_MEAN")
        out["FR11"] = out.get("SLAB_NORMAL_P90")
        out["FR12"] = (max(slab) if slab else None)

    # ---- the geometric-feature survey
    survey_path = results / "disagreement_survey.json"
    if survey_path.exists():
        survey = json.loads(survey_path.read_text())
        out["_survey"] = survey
        per = survey["per_protein"]
        n = len(per)
        exact = survey.get("n_exactly_curated")
        out["SURVEY_EXACT"] = f"{exact} / {n}"
        out["SURVEY_EXACT_PCT"] = f"{100.0 * exact / n:.1f} %" if n else None
        out["SURVEY_DIFF"] = f"{n - exact} / {n}" if n else None
        f1 = [p["f1"] for p in per if p.get("f1") is not None]
        out["SURVEY_F1"] = float(np.mean(f1)) if f1 else None
        for key, field in (("SURVEY_BROKEN", "broken_helix"),
                           ("SURVEY_UNWOUND", "unwound_crossing"),
                           ("SURVEY_INTERFACIAL", "interfacial_helix"),
                           ("SURVEY_REENTRANT", "reentrant")):
            out[key] = sum(1 for p in per if p.get(field))
    # ---- the two user stories
    story_path = results / "user_story_benchmark.json"
    if story_path.exists():
        story = json.loads(story_path.read_text())
        a = story.get("story_a", [])
        b = story.get("story_b", [])
        if a:
            f1 = [r["f1_vs_curation"] for r in a if r.get("f1_vs_curation") is not None]
            bnd = [r["boundary_median"] for r in a if r.get("boundary_median") is not None]
            times = [r["seconds"] for r in a]
            extra = sum(r.get("n_broken_helix", 0) + r.get("n_interfacial", 0)
                        + r.get("n_reentrant", 0) + r.get("n_unwound", 0) for r in a)
            out["STORY_A_N"] = len(a)
            out["STORY_A_F1"] = round(float(np.mean(f1)), 3) if f1 else None
            out["STORY_A_BND"] = round(float(np.median(bnd)), 1) if bnd else None
            out["STORY_A_T"] = round(float(np.median(times)), 1) if times else None
            out["STORY_A_EXTRA"] = extra
            for key, field in (("STORYA_BROKEN", "n_broken_helix"),
                               ("STORYA_UNWOUND", "n_unwound"),
                               ("STORYA_INTERFACIAL", "n_interfacial"),
                               ("STORYA_REENTRANT", "n_reentrant")):
                out[key] = sum(int(r.get(field, 0)) for r in a)
            out["STORY_A_UNCERTAIN"] = sum(1 for r in a
                                           if r.get("placement_quality") != "ok")
        if b:
            times = [r["seconds"] for r in b]
            out["STORY_B_N"] = len(b)
            out["STORY_B_ANSWERED"] = sum(1 for r in b if r.get("answered"))
            out["STORY_B_T"] = round(float(np.median(times)), 1) if times else None
            out["STORY_B_UNCERTAIN"] = sum(1 for r in b
                                           if r.get("placement_quality") != "ok")

    # ---- the external validation of the vendored TMHMM port
    port_path = results / "tmhmm_port_validation.json"
    if port_path.exists():
        port = json.loads(port_path.read_text())
        model = port.get("model_file", {})
        parsed = port.get("parsed_model", {})
        out["PORT_MODEL_IDENTICAL"] = ("identical" if model.get(
            "identical_after_normalising_line_endings") else "different")
        out["PORT_STATES"] = parsed.get("n_states")
        out["PORT_TRANS_DIFF"] = parsed.get("max_transition_difference")

    return out


def normalise(key: str) -> str:
    """The key of a placeholder, as it is stored in the result files.

    In the LaTeX source the key is verbatim, but a pipe that came from a Markdown
    table is still escaped (``T\\|3d_energy\\|f1``).
    """
    return key.replace(r"\|", "|").replace(r"$\vert$", "|").replace("\\_", "_")


def render(key: str, values: dict) -> str | None:
    """The LaTeX for one placeholder, or None when it is not a number we know."""
    key = normalise(key)
    if key.startswith("T|"):
        # the cells of table 1: the number of decimals depends on the column
        parts = key.split("|")
        if len(parts) == 3 and values.get(key) is not None:
            value = values[key]
            if parts[2] == "nref":
                return fmt(value, decimals=0)
            if parts[2] == "time":
                return fmt(value, decimals=2)
            if parts[2] == "bnd":
                return fmt(value, decimals=2)
            return fmt(value, decimals=3)
    if key in values and values[key] is not None:
        return fmt(values[key])
    if key.startswith("D|"):
        parts = key.split("|")
        if len(parts) == 3:
            name, field = parts[1], parts[2]
            if field == "diff":
                value = values.get(f"{name}_DIFF")
                return f"{value:+.3f}".lstrip("+") if value is not None else None
            if field == "lo":
                value = values.get(f"{name}_LO")
                return f"{value:+.3f}".lstrip("+") if value is not None else None
            if field == "hi":
                value = values.get(f"{name}_HI")
                return f"{value:+.3f}".lstrip("+") if value is not None else None
            if field == "p":
                value = values.get(f"{name}_P")
                if value is None:
                    return None
                return "$< 0.001$" if value < 0.001 else f"${value:.3f}$"
            if field == "ties":
                value = values.get(f"{name}_TIES")
                return str(value) if value is not None else None
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="bench/results")
    ap.add_argument("--template", default="manuscript/ProteinLab_manuscript_submission.tex")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    results = Path(args.results)
    template = Path(args.template)
    out_path = Path(args.out) if args.out else template.with_name(
        template.stem.replace("ProteinLab_manuscript", "ProteinLab_manuscript_filled")
        + ".tex")

    values = collect(results)
    text = template.read_text(encoding="utf-8")
    unresolved: list[str] = []

    def replace(match):
        key = match.group(1) if match.group(1) is not None else match.group(2)
        value = render(key, values)
        if value is None:
            unresolved.append(key)
            return match.group(0)
        return value

    text = PLACEHOLDER.sub(replace, text)
    out_path.write_text(text, encoding="utf-8")
    print(f"written to {out_path}")
    if unresolved:
        # the reference list keeps its markers on purpose: every reference that still
        # has to be checked stays visible
        print(f"\n{len(unresolved)} placeholder(s) left open (author metadata, the "
              f"declarations and the references that still have to be checked):")
        for key in sorted(set(unresolved)):
            print("  ", key)
    else:
        print("every placeholder was filled from the result files")


if __name__ == "__main__":
    main()
