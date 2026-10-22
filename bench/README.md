# Benchmarks — structure-derived topology and membrane placement

This directory is the validation harness behind the manuscript. Every number quoted
in the paper is produced by a script here and stored under `results/`; nothing is
typed by hand. The scripts are deterministic: the cohorts are drawn with fixed
seeds, the thresholds are constants in the scripts, and the UniProtKB release is
recorded with every run.

## The scripts

| Script | Purpose |
|---|---|
| `benchmark_tm_topology.py` | The main harness: cohort selection, prediction, scoring, soluble control |
| `compare_methods.py` | Paired Wilcoxon signed-rank tests and bootstrap confidence intervals between providers |
| `frame_cohort.py` | Selects the membrane-frame cohort (structures with an OPM/PPM reference) |
| `membrane_frame_benchmark.py` | Fits the bilayer to each frame-cohort structure and compares the normal, offset and thickness with OPM/PPM |
| `disagreement_survey.py` | Where the providers disagree, and whether the structure resolves the disagreement |
| `unroll_benchmark.py` | The two-dimensional unrolling audit: loop closure, path deviation, clashes, bonds, disulfides |
| `user_story_benchmark.py` | The two end-to-end user stories, timed through the public service |
| `curvature_survey.py` | How often the spherical mid-surface is accepted over the 150-protein cohort |
| `make_figures.py` | Article Figures 1--4 and the depth profile of Figure 1C |
| `validate_tmhmm_port.py` | Checks the bundled TMHMM port against the upstream implementation |
| `competitive_matrix.py` | Reads the competitors' repositories and writes the comparison matrix (Table 7) |
| `case_hsert.py` | Case study: computed vs curated topology for the serotonin transporter 5I6X |
| `case_alphafold.py` | Case study: an AlphaFold DB model, with and without the pLDDT filter |
| `case_plddt.py` | The pLDDT-filter effect on a named PDB structure |
| `case_piezo1.py` | Case study: PIEZO1 8YEZ and 8ZU3, plane vs sphere, and the hSERT control |
| `make_source_archive.py` | Builds the source archive that accompanies the article (Additional file 5) |

## What the benchmarks measure

| Question | Result file | Where it appears |
|---|---|---|
| Does the computed topology reproduce curated topology? | `benchmark.json`, `summary.csv` | Table 1 |
| Are the differences between providers significant? | `paired_comparisons.json` | Table 1 and the text |
| Is the fitted bilayer in the right place? | `membrane_frame_benchmark.json` | Table 3, Figure 2 |
| Is a soluble protein mistaken for a membrane protein? | `soluble_control.json`, `soluble_control_summary.json` | Table 2 |
| What does the two-dimensional unrolling cost? | `unroll_benchmark.json` | Table 6, Figure 4 |
| Do the two user stories work end to end? | `user_story_benchmark.json` | Table 5, Figure 3 |
| When does the curved mid-surface engage? | `curvature_survey.json` | Results, curvature subsection |
| What are the PIEZO1 plane and sphere energies? | `piezo1_case.json` | Results, curvature subsection |
| Does the TMHMM port agree with upstream? | `tmhmm_port_validation.json` | Results, port check |
| What can each competitor actually do? | `competitive_matrix.json` | Table 7 |
| Where do the providers disagree? | `disagreement_survey.json` | Results, disagreement subsection |
| What does the depth profile of one structure look like? | `figure1_depth_profile.json` | Figure 1C |

## Cohorts and determinism

| Benchmark | Cohort | Seed / release |
|---|---|---|
| Topology, user stories, curvature | 150 human Swiss-Prot entries with curated `Transmembrane` features and a PDB cross-reference | `SEED = 20261021`; UniProtKB 2026\_03 |
| Membrane frame | 211 structures with an OPM/PPM reference bilayer | frame seed 20261022 |
| Unrolling | 16 structures with at least three crossings | seeded `mulberry32` layout |
| Soluble control | 13 soluble Swiss-Prot entries | same release |

The exact UniProt query, the thresholds, the timestamps and the release are written
to `results/meta.json` with every run. The frame cohort is written to
`frame/frame_cohort.json`.

## Running it

The commands below are the ones in Additional file 2, S8.1. On Windows, replace
`./.venv/bin/python` with `.\.venv\Scripts\python.exe`.

```bash
python -m venv .venv
./.venv/bin/pip install -r backend/requirements.txt
./.venv/bin/pip install matplotlib scipy pytest
./.venv/bin/pip install ./backend          # the `proteinlab` command

./.venv/bin/python -m pytest backend/tests -q

./.venv/bin/python bench/benchmark_tm_topology.py \
    --out bench/results --limit 150 --soluble --workers 8
./.venv/bin/python bench/compare_methods.py --results bench/results
./.venv/bin/python bench/frame_cohort.py --out bench/frame --limit 220
./.venv/bin/python bench/membrane_frame_benchmark.py \
    --results bench/results --cohort bench/frame/frame_cohort.json \
    --out bench/results
./.venv/bin/python bench/disagreement_survey.py --results bench/results
./.venv/bin/python bench/unroll_benchmark.py --results bench/results --limit 16
./.venv/bin/python bench/make_figures.py \
    --results bench/results --out bench/figures
./.venv/bin/python bench/validate_tmhmm_port.py --upstream <pyTMHMM checkout>
./.venv/bin/python bench/competitive_matrix.py --clones <path to the clones>
./.venv/bin/python bench/user_story_benchmark.py \
    --results bench/results --story-a 150 --story-b 30
./.venv/bin/python bench/curvature_survey.py --out bench/results --workers 8
./.venv/bin/python bench/case_piezo1.py
./.venv/bin/python bench/make_source_archive.py \
    --out Additional_file_5_source_archive.zip
```

The TMHMM comparison is the only step that needs a file the release does not ship:
the TMHMM 2.0 model is not redistributed. Without it the provider answers HTTP 503,
every other provider runs, and the reported TMHMM predictions remain in
`bench/results/` for inspection.

## Outputs

| File | Content |
|---|---|
| `benchmark.json` | one record per protein: accession, PDB ID, chain, name, keywords, reference helix count, and the full metric set for each method |
| `summary.csv` / `summary.json` | the pooled table — the columns of Table 1 |
| `paired_comparisons.json` | the paired Wilcoxon tests and bootstrap intervals |
| `membrane_frame_benchmark.json` | the fitted frame per structure: normal error, offset, thickness, ΔG, curvature |
| `soluble_control.json` | per-protein transfer free energy of the soluble controls |
| `soluble_control_summary.json` | the pooled soluble statistics — Table 2 |
| `unroll_benchmark.json` | per-structure unrolling audit: crossings, loops, closure, clashes, bonds, disulfides |
| `user_story_benchmark.json` | the two user stories, per protein, with timings and exports |
| `curvature_survey.json` | per-protein plane/sphere decision over the 150-protein cohort |
| `piezo1_case.json` | PIEZO1 plane and sphere energies, radii and topology counts |
| `tmhmm_port_validation.json` | the port check against upstream |
| `disagreement_survey.json` | the disagreement survey |
| `figure1_depth_profile.json` | the per-residue depth profile drawn in Figure 1C |
| `competitive_matrix.json` | the competitor matrix behind Table 7 |
| `meta.json` | query, seed, thresholds, timestamps, match parameters, UniProtKB release |
| `structures/` | downloaded coordinate files (not committed; re-downloaded) |
| `alphafold/` | AlphaFold DB models used in the case studies (not committed) |

## Reading the results honestly

- The `uniprot_api` row is **not** an independent method: it paints the reference
  annotation onto the structure by the same code path that builds the reference set.
  It scores close to 1 and would score exactly 1 if the rules were identical. It is a
  consistency check of the numbering machinery, and it failed for 11 of the 150
  proteins (8 sequence mismatches, 2 missing accessions, 1 missing sequence) —
  failures of the reference, not of the prediction.
- The cohort contains **curated** proteins. Hydropathy- and profile-based predictors
  were trained on curated annotations of the same kind, so the comparison flatters
  the sequence methods relative to what a user meets on an uncurated protein. The
  value of the structure-based providers lies in the uncurated and predicted-model
  cases, which this benchmark does not cover.
- `precision` below 1 reflects **over-segmentation**: a provider that splits one
  curated helix into two predicted segments loses precision while keeping recall.
  Read `recall` and `boundary` together.
- **The two structure-based providers are statistically indistinguishable** on
  topology (paired Wilcoxon *P* = 0.10, 138 of 150 proteins tied). Do not report the
  energy model as "more accurate" than the slab fit. Its value is elsewhere: it
  produces a free energy, a fitted thickness and a curved mid-surface, and on the
  frame benchmark it is much closer to the OPM/PPM normal than the geometric slab
  control (median 3.25° vs 11.94°; 174 of 211 structures closer).
- The **soluble control is not a large sample**. Thirteen proteins show that the
  separation is clean, but they do not characterise the failure rate of the
  membrane-versus-soluble decision.
- The **unrolling is a rendering capability with a measured limitation**, not a
  validated reconstruction: 27 of 101 loops cannot be closed and are drawn
  schematically, and the layout introduces steric clashes whose median is 60.5
  heavy-atom pairs per structure.
- The **curved mid-surface is exercised by 7 of 150 structures**, at radii of
  60--142 Å. It is accepted only when it beats the plane by a wide margin, so it is
  conservative by construction; its positive behaviour rests on a small number of
  cases.

## What is not measured

The benchmarks are human, curated and small enough to audit. They say nothing
directly about bacterial or archaeal membranes, about proteome-scale uncurated
proteins, or about structures whose bilayer is not approximately planar or spherical.
The AlphaFold models in `case_alphafold.py` are case studies, not a benchmark.

## Rebuilding the source archive

The archive is built from the repository, not from `bench/` alone:

```bash
./.venv/bin/python bench/make_source_archive.py \
    --out Additional_file_5_source_archive.zip
```

It contains the backend, the bench harness, the tests and the committed result
files, and it excludes the manuscript, the supplementary information and the
downloaded coordinate files.
