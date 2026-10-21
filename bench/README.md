# Benchmark — transmembrane topology against curated UniProtKB annotation

`benchmark_tm_topology.py` is the validation harness for the ProteinLab topology
engine. It is the script that produced Tables 1 and 2 of the manuscript, and it is
committed so that the numbers in the paper can be regenerated.

## The scripts

| Script | Purpose |
|---|---|
| `benchmark_tm_topology.py` | The main harness: cohort selection, prediction, scoring, soluble control |
| `compare_methods.py` | Paired Wilcoxon signed-rank tests and bootstrap confidence intervals between providers |
| `membrane_frame.py` | Fitted membrane frame for every protein in the cohort: transfer free energy, hydrophobic thickness, curvature, chain count |
| `make_figures.py` | Figure 3 and Supplementary Figure S1 |
| `case_hsert.py` | Case study: computed vs curated topology for the serotonin transporter 5I6X, residue by residue |
| `case_alphafold.py` | Case study: an AlphaFold DB model, with and without the pLDDT filter |
| `case_plddt.py` | The pLDDT-filter effect on a named PDB structure |
| `case_piezo1.py` | Case study: the PIEZO1 structures 8YEZ and 8ZU3 |

## What it measures

| Question | Answer |
|---|---|
| Does the computed topology reproduce curated topology? | recall / precision / F1 against UniProtKB `Transmembrane` features, per protein and pooled |
| How far off are the boundaries? | mean and median boundary error in residues; fraction within 5 residues |
| Is a soluble protein mistaken for a membrane protein? | transfer free energy of a fixed soluble control set against the −10 kcal/mol threshold |
| Does the curved model engage where it should? | the `curved R=` field of the labeler string and the curvature warning |

## Protocol

**Cohort.** Swiss-Prot (reviewed UniProtKB) entries with curated `Transmembrane`
features and a PDB cross-reference:

```
(reviewed:true) AND (ft_transmem:*) AND (organism_id:9606)
```

The sample is drawn with a fixed seed (`SEED = 20261021`) after a deterministic
shuffle, so the cohort is reproducible. The UniProtKB release, the query, the seed,
the thresholds and the timestamps are written to `meta.json` with every run.

**Ground truth.** The curated `Transmembrane` features of the UniProtKB entry, mapped
onto the analysed chain by `app.services.topology.residues.map_to_reference` — the same
numbering-verified alignment the application uses at run time. A curated helix enters
the reference set only if at least `MIN_RESOLVED_FRACTION = 0.60` of its residues are
resolved in the chain, so the benchmark measures prediction error rather than
crystallographic incompleteness.

**Methods.** `3d_energy`, `3d_slab_geom`, `kyte_doolittle_seq`, `tmhmm` and
`uniprot_api`. The last one paints the curated annotation onto the same structure and
is a consistency check of the numbering machinery, not a competing predictor.

**Matching.** Greedy, one-to-one, by overlap: a predicted segment and a reference
segment are matched when their overlap covers at least `MATCH_IOU = 0.50` of the
smaller one.

**Control cohort.** `SOLUBLE_ACCESSIONS` — twelve soluble Swiss-Prot entries
(lysozyme, polyubiquitin, myoglobin, cytochrome *c*, GFP, fibronectin, transthyretin,
catalase, adenylate kinase, crambin, trypsin). For each, the PDB entry and chain are
chosen by the script (best resolution); the structure is called a membrane protein
when the transfer free energy is at or below −10 kcal/mol.

## Running it

```powershell
cd D:\ngan-work\protein-visualize
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r backend\requirements.txt matplotlib scipy
.\.venv\Scripts\python.exe bench\benchmark_tm_topology.py --out bench\results --limit 150 --soluble --workers 8
.\.venv\Scripts\python.exe bench\compare_methods.py --results bench\results
.\.venv\Scripts\python.exe bench\membrane_frame.py --results bench\results --workers 8
.\.venv\Scripts\python.exe bench\make_figures.py --results bench\results --out bench\figures
```

Useful flags:

| Flag | Meaning |
|---|---|
| `--limit N` | number of proteins to analyse (default 60) |
| `--query Q` | UniProt query for the cohort |
| `--soluble` | also run the soluble control cohort |
| `--workers N` | number of parallel worker processes (the predictors are CPU-bound; 8 workers cut a 150-protein run from ~90 min to ~8 min on 24 logical cores) |
| `--out DIR` | output directory (default `bench/results`) |

The run committed in `bench/results/` used `--limit 150 --soluble --workers 8` and
took about eight minutes on 24 logical cores.

## Outputs

| File | Content |
|---|---|
| `benchmark.json` | one record per protein: accession, PDB ID, chain, protein name, keywords, reference helix count, and the full metric set for each method |
| `summary.csv` / `summary.json` | the pooled table — the columns of manuscript Table 1 |
| `paired_comparisons.json` | the paired Wilcoxon tests and bootstrap intervals |
| `membrane_frame.json` | the fitted membrane frame per protein: ΔG, thickness, curvature radius, chain count |
| `soluble_control.json` | per-protein transfer free energy of the soluble controls |
| `soluble_control_summary.json` | the pooled soluble statistics — manuscript Table 2 |
| `meta.json` | query, seed, thresholds, timestamps, match parameters, UniProtKB release |
| `structures/` | the downloaded coordinate files (not committed; they are re-downloaded) |
| `alphafold/` | the AlphaFold DB models used in the case studies (not committed) |

## Reading the results honestly

- The `uniprot_api` row is **not** an independent method: it is the reference
  annotation painted on the structure by the same code path that builds the reference
  set, with a slightly different inclusion rule. It scores close to 1 and would score
  exactly 1 if the rules were identical. It is reported as a consistency check, and it
  failed for 11 of the 150 proteins (8 sequence mismatches, 2 missing accessions, 1
  missing sequence) — failures of the reference, not of the prediction.
- The cohort contains **curated** proteins. Hydropathy- and profile-based predictors
  were trained on curated annotations of the same kind, so the comparison is
  favourable to the sequence methods relative to what a user would meet on an
  uncurated protein. The value of the structure-based providers lies precisely in the
  uncurated and predicted-model cases, which this benchmark does not cover — see the
  case study in §3.5 of the manuscript.
- `precision` below 1 reflects **over-segmentation**: a provider that splits one
  curated helix into two predicted segments loses precision while keeping recall. Read
  `recall` and `boundary` together.
- **The two structure-based providers are statistically indistinguishable** on this
  cohort (paired Wilcoxon *P* = 0.10, 138 of 150 proteins tied). Do not report the
  energy model as "more accurate" than the slab fit. Its value is elsewhere: it
  produces a free energy, a fitted thickness and a curved mid-surface.
- The **soluble control is not a large sample**. Thirteen proteins is enough to show
  that the separation is clean, but not enough to characterise the failure rate of the
  membrane-versus-soluble decision; a larger control set, including elongated
  soluble α-helical proteins, is the obvious next step.

## A second benchmark worth writing

The current harness measures **topology**, not **membrane frame**. A companion script
should compare the fitted membrane frame against OPM/PPM for structures that carry a
reference bilayer: the angular error of the normal, the offset of the mid-plane in Å,
and the hydrophobic thickness in Å. The numbers claimed in the source-code docstring
of `membrane_energy.py` are not reproducible and must not appear in the paper until
such a script exists.
