# ProteinLab

**Structure-derived transmembrane topology annotation and interactive visualization of
membrane proteins.**

ProteinLab derives a residue-level, side-resolved transmembrane topology directly from
atomic coordinates. It places a lipid bilayer on a structure with an implicit-membrane
transfer-energy model, merges the resulting topology with an independently computed
secondary-structure assignment, resolves the geometric cases that sequence-based
methods cannot see — hairpins, unwound helices, interfacial helices, re-entrant
segments — and presents the result as an interactive, publication-quality figure with
machine-readable export.

- **Web application:** <https://protein-visualize.vercel.app>
- **API:** <https://protein-visualize.onrender.com> (OpenAPI at `/docs`)
- **Source:** <https://github.com/soratori1025/Protein-Visualize>

---

## What it does

1. **Membrane placement.** The bilayer is placed by minimising an implicit transfer
   free energy built from solvent-accessible surface area, per-residue lipid-transfer
   scales (Goldman–Engelman–Steitz for the hydrocarbon core, Wimley–White for the
   head-group interface) and a correction for residues that line an aqueous pore. The
   global search over the membrane normal, the mid-plane offset and the hydrophobic
   thickness is evaluated for all three at once by a one-dimensional fast Fourier
   transform. A curved (spherical) mid-surface is fitted when it is clearly preferred.
2. **Topology annotation.** The topology block and the secondary-structure block are
   merged on a single shared residue frame, and crossings, sides, hairpins, kinks,
   unwound helices, interfacial helices and membrane transitions are resolved by
   explicitly documented geometric rules.
3. **Figures.** The topology is drawn as an interactive SVG figure and can be exported
   as a raster image, as a `.3line` topology string, or as GFF3. Large assemblies can
   additionally be re-rendered as a flat topology map in which the transmembrane bodies
   keep their native conformation while the connecting loops are rebuilt kinematically —
   with the geometric distortion that the layout introduces measured and displayed.
4. **Interoperability.** A documented JSON API, a command-line interface, and machine
   -readable exports in two standard formats.

## Install

```bash
python -m venv .venv
./.venv/Scripts/python.exe -m pip install -r backend/requirements.txt   # Windows
./.venv/bin/pip install -r backend/requirements.txt                       # Linux/macOS

# as a package, with the command-line interface
./.venv/Scripts/python.exe -m pip install ./backend
```

### Command line

```bash
proteinlab topology 5i6x.pdb --algorithm 3d_energy --chain A
proteinlab topology 5i6x.pdb --algorithm 3d_energy --json result.json
proteinlab secondary-structure 5i6x.pdb --method dssp
proteinlab capabilities
```

### Backend

```bash
cd backend
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

### Frontend

```bash
cd frontend
npm install
npm run dev        # http://localhost:5173, talks to http://127.0.0.1:8000
```

The frontend reads the API location from `VITE_API_URL` (default
`http://127.0.0.1:8000`).

### External programs

Secondary structure is computed by **DSSP** and **STRIDE**, which are not bundled with
the Python packages:

| system | how to get it |
|---|---|
| Linux | `apt-get install dssp` (DSSP 4.x) and `conda install -c bioconda stride` |
| Windows | the executables in `tools/bin/` (`mkdssp.exe`, `stride.exe`) are added to `PATH` automatically by `app.main` |

`proteinlab capabilities` reports which of them the installation can actually use.
The DSSP adapter calls the program on the coordinate file itself, so the mmCIF path
works with the mmCIF-only builds of DSSP 4.4 and later; an input that the program
cannot read raises an error rather than returning an empty assignment.

## Tests

```bash
./.venv/Scripts/python.exe -m pytest backend/tests -q
```

The suite is the numerical contract of the method: the solvent-accessible surface is
checked against its closed form, the FFT search against a direct evaluation of the
objective, and the residue frame, the numbering-verified alignment, the consensus rules
and the full flow against their invariants.

## Benchmarks

Every accuracy number quoted in the publication is produced by a committed script and
committed with its output:

```bash
./.venv/Scripts/python.exe bench/benchmark_tm_topology.py --out bench/results --limit 150 --soluble --workers 8
./.venv/Scripts/python.exe bench/compare_methods.py --results bench/results
./.venv/Scripts/python.exe bench/membrane_frame_benchmark.py --results bench/results
./.venv/Scripts/python.exe bench/unroll_benchmark.py --results bench/results --limit 40
./.venv/Scripts/python.exe bench/make_figures.py --results bench/results --out bench/figures
```

| script | question |
|---|---|
| `benchmark_tm_topology.py` | does the computed topology reproduce curated UniProtKB topology, and does the transfer free energy separate membrane from soluble proteins? |
| `compare_methods.py` | how do the providers compare, paired on the same proteins? |
| `membrane_frame_benchmark.py` | how accurately is the bilayer placed, against OPM/PPM reference frames? |
| `unroll_benchmark.py` | what does the two-dimensional layout do to the structure? |
| `case_hsert.py`, `case_alphafold.py`, `case_piezo1.py`, `case_plddt.py` | the case studies |

See `bench/README.md` for the protocol and for how to read the results honestly.

## Repository layout

```
backend/          the FastAPI service, the topology engine, the unit tests
protein_engine/   the shared structure-parsing and secondary-structure library
frontend/         the React + TypeScript single-page application
bench/            the validation harness, its results and the figures
docs/             the publication analysis and the strengthening plan
manuscript/       the manuscript, the placeholder register, the submission checklist
```

## Citation

See `CITATION.cff`. The paper is in preparation; please cite the repository until the
DOI is assigned.

## Manuscript

The manuscript and the supplementary information are LaTeX and are kept **outside the
repository** (`manuscript/` is in `.gitignore`): several journals require that a
submitted paper is not publicly available before publication. See `manuscript/README.md`
for the compilation instructions and for how every number in the paper is filled
automatically from the committed benchmark results.

The paper is being written for ***BMC Bioinformatics*** as a **Software article**, and
the source archive that the journal asks to be attached to the submission is built
with:

```bash
python bench/make_source_archive.py --out manuscript/ProteinLab-source-archive.zip
```

The archive contains the code, the tests, the benchmark harness, the result files and
the figures — and no part of the paper.

## Licence

MIT, with the exceptions listed in `LICENSE`: the vendored TMHMM 2.0 model file and the
Windows builds of DSSP and STRIDE keep their own terms.
