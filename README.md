# Protein-Visualize

ProteinLab MVP for protein structure visualization and analysis.

## Current capabilities

- Upload PDB and mmCIF files from the React frontend.
- Parse chains, residues, sequences, and atom coordinates with Biopython.
- Inspect a lightweight SVG 3D backbone projection and 2D sequence/topology views.
- View the complete biological assembly as a colored 3D cartoon overview.
- Explore a separate interactive 3D viewer below the overview with residue selection.
- Expand chains/compounds into a horizontal 2D map with internal strips and external assembly connectors.
- Synchronize residue selection between the 3D projection, sequence, and topology map.
- Select DSSP, STRIDE, comparison, or manual annotation modes.
- Run DSSP/STRIDE through normalized adapters when `mkdssp` or `stride` is installed on PATH.
- Visualize DSSP/STRIDE output as a residue-level colored assignment track with semantic legends.
- Run a chain research toolbox with composition, molecular weight, pI, aromaticity, instability, hydropathy, CA geometry, and contact pairs.

## Run locally

Frontend:

```powershell
cd frontend
npm install
npm run dev
```

Backend:

```powershell
..\.venv\Scripts\python.exe -m uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

API endpoints:

- `GET /api/health`
- `POST /api/structure/upload`
- `POST /api/secondary-structure/{filename}?method=DSSP|STRIDE`
- `GET /api/analysis/{filename}/chain/{chain_id}?contact_cutoff=8`

### DSSP and STRIDE prerequisites

The adapters call the native executables `mkdssp` and `stride`; they are not bundled by pip. The app exposes `GET /api/secondary-structure/capabilities` and labels each method as `ready` or `not installed`.

On Windows:

1. Build or install `mkdssp` from [PDB-REDO/dssp](https://github.com/PDB-REDO/dssp).
2. Build or install `stride` from [MDAnalysis/stride](https://github.com/MDAnalysis/stride).
3. Add the folders containing `mkdssp.exe` and `stride.exe` to `PATH`.
4. Restart the backend and click `Check API` in ProteinLab.

DSSP receives PDB/mmCIF directly. STRIDE receives PDB; when the user uploads mmCIF, ProteinLab converts the selected model to a temporary PDB before invoking STRIDE. The temporary file is deleted after the run.

The annotation panel reports executable readiness, native errors, assignment counts, and a residue-level color track. Until the tools are installed, the API returns an explicit `503` instead of silently returning an empty result.

The analysis response records the input filename, chain, algorithm name, version, and parameters so results can be reproduced later. Current calculations use Biopython's `ProteinAnalysis` plus coordinate-based C-alpha contact geometry; SciPy/NetworkX remain available for the next contact graph and clustering modules.

The interactive viewer uses 3Dmol.js. The overview keeps the full assembly together, while the spread map is a 2D architecture view intended to make large multi-chain structures easier to read.