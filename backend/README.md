# Protein-Visualize Backend

This is the FastAPI backend for the Protein-Visualize application. It provides powerful structural biology tools for computing Secondary Structure (SS) assignments and predicting Transmembrane (TM) topologies using both geometric and physics-based models.

## 🏗️ Architecture

The backend follows the **Strategy Pattern** and a robust **Orchestrator** model. It decouples the specific prediction algorithms (Providers) from the consensus builder (Orchestrator).

```text
backend/app/
├── api/                  # Controllers: Route definitions
├── core/                 # Config & Constants: Hydropathy scales (GES, Kyte-Doolittle)
├── schemas/              # Data Models: Pydantic DTOs for request/response
└── services/
    └── topology/
        ├── service.py             # Factory registering TM & SS Algorithms
        ├── topology_predictor.py  # TopologyOrchestrator: Merges SS and TM data
        ├── membrane_energy.py     # Core physics engine for the 3D Energy model (SASA, FFT)
        └── providers/             # Strategy implementations
            ├── ss/                # Secondary Structure (DSSP, STRIDE)
            └── tm/                # Transmembrane Predictors (Sequence, Geometry, Energy, UniProt)
```

## 🧬 Transmembrane (TM) Algorithms

The application supports multiple algorithms for identifying transmembrane regions, allowing users to choose the most suitable method for their data (from sequence-only to full-atom physics).

### 1. 3D Energy (Implicit Membrane Model)
This is our most advanced, **completely self-implemented** biophysics algorithm. It determines the globally optimal membrane placement by minimizing the transfer free energy of the protein into a lipid bilayer.
* **SASA Calculation:** Implements the Shrake-Rupley algorithm via vectorised NumPy operations to estimate the Solvent Accessible Surface Area for each atom.
* **Energy Scale:** Computes implicit transfer free energy using the Goldman-Engelman-Steitz (GES) hydrophobicity scale and a Wimley-White interface penalty.
* **FFT Optimization:** It samples rotational states using a Fibonacci sphere and utilizes 1D Fast Fourier Transform (`np.fft.rfft`) for a lightning-fast exhaustive grid search along the Z-axis (translation and thickness).

### 2. 3D Slab Geometry
A fast geometric heuristic that identifies the membrane by clustering hydrophobic residues in 3D space.
* **Fibonacci Sphere Sampling:** Tests hundreds of possible membrane normal vectors.
* **Slab Scoring:** Scores each orientation based on the density of hydrophobic C-alpha atoms inside the defined slab minus the penalty for charged residues.
* **Consensus Mapping:** After finding the optimal plane, it slices the sequence and ensures runs of residues spanning across the slab are classified as TM helices/strands.

### 3. Kyte-Doolittle (Sequence-based)
A classic 1D sliding-window hydropathy scan. 
* Operates solely on the primary amino acid sequence.
* Uses a **hysteresis** thresholding system (upper and lower bounds) to cleanly separate hydrophobic core segments without being interrupted by single hydrophilic mutations.

### 4. UniProt API (Database Lookup)
Fetches curated annotations directly from the UniProt database using the protein's accession ID. Reliable for established proteins, returning exact boundaries for Extracellular, Transmembrane, and Cytoplasmic domains.

## 🛠️ How to Add a New Algorithm

The architecture relies heavily on **Providers**. If you want to add a new algorithm (e.g., an ML-based predictor), you do not need to rewrite the orchestration logic.

1. **Create a Provider:** 
   In `app/services/topology/providers/tm/`, create a new class implementing `TMProvider`.
   ```python
   from app.services.topology.providers.tm.base import TMProvider
   
   class MyMLProvider(TMProvider):
       def extract_tm_segments(self, file_path, params):
           # Your AI/ML logic here
           return [ ... list of segments ... ]
   ```

2. **Register it:** 
   Open `app/services/topology/service.py` and map your new algorithm key in the `TM_ALGORITHMS` dictionary.
   ```python
   TM_ALGORITHMS = {
       "kyte_doolittle_seq": SequenceTMProvider,
       "3d_slab_geom": _geometry_provider,
       "3d_energy": _energy_provider,
       "my_ml_algo": MyMLProvider,  # <--- Registered here
   }
   ```

3. **Call from Frontend:** 
   The frontend can now pass `algorithm=my_ml_algo` in the API request, and the `TopologyOrchestrator` will automatically route the request, merge the TM predictions with the chosen SS provider (DSSP/STRIDE), and return the standard JSON output.

## 🚀 Development Setup

To run the backend locally:

1. Create a virtual environment and install dependencies:
   ```bash
   python -m venv .venv
   .venv\Scripts\activate
   pip install -r requirements.txt
   ```
2. Start the FastAPI server:
   ```bash
   uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
   ```
3. Open `http://localhost:8000/docs` to view the interactive Swagger UI.
