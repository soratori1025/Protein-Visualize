# Protein-Visualize Backend

This is the FastAPI backend for the Protein-Visualize application. It is responsible for parsing protein structures, calculating secondary structures (using DSSP/STRIDE), and predicting transmembrane (TM) topologies.

## 🏗 Architecture

The backend follows **SOLID** and **Object-Oriented Design (OOD)** principles, separating concerns into distinct layers:

```text
backend/app/
├── api/                  # Controllers: Route definitions and Dependency Injection
├── core/                 # Config & Constants: Hydropathy scales, default parameters
├── schemas/              # Data Models: Pydantic schemas (DTOs) for request/response validation
└── services/             # Business Logic: Topology predictors, core algorithms
    └── topology/
        ├── base.py                 # Abstract Base Class for predictors
        ├── sequence_predictor.py   # Kyte-Doolittle TM scan + DSSP rescue (Sequence-based)
        └── structure_predictor.py  # Slab-geometry + SS-element parsing (Structure-based)
```

## 🚀 How to Extend and Add Your Own Algorithm

The application uses the **Strategy Pattern** to swap out topology prediction algorithms dynamically. If you want to add your own algorithm (for instance, an ML-based predictor or AlphaFold integration), you do not need to modify existing complex code. You simply extend the `BaseTopologyPredictor`.

### Step 1: Create Your Predictor Class

Create a new file in `app/services/topology/` (e.g., `my_custom_predictor.py`) and inherit from `BaseTopologyPredictor`.

```python
from pathlib import Path
from typing import Optional
from app.services.topology.base import BaseTopologyPredictor
from app.schemas.topology import TMParams, TopologyResponse, TopologyRegion

class MyCustomPredictor(BaseTopologyPredictor):
    def predict(self, file_path: Path, labeler: str = "DSSP", params: Optional[TMParams] = None, **kwargs) -> TopologyResponse:
        # 1. Parse the structure at file_path
        # 2. Run your custom algorithm
        
        # Example dummy regions
        regions = [
            TopologyRegion(type="Topological domain", start=1, end=20, description="Extracellular"),
            TopologyRegion(type="Transmembrane", start=21, end=41, description="Transmembrane Alpha Helix"),
            TopologyRegion(type="Topological domain", start=42, end=60, description="Cytoplasmic"),
        ]
        
        # 3. Return the standard Pydantic response
        return TopologyResponse(
            uniprot_id="CUSTOM",
            protein_name="Custom Prediction",
            gene_name="",
            organism="Computed",
            membrane_score=1.0,
            membrane_normal=[0.0, 0.0, 1.0],  # Optional
            labeler=labeler,
            parameters_used=params.to_response_dict() if params else {},
            regions=regions
        )
```

### Step 2: Register Your Predictor

Open `app/api/dependencies.py`. Modify the `get_topology_predictor` function to inject your custom class when a specific `algorithm` string is provided by the frontend.

```python
def get_topology_predictor(algorithm: str = "dssp_slab") -> BaseTopologyPredictor:
    if algorithm == "my_custom_algo":
        from app.services.topology.my_custom_predictor import MyCustomPredictor
        return MyCustomPredictor()
    
    # ... existing logic ...
```

### Step 3: Call Your Algorithm from the Frontend

Your new algorithm is now ready to use! In the frontend, just make an API call to `/api/secondary-structure/predict-topology/YOUR_FILE.pdb?algorithm=my_custom_algo`.

## ⚙️ Development Setup

To run the backend locally:

1. Install requirements:
   ```bash
   pip install -r requirements.txt
   ```
2. Start the FastAPI server (using Uvicorn):
   ```bash
   npm run dev 
   # or natively: uvicorn app.main:app --reload --port 8000
   ```
3. Open `http://localhost:8000/docs` to view the Swagger UI and test the endpoints.
