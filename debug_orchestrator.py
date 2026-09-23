import sys
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path("d:/Protein-Visualize/backend")))

from app.api.dependencies import get_topology_orchestrator
from app.schemas.topology import TMParams

print("Loading orchestrator...")
orchestrator = get_topology_orchestrator(tm_algo="3d_slab_geom", ss_algo="dssp", flow_type="ss_then_tm")
params = TMParams()  # Use actual default values, not FastAPI Query objects
pdb_path = Path("d:/Protein-Visualize/backend/data/uploads/hSERT_5I6X_ions.pdb")

print("Executing orchestrator...")
try:
    response = orchestrator.execute(pdb_path, params=params)
    print("Success!")
    print(response.dict())
except Exception as e:
    import traceback
    traceback.print_exc()
