import json
from pathlib import Path
from app.api.topology_predictor import predict_topology

pdb_file = Path(r"D:\Protein-Visualize\tests\hNET_8WTW_ions.pdb")
try:
    print("Testing DSSP Slab Algorithm...")
    result = predict_topology(pdb_file, "dssp_slab")
    for r in result.get('regions', []):
        if r['type'] == 'Transmembrane':
            print(f"Transmembrane: {r['start']}-{r['end']}")
    
    print("\nTesting KD Slab Algorithm...")
    result_kd = predict_topology(pdb_file, "kd_slab")
    for r in result_kd.get('regions', []):
        if r['type'] == 'Transmembrane':
            print(f"Transmembrane: {r['start']}-{r['end']}")
except Exception as e:
    print(f"Error: {e}")
