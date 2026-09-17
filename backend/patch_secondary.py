import sys

with open('app/api/secondary_structure.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Add import at the top
if 'from .tm_sequence_predictor import predict_topology_from_sequence' not in content:
    content = content.replace('from .topology_predictor import predict_topology', 'from .topology_predictor import predict_topology\nfrom .tm_sequence_predictor import predict_topology_from_sequence')

# Replace the endpoint
old_endpoint = '''@router.get("/predict-topology/{filename}")
def predict_topology_from_structure(filename: str) -> dict:
    path = UPLOAD_DIR / Path(filename).name
    if not path.is_file():
        raise HTTPException(status_code=404, detail="Uploaded structure not found")
    
    try:
        return predict_topology(path)
    except Exception as error:
        raise HTTPException(status_code=500, detail=str(error)) from error'''

new_endpoint = '''@router.get("/predict-topology/{filename}")
def predict_topology_from_structure(filename: str, algorithm: str = "tmhmm_seq") -> dict:
    path = UPLOAD_DIR / Path(filename).name
    if not path.is_file():
        raise HTTPException(status_code=404, detail="Uploaded structure not found")
    
    try:
        if algorithm == "tmhmm_seq":
            return predict_topology_from_sequence(path)
        return predict_topology(path, algorithm)
    except Exception as error:
        raise HTTPException(status_code=500, detail=str(error)) from error'''

if old_endpoint in content:
    content = content.replace(old_endpoint, new_endpoint)

with open('app/api/secondary_structure.py', 'w', encoding='utf-8') as f:
    f.write(content)
print("Updated secondary_structure.py successfully")
