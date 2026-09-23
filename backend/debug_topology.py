from pathlib import Path
from Bio.PDB import MMCIFParser, PDBParser
from app.services.topology.orchestrator import TopologyOrchestrator
from app.services.topology.providers.tm.geometry import GeometryTMProvider
from app.services.topology.providers.ss.dssp import DSSPProvider
from app.schemas.topology import TMParams

upload_dir = Path("data/uploads")
files = list(upload_dir.glob("*.pdb")) + list(upload_dir.glob("*.cif"))
if not files:
    print("No files found")
else:
    file = files[0]
    print(f"Testing {file}")
    
    parser = MMCIFParser() if file.suffix.lower() in ('.cif', '.mmcif') else PDBParser()
    structure = parser.get_structure('protein', str(file))
    chain = next(iter(next(iter(structure))))
    
    res_ids = [res.get_id()[1] for res in chain if 'CA' in res]
    print(f"PDB Residue IDs (first 20): {res_ids[:20]}")
    
    from app.services.topology.providers.tm.uniprot import UniprotTMProvider

    tm_prov = UniprotTMProvider()
    ss_prov = DSSPProvider()
    
    tm_pred = tm_prov.predict_tm(file, uniprot_id="P23975")
    print("TM Boundaries from Geometry:")
    for b in tm_pred.regions:
        print(f"  {b.start} - {b.end} ({b.type})")
        
    ss_map = ss_prov.get_secondary_structure(file)
    print("SS Map first 20 keys:")
    print(list(ss_map.keys())[:20])
    
    orch = TopologyOrchestrator(tm_prov, ss_prov, "ss_then_tm")
    res = orch.execute(file, params=TMParams(), uniprot_id="P23975")
    
    print("\nOrchestrator Regions:")
    for r in res.regions:
        print(f"  {r.type}: {r.start} - {r.end} ({r.description})")
