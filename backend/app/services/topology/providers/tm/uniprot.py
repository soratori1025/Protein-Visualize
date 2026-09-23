import json
import urllib.request
from urllib.error import HTTPError, URLError
from pathlib import Path
from fastapi import HTTPException
import warnings
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from app.schemas.topology import TMParams
from app.services.topology.providers.tm.base import TMProvider, TMPrediction, TMBoundary

class UniprotTMProvider(TMProvider):
    def _extract_uniprot_id_from_pdb(self, file_path: Path) -> str:
        """Attempt to extract UniProt ID from PDB or mmCIF headers."""
        # Simple extraction logic (can be expanded)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            try:
                if file_path.suffix.lower() in ('.cif', '.mmcif'):
                    parser = MMCIFParser()
                    structure = parser.get_structure('protein', str(file_path))
                    # mmCIF headers might have _struct_ref.pdbx_db_accession
                    header = parser.get_dict()
                    if '_struct_ref.pdbx_db_accession' in header:
                        accessions = header['_struct_ref.pdbx_db_accession']
                        if accessions:
                            return accessions[0] if isinstance(accessions, list) else accessions
                else:
                    parser = PDBParser()
                    structure = parser.get_structure('protein', str(file_path))
                    header = structure.header
                    if 'dbref' in header and header['dbref']:
                        for dbref in header['dbref']:
                            if dbref.get('db') == 'UNP': # UniProt
                                return dbref.get('accession')
            except Exception:
                pass
        return ""

    def predict_tm(self, file_path: Path, params: TMParams = None, **kwargs) -> TMPrediction:
        uniprot_id = kwargs.get("uniprot_id")
        
        if not uniprot_id:
            uniprot_id = self._extract_uniprot_id_from_pdb(file_path)
            
        if not uniprot_id:
            raise HTTPException(
                status_code=400, 
                detail="UniProt ID is required for the UniProt algorithm. Could not extract it from the PDB file. Please provide it manually."
            )
            
        clean_id = uniprot_id.strip().upper()
        url = f"https://rest.uniprot.org/uniprotkb/{clean_id}.json"
        req = urllib.request.Request(url, headers={"User-Agent": "ProteinVisualizeApp/1.0"})
        
        try:
            with urllib.request.urlopen(req, timeout=10) as response:
                if response.status != 200:
                    raise HTTPException(status_code=response.status, detail="Failed to fetch data from UniProt API")
                data = json.loads(response.read().decode("utf-8"))
        except HTTPError as err:
            if err.code == 404:
                raise HTTPException(status_code=404, detail=f"UniProt ID {clean_id} not found") from err
            raise HTTPException(status_code=502, detail=f"UniProt API error: {err.reason}") from err
        except URLError as err:
            raise HTTPException(status_code=504, detail=f"Network error connecting to UniProt API: {err.reason}") from err
            
        features = data.get("features", [])
        boundaries = []
        
        for f in features:
            ftype = f.get("type")
            if ftype in ("Transmembrane", "Intramembrane"):
                loc = f.get("location", {})
                start = loc.get("start", {}).get("value")
                end = loc.get("end", {}).get("value")
                if start is not None and end is not None:
                    boundaries.append(TMBoundary(start=start, end=end))
                    
        return TMPrediction(
            boundaries=boundaries,
            labeler=f"UniProt_{clean_id}"
        )
