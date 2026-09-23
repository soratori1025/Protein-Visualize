import numpy as np
from pathlib import Path
import warnings
from typing import Optional
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from app.core.constants import KYTE_DOOLITTLE, SMOOTH_WINDOW, TM_HYDRO_THRESHOLD, TM_MIN_LENGTH, TM_MERGE_GAP
from app.schemas.topology import TMParams
from app.services.topology.providers.tm.base import TMProvider, TMPrediction, TMBoundary

def _smooth_hydrophobicity(sequence, window: int = SMOOTH_WINDOW) -> np.ndarray:
    scores = [KYTE_DOOLITTLE.get(res, 0.0) for res in sequence]
    half = window // 2
    out = []
    for i in range(len(scores)):
        lo, hi = max(0, i - half), min(len(scores), i + half + 1)
        out.append(float(np.mean(scores[lo:hi])))
    return np.array(out)

def _find_hydrophobic_segments(hydro, threshold, min_length, merge_gap=TM_MERGE_GAP):
    raw_segments = []
    in_seg = False
    start = 0
    for i, h in enumerate(hydro):
        if h > threshold:
            if not in_seg:
                in_seg = True
                start = i
        else:
            if in_seg:
                in_seg = False
                raw_segments.append((start, i - 1))
    if in_seg:
        raw_segments.append((start, len(hydro) - 1))

    merged = []
    half = SMOOTH_WINDOW // 2
    for seg in raw_segments:
        s = max(0, seg[0] - half)
        e = min(len(hydro) - 1, seg[1] + half)
        if merged and s - merged[-1][1] - 1 <= merge_gap:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))

    return [(s, e) for s, e in merged if (e - s + 1) >= min_length]

class SequenceTMProvider(TMProvider):
    def predict_tm(self, file_path: Path, params: Optional[TMParams] = None, **kwargs) -> TMPrediction:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            parser = MMCIFParser() if file_path.suffix.lower() in ('.cif', '.mmcif') else PDBParser()
            structure = parser.get_structure('protein', str(file_path))

        model = next(iter(structure))
        chain = next(iter(model))

        residues_data = []
        for residue in chain:
            if 'CA' in residue and residue.get_resname() in KYTE_DOOLITTLE:
                residues_data.append({'id': residue.get_id()[1], 'name': residue.get_resname()})
                
        sequence = [r['name'] for r in residues_data]
        hydro = _smooth_hydrophobicity(sequence)
        
        segments = _find_hydrophobic_segments(hydro, TM_HYDRO_THRESHOLD, TM_MIN_LENGTH)
        
        boundaries = []
        for s, e in segments:
            boundaries.append(TMBoundary(start=residues_data[s]['id'], end=residues_data[e]['id']))
            
        return TMPrediction(
            boundaries=boundaries,
            labeler="Kyte-Doolittle_Sequence"
        )
