from typing import Optional
from pathlib import Path
from pydantic import BaseModel

from app.schemas.topology import TMParams, TopologyResponse, TopologyRegion
from app.services.topology.providers.tm.base import TMProvider, TMPrediction
from app.services.topology.providers.ss.base import SSProvider
from app.core.constants import MAX_SNAP

class TopologyOrchestrator:
    def __init__(self, tm_provider: TMProvider, ss_provider: Optional[SSProvider], flow_type: str):
        """
        flow_type can be:
        - "tm_then_ss": TM runs first to find boundaries. SS is then applied only within those boundaries.
        - "ss_then_tm": SS runs first to find elements. TM finds boundaries, then snaps to nearest SS element.
        - "parallel_merge": Both run independently. Residues that are both in TM boundary and SS element are labeled.
        """
        self.tm_provider = tm_provider
        self.ss_provider = ss_provider
        self.flow_type = flow_type

    def execute(self, file_path: Path, params: Optional[TMParams] = None, **kwargs) -> TopologyResponse:
        if params is None:
            params = TMParams()

        # 1. Run TM Prediction
        tm_pred = self.tm_provider.predict_tm(file_path, params=params, **kwargs)
        
        # Fast exit if no TM found
        if not tm_pred.boundaries:
            return TopologyResponse(
                uniprot_id="CALCULATED",
                protein_name=f"Predicted ({tm_pred.labeler})",
                gene_name="",
                organism="Computed",
                membrane_score=tm_pred.membrane_score,
                membrane_normal=tm_pred.membrane_normal,
                labeler=tm_pred.labeler,
                parameters_used=params.to_response_dict(),
                regions=[]
            )

        # 2. Run SS Prediction (if provider exists)
        ss_map = {}
        ss_labeler = "None"
        if self.ss_provider:
            ss_map = self.ss_provider.get_secondary_structure(file_path)
            ss_labeler = self.ss_provider.__class__.__name__.replace("Provider", "")
            
        final_labeler = f"{tm_pred.labeler} + {ss_labeler} ({self.flow_type})"

        # 3. Apply the specific flow
        regions = []
        if self.flow_type == "tm_then_ss":
            regions = self._flow_tm_then_ss(tm_pred, ss_map)
        elif self.flow_type == "ss_then_tm":
            regions = self._flow_ss_then_tm(tm_pred, ss_map, max_snap=MAX_SNAP)
        elif self.flow_type == "parallel_merge":
            regions = self._flow_parallel_merge(tm_pred, ss_map)
        else:
            # Default to basic boundaries if flow is unknown or None
            regions = [TopologyRegion(type="Transmembrane", start=b.start, end=b.end, description="Transmembrane region") for b in tm_pred.boundaries]

        return TopologyResponse(
            uniprot_id="CALCULATED",
            protein_name=f"Predicted Topology",
            gene_name="",
            organism="Computed",
            membrane_score=tm_pred.membrane_score,
            membrane_normal=tm_pred.membrane_normal,
            labeler=final_labeler,
            parameters_used=params.to_response_dict(),
            regions=regions
        )

    def _flow_tm_then_ss(self, tm_pred: TMPrediction, ss_map: dict) -> list[TopologyRegion]:
        """
        Strict TM boundaries. For each TM boundary, determine the majority SS within it.
        """
        regions = []
        for b in tm_pred.boundaries:
            # Count SS types within this boundary
            counts = {'H': 0, 'E': 0, 'C': 0}
            for i in range(b.start, b.end + 1):
                c = ss_map.get(i, 'C')
                if c in ['H', 'G', 'I']: counts['H'] += 1
                elif c in ['E', 'B']: counts['E'] += 1
                else: counts['C'] += 1
                
            majority = max(counts, key=counts.get)
            desc = "Transmembrane Alpha Helix" if majority == 'H' else "Transmembrane Beta Strand" if majority == 'E' else "Transmembrane Region"
            
            regions.append(TopologyRegion(type="Transmembrane", start=b.start, end=b.end, description=desc))
        return regions

    def _flow_ss_then_tm(self, tm_pred: TMPrediction, ss_map: dict, max_snap: int = 4) -> list[TopologyRegion]:
        """
        Snap TM boundaries to nearest complete SS element.
        """
        # First group SS into contiguous blocks
        ss_blocks = []
        current_block = None
        sorted_ids = sorted(ss_map.keys())
        if not sorted_ids:
            return self._flow_tm_then_ss(tm_pred, ss_map)
            
        for rid in sorted_ids:
            code = ss_map[rid]
            c_type = 'H' if code in ['H', 'G', 'I'] else 'E' if code in ['E', 'B'] else 'C'
            
            if current_block is None:
                current_block = {'type': c_type, 'start': rid, 'end': rid}
            elif current_block['type'] == c_type and rid == current_block['end'] + 1:
                current_block['end'] = rid
            else:
                ss_blocks.append(current_block)
                current_block = {'type': c_type, 'start': rid, 'end': rid}
        if current_block:
            ss_blocks.append(current_block)
            
        regions = []
        for b in tm_pred.boundaries:
            # Find all SS blocks that overlap with this TM boundary
            overlapping = [bk for bk in ss_blocks if bk['start'] <= b.end and bk['end'] >= b.start and bk['type'] in ['H', 'E']]
            
            if not overlapping:
                regions.append(TopologyRegion(type="Transmembrane", start=b.start, end=b.end, description="Transmembrane Region"))
                continue
                
            # Merge overlapping blocks that are close
            merged_start = min([bk['start'] for bk in overlapping])
            merged_end = max([bk['end'] for bk in overlapping])
            
            # Apply snapping constraints
            final_start = max(b.start - max_snap, merged_start)
            final_end = min(b.end + max_snap, merged_end)
            
            # Majority type
            h_len = sum([min(final_end, bk['end']) - max(final_start, bk['start']) + 1 for bk in overlapping if bk['type'] == 'H'])
            e_len = sum([min(final_end, bk['end']) - max(final_start, bk['start']) + 1 for bk in overlapping if bk['type'] == 'E'])
            desc = "Transmembrane Alpha Helix" if h_len >= e_len else "Transmembrane Beta Strand"
            
            regions.append(TopologyRegion(type="Transmembrane", start=final_start, end=final_end, description=desc))
            
        return regions

    def _flow_parallel_merge(self, tm_pred: TMPrediction, ss_map: dict) -> list[TopologyRegion]:
        """
        Only residues that are BOTH in TM and in an SS element are classified as Transmembrane Helix/Strand.
        """
        # Create a boolean mask of TM residues
        tm_set = set()
        for b in tm_pred.boundaries:
            tm_set.update(range(b.start, b.end + 1))
            
        regions = []
        in_region = False
        start = -1
        current_type = None
        
        sorted_ids = sorted(ss_map.keys()) if ss_map else sorted(list(tm_set))
        if not sorted_ids:
            return []
            
        max_id = max(sorted_ids)
        for rid in range(min(sorted_ids), max_id + 1):
            is_tm = rid in tm_set
            code = ss_map.get(rid, 'C')
            is_ss = code in ['H', 'G', 'I', 'E', 'B']
            
            if is_tm and is_ss:
                c_type = 'H' if code in ['H', 'G', 'I'] else 'E'
                if not in_region:
                    in_region = True
                    start = rid
                    current_type = c_type
                elif current_type != c_type:
                    # Type changed mid-region (e.g. H -> E)
                    desc = "Transmembrane Alpha Helix" if current_type == 'H' else "Transmembrane Beta Strand"
                    regions.append(TopologyRegion(type="Transmembrane", start=start, end=rid-1, description=desc))
                    start = rid
                    current_type = c_type
            else:
                if in_region:
                    desc = "Transmembrane Alpha Helix" if current_type == 'H' else "Transmembrane Beta Strand"
                    regions.append(TopologyRegion(type="Transmembrane", start=start, end=rid-1, description=desc))
                    in_region = False
                    
        if in_region:
            desc = "Transmembrane Alpha Helix" if current_type == 'H' else "Transmembrane Beta Strand"
            regions.append(TopologyRegion(type="Transmembrane", start=start, end=max_id, description=desc))
            
        return regions
