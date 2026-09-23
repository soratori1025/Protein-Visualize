from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Mapping, Optional

from app.services.topology.residues import (
    ResidueFrame, SSAssignment, load_residue_frame, map_ss_records, parse_ss_records,
)


class SSProvider(ABC):
    """Secondary-structure provider.

    Subclasses only return the tool's raw per-residue records. Turning those into
    per-position codes (chain filter, insertion codes, numbering verification, offset /
    alignment fallback) is done once, here, for every tool."""

    labeler_name: str = ""

    @property
    def name(self) -> str:
        return self.labeler_name or self.__class__.__name__.replace("Provider", "")

    @abstractmethod
    def raw_residues(self, file_path: Path) -> List[Mapping]:
        """The tool's per-residue records in its own numbering, e.g.
        ``[{'residue_number': '100A', 'code': 'H', 'chain': 'A', 'amino_acid': 'L'}, ...]``."""

    def assign(self, file_path: Path, frame: Optional[ResidueFrame] = None) -> SSAssignment:
        """Per-position SS codes on ``frame`` (None where the tool gave nothing)."""
        frame = frame if frame is not None else load_residue_frame(file_path)
        records, skipped = parse_ss_records(self.raw_residues(file_path))
        result = map_ss_records(frame, records, labeler=self.name)
        if skipped:
            result.warnings.append(f"{self.name}: {skipped} records without a usable residue "
                                   "number were ignored (chain-break markers etc.)")
        return result

    def get_secondary_structure(self, file_path: Path,
                                frame: Optional[ResidueFrame] = None) -> Dict[int, str]:
        """Legacy view {author resseq: code} for the analysed chain. Lossy for insertion
        codes (100 and 100A share a key) - new code should call ``assign``."""
        frame = frame if frame is not None else load_residue_frame(file_path)
        result = self.assign(file_path, frame)
        out: Dict[int, str] = {}
        for residue, code in zip(frame.residues, result.codes):
            if code is not None:
                out.setdefault(residue.resseq, code)
        return out
