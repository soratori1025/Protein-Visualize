from dataclasses import asdict, dataclass
from typing import Optional

@dataclass
class ResidueSecondaryStructure:
    chain_id: str
    residue_number: int
    residue_name: str
    code: str
    phi: Optional[float] = None
    psi: Optional[float] = None
    asa: Optional[float] = None

@dataclass
class SecondaryStructureResult:
    method: str
    residues: list[ResidueSecondaryStructure]
    def to_dict(self) -> dict:
        return {
            "method": self.method,
            "residues": [asdict(residue) for residue in self.residues],
        }