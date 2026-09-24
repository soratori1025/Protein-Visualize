import math
from dataclasses import asdict, dataclass
from typing import Any, Optional


@dataclass
class ResidueSecondaryStructure:
    chain_id: str
    residue_number: int
    residue_name: str
    code: str
    phi: Optional[float] = None
    psi: Optional[float] = None
    asa: Optional[float] = None          # RELATIVE accessibility 0..1 (same scale for DSSP and STRIDE)
    # PDB insertion code ("A" in residue 100A); "" when blank. Without it 100 and
    # 100A collapse onto the same residue_number and cannot be told apart downstream.
    insertion_code: str = ""


@dataclass
class SecondaryStructureResult:
    method: str
    residues: list[ResidueSecondaryStructure]

    def to_dict(self) -> dict:
        return {
            "method": self.method,
            "residues": [asdict(residue) for residue in self.residues],
        }


def clean_float(value: Any) -> Optional[float]:
    """None for missing / non-numeric values ('NA', '', NaN) - keeps the JSON typed."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def clean_angle(value: Any) -> Optional[float]:
    """DSSP and STRIDE print 360.0 for an undefined phi/psi (chain ends, breaks)."""
    number = clean_float(value)
    return None if number is None or abs(number) >= 360.0 else number
