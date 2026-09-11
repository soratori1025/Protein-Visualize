import os
import subprocess
import tempfile
from pathlib import Path
from shutil import which

from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from protein_engine.secondary_structure.base import SecondaryStructureMethod
from protein_engine.secondary_structure.result import (
    ResidueSecondaryStructure,
    SecondaryStructureResult,
)

ROOT = Path(__file__).resolve().parents[2]
TOOLS_BIN = ROOT / "tools" / "bin"

class STRIDEMethod(SecondaryStructureMethod):
    name = "STRIDE"

    def __init__(self, executable: str = "stride"):
        if TOOLS_BIN.is_dir() and str(TOOLS_BIN) not in os.environ.get("PATH", ""):
            os.environ["PATH"] = f"{TOOLS_BIN}{os.pathsep}{os.environ.get('PATH', '')}"

        resolved = which(executable)
        if not resolved and (TOOLS_BIN / f"{executable}.exe").is_file():
            resolved = str(TOOLS_BIN / f"{executable}.exe")
        elif not resolved and (TOOLS_BIN / executable).is_file():
            resolved = str(TOOLS_BIN / executable)

        self.executable = resolved or executable

    def assign(self, structure_path: Path, model_id: int = 0) -> SecondaryStructureResult:

        if structure_path.suffix.lower() in {".cif", ".mmcif"}:
            with tempfile.TemporaryDirectory(prefix="proteinlab-stride-") as temporary_dir:
                pdb_path = Path(temporary_dir) / f"{structure_path.stem}.pdb"
                self._convert_to_pdb(structure_path, pdb_path, model_id)
                completed = self._run_stride(pdb_path)
        else:
            completed = self._run_stride(structure_path)

        residues = []
        for line in completed.stdout.splitlines():
            if not line.startswith("ASG"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                res_num = int(parts[3])
            except ValueError:
                res_num = int(parts[4]) if parts[4].isdigit() else 0

            residues.append(ResidueSecondaryStructure(
                chain_id=parts[2],
                residue_number=res_num,
                residue_name=parts[1],
                code=parts[5],
            ))
        return SecondaryStructureResult(method=self.name, residues=residues)

    def _run_stride(self, structure_path: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [self.executable, str(structure_path)],
            check=True,
            capture_output=True,
            text=True,
        )

    @staticmethod
    def _convert_to_pdb(structure_path: Path, output_path: Path, model_id: int) -> None:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(structure_path.stem, str(structure_path))
        io = PDBIO()
        io.set_structure(structure[model_id])
        io.save(str(output_path))