import os
import tempfile
from pathlib import Path
from shutil import which
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.DSSP import DSSP
from protein_engine.secondary_structure.base import SecondaryStructureMethod
from protein_engine.secondary_structure.result import (
    ResidueSecondaryStructure,
    SecondaryStructureResult,
)

ROOT = Path(__file__).resolve().parents[2]
TOOLS_BIN = ROOT / "tools" / "bin"
TOOLS_SHARE = ROOT / "tools" / "share" / "libcifpp"

class DSSPMethod(SecondaryStructureMethod):
    name = "DSSP"
    def __init__(self, executable: str = "mkdssp", asa_scale: str = "Sander"):
        if TOOLS_BIN.is_dir() and str(TOOLS_BIN) not in os.environ.get("PATH", ""):
            os.environ["PATH"] = f"{TOOLS_BIN}{os.pathsep}{os.environ.get('PATH', '')}"
        if TOOLS_SHARE.is_dir() and "LIBCIFPP_DATA_DIR" not in os.environ:
            os.environ["LIBCIFPP_DATA_DIR"] = str(TOOLS_SHARE)
        resolved = which(executable)
        if not resolved and (TOOLS_BIN / f"{executable}.exe").is_file():
            resolved = str(TOOLS_BIN / f"{executable}.exe")
        elif not resolved and (TOOLS_BIN / executable).is_file():
            resolved = str(TOOLS_BIN / executable)
        self.executable = resolved or executable
        self.asa_scale = asa_scale

    def assign(self, structure_path: Path, model_id: int = 0) -> SecondaryStructureResult:
        working_path = structure_path
        temp_dir_obj = None
        if structure_path.suffix.lower() == ".pdb":
            content = structure_path.read_text(encoding="utf-8", errors="ignore")
            lines = content.splitlines()
            if any(line.startswith("REMARK") and not line[6:10].strip().isdigit() for line in lines):
                temp_dir_obj = tempfile.TemporaryDirectory(prefix="proteinlab-dssp-")
                clean_lines = [
                    line for line in lines
                    if not (line.startswith("REMARK") and not line[6:10].strip().isdigit())
                ]
                clean_file = Path(temp_dir_obj.name) / structure_path.name
                clean_file.write_text("\n".join(clean_lines), encoding="utf-8")
                working_path = clean_file
        try:
            parser = MMCIFParser(QUIET=True) if working_path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
            structure = parser.get_structure(working_path.stem, working_path)
            dssp = DSSP(structure[model_id], str(working_path), dssp=self.executable, acc_array=self.asa_scale)
            residues = []
            for key in dssp.keys():
                value = dssp[key]
                residues.append(ResidueSecondaryStructure(
                    chain_id=key[0],
                    residue_number=key[1][1],
                    residue_name=value[1],
                    code=value[2],
                    phi=value[4],
                    psi=value[5],
                    asa=value[3],
                ))
            return SecondaryStructureResult(method=self.name, residues=residues)
        finally:
            if temp_dir_obj:
                temp_dir_obj.cleanup()
