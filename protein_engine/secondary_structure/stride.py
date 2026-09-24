import os
import re
import subprocess
import tempfile
from pathlib import Path
from shutil import which

from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from Bio.PDB.DSSP import residue_max_acc

from protein_engine.secondary_structure.base import SecondaryStructureMethod
from protein_engine.secondary_structure.result import (
    ResidueSecondaryStructure,
    SecondaryStructureResult,
    clean_angle,
    clean_float,
)

ROOT = Path(__file__).resolve().parents[2]
TOOLS_BIN = ROOT / "tools" / "bin"

CIF_SUFFIXES = {".cif", ".mmcif"}
STRIDE_TIMEOUT_S = 300
# PDB residue field as STRIDE prints it: number plus optional insertion code ("100A")
_RESNUM_RE = re.compile(r"^(-?\d+)([A-Za-z]?)$")
_SINGLE_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


class STRIDEMethod(SecondaryStructureMethod):
    name = "STRIDE"

    def __init__(self, executable: str = "stride", asa_scale: str = "Sander"):
        if TOOLS_BIN.is_dir() and str(TOOLS_BIN) not in os.environ.get("PATH", ""):
            os.environ["PATH"] = f"{TOOLS_BIN}{os.pathsep}{os.environ.get('PATH', '')}"

        resolved = which(executable)
        if not resolved and (TOOLS_BIN / f"{executable}.exe").is_file():
            resolved = str(TOOLS_BIN / f"{executable}.exe")
        elif not resolved and (TOOLS_BIN / executable).is_file():
            resolved = str(TOOLS_BIN / executable)

        self.executable = resolved or executable
        self.asa_scale = asa_scale

    def assign(self, structure_path: Path, model_id: int = 0) -> SecondaryStructureResult:
        structure_path = Path(structure_path)
        chain_names: dict[str, str] = {}
        is_cif = structure_path.suffix.lower() in CIF_SUFFIXES
        # STRIDE reads PDB only. mmCIF (and a non-first model) is written out as PDB first.
        if is_cif or model_id != 0:
            with tempfile.TemporaryDirectory(prefix="proteinlab-stride-") as temporary_dir:
                pdb_path = Path(temporary_dir) / f"{structure_path.stem}.pdb"
                chain_names = self._convert_to_pdb(structure_path, pdb_path, model_id, is_cif)
                completed = self._run_stride(pdb_path)
        else:
            completed = self._run_stride(structure_path)
        return SecondaryStructureResult(
            method=self.name, residues=self._parse_asg(completed.stdout, chain_names))

    # ----------------------------------------------------------------- parsing
    def _parse_asg(self, stdout: str, chain_names: dict[str, str]):
        """ASG  ALA A  100A   57    H    AlphaHelix    -62.1    -41.3      12.0      ~~~~
        fields: resname, chain ('-' = blank), PDB number (+ insertion code),
        STRIDE ordinal, code, name, phi, psi, area (A^2)."""
        max_acc = residue_max_acc[self.asa_scale]
        residues = []
        for line in stdout.splitlines():
            if not line.startswith("ASG"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            match = _RESNUM_RE.match(parts[3])
            if not match:
                # never fall back to the STRIDE ordinal (parts[4]): it is a different
                # numbering and silently shifts the residue onto another one
                continue
            chain = " " if parts[2] == "-" else parts[2]          # STRIDE prints blank chain as '-'
            resname = parts[1].upper()
            area = clean_float(parts[9]) if len(parts) > 9 else None
            rel = min(1.0, area / max_acc[resname]) if area is not None and resname in max_acc else None
            residues.append(ResidueSecondaryStructure(
                chain_id=chain_names.get(chain, chain),
                residue_number=int(match.group(1)),
                residue_name=resname,
                code=parts[5],
                phi=clean_angle(parts[7]) if len(parts) > 7 else None,
                psi=clean_angle(parts[8]) if len(parts) > 8 else None,
                asa=clean_float(rel),
                insertion_code=match.group(2),
            ))
        return residues

    # --------------------------------------------------------------- execution
    def _run_stride(self, structure_path: Path) -> subprocess.CompletedProcess[str]:
        try:
            return subprocess.run(
                [self.executable, str(structure_path)],
                check=True,
                capture_output=True,
                text=True,
                timeout=STRIDE_TIMEOUT_S,
            )
        except subprocess.CalledProcessError as error:
            tail = (error.stderr or error.stdout or "").strip().splitlines()[-3:]
            raise RuntimeError(f"STRIDE failed (exit {error.returncode}): {' | '.join(tail)}") from error

    @staticmethod
    def _convert_to_pdb(structure_path: Path, output_path: Path, model_id: int,
                        is_cif: bool) -> dict[str, str]:
        """Write one model as PDB. PDB allows 1-character chain IDs only, so longer
        mmCIF chain IDs are renamed to free single characters; the returned map
        {temporary_id: original_id} is used to put the original IDs back."""
        parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)
        structure = parser.get_structure(structure_path.stem, str(structure_path))
        model = structure[model_id]
        used = {chain.id for chain in model if len(chain.id) == 1}
        free = (c for c in _SINGLE_CHAIN_IDS if c not in used)
        renamed: dict[str, str] = {}
        for chain in list(model):
            if len(chain.id) != 1:
                new_id = next(free, None)
                if new_id is None:
                    raise RuntimeError("Too many chains to convert to PDB format for STRIDE")
                renamed[new_id] = chain.id
                chain.id = new_id
        io = PDBIO()
        io.set_structure(model)
        io.save(str(output_path))
        return renamed
