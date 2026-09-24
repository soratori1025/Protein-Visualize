import os
import re
import subprocess
import tempfile
from pathlib import Path
from shutil import which

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.DSSP import DSSP, dssp_dict_from_pdb_file, residue_max_acc
from Bio.PDB.PDBExceptions import PDBException

from protein_engine.secondary_structure.base import SecondaryStructureMethod
from protein_engine.secondary_structure.result import (
    ResidueSecondaryStructure,
    SecondaryStructureResult,
    clean_angle,
    clean_float,
)

ROOT = Path(__file__).resolve().parents[2]
TOOLS_BIN = ROOT / "tools" / "bin"
TOOLS_SHARE = ROOT / "tools" / "share" / "libcifpp"

CIF_SUFFIXES = {".cif", ".mmcif"}
# mkdssp 4.x refuses some PDB files that have no CRYST1 record; a unit cell of 1 A
# is the conventional placeholder for non-crystallographic models.
DUMMY_CRYST1 = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1"

_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


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

    # ------------------------------------------------------------------ public
    def assign(self, structure_path: Path, model_id: int = 0) -> SecondaryStructureResult:
        structure_path = Path(structure_path)
        is_cif = structure_path.suffix.lower() in CIF_SUFFIXES
        temp_dir_obj = None
        working_path = structure_path
        if not is_cif:
            cleaned = self._clean_pdb_lines(structure_path)
            if cleaned is not None:
                temp_dir_obj = tempfile.TemporaryDirectory(prefix="proteinlab-dssp-")
                # keep a .pdb suffix whatever the original was (.ent, .pdb1 ...)
                working_path = Path(temp_dir_obj.name) / f"{structure_path.stem}.pdb"
                working_path.write_text("\n".join(cleaned) + "\n", encoding="utf-8")
        try:
            try:
                residues = self._via_biopython(working_path, is_cif, model_id)
            except (PDBException, KeyError) as error:
                # Bio.PDB.DSSP cross-checks every residue against the structure and aborts
                # the WHOLE run on one disagreement (e.g. MSE printed as 'M' by mkdssp).
                # The raw DSSP table does not need that check.
                print(f"[DSSP] Bio.PDB.DSSP mapping failed ({error}); using raw DSSP output")
                residues = self._via_raw_output(working_path)
            return SecondaryStructureResult(method=self.name, residues=residues)
        finally:
            if temp_dir_obj:
                temp_dir_obj.cleanup()

    # ----------------------------------------------------------------- helpers
    @staticmethod
    def _clean_pdb_lines(path: Path):
        """Lines of a cleaned copy, or None when the file can be used as is.
        Drops malformed REMARK lines and adds a CRYST1 record when missing."""
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
        bad_remark = [line.startswith("REMARK") and not line[6:10].strip().isdigit() for line in lines]
        has_cryst1 = any(line.startswith("CRYST1") for line in lines)
        needs_suffix = path.suffix.lower() != ".pdb"
        if not any(bad_remark) and has_cryst1 and not needs_suffix:
            return None
        cleaned = [line for line, bad in zip(lines, bad_remark) if not bad]
        if not has_cryst1:
            insert_at = next((i for i, line in enumerate(cleaned)
                              if line.startswith(("ATOM", "HETATM", "MODEL"))), 0)
            cleaned.insert(insert_at, DUMMY_CRYST1)
        return cleaned

    def _via_biopython(self, working_path: Path, is_cif: bool, model_id: int):
        parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)
        structure = parser.get_structure(working_path.stem, str(working_path))
        dssp = DSSP(structure[model_id], str(working_path), dssp=self.executable,
                    acc_array=self.asa_scale, file_type="MMCIF" if is_cif else "PDB")
        residues = []
        for key in dssp.keys():
            chain_id, (_, resseq, icode) = key
            value = dssp[key]
            residues.append(ResidueSecondaryStructure(
                chain_id=chain_id,
                residue_number=int(resseq),
                residue_name=value[1],
                code=value[2],
                phi=clean_angle(value[4]),
                psi=clean_angle(value[5]),
                asa=clean_float(value[3]),              # 'NA' for non-standard residues
                insertion_code=(icode or "").strip(),
            ))
        return residues

    def _via_raw_output(self, working_path: Path):
        version_text = subprocess.check_output([self.executable, "--version"], text=True)
        match = re.search(r"\s*([\d.]+)", version_text)
        table, keys = dssp_dict_from_pdb_file(str(working_path), self.executable,
                                              match.group(1) if match else "4.0.0")
        max_acc = residue_max_acc[self.asa_scale]
        residues = []
        for key in keys:
            chain_id, (_, resseq, icode) = key
            aa, ss, acc, phi, psi = table[key][:5]
            aa = "C" if aa.islower() else aa                     # a..z = bridged Cys
            three = _ONE_TO_THREE.get(aa)
            rel = min(1.0, acc / max_acc[three]) if three in max_acc else None
            residues.append(ResidueSecondaryStructure(
                chain_id=chain_id,
                residue_number=int(resseq),
                residue_name=aa,
                code=ss,
                phi=clean_angle(phi),
                psi=clean_angle(psi),
                asa=clean_float(rel),
                insertion_code=(icode or "").strip(),
            ))
        return residues
