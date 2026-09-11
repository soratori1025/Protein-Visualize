from pathlib import Path

from Bio.PDB import MMCIFParser, PDBParser


def parse_pdb(path: str):
    """Parse a PDB or mmCIF structure file into a Biopython Structure object."""
    structure_path = Path(path)
    structure_id = structure_path.stem
    if structure_path.suffix.lower() in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
    return parser.get_structure(structure_id, str(structure_path))
