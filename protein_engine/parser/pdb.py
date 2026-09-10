from pathlib import Path

from Bio.PDB import PDBParser


def parse_pdb(path: str):
    """Parse a PDB or mmCIF-like structure file into a Biopython Structure object."""
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure_id = Path(path).stem
    return parser.get_structure(structure_id, path)
