import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "backend"))

from app.api.secondary_structure import secondary_structure_capabilities
from protein_engine.secondary_structure.dssp import DSSPMethod
from protein_engine.secondary_structure.stride import STRIDEMethod

TEST_CIF = ROOT / "tests" / "12BN.cif"


class TestSecondaryStructure(unittest.TestCase):
    def test_capabilities(self):
        caps = secondary_structure_capabilities()
        self.assertTrue(caps["DSSP"]["available"], "DSSP should be available")
        self.assertTrue(caps["STRIDE"]["available"], "STRIDE should be available")

    def test_dssp(self):
        method = DSSPMethod()
        result = method.assign(TEST_CIF)
        self.assertEqual(result.method, "DSSP")
        self.assertGreater(len(result.residues), 0)
        first_res = result.residues[0]
        self.assertEqual(first_res.chain_id, "E")
        self.assertEqual(first_res.residue_name, "Y")

    def test_stride(self):
        method = STRIDEMethod()
        result = method.assign(TEST_CIF)
        self.assertEqual(result.method, "STRIDE")
        self.assertGreater(len(result.residues), 0)
        first_res = result.residues[0]
        self.assertEqual(first_res.chain_id, "E")
        self.assertEqual(first_res.residue_name, "TYR")


if __name__ == "__main__":
    unittest.main()
