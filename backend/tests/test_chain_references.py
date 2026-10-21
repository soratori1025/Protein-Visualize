"""Tests for the per-chain UniProt mapping read from a structure file.

A file may hold several proteins, each with its own entry: 7YK7 is a fusion of a receptor,
a G-protein subunit and a soluble partner.  The accession used for the curated topology has
to be the one of the chain being analysed, not the one of whichever chain comes first.
"""
from __future__ import annotations

from pathlib import Path

from app.services.topology.providers.tm.uniprot import chain_references

CIF = """data_TEST
#
loop_
_struct_ref.id
_struct_ref.db_name
_struct_ref.db_code
_struct_ref.pdbx_db_accession
1 UNP A_HUMAN P00001
2 UNP B_HUMAN Q00002
3 PDB . TEST
#
loop_
_struct_ref_seq.ref_id
_struct_ref_seq.pdbx_strand_id
1 A,C
2 B
#
"""


def _dbref(chain: str, accession: str, db_id: str, seq_begin: int, seq_end: int) -> str:
    """A PDB DBREF record, at the columns the format defines.

    Columns (1-based): 1-6 record, 8-11 idCode, 13 chainID, 15-18 seqBegin, 20-23 seqEnd,
    26-28 database, 34-41 accession, 43-54 dbIdCode, 56-60 dbseqBegin, 63-67 dbseqEnd.
    """
    line = [" "] * 67

    def put(start: int, text: str) -> None:
        line[start - 1:start - 1 + len(text)] = text

    put(1, "DBREF ")
    put(8, "1ABC")
    put(13, chain)
    put(15, f"{seq_begin:>4}")
    put(20, f"{seq_end:>4}")
    put(26, "UNP")
    put(34, f"{accession:<8}")
    put(43, f"{db_id:<12}")
    put(56, f"{seq_begin:>5}")
    put(63, f"{seq_end:>5}")
    return "".join(line)


PDB = "\n".join([
    _dbref("A", "P00001", "A_HUMAN", 1, 100),
    _dbref("B", "Q00002", "B_HUMAN", 1, 50),
    "ATOM      1  N   MET A   1      11.104  13.207   2.100  1.00 20.00           N",
]) + "\n"


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_each_chain_gets_its_own_accession(tmp_path: Path) -> None:
    refs = chain_references(_write(tmp_path, "test.cif", CIF))
    assert refs == {"A": "P00001", "C": "P00001", "B": "Q00002"}


def test_pdb_dbref_records_are_mapped_per_chain(tmp_path: Path) -> None:
    refs = chain_references(_write(tmp_path, "test.pdb", PDB))
    assert refs == {"A": "P00001", "B": "Q00002"}


def test_each_chain_of_a_fusion_gets_its_own_entry(tmp_path: Path) -> None:
    """A fusion construct: one chain is a soluble partner, the other the receptor."""
    cif = """data_FUSION
#
loop_
_struct_ref.id
_struct_ref.db_name
_struct_ref.db_code
_struct_ref.pdbx_db_accession
1 UNP A_HUMAN P00001
2 UNP B_HUMAN Q00002
#
loop_
_struct_ref_seq.ref_id
_struct_ref_seq.pdbx_strand_id
1 A
2 B,C
#
"""
    refs = chain_references(_write(tmp_path, "fusion.cif", cif))
    assert refs["A"] == "P00001"
    assert refs["B"] == refs["C"] == "Q00002"


def test_a_file_without_records_gives_an_empty_map(tmp_path: Path) -> None:
    assert chain_references(_write(tmp_path, "empty.pdb", "ATOM\n")) == {}


def test_a_file_that_cannot_be_parsed_does_not_raise(tmp_path: Path) -> None:
    broken = _write(tmp_path, "broken.cif", "data_BROKEN\n#\nloop_\n_struct_ref.id\n1 2\n")
    assert chain_references(broken) == {}
    assert chain_references(tmp_path / "missing.cif") == {}
