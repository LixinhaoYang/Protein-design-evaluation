"""Synthetic mmCIF conversion checks; never read real model files."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from stage4_md_stability import cif_to_pdb

HEADER = """data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
"""


def row(serial, element, name, residue="ALA", chain="A", seq=1, model=1, group="ATOM", alt="."):
    return (f"{group} {serial} {element} {name} {alt} {residue} A 1 {seq} ? "
            f"1.000 2.000 3.000 1.00 80.00 {seq} {chain} {model}\n")


def cif(*rows: str) -> str:
    return HEADER + "".join(rows) + "#\n"


class TestConvert(unittest.TestCase):
    def convert(self, text: str, chain: str = "A"):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "model.cif"
            source.write_text(text)
            output = Path(tmp) / "out.pdb"
            count = cif_to_pdb.convert(source, output, chain)
            return count, output.read_text()

    def test_writes_fixed_pdb_columns(self):
        count, text = self.convert(cif(row(1, "N", "N"), row(2, "C", "CA")))
        self.assertEqual(count, 2)
        line = text.splitlines()[0]
        self.assertEqual(line[:6], "ATOM  ")
        self.assertEqual(line[12:16], " N  ")      # name, element letter in column 14
        self.assertEqual(line[16], " ")            # altLoc stays blank
        self.assertEqual(line[17:20], "ALA")       # resName
        self.assertEqual(line[21], "A")            # chainID
        self.assertEqual(line[22:26], "   1")      # resSeq
        self.assertEqual(line[76:78], " N")        # element

    def test_four_character_name_starts_at_column_13(self):
        _, text = self.convert(cif(row(1, "C", "CD11")))
        self.assertEqual(text.splitlines()[0][12:16], "CD11")

    def test_drops_hydrogens_and_oxt_and_renumbers(self):
        count, text = self.convert(
            cif(row(1, "N", "N"), row(2, "H", "HA"), row(3, "O", "OXT"), row(4, "C", "CB")))
        self.assertEqual(count, 2)
        serials = [line[6:11].strip() for line in text.splitlines() if line.startswith("ATOM")]
        self.assertEqual(serials, ["1", "2"])

    def test_terminates_the_chain(self):
        _, text = self.convert(cif(row(1, "N", "N")))
        self.assertEqual(text.splitlines()[-2:], ["TER", "END"])

    def test_selects_only_the_requested_chain(self):
        count, text = self.convert(cif(row(1, "N", "N", chain="A"), row(2, "N", "N", chain="B")))
        self.assertEqual(count, 1)
        self.assertNotIn("  B", text)

    def test_rejects_absent_chain(self):
        with self.assertRaisesRegex(ValueError, "not present"):
            self.convert(cif(row(1, "N", "N")), chain="Z")

    def test_rejects_hetatm_ligand_in_chain(self):
        with self.assertRaisesRegex(ValueError, "HETATM"):
            self.convert(cif(row(1, "N", "N"), row(2, "ZN", "ZN", residue="ZN", group="HETATM")))

    def test_rejects_second_model(self):
        with self.assertRaisesRegex(ValueError, "more than one model"):
            self.convert(cif(row(1, "N", "N"), row(2, "N", "N", model=2)))

    def test_rejects_alternate_conformer(self):
        with self.assertRaisesRegex(ValueError, "alternate conformers"):
            self.convert(cif(row(1, "N", "N"), row(2, "N", "N", alt="B")))

    def test_rejects_nonstandard_residue(self):
        with self.assertRaisesRegex(ValueError, "non-standard residue"):
            self.convert(cif(row(1, "N", "N", residue="MSE")))

    def test_rejects_missing_atom_site_loop(self):
        with self.assertRaisesRegex(ValueError, "no _atom_site loop"):
            self.convert("data_test\n#\n")

    def test_rejects_missing_required_column(self):
        # Drop the column and its field, so the field-count check passes and the
        # missing-column check is what rejects the file.
        text = cif(row(1, "N", "N")).replace("_atom_site.auth_asym_id\n", "")
        text = text.replace("1.00 80.00 1 A 1\n", "1.00 80.00 1 1\n")
        with self.assertRaisesRegex(ValueError, "missing auth_asym_id"):
            self.convert(text)

    def test_rejects_row_with_wrong_field_count(self):
        with self.assertRaisesRegex(ValueError, "column count"):
            self.convert(cif("ATOM 1 N N . ALA A 1 1 ?\n"))

    def test_rejects_empty_atom_site_loop(self):
        with self.assertRaisesRegex(ValueError, "no atoms"):
            self.convert(cif())

    def test_rejects_multi_character_chain(self):
        with self.assertRaisesRegex(ValueError, "one character"):
            self.convert(cif(row(1, "N", "N")), chain="AB")


class TestCli(unittest.TestCase):
    def test_refuses_to_overwrite_existing_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "cand_model.cif"
            source.write_text(cif(row(1, "N", "N")))
            outdir = Path(tmp) / "out"
            self.assertEqual(cif_to_pdb.main([str(source), "--outdir", str(outdir)]), 0)
            self.assertTrue((outdir / "cand.pdb").exists())
            self.assertEqual(cif_to_pdb.main([str(source), "--outdir", str(outdir)]), 2)


if __name__ == "__main__":
    unittest.main()
