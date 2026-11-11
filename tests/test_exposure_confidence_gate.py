"""The exposed-hydrophobicity flag must not fire on a low-confidence predicted model.

Exposure is measured on a structure, so on a prediction it inherits the prediction's
errors: a model the predictor got wrong is loose and reads as exposed whether or not the
design is. Measured case behind this gate -- a design whose ESMFold model scored pLDDT 55
sat 10.6 A from its NMR structure, while its design model was 1.3 A from it. The design
was right, the prediction was wrong, and the exposure flag saw only the wrong structure.

The gate is opt-in because the B-factor column means pLDDT only for predicted structures;
in an experimental structure it is a temperature factor.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from stage2_fast_physics.exposed_hydrophobicity import (  # noqa: E402
    DEFAULT_MIN_PLDDT,
    mean_plddt_from_bfactor,
    summarize_pdbs,
)

# One leucine per residue, residues placed far apart so every side chain is fully exposed
# and the flag fires whenever it is allowed to. Only the B-factor column varies.
LEUCINE_ATOMS = [
    ("N", "N", 0.000, 0.000, 0.000),
    ("CA", "C", 1.458, 0.000, 0.000),
    ("C", "C", 2.009, 1.420, 0.000),
    ("O", "O", 1.251, 2.390, 0.000),
    ("CB", "C", 1.988, -0.773, 1.206),
    ("CG", "C", 3.500, -0.900, 1.300),
    ("CD1", "C", 4.000, -1.700, 2.500),
    ("CD2", "C", 4.100, 0.470, 1.400),
]


def pdb_line(serial: int, name: str, element: str, resseq: int,
             x: float, y: float, z: float, bfactor: float) -> str:
    """One fixed-width ATOM record. Columns follow the PDB spec exactly."""
    atom_name = f" {name:<3}" if len(name) < 4 else name
    return (
        f"ATOM  {serial:>5} {atom_name}{'':1}LEU {'A'}{resseq:>4}{'':1}   "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{bfactor:>6.2f}"
        f"{'':10}{element:>2}"
    )


def write_model(directory: Path, name: str, plddt: float, n_residues: int = 3) -> Path:
    path = directory / f"{name}.pdb"
    lines, serial = [], 1
    for index in range(n_residues):
        shift = index * 25.0  # far enough apart that residues never contact each other
        for atom_name, element, x, y, z in LEUCINE_ATOMS:
            lines.append(pdb_line(serial, atom_name, element, index + 1,
                                  x + shift, y, z, plddt))
            serial += 1
    path.write_text("\n".join(lines) + "\n")
    return path


class TestMeanPlddtFromBfactor(unittest.TestCase):
    def test_reads_the_bfactor_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_model(Path(tmp), "confident", 88.5)
            self.assertAlmostEqual(mean_plddt_from_bfactor(path), 88.5, places=2)

    def test_returns_none_without_ca_atoms(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.pdb"
            path.write_text("REMARK no atoms here\n")
            self.assertIsNone(mean_plddt_from_bfactor(path))


class TestConfidenceGate(unittest.TestCase):
    def rows_for(self, plddt: float, **kwargs):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_model(Path(tmp), "model", plddt)
            return summarize_pdbs([path], max_exposed_hydro_per_res=0.1, **kwargs)[0]

    def test_gate_off_by_default(self):
        row = self.rows_for(40.0)
        self.assertEqual(row["confidence_status"], "not_checked")
        self.assertEqual(row["mean_plddt"], "")
        self.assertEqual(row["flag_exposed_hydrophobicity"], 1,
                         "with the gate off the flag must behave as before")

    def test_low_confidence_withholds_the_flag(self):
        row = self.rows_for(55.0, min_plddt=DEFAULT_MIN_PLDDT)
        self.assertEqual(row["confidence_status"], "prediction_not_trusted")
        self.assertEqual(row["flag_exposed_hydrophobicity"], 0)
        self.assertEqual(row["total_flags"], 0)

    def test_metrics_are_still_reported_when_not_trusted(self):
        """The numbers stay in the table as features; only the verdict is withheld."""
        row = self.rows_for(55.0, min_plddt=DEFAULT_MIN_PLDDT)
        self.assertGreater(row["exposed_hydrophobic_sasa"], 0.0)
        self.assertEqual(row["mean_plddt"], 55.0)

    def test_confident_model_still_flags(self):
        row = self.rows_for(88.0, min_plddt=DEFAULT_MIN_PLDDT)
        self.assertEqual(row["confidence_status"], "trusted")
        self.assertEqual(row["flag_exposed_hydrophobicity"], 1)

    def test_threshold_is_inclusive(self):
        self.assertEqual(
            self.rows_for(DEFAULT_MIN_PLDDT, min_plddt=DEFAULT_MIN_PLDDT)["confidence_status"],
            "trusted")
        self.assertEqual(
            self.rows_for(DEFAULT_MIN_PLDDT - 0.1, min_plddt=DEFAULT_MIN_PLDDT)["confidence_status"],
            "prediction_not_trusted")


if __name__ == "__main__":
    unittest.main()
