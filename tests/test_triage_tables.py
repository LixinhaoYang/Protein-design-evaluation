"""Regression checks for missing results and candidate matching.

All values are synthetic; these tests do not require structures or PyRosetta.
"""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from stage2_fast_physics.triage_features import summarize
from triage_table import build_rows, read_csv_by_id


ROOT = Path(__file__).resolve().parents[1]


def relaxation_row() -> dict:
    return {
        "structure_id": "model_001",
        "relaxed_pdb": "relaxed/model_001_restrained_relax.pdb",
        "length": "100",
        "pre_relax_clashscore": "0",
        "post_relax_clashscore": "0",
        "relax_CA_RMSD": "0.2",
        "relax_motif_RMSD": "",
        "max_CA_displacement": "0.5",
        "total_flags": "0",
        "warn_pre_relax_clash": "0",
    }


def rosetta_row() -> dict:
    return {
        "structure_id": "model_001_restrained_relax",
        "source": "relaxed/model_001_restrained_relax.pdb",
        "length": "100",
        "rosetta_score_per_residue": "-2.0",
        "fa_rep_weighted": "0",
        "rosetta_buried_unsat_hbonds": "0",
        "rosetta_sap_score": "0",
        "rosetta_packstat": "0.7",
    }


def stage_row(status: str = "PASS", n_flags: int = 0) -> dict:
    return {
        "sequence_id": "candidate_001",
        "length": "100",
        "n_flags": str(n_flags),
        "triage_status": status,
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class Stage2CompletenessTests(unittest.TestCase):
    def test_complete_row_allows_zero_metrics_and_optional_blank_motif(self):
        row = summarize(relax_rows=[relaxation_row()], rosetta_rows=[rosetta_row()])[0]
        self.assertEqual(row["triage_status"], "PASS")
        self.assertEqual(row["triage_notes"], "")
        self.assertEqual(row["n_flags"], 0)
        self.assertEqual(row["fa_rep"], "0")
        self.assertEqual(row["relax_motif_RMSD"], "")

    def test_structure_suffix_fallback_and_explicit_id_map(self):
        rosetta = rosetta_row()
        rosetta["source"] = "/other/path/model_001_restrained_relax.pdb"
        row = summarize(
            relax_rows=[relaxation_row()], rosetta_rows=[rosetta],
            id_map={"model_001": "candidate_001"},
        )[0]
        self.assertEqual(row["sequence_id"], "candidate_001")
        self.assertEqual(row["triage_status"], "PASS")
        self.assertEqual(row["rosetta_score_per_residue"], "-2.0")

    def test_missing_or_unmatched_rosetta_row_is_incomplete(self):
        unrelated = rosetta_row()
        unrelated.update(structure_id="other", source="other.pdb")
        for rows in ([], [unrelated]):
            with self.subTest(rows=rows):
                row = summarize(relax_rows=[relaxation_row()], rosetta_rows=rows)[0]
                self.assertEqual(row["triage_status"], "INCOMPLETE")
                self.assertIn("missing Rosetta row", row["triage_notes"])
                self.assertEqual(row["n_flags"], 0)

    def test_unavailable_rosetta_metrics_are_incomplete(self):
        fields = (
            "rosetta_score_per_residue", "fa_rep_weighted",
            "rosetta_buried_unsat_hbonds", "rosetta_sap_score", "rosetta_packstat",
        )
        for field in fields:
            for value in ("", None, "nan", "inf", "-inf", "not-a-number"):
                with self.subTest(field=field, value=value):
                    rosetta = rosetta_row()
                    rosetta[field] = value
                    row = summarize(relax_rows=[relaxation_row()], rosetta_rows=[rosetta])[0]
                    self.assertEqual(row["triage_status"], "INCOMPLETE")
                    self.assertIn(field, row["triage_notes"])

    def test_missing_relaxation_results_are_incomplete(self):
        for field in ("post_relax_clashscore", "total_flags", "warn_pre_relax_clash"):
            with self.subTest(field=field):
                relax = relaxation_row()
                del relax[field]
                row = summarize(relax_rows=[relax], rosetta_rows=[rosetta_row()])[0]
                self.assertEqual(row["triage_status"], "INCOMPLETE")
                self.assertIn(field, row["triage_notes"])

    def test_known_failure_is_preserved_when_rosetta_is_missing(self):
        relax = relaxation_row()
        relax["total_flags"] = "1"
        row = summarize(relax_rows=[relax], rosetta_rows=[])[0]
        self.assertEqual(row["triage_status"], "FAIL")
        self.assertEqual(row["n_flags"], 1)
        self.assertIn("missing Rosetta row", row["triage_notes"])

    def test_pre_relax_warning_requires_complete_metrics(self):
        relax = relaxation_row()
        relax["warn_pre_relax_clash"] = "1"
        complete = summarize(relax_rows=[relax], rosetta_rows=[rosetta_row()])[0]
        incomplete = summarize(relax_rows=[relax], rosetta_rows=[])[0]
        self.assertEqual(complete["triage_status"], "WARN")
        self.assertEqual(incomplete["triage_status"], "INCOMPLETE")
        self.assertEqual(complete["n_flags"], 1)

    def test_duplicate_rosetta_matches_are_rejected(self):
        original = rosetta_row()
        same_source = {**original, "structure_id": "other"}
        same_id = {**original, "source": "other.pdb"}
        for duplicate in (same_source, same_id):
            with self.subTest(duplicate=duplicate):
                with self.assertRaisesRegex(ValueError, "duplicate Rosetta"):
                    summarize(relax_rows=[relaxation_row()], rosetta_rows=[original, duplicate])

    def test_empty_structure_id_is_rejected(self):
        relax = relaxation_row()
        relax["structure_id"] = ""
        with self.assertRaisesRegex(ValueError, "empty structure_id"):
            summarize(relax_rows=[relax], rosetta_rows=[rosetta_row()])


class AggregateCompletenessTests(unittest.TestCase):
    def test_pre_af3_table_does_not_require_stage3(self):
        stage = {"candidate_001": stage_row()}
        row = build_rows(stage, stage, None)[0]
        self.assertEqual(row["triage_status"], "PASS")
        self.assertEqual(row["triage_notes"], "")

    def test_missing_required_stage_is_incomplete(self):
        stage = {"candidate_001": stage_row()}
        cases = (({}, stage, None, "Stage 1"), (stage, {}, None, "Stage 2"),
                 (stage, stage, {}, "Stage 3"))
        for sequence, physics, fold, missing in cases:
            with self.subTest(missing=missing):
                row = build_rows(sequence, physics, fold)[0]
                self.assertEqual(row["triage_status"], "INCOMPLETE")
                self.assertIn(f"missing {missing} row", row["triage_notes"])

    def test_supplied_complete_stage3_table_passes(self):
        stage = {"candidate_001": stage_row()}
        self.assertEqual(build_rows(stage, stage, stage)[0]["triage_status"], "PASS")

    def test_candidate_union_retains_unmatched_ids(self):
        rows = build_rows({"sequence_only": stage_row()}, {"physics_only": stage_row()}, None)
        self.assertEqual({row["sequence_id"] for row in rows}, {"sequence_only", "physics_only"})
        self.assertTrue(all(row["triage_status"] == "INCOMPLETE" for row in rows))

    def test_failure_takes_precedence_over_incomplete(self):
        row = build_rows({"candidate_001": stage_row("FAIL", 1)}, {}, None)[0]
        self.assertEqual(row["triage_status"], "FAIL")
        self.assertEqual(row["n_flags"], 1)
        self.assertIn("missing Stage 2 row", row["triage_notes"])

    def test_stage2_incompleteness_and_notes_propagate(self):
        physics = stage_row("INCOMPLETE")
        physics["triage_notes"] = "missing Rosetta row"
        row = build_rows({"candidate_001": stage_row("WARN", 2)},
                         {"candidate_001": physics}, None)[0]
        self.assertEqual(row["triage_status"], "INCOMPLETE")
        self.assertEqual(row["triage_notes"], "Stage 2: missing Rosetta row")
        self.assertEqual(row["n_flags"], 2)

    def test_duplicate_or_empty_sequence_ids_are_rejected(self):
        for rows in ([stage_row(), stage_row()], [{**stage_row(), "sequence_id": ""}]):
            with self.subTest(rows=rows), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "stage.csv"
                write_csv(path, rows)
                with self.assertRaisesRegex(ValueError, "duplicate sequence_id|empty sequence_id"):
                    read_csv_by_id(path)

    def test_wrong_summary_schema_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "summary.csv"
            write_csv(path, [{"sequence_id": "candidate_001", "total_flags": "0"}])
            with self.assertRaisesRegex(ValueError, "missing columns"):
                read_csv_by_id(path)

    def test_blank_or_unknown_status_is_rejected(self):
        for status in ("", "PENDING", "pass"):
            with self.subTest(status=status), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "stage.csv"
                write_csv(path, [stage_row(status)])
                with self.assertRaisesRegex(ValueError, "invalid triage_status"):
                    read_csv_by_id(path)

    def test_missing_or_invalid_flag_count_is_rejected(self):
        for value in ("", "nan", "inf", "-1", "0.5", "not-a-number"):
            with self.subTest(value=value), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "stage.csv"
                write_csv(path, [{**stage_row(), "n_flags": value}])
                with self.assertRaisesRegex(ValueError, "nonnegative integer n_flags"):
                    read_csv_by_id(path)


class CommandLineTests(unittest.TestCase):
    def run_script(self, script: str, *args: object) -> subprocess.CompletedProcess:
        return subprocess.run(
            [sys.executable, str(ROOT / script), *(str(arg) for arg in args)],
            cwd=ROOT, capture_output=True, text=True, check=False,
        )

    def test_stage2_and_aggregate_csv_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            relax_path = directory / "relax.csv"
            rosetta_path = directory / "rosetta.csv"
            map_path = directory / "ids.csv"
            sequence_path = directory / "sequence.csv"
            physics_path = directory / "stage2.csv"
            output_path = directory / "nested" / "triage.csv"
            write_csv(relax_path, [relaxation_row()])
            write_csv(map_path, [{"structure_id": "model_001", "sequence_id": "candidate_001"}])
            write_csv(sequence_path, [stage_row()])

            for sap, expected in (("0", "PASS"), ("", "INCOMPLETE")):
                with self.subTest(expected=expected):
                    write_csv(rosetta_path, [{**rosetta_row(), "rosetta_sap_score": sap}])
                    physics_run = self.run_script(
                        "stage2_fast_physics/triage_features.py",
                        "--relax-csv", relax_path, "--rosetta-csv", rosetta_path,
                        "--id-map-csv", map_path, "--csv", physics_path,
                    )
                    self.assertEqual(physics_run.returncode, 0, physics_run.stderr)
                    aggregate_run = self.run_script(
                        "triage_table.py", "--sequence-csv", sequence_path,
                        "--physics-csv", physics_path, "--csv", output_path,
                    )
                    self.assertEqual(aggregate_run.returncode, 0, aggregate_run.stderr)
                    row = read_csv_by_id(output_path)["candidate_001"]
                    self.assertEqual(row["triage_status"], expected)
                    self.assertEqual(row["SAP"], sap)
                    if expected == "INCOMPLETE":
                        self.assertIn("rosetta_sap_score", row["triage_notes"])

    def test_duplicate_input_exits_with_error_before_writing_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            duplicate_path = directory / "duplicate.csv"
            physics_path = directory / "physics.csv"
            output_path = directory / "triage.csv"
            write_csv(duplicate_path, [stage_row(), stage_row()])
            write_csv(physics_path, [stage_row()])
            result = self.run_script(
                "triage_table.py", "--sequence-csv", duplicate_path,
                "--physics-csv", physics_path, "--csv", output_path,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("duplicate sequence_id", result.stderr)
            self.assertFalse(output_path.exists())


if __name__ == "__main__":
    unittest.main()
