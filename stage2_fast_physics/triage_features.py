#!/usr/bin/env python3
"""Stage 2 triage feature table.

This joins `restrained_relax.py` output and `rosetta_official_metrics.py`
output into the Stage 2 subset used by the cross-step triage table.
Unmatched rows or unavailable required metrics produce INCOMPLETE, with details
in triage_notes. Existing relaxation failures retain FAIL status.

Usage
    python stage2_fast_physics/triage_features.py --relax-csv relax.csv --rosetta-csv rosetta_official_metrics.csv --csv stage2_triage.csv
    python stage2_fast_physics/triage_features.py --relax-csv relax.csv --rosetta-csv rosetta_official_metrics.csv --id-map-csv candidate_metrics.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import TextIO


CSV_COLUMNS = [
    "sequence_id",
    "length",
    "pre_relax_clashscore",
    "post_relax_clashscore",
    "relax_CA_RMSD",
    "relax_motif_RMSD",
    "max_CA_displacement",
    "rosetta_score_per_residue",
    "fa_rep",
    "buried_unsats",
    "SAP",
    "PackStat",
    "n_flags",
    "triage_status",
    "triage_notes",
]


RELAX_METRICS = (
    "length",
    "pre_relax_clashscore",
    "post_relax_clashscore",
    "relax_CA_RMSD",
    "max_CA_displacement",
    "total_flags",
    "warn_pre_relax_clash",
)

ROSETTA_METRICS = (
    "rosetta_score_per_residue",
    "fa_rep_weighted",
    "rosetta_buried_unsat_hbonds",
    "rosetta_sap_score",
    "rosetta_packstat",
)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def read_id_maps(paths: list[Path]) -> dict[str, str]:
    """Return structure_id -> sequence_id.

    Supported schemas:
      * structure_id,sequence_id
      * gen_pdb_name,entry from example `candidate_metrics.csv`
    """
    mapping: dict[str, str] = {}
    for path in paths:
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            fields = set(reader.fieldnames or [])
            if {"structure_id", "sequence_id"} <= fields:
                for row in reader:
                    if row.get("structure_id") and row.get("sequence_id"):
                        mapping[row["structure_id"]] = row["sequence_id"]
            elif {"gen_pdb_name", "entry"} <= fields:
                for row in reader:
                    if row.get("gen_pdb_name") and row.get("entry"):
                        mapping[row["gen_pdb_name"]] = row["entry"]
            else:
                raise ValueError(
                    f"{path} must contain either structure_id,sequence_id or "
                    "gen_pdb_name,entry columns"
                )
    return mapping


def strip_restrained_relax_suffix(structure_id: str) -> str:
    suffix = "_restrained_relax"
    return structure_id[:-len(suffix)] if structure_id.endswith(suffix) else structure_id


def normalize_source(path_text: str) -> str:
    if not path_text:
        return ""
    return str(Path(path_text))


def rosetta_rows_by_relaxed_pdb(rows: list[dict]) -> dict[str, dict]:
    by_source = {}
    for row in rows:
        source = normalize_source(row.get("source", ""))
        if source:
            if source in by_source:
                raise ValueError(f"duplicate Rosetta source {source!r}")
            by_source[source] = row
    return by_source


def rosetta_rows_by_structure_id(rows: list[dict]) -> dict[str, dict]:
    by_id = {}
    for row in rows:
        structure_id = strip_restrained_relax_suffix(row.get("structure_id", ""))
        if structure_id:
            if structure_id in by_id:
                raise ValueError(f"duplicate Rosetta structure_id {structure_id!r}")
            by_id[structure_id] = row
    return by_id


def missing_numeric_metrics(row: dict, names: tuple[str, ...]) -> list[str]:
    missing = []
    for name in names:
        try:
            available = math.isfinite(float(row.get(name, "")))
        except (TypeError, ValueError):
            available = False
        if not available:
            missing.append(name)
    return missing


def int_value(value: str | int | None) -> int:
    if value in (None, ""):
        return 0
    return int(float(value))


def status_from_flags(
    *, n_fail_flags: int, warn_pre_relax_clash: int, incomplete: bool = False,
) -> str:
    if n_fail_flags:
        return "FAIL"
    if incomplete:
        return "INCOMPLETE"
    if warn_pre_relax_clash:
        return "WARN"
    return "PASS"


def summarize(
    *,
    relax_rows: list[dict],
    rosetta_rows: list[dict],
    id_map: dict[str, str] | None = None,
) -> list[dict]:
    id_map = id_map or {}
    rosetta_by_source = rosetta_rows_by_relaxed_pdb(rosetta_rows)
    rosetta_by_id = rosetta_rows_by_structure_id(rosetta_rows)
    rows: list[dict] = []

    for relax in relax_rows:
        structure_id = relax.get("structure_id", "")
        if not structure_id or not structure_id.strip():
            raise ValueError("relaxation row has an empty structure_id")
        sequence_id = id_map.get(structure_id, structure_id)
        relaxed_pdb = normalize_source(relax.get("relaxed_pdb", ""))
        rosetta = rosetta_by_source.get(relaxed_pdb) or rosetta_by_id.get(structure_id, {})
        notes: list[str] = []
        missing_relax = missing_numeric_metrics(relax, RELAX_METRICS)
        if missing_relax:
            notes.append(f"missing or invalid relaxation metrics: {', '.join(missing_relax)}")
        if not rosetta:
            notes.append("missing Rosetta row")
        else:
            missing_rosetta = missing_numeric_metrics(rosetta, ROSETTA_METRICS)
            if missing_rosetta:
                notes.append(f"missing or invalid Rosetta metrics: {', '.join(missing_rosetta)}")
        n_fail_flags = (
            int_value(relax.get("total_flags"))
            if "total_flags" not in missing_relax else 0
        )
        warn_pre = (
            int_value(relax.get("warn_pre_relax_clash"))
            if "warn_pre_relax_clash" not in missing_relax else 0
        )
        n_flags = n_fail_flags + warn_pre
        rows.append({
            "sequence_id": sequence_id,
            "length": relax.get("length") or rosetta.get("length", ""),
            "pre_relax_clashscore": relax.get("pre_relax_clashscore", ""),
            "post_relax_clashscore": relax.get("post_relax_clashscore", ""),
            "relax_CA_RMSD": relax.get("relax_CA_RMSD", ""),
            "relax_motif_RMSD": relax.get("relax_motif_RMSD", ""),
            "max_CA_displacement": relax.get("max_CA_displacement", ""),
            "rosetta_score_per_residue": rosetta.get("rosetta_score_per_residue", ""),
            "fa_rep": rosetta.get("fa_rep_weighted", ""),
            "buried_unsats": rosetta.get("rosetta_buried_unsat_hbonds", ""),
            "SAP": rosetta.get("rosetta_sap_score", ""),
            "PackStat": rosetta.get("rosetta_packstat", ""),
            "n_flags": n_flags,
            "triage_status": status_from_flags(
                n_fail_flags=n_fail_flags,
                warn_pre_relax_clash=warn_pre,
                incomplete=bool(notes),
            ),
            "triage_notes": "; ".join(notes),
        })
    return rows


def write_rows(rows: list[dict], fh: TextIO) -> None:
    writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--relax-csv", type=Path, required=True)
    parser.add_argument("--rosetta-csv", type=Path, required=True)
    parser.add_argument("--id-map-csv", type=Path, action="append", default=[],
                        help="CSV mapping structure IDs to sequence IDs; may be repeated")
    parser.add_argument("--csv", type=Path, help="write CSV here instead of stdout")
    args = parser.parse_args(argv)

    try:
        rows = summarize(
            relax_rows=read_csv(args.relax_csv),
            rosetta_rows=read_csv(args.rosetta_csv),
            id_map=read_id_maps(args.id_map_csv),
        )
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="") as fh:
            write_rows(rows, fh)
    else:
        write_rows(rows, sys.stdout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
