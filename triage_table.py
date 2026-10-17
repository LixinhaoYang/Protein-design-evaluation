#!/usr/bin/env python3
"""Join sequential triage features into a final table.

The execution order is sequence screening (Stage 1), fast physical screening
(Stage 2), then optional AF3 fold/motif recovery (Stage 3). Stage 3 is omitted
when building the pre-AF3 triage table.

Usage
    python triage_table.py --sequence-csv stage1.csv --physics-csv stage2.csv \\
        --fold-csv stage3.csv --csv triage.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import TextIO


FINAL_COLUMNS = [
    "sequence_id",
    "length",
    "min_entropy_12",
    "repeat_score",
    "max_KD19",
    "min_KD19",
    "q_norm",
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
    "mean_pLDDT",
    "motif_pLDDT",
    "pTM",
    "global_CA_RMSD",
    "TM_score",
    "motif_CA_RMSD",
    "motif_sidechain_RMSD",
    "n_flags",
    "triage_status",
]


def read_csv_by_id(path: Path) -> dict[str, dict]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if "sequence_id" not in (reader.fieldnames or []):
            raise ValueError(f"{path} must contain a sequence_id column")
        return {
            row["sequence_id"]: row
            for row in reader
            if row.get("sequence_id")
        }


def numeric_int(value: str | int | None) -> int:
    if value in (None, ""):
        return 0
    return int(float(value))


def combine_status(*statuses: str) -> str:
    statuses = set(statuses)
    if "FAIL" in statuses:
        return "FAIL"
    if "WARN" in statuses:
        return "WARN"
    return "PASS"


def build_rows(
    sequence: dict[str, dict],
    physics: dict[str, dict],
    fold: dict[str, dict] | None,
) -> list[dict]:
    fold = fold or {}
    ids = sorted(set(sequence) | set(physics) | set(fold))
    rows: list[dict] = []
    for seq_id in ids:
        stage1 = sequence.get(seq_id, {})
        stage2 = physics.get(seq_id, {})
        stage3 = fold.get(seq_id, {})
        n_flags = sum(
            numeric_int(stage.get("n_flags"))
            for stage in (stage1, stage2, stage3)
        )
        row = {column: "" for column in FINAL_COLUMNS}
        row.update({
            "sequence_id": seq_id,
            "length": (
                stage1.get("length")
                or stage2.get("length")
                or stage3.get("length", "")
            ),
            "min_entropy_12": stage1.get("min_entropy_12", ""),
            "repeat_score": stage1.get("repeat_score", ""),
            "max_KD19": stage1.get("max_KD19", ""),
            "min_KD19": stage1.get("min_KD19", ""),
            "q_norm": stage1.get("q_norm", ""),
            "pre_relax_clashscore": stage2.get("pre_relax_clashscore", ""),
            "post_relax_clashscore": stage2.get("post_relax_clashscore", ""),
            "relax_CA_RMSD": stage2.get("relax_CA_RMSD", ""),
            "relax_motif_RMSD": stage2.get("relax_motif_RMSD", ""),
            "max_CA_displacement": stage2.get("max_CA_displacement", ""),
            "rosetta_score_per_residue": stage2.get("rosetta_score_per_residue", ""),
            "fa_rep": stage2.get("fa_rep", ""),
            "buried_unsats": stage2.get("buried_unsats", ""),
            "SAP": stage2.get("SAP", ""),
            "PackStat": stage2.get("PackStat", ""),
            "mean_pLDDT": stage3.get("mean_pLDDT", ""),
            "motif_pLDDT": stage3.get("motif_pLDDT", ""),
            "pTM": stage3.get("pTM", ""),
            "global_CA_RMSD": stage3.get("global_CA_RMSD", ""),
            "TM_score": stage3.get("TM_score", ""),
            "motif_CA_RMSD": stage3.get("motif_CA_RMSD", ""),
            "motif_sidechain_RMSD": stage3.get("motif_sidechain_RMSD", ""),
            "n_flags": n_flags,
            "triage_status": combine_status(
                stage1.get("triage_status", "PASS"),
                stage2.get("triage_status", "PASS"),
                stage3.get("triage_status", "PASS"),
            ),
        })
        rows.append(row)
    return rows


def write_rows(rows: list[dict], fh: TextIO) -> None:
    writer = csv.DictWriter(fh, fieldnames=FINAL_COLUMNS, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--sequence-csv", type=Path, required=True)
    parser.add_argument("--physics-csv", type=Path, required=True)
    parser.add_argument("--fold-csv", type=Path)
    parser.add_argument("--csv", type=Path, help="write CSV here instead of stdout")
    args = parser.parse_args(argv)

    try:
        rows = build_rows(
            read_csv_by_id(args.sequence_csv),
            read_csv_by_id(args.physics_csv),
            read_csv_by_id(args.fold_csv) if args.fold_csv else None,
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
