#!/usr/bin/env python3
"""Stage 2 subcheck -- Rosetta score features via PyRosetta.

This is the runnable scorer counterpart to `rosetta_score.py`, which only
parses existing Rosetta output. It reads PDB files, scores them with the
PyRosetta full-atom scorefunction, and emits one CSV row per structure. It can
also write residue-level scores and rolling local score windows.

Usage
    conda run -n pyrosetta python stage2_fast_physics/pyrosetta_score.py model.pdb
    conda run -n pyrosetta python stage2_fast_physics/pyrosetta_score.py *.pdb --csv scores.csv
    conda run -n pyrosetta python stage2_fast_physics/pyrosetta_score.py *.pdb --residue-csv residue_scores.csv --window-csv local_scores.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import TextIO

try:
    from .pyrosetta_common import rolling_score_windows, score_pdb_residues
except ImportError:  # pragma: no cover - direct script execution
    from pyrosetta_common import rolling_score_windows, score_pdb_residues


DEFAULT_WINDOW_SIZE = 9

CSV_COLUMNS = [
    "structure_id",
    "source",
    "length",
    "rosetta_total_score",
    "rosetta_score_per_residue",
    "max_residue_score",
    "max_residue_score_position",
    "rolling_window_size",
    "max_window_score_per_residue",
    "max_window_start",
    "max_window_end",
    "flag_rosetta_score",
    "total_flags",
]

RESIDUE_CSV_COLUMNS = [
    "structure_id",
    "source",
    "residue_index",
    "chain",
    "pdb_residue",
    "aa",
    "residue_total_score",
]

WINDOW_CSV_COLUMNS = [
    "structure_id",
    "source",
    "window_start",
    "window_end",
    "window_size",
    "window_total_score",
    "window_score_per_residue",
]


def collect_score_tables(
    paths: list[Path],
    *,
    scorefxn_name: str = "ref2015",
    max_score_per_residue: float | None = None,
    window_size: int = DEFAULT_WINDOW_SIZE,
) -> tuple[list[dict], list[dict], list[dict]]:
    summary_rows: list[dict] = []
    residue_rows: list[dict] = []
    window_rows: list[dict] = []
    for path in paths:
        score, residues = score_pdb_residues(path, scorefxn_name=scorefxn_name)
        windows = rolling_score_windows(residues, window_size=window_size)
        score_per_residue = score.score_per_residue
        flag = int(
            max_score_per_residue is not None
            and score_per_residue > max_score_per_residue
        )
        max_residue = max(residues, key=lambda row: row.residue_total_score)
        max_window = max(windows, key=lambda row: row.window_score_per_residue)
        summary_rows.append({
            "structure_id": score.structure_id,
            "source": str(score.source),
            "length": score.length,
            "rosetta_total_score": round(score.total_score, 3),
            "rosetta_score_per_residue": round(score_per_residue, 3),
            "max_residue_score": round(max_residue.residue_total_score, 3),
            "max_residue_score_position": max_residue.position,
            "rolling_window_size": max_window.window_size,
            "max_window_score_per_residue": round(max_window.window_score_per_residue, 3),
            "max_window_start": max_window.window_start,
            "max_window_end": max_window.window_end,
            "flag_rosetta_score": flag,
            "total_flags": flag,
        })
        for residue in residues:
            residue_rows.append({
                "structure_id": residue.structure_id,
                "source": str(residue.source),
                "residue_index": residue.residue_index,
                "chain": residue.chain,
                "pdb_residue": residue.pdb_residue,
                "aa": residue.aa,
                "residue_total_score": round(residue.residue_total_score, 3),
            })
        for window in windows:
            window_rows.append({
                "structure_id": window.structure_id,
                "source": str(window.source),
                "window_start": window.window_start,
                "window_end": window.window_end,
                "window_size": window.window_size,
                "window_total_score": round(window.window_total_score, 3),
                "window_score_per_residue": round(window.window_score_per_residue, 3),
            })
    return summary_rows, residue_rows, window_rows


def summarize_pdbs(
    paths: list[Path],
    *,
    scorefxn_name: str = "ref2015",
    max_score_per_residue: float | None = None,
    window_size: int = DEFAULT_WINDOW_SIZE,
) -> list[dict]:
    rows, _residue_rows, _window_rows = collect_score_tables(
        paths,
        scorefxn_name=scorefxn_name,
        max_score_per_residue=max_score_per_residue,
        window_size=window_size,
    )
    return rows


def write_rows(rows: list[dict], fh: TextIO, *, columns: list[str] = CSV_COLUMNS) -> None:
    writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("pdbs", nargs="+", type=Path, help="PDB file(s) to score")
    parser.add_argument("--scorefxn", default="ref2015",
                        help="PyRosetta scorefunction name")
    parser.add_argument("--max-score-per-residue", type=float,
                        default=None,
                        help="optional WARN cutoff; disabled by default")
    parser.add_argument("--window-size", type=int, default=DEFAULT_WINDOW_SIZE,
                        help="rolling local-score window size in residues")
    parser.add_argument("--csv", type=Path, help="write CSV here instead of stdout")
    parser.add_argument("--residue-csv", type=Path,
                        help="write per-residue Rosetta scores here")
    parser.add_argument("--window-csv", type=Path,
                        help="write rolling-window Rosetta scores here")
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 if any row has one or more flags")
    args = parser.parse_args(argv)

    try:
        rows, residue_rows, window_rows = collect_score_tables(
            args.pdbs,
            scorefxn_name=args.scorefxn,
            max_score_per_residue=args.max_score_per_residue,
            window_size=args.window_size,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="") as fh:
            write_rows(rows, fh)
    else:
        write_rows(rows, sys.stdout)

    if args.residue_csv:
        args.residue_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.residue_csv.open("w", newline="") as fh:
            write_rows(residue_rows, fh, columns=RESIDUE_CSV_COLUMNS)

    if args.window_csv:
        args.window_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.window_csv.open("w", newline="") as fh:
            write_rows(window_rows, fh, columns=WINDOW_CSV_COLUMNS)

    if args.strict:
        return 1 if any(row["total_flags"] for row in rows) else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
