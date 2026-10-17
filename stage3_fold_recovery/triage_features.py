#!/usr/bin/env python3
"""Stage 3 triage feature table.

This converts `structure_recovery.py` output into the Stage 3 subset used by the
cross-step triage table.

Usage
    conda run -n petase-zeroshot python stage3_fold_recovery/triage_features.py --metrics-csv structure_recovery_metrics.csv --csv stage3_triage.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import TextIO


CSV_COLUMNS = [
    "sequence_id",
    "length",
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


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        required = {
            "sequence_id",
            "length",
            "pLDDT",
            "pLDDT_motif",
            "pTM",
            "global_CA_RMSD",
            "motif_CA_RMSD_global_align",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")
        return list(reader)


def summarize(rows: list[dict]) -> list[dict]:
    summarized = []
    for row in rows:
        summarized.append({
            "sequence_id": row.get("sequence_id", ""),
            "length": row.get("length", ""),
            "mean_pLDDT": row.get("pLDDT", ""),
            "motif_pLDDT": row.get("pLDDT_motif", ""),
            "pTM": row.get("pTM", ""),
            "global_CA_RMSD": row.get("global_CA_RMSD", ""),
            "TM_score": "",
            "motif_CA_RMSD": row.get("motif_CA_RMSD_global_align", ""),
            "motif_sidechain_RMSD": "",
            # Stage 3 metrics are reference ranges for now, not calibrated flags.
            "n_flags": 0,
            "triage_status": "PASS",
        })
    return summarized


def write_rows(rows: list[dict], fh: TextIO) -> None:
    writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--metrics-csv", type=Path, required=True)
    parser.add_argument("--csv", type=Path, help="write CSV here instead of stdout")
    args = parser.parse_args(argv)

    try:
        rows = summarize(read_rows(args.metrics_csv))
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
