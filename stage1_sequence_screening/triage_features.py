#!/usr/bin/env python3
"""Step 1 triage feature table.

This emits the sequence-sanity subset used by the cross-step triage table:

`sequence_id`, `length`, `min_entropy_12`, `repeat_score`, `max_KD19`,
`min_KD19`, `q_norm`, `n_flags`, `triage_status`.

Usage
    python stage1_sequence_screening/triage_features.py designs.fasta --csv step1_triage.csv
    python stage1_sequence_screening/triage_features.py --seq MKVLAA...
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import TextIO

try:
    from . import basic_validity, charge_pi, hydrophobicity, low_complexity
    from .basic_validity import Record, STANDARD_AA, _WHITESPACE, read_fasta
except ImportError:  # pragma: no cover - direct script execution
    import basic_validity
    import charge_pi
    import hydrophobicity
    import low_complexity
    from basic_validity import Record, STANDARD_AA, _WHITESPACE, read_fasta


CSV_COLUMNS = [
    "sequence_id",
    "length",
    "min_entropy_12",
    "repeat_score",
    "max_KD19",
    "min_KD19",
    "q_norm",
    "n_flags",
    "triage_status",
]


def triage_sequence_id(seq_id: str) -> str:
    return seq_id.split("|", 1)[0]


def metric_valid(record: Record, result: basic_validity.Result) -> bool:
    return result.ok and set(record.sequence) <= STANDARD_AA


def status_from_flags(*, invalid: bool, n_flags: int) -> str:
    if invalid:
        return "FAIL"
    if n_flags:
        return "WARN"
    return "PASS"


def summarize_records(records: list[Record]) -> list[dict]:
    validity_results = basic_validity.validate_all(records)
    rows: list[dict] = []
    for record, valid_result in zip(records, validity_results):
        invalid = not metric_valid(record, valid_result)
        if invalid:
            rows.append({
                "sequence_id": triage_sequence_id(record.seq_id),
                "length": len(record.sequence),
                "min_entropy_12": "",
                "repeat_score": "",
                "max_KD19": "",
                "min_KD19": "",
                "q_norm": "",
                "n_flags": 1,
                "triage_status": "FAIL",
            })
            continue

        lc = low_complexity.analyse(record)
        hp = hydrophobicity.analyse(record)
        ch = charge_pi.analyse(record)
        n_flags = int(bool(lc.issues)) + int(bool(hp.issues)) + int(bool(ch.issues))
        rows.append({
            "sequence_id": triage_sequence_id(record.seq_id),
            "length": len(record.sequence),
            "min_entropy_12": lc.metrics[f"min_entropy_w{low_complexity.DEFAULT_WINDOW}"],
            "repeat_score": lc.metrics["max_repeat_fraction"],
            "max_KD19": hp.metrics[f"max_hydropathy_w{hydrophobicity.DEFAULT_WINDOW}"],
            "min_KD19": hp.metrics[f"min_hydropathy_w{hydrophobicity.DEFAULT_WINDOW}"],
            "q_norm": ch.metrics["charge_density"],
            "n_flags": n_flags,
            "triage_status": status_from_flags(invalid=False, n_flags=n_flags),
        })
    return rows


def write_rows(rows: list[dict], fh: TextIO) -> None:
    writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)


def read_inputs(paths: list[Path], seq: str | None) -> list[Record]:
    records: list[Record] = []
    if seq:
        records.append(Record(seq_id="<--seq>", sequence=_WHITESPACE.sub("", seq).upper()))
    for path in paths:
        records.extend(read_fasta(path))
    return records


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("fasta", nargs="*", type=Path, help="FASTA file(s)")
    parser.add_argument("--seq", help="one sequence on the command line")
    parser.add_argument("--csv", type=Path, help="write CSV here instead of stdout")
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 if any row is WARN or FAIL")
    args = parser.parse_args(argv)

    if not args.fasta and not args.seq:
        parser.error("give at least one FASTA file, or --seq")

    try:
        rows = summarize_records(read_inputs(args.fasta, args.seq))
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="") as fh:
            write_rows(rows, fh)
    else:
        write_rows(rows, sys.stdout)

    if args.strict:
        return 1 if any(row["triage_status"] != "PASS" for row in rows) else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
