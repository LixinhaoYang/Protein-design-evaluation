#!/usr/bin/env python3
"""Step 1 wrapper -- one-row sequence sanity summary.

This is the thin orchestration layer over the current sequence-only modules:

  * basic_validity.py
  * aa_composition.py
  * low_complexity.py
  * hydrophobicity.py
  * charge_pi.py
  * cys_pattern.py

It emits one CSV row per sequence, suitable for quick ranking/filtering. The
underlying modules keep their detailed reports; this wrapper keeps only the
small set of metrics and conservative flags needed for a first-pass table.

Usage
    python stage1_sequence_screening/run.py designs.fasta
    python stage1_sequence_screening/run.py --seq MKVLAA...
    python stage1_sequence_screening/run.py designs.fasta --csv sanity.csv
    python stage1_sequence_screening/run.py designs.fasta --summary
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import TextIO

try:
    from . import (
        aa_composition,
        basic_validity,
        charge_pi,
        cys_pattern,
        hydrophobicity,
        low_complexity,
    )
    from .basic_validity import Record, STANDARD_AA, _WHITESPACE, read_fasta
except ImportError:  # pragma: no cover - direct script execution path
    import aa_composition
    import basic_validity
    import charge_pi
    import cys_pattern
    import hydrophobicity
    import low_complexity
    from basic_validity import Record, STANDARD_AA, _WHITESPACE, read_fasta


DETAIL_COLUMNS = [
    "sequence_id",
    "length",
    "valid",
    "n_invalid",
    "global_entropy",
    "min_local_entropy",
    "longest_homopolymer",
    "gravy",
    "max_local_hydropathy",
    "net_charge_simple",
    "normalized_charge",
    "pI",
    "n_cys",
    "cys_fraction",
    "hydrophobic_fraction",
    "charged_fraction",
    "flag_invalid",
    "flag_low_complexity",
    "flag_hydrophobic_stretch",
    "flag_extreme_charge",
    "flag_cys_pattern",
    "total_flags",
]

SUMMARY_COLUMNS = ["sequence_id", "total_flags"]


def _invalid_position_count(result: basic_validity.Result) -> int:
    """Count invalid residue positions.

    Non-positional hard failures such as length mismatch or duplicate ID still
    make `valid` false through `flag_invalid`, but they are not invalid residues
    and therefore do not add to `n_invalid`.
    """
    return sum(
        len(issue.positions)
        for issue in result.issues
        if issue.severity == "FAIL"
    )


def _metric_valid(record: Record, result: basic_validity.Result) -> bool:
    """True when downstream sequence-only metrics can be computed honestly."""
    return result.ok and set(record.sequence) <= STANDARD_AA


def _blank_metrics(row: dict) -> None:
    for key in (
        "global_entropy",
        "min_local_entropy",
        "longest_homopolymer",
        "gravy",
        "max_local_hydropathy",
        "net_charge_simple",
        "normalized_charge",
        "pI",
        "n_cys",
        "cys_fraction",
        "hydrophobic_fraction",
        "charged_fraction",
    ):
        row[key] = ""


def summarize_records(records: list[Record]) -> list[dict]:
    """Return one flat CSV-ready dict per sequence."""
    validity_results = basic_validity.validate_all(records)
    rows: list[dict] = []

    for record, valid_result in zip(records, validity_results):
        metric_valid = _metric_valid(record, valid_result)
        row = {
            "sequence_id": record.seq_id,
            "length": len(record.sequence),
            "valid": int(metric_valid),
            "n_invalid": _invalid_position_count(valid_result),
            "flag_invalid": int(not metric_valid),
            "flag_low_complexity": 0,
            "flag_hydrophobic_stretch": 0,
            "flag_extreme_charge": 0,
            "flag_cys_pattern": 0,
        }

        if not metric_valid:
            _blank_metrics(row)
            row["total_flags"] = row["flag_invalid"]
            rows.append(row)
            continue

        comp = aa_composition.analyse(record)
        lc = low_complexity.analyse(record)
        hp = hydrophobicity.analyse(record)
        ch = charge_pi.analyse(record)
        cy = cys_pattern.analyse(record)

        row.update({
            "global_entropy": lc.metrics["entropy_whole"],
            "min_local_entropy": lc.metrics[f"min_entropy_w{low_complexity.DEFAULT_WINDOW}"],
            "longest_homopolymer": lc.metrics["longest_homopolymer"],
            "gravy": hp.metrics["gravy"],
            "max_local_hydropathy": hp.metrics[
                f"max_hydropathy_w{hydrophobicity.DEFAULT_WINDOW}"
            ],
            "net_charge_simple": ch.metrics["net_charge_simple"],
            "normalized_charge": ch.metrics["charge_density"],
            "pI": "" if ch.metrics["isoelectric_point"] is None
                  else ch.metrics["isoelectric_point"],
            "n_cys": cy.metrics["n_cys"],
            "cys_fraction": cy.metrics["cys_fraction"],
            "hydrophobic_fraction": round(comp.groups["hydrophobic"], 4),
            "charged_fraction": round(comp.groups["charged"], 4),
            "flag_low_complexity": int(bool(lc.issues)),
            "flag_hydrophobic_stretch": int(
                any(i.check == "tm_like_hydrophobic_stretch" for i in hp.issues)
            ),
            "flag_extreme_charge": int(
                any(i.check == "charge_density_outlier" for i in ch.issues)
            ),
            "flag_cys_pattern": int(bool(cy.issues)),
        })
        row["total_flags"] = sum(
            row[key] for key in (
                "flag_invalid",
                "flag_low_complexity",
                "flag_hydrophobic_stretch",
                "flag_extreme_charge",
                "flag_cys_pattern",
            )
        )
        rows.append(row)

    return rows


def write_rows(rows: list[dict], fh: TextIO, summary: bool = False) -> None:
    columns = SUMMARY_COLUMNS if summary else DETAIL_COLUMNS
    writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)


def _read_inputs(paths: list[Path], seq: str | None) -> list[Record]:
    records: list[Record] = []
    if seq:
        records.append(
            Record(seq_id="<--seq>", sequence=_WHITESPACE.sub("", seq).upper())
        )
    for path in paths:
        records.extend(read_fasta(path))
    return records


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("fasta", nargs="*", type=Path, help="FASTA file(s)")
    parser.add_argument("--seq", help="check one sequence given on the command line")
    parser.add_argument("--csv", type=Path,
                        help="write CSV here instead of stdout")
    parser.add_argument("--summary", action="store_true",
                        help="emit only sequence_id,total_flags")
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 if any row has one or more flags")
    args = parser.parse_args(argv)

    if not args.fasta and not args.seq:
        parser.error("give at least one FASTA file, or --seq")

    try:
        records = _read_inputs(args.fasta, args.seq)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    rows = summarize_records(records)

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with open(args.csv, "w", newline="") as fh:
            write_rows(rows, fh, summary=args.summary)
    else:
        write_rows(rows, sys.stdout, summary=args.summary)

    if args.strict:
        return 1 if any(row["total_flags"] for row in rows) else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
