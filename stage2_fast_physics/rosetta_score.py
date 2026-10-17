#!/usr/bin/env python3
"""Stage 2 subcheck -- Rosetta score per residue.

This reads Rosetta score outputs that already exist. It does not run Rosetta;
that runner belongs in a later pass once the Rosetta executable and relaxation
protocol are fixed for the project.

Supported inputs:
  * Rosetta scorefiles with a `SCORE:` header containing `total_score`
  * Rosetta-scored PDBs containing a `#BEGIN_POSE_ENERGIES_TABLE` block

The output is one CSV row per scored structure. Rosetta score per residue is
treated as a continuous feature by default. Add `--max-score-per-residue` only
after calibrating a cutoff on a larger scored set.

Usage
    conda run -n base python stage2_fast_physics/rosetta_score.py scores.sc --length 120
    conda run -n base python stage2_fast_physics/rosetta_score.py scored_model.pdb
    conda run -n base python stage2_fast_physics/rosetta_score.py scores.sc --lengths-csv lengths.csv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


CSV_COLUMNS = [
    "structure_id",
    "source",
    "length",
    "rosetta_total_score",
    "rosetta_score_per_residue",
    "flag_rosetta_score",
    "total_flags",
]

AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}


@dataclass(frozen=True)
class ScoreRecord:
    structure_id: str
    source: Path
    total_score: float
    length: int | None = None


def count_protein_residues(pdb_path: Path) -> int:
    """Count unique standard amino-acid residues in ATOM records."""
    residues: set[tuple[str, str, str]] = set()
    with pdb_path.open() as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            resname = line[17:20].strip().upper()
            if resname not in AA3:
                continue
            chain = line[21].strip()
            resseq = line[22:26].strip()
            icode = line[26].strip()
            residues.add((chain, resseq, icode))
    return len(residues)


def parse_scorefile(path: Path) -> list[ScoreRecord]:
    """Parse a Rosetta `.sc` scorefile."""
    header: list[str] | None = None
    records: list[ScoreRecord] = []

    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line.startswith("SCORE:"):
                continue
            fields = line.split()[1:]
            if not fields:
                continue
            if "total_score" in fields and "description" in fields:
                header = fields
                continue
            if header is None:
                continue
            if len(fields) < len(header):
                continue
            row = dict(zip(header, fields))
            records.append(
                ScoreRecord(
                    structure_id=row["description"],
                    source=path,
                    total_score=float(row["total_score"]),
                )
            )

    if header is None:
        raise ValueError(f"{path} is not a Rosetta scorefile with total_score")
    if not records:
        raise ValueError(f"{path} has a Rosetta score header but no score rows")
    return records


def parse_rosetta_pdb_total_score(path: Path) -> float:
    """Read total score from a Rosetta PDB pose energies table."""
    in_table = False
    header: list[str] | None = None

    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if line == "#BEGIN_POSE_ENERGIES_TABLE":
                in_table = True
                continue
            if line == "#END_POSE_ENERGIES_TABLE":
                break
            if not in_table:
                continue
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "label":
                header = fields
                continue
            if fields[0] != "pose" or header is None:
                continue
            for total_key in ("total", "total_score"):
                if total_key in header:
                    idx = header.index(total_key)
                    if idx >= len(fields):
                        raise ValueError(f"{path} pose energy row lacks {total_key}")
                    return float(fields[idx])

    remark_score = _parse_remark_total_score(path)
    if remark_score is not None:
        return remark_score

    raise ValueError(f"{path} does not contain a Rosetta total score")


def _parse_remark_total_score(path: Path) -> float | None:
    pattern = re.compile(r"\btotal_score\b\s*:?\s*(-?\d+(?:\.\d+)?)")
    with path.open() as fh:
        for line in fh:
            if not line.startswith("REMARK"):
                continue
            match = pattern.search(line)
            if match:
                return float(match.group(1))
    return None


def parse_scored_pdb(path: Path) -> ScoreRecord:
    length = count_protein_residues(path)
    if length <= 0:
        raise ValueError(f"{path} has no standard protein residues")
    return ScoreRecord(
        structure_id=path.stem,
        source=path,
        total_score=parse_rosetta_pdb_total_score(path),
        length=length,
    )


def read_lengths_csv(path: Path) -> dict[str, int]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if "structure_id" not in (reader.fieldnames or []):
            raise ValueError(f"{path} must contain a structure_id column")
        if "length" not in (reader.fieldnames or []):
            raise ValueError(f"{path} must contain a length column")
        lengths = {}
        for row in reader:
            lengths[row["structure_id"]] = int(row["length"])
    return lengths


def read_inputs(paths: list[Path]) -> list[ScoreRecord]:
    records: list[ScoreRecord] = []
    for path in paths:
        suffix = path.suffix.lower()
        if suffix == ".pdb":
            records.append(parse_scored_pdb(path))
        else:
            records.extend(parse_scorefile(path))
    return records


def summarize_records(
    records: list[ScoreRecord],
    *,
    default_length: int | None = None,
    lengths: dict[str, int] | None = None,
    max_score_per_residue: float | None = None,
) -> list[dict]:
    rows: list[dict] = []
    length_map = lengths or {}
    for record in records:
        length = record.length or length_map.get(record.structure_id) or default_length
        if length is None:
            raise ValueError(
                f"no length available for {record.structure_id}; use --length "
                "or --lengths-csv"
            )
        if length <= 0:
            raise ValueError(f"length must be positive for {record.structure_id}")

        score_per_residue = record.total_score / length
        flag = int(
            max_score_per_residue is not None
            and score_per_residue > max_score_per_residue
        )
        rows.append({
            "structure_id": record.structure_id,
            "source": str(record.source),
            "length": length,
            "rosetta_total_score": round(record.total_score, 3),
            "rosetta_score_per_residue": round(score_per_residue, 3),
            "flag_rosetta_score": flag,
            "total_flags": flag,
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
    parser.add_argument("inputs", nargs="+", type=Path,
                        help="Rosetta scorefile(s) or Rosetta-scored PDB(s)")
    parser.add_argument("--length", type=int,
                        help="use this length for scorefile rows lacking length")
    parser.add_argument("--lengths-csv", type=Path,
                        help="CSV with structure_id,length columns")
    parser.add_argument("--max-score-per-residue", type=float,
                        default=None,
                        help="optional WARN cutoff; disabled by default")
    parser.add_argument("--csv", type=Path,
                        help="write CSV here instead of stdout")
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 if any row has one or more flags")
    args = parser.parse_args(argv)

    try:
        lengths = read_lengths_csv(args.lengths_csv) if args.lengths_csv else None
        rows = summarize_records(
            read_inputs(args.inputs),
            default_length=args.length,
            lengths=lengths,
            max_score_per_residue=args.max_score_per_residue,
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

    if args.strict:
        return 1 if any(row["total_flags"] for row in rows) else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
