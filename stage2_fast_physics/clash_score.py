#!/usr/bin/env python3
"""Stage 2 subcheck -- simple heavy-atom clash score.

The score is severe nonlocal heavy-atom overlaps per 1000 atoms, similar in
spirit to a MolProbity clashscore but deliberately simpler. Adjacent residues
are ignored so covalent peptide geometry is not counted as a clash.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import TextIO

try:
    from .pdb_geometry import nonlocal_residue_pair, nearby_atom_pairs, parse_pdb_atoms, residue_count
except ImportError:  # pragma: no cover
    from pdb_geometry import nonlocal_residue_pair, nearby_atom_pairs, parse_pdb_atoms, residue_count


DEFAULT_MIN_OVERLAP = 0.4
DEFAULT_MAX_CLASHSCORE = 30.0

CSV_COLUMNS = [
    "structure_id",
    "source",
    "length",
    "n_atoms",
    "n_clashes",
    "clashscore",
    "max_overlap",
    "flag_clash",
    "total_flags",
]


def clash_metrics(path: Path, *, min_overlap: float = DEFAULT_MIN_OVERLAP) -> dict:
    atoms = parse_pdb_atoms(path)
    max_pair_radius = max(a.radius for a in atoms) * 2.0
    cutoff = max_pair_radius - min_overlap
    n_clashes = 0
    max_overlap = 0.0
    for i, j, dist in nearby_atom_pairs(atoms, cutoff):
        a = atoms[i]
        b = atoms[j]
        if not nonlocal_residue_pair(a, b):
            continue
        overlap = a.radius + b.radius - dist
        if overlap >= min_overlap:
            n_clashes += 1
            max_overlap = max(max_overlap, overlap)
    n_atoms = len(atoms)
    return {
        "structure_id": path.stem,
        "source": str(path),
        "length": residue_count(atoms),
        "n_atoms": n_atoms,
        "n_clashes": n_clashes,
        "clashscore": 1000.0 * n_clashes / n_atoms if n_atoms else 0.0,
        "max_overlap": max_overlap,
    }


def summarize_pdbs(
    paths: list[Path],
    *,
    min_overlap: float = DEFAULT_MIN_OVERLAP,
    max_clashscore: float = DEFAULT_MAX_CLASHSCORE,
) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        row = clash_metrics(path, min_overlap=min_overlap)
        flag = int(row["clashscore"] > max_clashscore)
        row["clashscore"] = round(row["clashscore"], 3)
        row["max_overlap"] = round(row["max_overlap"], 3)
        row["flag_clash"] = flag
        row["total_flags"] = flag
        rows.append(row)
    return rows


def write_rows(rows: list[dict], fh: TextIO) -> None:
    writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdbs", nargs="+", type=Path)
    parser.add_argument("--min-overlap", type=float, default=DEFAULT_MIN_OVERLAP)
    parser.add_argument("--max-clashscore", type=float, default=DEFAULT_MAX_CLASHSCORE)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)

    try:
        rows = summarize_pdbs(
            args.pdbs,
            min_overlap=args.min_overlap,
            max_clashscore=args.max_clashscore,
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
