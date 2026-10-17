#!/usr/bin/env python3
"""Stage 2 fallback -- geometric buried-unsat proxy.

This is a fast geometric proxy: a polar heavy atom is counted as buried when
its atom SASA is below a small cutoff, and unsatisfied when no other nonlocal
polar heavy atom sits within the hydrogen-bond distance cutoff.

This is not Rosetta BuriedUnsatHbondFilter. Use rosetta_official_metrics.py for
the Rosetta-native buried unsatisfied H-bond metric.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import TextIO

try:
    from .pdb_geometry import atom_sasa, distance, nonlocal_residue_pair, parse_pdb_atoms, residue_count
except ImportError:  # pragma: no cover
    from pdb_geometry import atom_sasa, distance, nonlocal_residue_pair, parse_pdb_atoms, residue_count


DEFAULT_BURIED_SASA = 1.0
DEFAULT_HBOND_DISTANCE = 3.5

CSV_COLUMNS = [
    "structure_id",
    "source",
    "length",
    "n_polar_atoms",
    "n_buried_polar_atoms",
    "n_buried_unsat_atoms",
    "buried_unsats_per_100_res",
    "flag_buried_unsats",
    "total_flags",
]


def has_polar_partner(atoms, idx: int, *, hbond_distance: float) -> bool:
    atom = atoms[idx]
    for other_idx, other in enumerate(atoms):
        if other_idx == idx or not other.is_polar_heavy:
            continue
        if not nonlocal_residue_pair(atom, other, skip_adjacent=False):
            continue
        if distance(atom, other) <= hbond_distance:
            return True
    return False


def unsat_metrics(
    path: Path,
    *,
    buried_sasa: float = DEFAULT_BURIED_SASA,
    hbond_distance: float = DEFAULT_HBOND_DISTANCE,
    sasa_points: int = 64,
) -> dict:
    atoms = parse_pdb_atoms(path)
    sasa = atom_sasa(atoms, n_points=sasa_points)
    polar_indices = [idx for idx, atom in enumerate(atoms) if atom.is_polar_heavy]
    buried = [idx for idx in polar_indices if sasa[idx] <= buried_sasa]
    unsat = [
        idx for idx in buried
        if not has_polar_partner(atoms, idx, hbond_distance=hbond_distance)
    ]
    length = residue_count(atoms)
    return {
        "structure_id": path.stem,
        "source": str(path),
        "length": length,
        "n_polar_atoms": len(polar_indices),
        "n_buried_polar_atoms": len(buried),
        "n_buried_unsat_atoms": len(unsat),
        "buried_unsats_per_100_res": 100.0 * len(unsat) / length if length else 0.0,
    }


def summarize_pdbs(
    paths: list[Path],
    *,
    buried_sasa: float = DEFAULT_BURIED_SASA,
    hbond_distance: float = DEFAULT_HBOND_DISTANCE,
    max_unsats_per_100: float | None = None,
    sasa_points: int = 64,
) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        row = unsat_metrics(
            path,
            buried_sasa=buried_sasa,
            hbond_distance=hbond_distance,
            sasa_points=sasa_points,
        )
        flag = int(
            max_unsats_per_100 is not None
            and row["buried_unsats_per_100_res"] > max_unsats_per_100
        )
        row["buried_unsats_per_100_res"] = round(row["buried_unsats_per_100_res"], 3)
        row["flag_buried_unsats"] = flag
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
    parser.add_argument("--buried-sasa", type=float, default=DEFAULT_BURIED_SASA)
    parser.add_argument("--hbond-distance", type=float, default=DEFAULT_HBOND_DISTANCE)
    parser.add_argument("--max-unsats-per-100", type=float,
                        default=None,
                        help="optional WARN cutoff; disabled by default")
    parser.add_argument("--sasa-points", type=int, default=64)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)

    try:
        rows = summarize_pdbs(
            args.pdbs,
            buried_sasa=args.buried_sasa,
            hbond_distance=args.hbond_distance,
            max_unsats_per_100=args.max_unsats_per_100,
            sasa_points=args.sasa_points,
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
