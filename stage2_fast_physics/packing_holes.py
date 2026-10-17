#!/usr/bin/env python3
"""Stage 2 fallback -- packing/holes proxy.

This approximates poor core packing by looking for buried hydrophobic residues
whose side-chain atoms make too few nonlocal heavy-atom contacts.

This is not Rosetta PackStat. Use rosetta_official_metrics.py for the
Rosetta-native PackStat metric.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import TextIO

try:
    from .pdb_geometry import (
        HYDROPHOBIC_RESNAMES,
        atom_sasa,
        nearby_atom_pairs,
        nonlocal_residue_pair,
        parse_pdb_atoms,
        residue_count,
        residue_resnames,
        sidechain_capacity,
        sidechain_sasa,
    )
except ImportError:  # pragma: no cover
    from pdb_geometry import (
        HYDROPHOBIC_RESNAMES,
        atom_sasa,
        nearby_atom_pairs,
        nonlocal_residue_pair,
        parse_pdb_atoms,
        residue_count,
        residue_resnames,
        sidechain_capacity,
        sidechain_sasa,
    )


DEFAULT_BURIED_RELATIVE_SASA = 0.15
DEFAULT_CONTACT_CUTOFF = 4.5
DEFAULT_MIN_CONTACTS_PER_ATOM = 2.0

CSV_COLUMNS = [
    "structure_id",
    "source",
    "length",
    "n_buried_hydrophobic_residues",
    "n_underpacked_buried_hydrophobics",
    "underpacked_buried_hydrophobic_fraction",
    "buried_hydrophobic_contacts_per_atom",
    "flag_packing_holes",
    "total_flags",
]


def packing_metrics(
    path: Path,
    *,
    buried_relative_sasa: float = DEFAULT_BURIED_RELATIVE_SASA,
    contact_cutoff: float = DEFAULT_CONTACT_CUTOFF,
    min_contacts_per_atom: float = DEFAULT_MIN_CONTACTS_PER_ATOM,
    sasa_points: int = 64,
) -> dict:
    atoms = parse_pdb_atoms(path)
    sasa = atom_sasa(atoms, n_points=sasa_points)
    side_sasa = sidechain_sasa(atoms, sasa)
    side_capacity = sidechain_capacity(atoms)
    names = residue_resnames(atoms)
    length = residue_count(atoms)

    buried_hydrophobic = {
        rid for rid, name in names.items()
        if name in HYDROPHOBIC_RESNAMES
        and side_capacity.get(rid, 0.0) > 0.0
        and side_sasa.get(rid, 0.0) / side_capacity[rid] <= buried_relative_sasa
    }

    sidechain_atom_counts: dict[tuple[str, str, str], int] = defaultdict(int)
    contacts: dict[tuple[str, str, str], int] = defaultdict(int)
    for atom in atoms:
        if atom.is_sidechain and atom.residue_id in buried_hydrophobic:
            sidechain_atom_counts[atom.residue_id] += 1

    for i, j, _dist in nearby_atom_pairs(atoms, contact_cutoff):
        a = atoms[i]
        b = atoms[j]
        if not nonlocal_residue_pair(a, b):
            continue
        if a.residue_id in buried_hydrophobic and a.is_sidechain:
            contacts[a.residue_id] += 1
        if b.residue_id in buried_hydrophobic and b.is_sidechain:
            contacts[b.residue_id] += 1

    underpacked = 0
    total_contacts = 0
    total_atoms = 0
    for rid in buried_hydrophobic:
        n_atoms = sidechain_atom_counts.get(rid, 0)
        if n_atoms <= 0:
            continue
        contact_density = contacts.get(rid, 0) / n_atoms
        total_contacts += contacts.get(rid, 0)
        total_atoms += n_atoms
        if contact_density < min_contacts_per_atom:
            underpacked += 1

    n_buried = len(buried_hydrophobic)
    return {
        "structure_id": path.stem,
        "source": str(path),
        "length": length,
        "n_buried_hydrophobic_residues": n_buried,
        "n_underpacked_buried_hydrophobics": underpacked,
        "underpacked_buried_hydrophobic_fraction": underpacked / n_buried if n_buried else 0.0,
        "buried_hydrophobic_contacts_per_atom": total_contacts / total_atoms if total_atoms else 0.0,
    }


def summarize_pdbs(
    paths: list[Path],
    *,
    buried_relative_sasa: float = DEFAULT_BURIED_RELATIVE_SASA,
    contact_cutoff: float = DEFAULT_CONTACT_CUTOFF,
    min_contacts_per_atom: float = DEFAULT_MIN_CONTACTS_PER_ATOM,
    max_underpacked_fraction: float | None = None,
    sasa_points: int = 64,
) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        row = packing_metrics(
            path,
            buried_relative_sasa=buried_relative_sasa,
            contact_cutoff=contact_cutoff,
            min_contacts_per_atom=min_contacts_per_atom,
            sasa_points=sasa_points,
        )
        flag = int(
            max_underpacked_fraction is not None
            and
            row["n_buried_hydrophobic_residues"] >= 5
            and row["underpacked_buried_hydrophobic_fraction"] > max_underpacked_fraction
        )
        row["underpacked_buried_hydrophobic_fraction"] = round(
            row["underpacked_buried_hydrophobic_fraction"], 3
        )
        row["buried_hydrophobic_contacts_per_atom"] = round(
            row["buried_hydrophobic_contacts_per_atom"], 3
        )
        row["flag_packing_holes"] = flag
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
    parser.add_argument("--buried-relative-sasa", type=float,
                        default=DEFAULT_BURIED_RELATIVE_SASA)
    parser.add_argument("--contact-cutoff", type=float, default=DEFAULT_CONTACT_CUTOFF)
    parser.add_argument("--min-contacts-per-atom", type=float,
                        default=DEFAULT_MIN_CONTACTS_PER_ATOM)
    parser.add_argument("--max-underpacked-fraction", type=float,
                        default=None,
                        help="optional WARN cutoff; disabled by default")
    parser.add_argument("--sasa-points", type=int, default=64)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)

    try:
        rows = summarize_pdbs(
            args.pdbs,
            buried_relative_sasa=args.buried_relative_sasa,
            contact_cutoff=args.contact_cutoff,
            min_contacts_per_atom=args.min_contacts_per_atom,
            max_underpacked_fraction=args.max_underpacked_fraction,
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
