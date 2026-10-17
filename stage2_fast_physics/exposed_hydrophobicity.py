#!/usr/bin/env python3
"""Stage 2 fallback -- SAP-like exposed-hydrophobicity proxy.

This is a SAP-like static proxy. It sums exposed side-chain SASA on hydrophobic
residues and reports the largest connected exposed hydrophobic patch.

This is not the Rosetta SAP score. Use rosetta_official_metrics.py for the
Rosetta-native SAP metric.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import TextIO

try:
    from .pdb_geometry import (
        HYDROPHOBIC_RESNAMES,
        atom_sasa,
        hydrophobic_residue_components,
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
        hydrophobic_residue_components,
        parse_pdb_atoms,
        residue_count,
        residue_resnames,
        sidechain_capacity,
        sidechain_sasa,
    )


DEFAULT_MIN_EXPOSED_SIDECHAIN_SASA = 10.0
DEFAULT_PATCH_CONTACT_CUTOFF = 4.5

CSV_COLUMNS = [
    "structure_id",
    "source",
    "length",
    "exposed_hydrophobic_sasa",
    "exposed_hydrophobic_sasa_per_residue",
    "exposed_hydrophobic_fraction",
    "n_exposed_hydrophobic_residues",
    "max_hydrophobic_patch_sasa",
    "flag_exposed_hydrophobicity",
    "total_flags",
]


def exposed_hydrophobic_metrics(
    path: Path,
    *,
    min_exposed_sidechain_sasa: float = DEFAULT_MIN_EXPOSED_SIDECHAIN_SASA,
    patch_contact_cutoff: float = DEFAULT_PATCH_CONTACT_CUTOFF,
    sasa_points: int = 64,
) -> dict:
    atoms = parse_pdb_atoms(path)
    sasa = atom_sasa(atoms, n_points=sasa_points)
    side_sasa = sidechain_sasa(atoms, sasa)
    side_capacity = sidechain_capacity(atoms)
    names = residue_resnames(atoms)
    length = residue_count(atoms)

    exposed_hydrophobic = {
        rid for rid, value in side_sasa.items()
        if names.get(rid) in HYDROPHOBIC_RESNAMES
        and value >= min_exposed_sidechain_sasa
    }
    exposed_sasa = sum(side_sasa.get(rid, 0.0) for rid in exposed_hydrophobic)
    hydrophobic_capacity = sum(
        side_capacity.get(rid, 0.0)
        for rid, name in names.items()
        if name in HYDROPHOBIC_RESNAMES
    )

    components = hydrophobic_residue_components(
        atoms, exposed_hydrophobic, contact_cutoff=patch_contact_cutoff
    )
    max_patch = 0.0
    for component in components:
        max_patch = max(max_patch, sum(side_sasa.get(rid, 0.0) for rid in component))

    return {
        "structure_id": path.stem,
        "source": str(path),
        "length": length,
        "exposed_hydrophobic_sasa": exposed_sasa,
        "exposed_hydrophobic_sasa_per_residue": exposed_sasa / length if length else 0.0,
        "exposed_hydrophobic_fraction": (
            exposed_sasa / hydrophobic_capacity if hydrophobic_capacity else 0.0
        ),
        "n_exposed_hydrophobic_residues": len(exposed_hydrophobic),
        "max_hydrophobic_patch_sasa": max_patch,
    }


def summarize_pdbs(
    paths: list[Path],
    *,
    min_exposed_sidechain_sasa: float = DEFAULT_MIN_EXPOSED_SIDECHAIN_SASA,
    max_exposed_hydro_per_res: float | None = None,
    max_patch_sasa: float | None = None,
    patch_contact_cutoff: float = DEFAULT_PATCH_CONTACT_CUTOFF,
    sasa_points: int = 64,
) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        row = exposed_hydrophobic_metrics(
            path,
            min_exposed_sidechain_sasa=min_exposed_sidechain_sasa,
            patch_contact_cutoff=patch_contact_cutoff,
            sasa_points=sasa_points,
        )
        flag = int(
            (
                max_exposed_hydro_per_res is not None
                and row["exposed_hydrophobic_sasa_per_residue"] > max_exposed_hydro_per_res
            )
            or (
                max_patch_sasa is not None
                and row["max_hydrophobic_patch_sasa"] > max_patch_sasa
            )
        )
        for key in (
            "exposed_hydrophobic_sasa",
            "exposed_hydrophobic_sasa_per_residue",
            "exposed_hydrophobic_fraction",
            "max_hydrophobic_patch_sasa",
        ):
            row[key] = round(row[key], 3)
        row["flag_exposed_hydrophobicity"] = flag
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
    parser.add_argument("--min-exposed-sidechain-sasa", type=float,
                        default=DEFAULT_MIN_EXPOSED_SIDECHAIN_SASA)
    parser.add_argument("--max-exposed-hydro-per-res", type=float,
                        default=None,
                        help="optional WARN cutoff; disabled by default")
    parser.add_argument("--max-patch-sasa", type=float, default=None,
                        help="optional WARN cutoff; disabled by default")
    parser.add_argument("--patch-contact-cutoff", type=float,
                        default=DEFAULT_PATCH_CONTACT_CUTOFF)
    parser.add_argument("--sasa-points", type=int, default=64)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)

    try:
        rows = summarize_pdbs(
            args.pdbs,
            min_exposed_sidechain_sasa=args.min_exposed_sidechain_sasa,
            max_exposed_hydro_per_res=args.max_exposed_hydro_per_res,
            max_patch_sasa=args.max_patch_sasa,
            patch_contact_cutoff=args.patch_contact_cutoff,
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
