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

# Exposure is read off a structure, so on a predicted structure it is only as good as the
# prediction. A model the predictor got wrong is loose and reads as exposed whether or not
# the design is. Measured case: a design whose ESMFold model scored pLDDT 55 sat 10.6 A
# from its NMR structure while its own design model was 1.3 A from it -- the design was
# right and the prediction was wrong, but the exposure flag saw only the wrong structure.
# 70 is the conventional low-confidence line. pLDDT scales are not interchangeable between
# predictors or averaging conventions, so re-check this against whatever produced the models.
DEFAULT_MIN_PLDDT = 70.0

CSV_COLUMNS = [
    "structure_id",
    "source",
    "length",
    "exposed_hydrophobic_sasa",
    "exposed_hydrophobic_sasa_per_residue",
    "exposed_hydrophobic_fraction",
    "n_exposed_hydrophobic_residues",
    "max_hydrophobic_patch_sasa",
    "mean_plddt",
    "confidence_status",
    "flag_exposed_hydrophobicity",
    "total_flags",
]


def mean_plddt_from_bfactor(path: Path) -> float | None:
    """Mean of the B-factor column over CA atoms, or None when there are none.

    AlphaFold and ESMFold write per-residue pLDDT into the B-factor column, so confidence
    can be read straight from the model. **An experimental structure's B-factor column is
    a temperature factor, not a confidence**, and would be misread here; that is why the
    gate below is opt-in rather than automatic.
    """
    values: list[float] = []
    with path.open() as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                values.append(float(line[60:66]))
            except ValueError:
                continue
    return sum(values) / len(values) if values else None


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
    min_plddt: float | None = None,
) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        row = exposed_hydrophobic_metrics(
            path,
            min_exposed_sidechain_sasa=min_exposed_sidechain_sasa,
            patch_contact_cutoff=patch_contact_cutoff,
            sasa_points=sasa_points,
        )
        plddt = mean_plddt_from_bfactor(path) if min_plddt is not None else None
        if min_plddt is None:
            row["mean_plddt"] = ""
            row["confidence_status"] = "not_checked"
        elif plddt is None:
            row["mean_plddt"] = ""
            row["confidence_status"] = "no_confidence_in_file"
        else:
            row["mean_plddt"] = round(plddt, 2)
            row["confidence_status"] = (
                "trusted" if plddt >= min_plddt else "prediction_not_trusted")
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
        # A model the predictor was unsure about yields no verdict rather than a flag:
        # the metrics stay in the table as features, but they are not evidence here.
        if row["confidence_status"] == "prediction_not_trusted":
            flag = 0
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
    parser.add_argument("--min-plddt", type=float, nargs="?", const=DEFAULT_MIN_PLDDT,
                        default=None,
                        help="treat the B-factor column as pLDDT and withhold the flag "
                             f"below this value (bare flag uses {DEFAULT_MIN_PLDDT:g}); "
                             "only for predicted structures -- an experimental B-factor "
                             "is a temperature factor, not a confidence")
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
            min_plddt=args.min_plddt,
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
