#!/usr/bin/env python3
"""Stage 2 subcheck -- Rosetta/PyRosetta official structure metrics.

Use this on restrained-relaxed PDBs when PyRosetta is available. The geometric
proxy scripts remain useful as fallbacks and correlation checks, but this script
is the preferred source for Rosetta-native buried unsats, SAP and PackStat.

Usage
    conda run -n pyrosetta python stage2_fast_physics/rosetta_official_metrics.py relaxed_pdbs/*.pdb --csv rosetta_official_metrics.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Callable, TextIO

try:
    from .pyrosetta_common import get_scorefxn, load_pose, score_pose
except ImportError:  # pragma: no cover - direct script execution
    from pyrosetta_common import get_scorefxn, load_pose, score_pose


CSV_COLUMNS = [
    "structure_id",
    "source",
    "length",
    "rosetta_total_score",
    "rosetta_score_per_residue",
    "fa_rep_raw",
    "fa_rep_weighted",
    "fa_rep_weighted_per_residue",
    "rosetta_buried_unsat_hbonds",
    "rosetta_buried_unsat_hbonds_per_100_res",
    "rosetta_sap_score",
    "rosetta_sap_score_per_residue",
    "rosetta_max_residue_sap",
    "rosetta_max_residue_sap_position",
    "rosetta_packstat",
    "metric_errors",
    "total_flags",
]


def blank_or_round(value: float | None, digits: int = 3) -> str | float:
    return "" if value is None else round(value, digits)


def residue_position(pose, residue_index: int) -> str:
    pdb_info = pose.pdb_info()
    if pdb_info is None:
        return str(residue_index)
    chain = str(pdb_info.chain(residue_index)).strip() or "_"
    number = str(pdb_info.number(residue_index)).strip()
    insertion_code = ""
    if hasattr(pdb_info, "icode"):
        insertion_code = str(pdb_info.icode(residue_index)).strip()
    return f"{chain}:{number}{insertion_code}"


def safe_metric(name: str, fn: Callable[[], float | tuple]) -> tuple[float | tuple | None, str | None]:
    try:
        return fn(), None
    except Exception as exc:  # pragma: no cover - depends on Rosetta internals
        return None, f"{name}:{type(exc).__name__}:{exc}"


def fa_rep_terms(pose, scorefxn) -> tuple[float, float]:
    pyrosetta = sys.modules["pyrosetta"]
    score_type = pyrosetta.rosetta.core.scoring.fa_rep
    raw = float(pose.energies().total_energies()[score_type])
    weighted = raw * float(scorefxn.get_weight(score_type))
    return raw, weighted


def buried_unsat_hbonds(pose) -> float:
    pyrosetta = sys.modules["pyrosetta"]
    filt = pyrosetta.rosetta.protocols.simple_filters.BuriedUnsatHbondFilter()
    filt.set_report_all_heavy_atom_unsats(True)
    return float(filt.compute(pose))


def sap_metrics(pose) -> tuple[float, float, str]:
    pyrosetta = sys.modules["pyrosetta"]
    sap = pyrosetta.rosetta.core.pack.guidance_scoreterms.sap
    selector = pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()
    options = sap.SapParameterOptions()
    total_sap = float(sap.calculate_sap(pose, selector, selector, selector, options))
    per_res_sap = sap.calculate_per_res_sap(pose, selector, selector, selector, options)
    values = [float(per_res_sap[idx]) for idx in range(1, len(per_res_sap) + 1)]
    if not values:
        return total_sap, 0.0, ""
    max_offset, max_value = max(enumerate(values, start=1), key=lambda item: item[1])
    return total_sap, max_value, residue_position(pose, max_offset)


def packstat(pose) -> float:
    pyrosetta = sys.modules["pyrosetta"]
    filt = pyrosetta.rosetta.protocols.simple_filters.PackStatFilter()
    return float(filt.compute(pose))


def summarize_pdbs(
    paths: list[Path],
    *,
    scorefxn_name: str = "ref2015",
    strict_metrics: bool = False,
) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        pose = load_pose(path)
        length = pose.total_residue()
        if length <= 0:
            raise ValueError(f"{path} has no residues after PyRosetta import")

        total_score = score_pose(pose, scorefxn_name=scorefxn_name)
        scorefxn = get_scorefxn(scorefxn_name)
        errors: list[str] = []

        fa_rep_raw, fa_rep_weighted = fa_rep_terms(pose, scorefxn)

        buried_unsats, err = safe_metric("buried_unsat_hbonds", lambda: buried_unsat_hbonds(pose))
        if err:
            errors.append(err)

        sap_result, err = safe_metric("sap", lambda: sap_metrics(pose))
        if err:
            errors.append(err)
            sap_score = max_residue_sap = None
            max_residue_sap_position = ""
        else:
            sap_score, max_residue_sap, max_residue_sap_position = sap_result

        packstat_result, err = safe_metric("packstat", lambda: packstat(pose))
        if err:
            errors.append(err)

        if strict_metrics and errors:
            raise RuntimeError("; ".join(errors))

        rows.append({
            "structure_id": path.stem,
            "source": str(path),
            "length": length,
            "rosetta_total_score": round(total_score, 3),
            "rosetta_score_per_residue": round(total_score / length, 3),
            "fa_rep_raw": round(fa_rep_raw, 3),
            "fa_rep_weighted": round(fa_rep_weighted, 3),
            "fa_rep_weighted_per_residue": round(fa_rep_weighted / length, 3),
            "rosetta_buried_unsat_hbonds": blank_or_round(buried_unsats),
            "rosetta_buried_unsat_hbonds_per_100_res": (
                blank_or_round(100.0 * buried_unsats / length)
                if buried_unsats is not None
                else ""
            ),
            "rosetta_sap_score": blank_or_round(sap_score),
            "rosetta_sap_score_per_residue": (
                blank_or_round(sap_score / length)
                if sap_score is not None
                else ""
            ),
            "rosetta_max_residue_sap": blank_or_round(max_residue_sap),
            "rosetta_max_residue_sap_position": max_residue_sap_position,
            "rosetta_packstat": blank_or_round(packstat_result),
            "metric_errors": "; ".join(errors),
            "total_flags": 0,
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
    parser.add_argument("pdbs", nargs="+", type=Path, help="PDB file(s) to score")
    parser.add_argument("--scorefxn", default="ref2015",
                        help="PyRosetta scorefunction name")
    parser.add_argument("--strict-metrics", action="store_true",
                        help="exit 2 if any official metric fails")
    parser.add_argument("--csv", type=Path, help="write CSV here instead of stdout")
    args = parser.parse_args(argv)

    try:
        rows = summarize_pdbs(
            args.pdbs,
            scorefxn_name=args.scorefxn,
            strict_metrics=args.strict_metrics,
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

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
