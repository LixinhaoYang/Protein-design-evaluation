#!/usr/bin/env python3
"""Stage 2 pre-processing -- restrained Rosetta relax via PyRosetta.

Raw ESMFold/AF2 structures can contain local strain that dominates Rosetta and
geometry metrics. This script performs a coordinate-constrained FastRelax pass
and writes relaxed PDBs for downstream scoring and static geometry screens.

Usage
    conda run -n pyrosetta python stage2_fast_physics/restrained_relax.py model.pdb --outdir relaxed_pdbs --csv relax.csv
    conda run -n pyrosetta python stage2_fast_physics/restrained_relax.py *.pdb --outdir relaxed_pdbs
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import TextIO

try:
    from .clash_score import DEFAULT_MAX_CLASHSCORE, DEFAULT_MIN_OVERLAP, clash_metrics
    from .pyrosetta_common import restrained_relax_pdb
except ImportError:  # pragma: no cover - direct script execution
    from clash_score import DEFAULT_MAX_CLASHSCORE, DEFAULT_MIN_OVERLAP, clash_metrics
    from pyrosetta_common import restrained_relax_pdb


DEFAULT_MAX_RELAX_MOTIF_RMSD = 1.0

CSV_COLUMNS = [
    "structure_id",
    "source",
    "relaxed_pdb",
    "length",
    "pre_relax_score",
    "post_relax_score",
    "pre_relax_score_per_residue",
    "post_relax_score_per_residue",
    "score_delta",
    "relax_CA_RMSD",
    "relax_all_atom_RMSD",
    "relax_motif_RMSD",
    "motif_residue_count",
    "motif_residues",
    "max_CA_displacement",
    "pre_relax_n_clashes",
    "pre_relax_clashscore",
    "pre_relax_max_overlap",
    "post_relax_n_clashes",
    "post_relax_clashscore",
    "post_relax_max_overlap",
    "warn_pre_relax_clash",
    "flag_post_relax_clash",
    "flag_relax_deformation",
    "total_flags",
    "coord_constraint_weight",
    "max_iter",
]


def output_path_for(path: Path, outdir: Path) -> Path:
    return outdir / f"{path.stem}_restrained_relax.pdb"


def read_motif_residue_csv(path: Path) -> dict[str, str]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        fields = set(reader.fieldnames or [])
        if "structure_id" not in fields or "motif_residues" not in fields:
            raise ValueError(f"{path} must contain structure_id,motif_residues columns")
        return {
            row["structure_id"]: row.get("motif_residues", "")
            for row in reader
            if row.get("structure_id")
        }


def metric_or_blank(value: float | None) -> str | float:
    return "" if value is None else round(value, 3)


def deformation_flag(
    *,
    relax_ca_rmsd: float,
    relax_motif_rmsd: float | None,
    max_ca_displacement: float,
    max_relax_ca_rmsd: float | None,
    max_relax_motif_rmsd: float | None,
    max_ca_displacement_cutoff: float | None,
) -> int:
    return int(
        (
            max_relax_ca_rmsd is not None
            and relax_ca_rmsd > max_relax_ca_rmsd
        )
        or (
            max_relax_motif_rmsd is not None
            and relax_motif_rmsd is not None
            and relax_motif_rmsd > max_relax_motif_rmsd
        )
        or (
            max_ca_displacement_cutoff is not None
            and max_ca_displacement > max_ca_displacement_cutoff
        )
    )


def relax_pdbs(
    paths: list[Path],
    *,
    outdir: Path,
    scorefxn_name: str = "ref2015",
    coord_constraint_weight: float = 1.0,
    max_iter: int = 200,
    motif_residues: str | None = None,
    motif_residue_map: dict[str, str] | None = None,
    min_overlap: float = DEFAULT_MIN_OVERLAP,
    max_post_relax_clashscore: float = DEFAULT_MAX_CLASHSCORE,
    max_relax_ca_rmsd: float | None = None,
    max_relax_motif_rmsd: float | None = DEFAULT_MAX_RELAX_MOTIF_RMSD,
    max_ca_displacement: float | None = None,
    overwrite: bool = False,
) -> list[dict]:
    rows: list[dict] = []
    motif_residue_map = motif_residue_map or {}
    for path in paths:
        relaxed_pdb = output_path_for(path, outdir)
        if relaxed_pdb.exists() and not overwrite:
            raise ValueError(
                f"{relaxed_pdb} already exists; use --overwrite to replace it"
            )
        pre_clash = clash_metrics(path, min_overlap=min_overlap)
        structure_motif_residues = motif_residue_map.get(path.stem, motif_residues)
        result = restrained_relax_pdb(
            path,
            relaxed_pdb,
            scorefxn_name=scorefxn_name,
            coord_constraint_weight=coord_constraint_weight,
            max_iter=max_iter,
            motif_residues=structure_motif_residues,
        )
        post_clash = clash_metrics(relaxed_pdb, min_overlap=min_overlap)
        warn_pre_relax_clash = int(pre_clash["clashscore"] > max_post_relax_clashscore)
        flag_post_relax_clash = int(post_clash["clashscore"] > max_post_relax_clashscore)
        flag_relax_deformation = deformation_flag(
            relax_ca_rmsd=result.relax_CA_RMSD,
            relax_motif_rmsd=result.relax_motif_RMSD,
            max_ca_displacement=result.max_CA_displacement,
            max_relax_ca_rmsd=max_relax_ca_rmsd,
            max_relax_motif_rmsd=max_relax_motif_rmsd,
            max_ca_displacement_cutoff=max_ca_displacement,
        )
        rows.append({
            "structure_id": result.structure_id,
            "source": str(result.source),
            "relaxed_pdb": str(result.relaxed_pdb),
            "length": result.length,
            "pre_relax_score": round(result.pre_relax_score, 3),
            "post_relax_score": round(result.post_relax_score, 3),
            "pre_relax_score_per_residue": round(result.pre_relax_score_per_residue, 3),
            "post_relax_score_per_residue": round(result.post_relax_score_per_residue, 3),
            "score_delta": round(result.score_delta, 3),
            "relax_CA_RMSD": round(result.relax_CA_RMSD, 3),
            "relax_all_atom_RMSD": round(result.relax_all_atom_RMSD, 3),
            "relax_motif_RMSD": metric_or_blank(result.relax_motif_RMSD),
            "motif_residue_count": result.motif_residue_count,
            "motif_residues": result.motif_residues,
            "max_CA_displacement": round(result.max_CA_displacement, 3),
            "pre_relax_n_clashes": pre_clash["n_clashes"],
            "pre_relax_clashscore": round(pre_clash["clashscore"], 3),
            "pre_relax_max_overlap": round(pre_clash["max_overlap"], 3),
            "post_relax_n_clashes": post_clash["n_clashes"],
            "post_relax_clashscore": round(post_clash["clashscore"], 3),
            "post_relax_max_overlap": round(post_clash["max_overlap"], 3),
            "warn_pre_relax_clash": warn_pre_relax_clash,
            "flag_post_relax_clash": flag_post_relax_clash,
            "flag_relax_deformation": flag_relax_deformation,
            "total_flags": flag_post_relax_clash + flag_relax_deformation,
            "coord_constraint_weight": result.coord_constraint_weight,
            "max_iter": result.max_iter,
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
    parser.add_argument("pdbs", nargs="+", type=Path, help="PDB file(s) to relax")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="directory for relaxed PDB files")
    parser.add_argument("--scorefxn", default="ref2015",
                        help="PyRosetta scorefunction name")
    parser.add_argument("--coord-cst-weight", type=float, default=1.0,
                        help="coordinate constraint weight during relax")
    parser.add_argument("--max-iter", type=int, default=200,
                        help="FastRelax max_iter value")
    parser.add_argument("--motif-residues",
                        help="common motif spec, e.g. A:160,A:206,A:237 or 10-15")
    parser.add_argument("--motif-residues-csv", type=Path,
                        help="CSV with structure_id,motif_residues columns")
    parser.add_argument("--min-overlap", type=float, default=DEFAULT_MIN_OVERLAP,
                        help="heavy-atom overlap cutoff used for pre/post clash diagnostics")
    parser.add_argument("--max-post-relax-clashscore", type=float,
                        default=DEFAULT_MAX_CLASHSCORE,
                        help="hard FAIL cutoff for post-relax clashscore")
    parser.add_argument("--max-relax-ca-rmsd", type=float, default=None,
                        help="optional deformation cutoff for global CA RMSD")
    parser.add_argument("--max-relax-motif-rmsd", type=float,
                        default=DEFAULT_MAX_RELAX_MOTIF_RMSD,
                        help="deformation cutoff for motif CA RMSD when motif is provided")
    parser.add_argument("--max-ca-displacement", type=float, default=None,
                        help="optional deformation cutoff for maximum single-residue CA displacement")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace existing relaxed PDBs")
    parser.add_argument("--csv", type=Path, help="write CSV here instead of stdout")
    args = parser.parse_args(argv)

    try:
        motif_residue_map = (
            read_motif_residue_csv(args.motif_residues_csv)
            if args.motif_residues_csv
            else None
        )
        rows = relax_pdbs(
            args.pdbs,
            outdir=args.outdir,
            scorefxn_name=args.scorefxn,
            coord_constraint_weight=args.coord_cst_weight,
            max_iter=args.max_iter,
            motif_residues=args.motif_residues,
            motif_residue_map=motif_residue_map,
            min_overlap=args.min_overlap,
            max_post_relax_clashscore=args.max_post_relax_clashscore,
            max_relax_ca_rmsd=args.max_relax_ca_rmsd,
            max_relax_motif_rmsd=args.max_relax_motif_rmsd,
            max_ca_displacement=args.max_ca_displacement,
            overwrite=args.overwrite,
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
