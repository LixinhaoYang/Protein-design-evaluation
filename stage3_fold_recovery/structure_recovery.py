#!/usr/bin/env python3
"""Fold and motif recovery metrics for predicted structures.

This script compares a predicted/refolded structure against the designed
structure and emits one CSV row per design.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


BACKBONE_ATOMS = frozenset({"N", "CA", "C"})


@dataclass(frozen=True)
class Atom:
    chain_id: str
    residue_id: int
    atom_name: str
    coord: np.ndarray
    plddt: float | None = None


@dataclass(frozen=True)
class Case:
    sequence_id: str
    design_structure: Path
    predicted_structure: Path
    summary_confidences: Path
    motif_mapping: Path


def _as_float(text: str) -> float:
    return float(text) if text not in {".", "?"} else math.nan


def parse_pdb_atoms(path: Path) -> dict[tuple[str, int, str], Atom]:
    """Parse ATOM records from a PDB file."""
    atoms: dict[tuple[str, int, str], Atom] = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            chain_id = (line[21].strip() or "A")
            residue_id = int(line[22:26])
            coord = np.array(
                [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ],
                dtype=float,
            )
            key = (chain_id, residue_id, atom_name)
            atoms[key] = Atom(chain_id, residue_id, atom_name, coord)
    return atoms


def parse_af3_cif_atoms(path: Path) -> dict[tuple[str, int, str], Atom]:
    """Parse simple AlphaFold 3 mmCIF ATOM rows.

    The AF3 output here is whitespace-tokenizable; this parser intentionally
    targets that output shape rather than a fully general mmCIF grammar.
    """
    atoms: dict[tuple[str, int, str], Atom] = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM ", "HETATM ")):
                continue
            parts = line.split()
            if len(parts) < 18:
                continue
            atom_name = parts[3]
            chain_id = parts[6]
            residue_id = int(parts[8])
            coord = np.array(
                [
                    _as_float(parts[10]),
                    _as_float(parts[11]),
                    _as_float(parts[12]),
                ],
                dtype=float,
            )
            plddt = _as_float(parts[14])
            key = (chain_id, residue_id, atom_name)
            atoms[key] = Atom(chain_id, residue_id, atom_name, coord, plddt)
    return atoms


def load_motif_positions(path: Path) -> list[int]:
    positions: set[int] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if "candidate_position_1based" not in reader.fieldnames:
            raise ValueError(
                f"{path} must contain candidate_position_1based column"
            )
        for row in reader:
            positions.add(int(row["candidate_position_1based"]))
    return sorted(positions)


def ca_by_position(
    atoms: dict[tuple[str, int, str], Atom], chain_id: str = "A"
) -> dict[int, Atom]:
    return {
        residue_id: atom
        for (chain, residue_id, atom_name), atom in atoms.items()
        if chain == chain_id and atom_name == "CA"
    }


def coords_for_positions(
    atoms_by_position: dict[int, Atom], positions: Iterable[int]
) -> np.ndarray:
    return np.array([atoms_by_position[pos].coord for pos in positions], dtype=float)


def plddts_for_positions(
    atoms_by_position: dict[int, Atom], positions: Iterable[int]
) -> list[float]:
    return [
        atoms_by_position[pos].plddt
        for pos in positions
        if atoms_by_position[pos].plddt is not None
        and not math.isnan(atoms_by_position[pos].plddt)
    ]


def mean(values: Iterable[float]) -> float:
    values = list(values)
    if not values:
        return math.nan
    return float(sum(values) / len(values))


def superpose_mobile_to_target(
    mobile: np.ndarray, target: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return aligned mobile coordinates plus rotation and translation."""
    if mobile.shape != target.shape:
        raise ValueError(f"shape mismatch: {mobile.shape} vs {target.shape}")
    if mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError(f"expected Nx3 coordinates, got {mobile.shape}")

    mobile_centroid = mobile.mean(axis=0)
    target_centroid = target.mean(axis=0)
    mobile_centered = mobile - mobile_centroid
    target_centered = target - target_centroid

    covariance = mobile_centered.T @ target_centered
    u, _, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt

    translation = target_centroid - mobile_centroid @ rotation
    aligned = mobile @ rotation + translation
    return aligned, rotation, translation


def rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    diff = coords_a - coords_b
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def transformed(coords: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return coords @ rotation + translation


def read_cases(path: Path) -> list[Case]:
    cases: list[Case] = []
    base = path.parent
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "sequence_id",
            "design_structure",
            "predicted_structure",
            "summary_confidences",
            "motif_mapping",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            cases.append(
                Case(
                    sequence_id=row["sequence_id"],
                    design_structure=(base / row["design_structure"]).resolve(),
                    predicted_structure=(base / row["predicted_structure"]).resolve(),
                    summary_confidences=(base / row["summary_confidences"]).resolve(),
                    motif_mapping=(base / row["motif_mapping"]).resolve(),
                )
            )
    return cases


def evaluate_case(case: Case) -> dict[str, object]:
    design_atoms = parse_pdb_atoms(case.design_structure)
    predicted_atoms = parse_af3_cif_atoms(case.predicted_structure)
    design_ca = ca_by_position(design_atoms)
    predicted_ca = ca_by_position(predicted_atoms)
    motif_positions_all = load_motif_positions(case.motif_mapping)

    common_positions = sorted(set(design_ca) & set(predicted_ca))
    if len(common_positions) < 3:
        raise ValueError(f"{case.sequence_id}: fewer than 3 matched CA atoms")

    motif_positions = [
        pos for pos in motif_positions_all if pos in design_ca and pos in predicted_ca
    ]
    if len(motif_positions) < 3:
        raise ValueError(f"{case.sequence_id}: fewer than 3 matched motif CA atoms")

    design_global = coords_for_positions(design_ca, common_positions)
    predicted_global = coords_for_positions(predicted_ca, common_positions)
    predicted_global_aligned, rotation, translation = superpose_mobile_to_target(
        predicted_global, design_global
    )

    design_motif = coords_for_positions(design_ca, motif_positions)
    predicted_motif = coords_for_positions(predicted_ca, motif_positions)
    predicted_motif_global = transformed(predicted_motif, rotation, translation)
    predicted_motif_local, _, _ = superpose_mobile_to_target(
        predicted_motif, design_motif
    )

    summary = json.loads(case.summary_confidences.read_text())
    ca_plddts = plddts_for_positions(predicted_ca, common_positions)
    motif_plddts = plddts_for_positions(predicted_ca, motif_positions)

    return {
        "sequence_id": case.sequence_id,
        "length": len(common_positions),
        "n_ca_aligned": len(common_positions),
        "n_motif_ca": len(motif_positions),
        "pLDDT": f"{mean(ca_plddts):.3f}",
        "pLDDT_motif": f"{mean(motif_plddts):.3f}",
        "pTM": f"{float(summary.get('ptm', math.nan)):.3f}",
        "ranking_score": f"{float(summary.get('ranking_score', math.nan)):.3f}",
        "global_CA_RMSD": f"{rmsd(predicted_global_aligned, design_global):.3f}",
        "motif_CA_RMSD_global_align": f"{rmsd(predicted_motif_global, design_motif):.3f}",
        "motif_CA_RMSD_local_align": f"{rmsd(predicted_motif_local, design_motif):.3f}",
        "design_structure": str(case.design_structure),
        "predicted_structure": str(case.predicted_structure),
        "motif_mapping": str(case.motif_mapping),
    }


def write_csv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sequence_id",
        "length",
        "n_ca_aligned",
        "n_motif_ca",
        "pLDDT",
        "pLDDT_motif",
        "pTM",
        "ranking_score",
        "global_CA_RMSD",
        "motif_CA_RMSD_global_align",
        "motif_CA_RMSD_local_align",
        "design_structure",
        "predicted_structure",
        "motif_mapping",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-csv", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path, help="Output metrics CSV")
    args = parser.parse_args(argv)

    rows = [evaluate_case(case) for case in read_cases(args.case_csv)]
    write_csv(rows, args.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
