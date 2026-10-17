"""PyRosetta helpers for Stage 2 scripts."""

from __future__ import annotations

import sys
from contextlib import redirect_stdout
from dataclasses import dataclass
from math import sqrt
from pathlib import Path
from typing import Sequence


@dataclass(frozen=True)
class PyRosettaScore:
    structure_id: str
    source: Path
    length: int
    total_score: float

    @property
    def score_per_residue(self) -> float:
        return self.total_score / self.length


@dataclass(frozen=True)
class ResidueScore:
    structure_id: str
    source: Path
    residue_index: int
    chain: str
    pdb_residue: str
    aa: str
    residue_total_score: float

    @property
    def position(self) -> str:
        chain = self.chain or "_"
        return f"{chain}:{self.pdb_residue}"


@dataclass(frozen=True)
class WindowScore:
    structure_id: str
    source: Path
    window_start: int
    window_end: int
    window_size: int
    window_total_score: float

    @property
    def window_score_per_residue(self) -> float:
        return self.window_total_score / self.window_size


@dataclass(frozen=True)
class RelaxResult:
    structure_id: str
    source: Path
    relaxed_pdb: Path
    length: int
    pre_relax_score: float
    post_relax_score: float
    relax_CA_RMSD: float
    relax_all_atom_RMSD: float
    max_CA_displacement: float
    relax_motif_RMSD: float | None
    motif_residue_count: int
    motif_residues: str
    coord_constraint_weight: float
    max_iter: int

    @property
    def pre_relax_score_per_residue(self) -> float:
        return self.pre_relax_score / self.length

    @property
    def post_relax_score_per_residue(self) -> float:
        return self.post_relax_score / self.length

    @property
    def score_delta(self) -> float:
        return self.post_relax_score - self.pre_relax_score


_PYROSETTA = None
_SCOREFXN = None


def require_pyrosetta():
    try:
        import pyrosetta
    except ImportError as exc:  # pragma: no cover - depends on local env
        raise RuntimeError(
            "PyRosetta is not importable. Run this command in an environment "
            "that provides a licensed, configured PyRosetta installation."
        ) from exc
    return pyrosetta


def init_pyrosetta(extra_options: str = "-mute all"):
    global _PYROSETTA
    pyrosetta = require_pyrosetta()
    if _PYROSETTA is None:
        # PyRosetta prints a license banner to stdout during init. CLI scripts
        # emit CSV on stdout, so keep the banner visible but move it to stderr.
        with redirect_stdout(sys.stderr):
            pyrosetta.init(extra_options)
        _PYROSETTA = pyrosetta
    return _PYROSETTA


def get_scorefxn(name: str = "ref2015"):
    global _SCOREFXN
    pyrosetta = init_pyrosetta()
    if name != "ref2015":
        return pyrosetta.create_score_function(name)
    if _SCOREFXN is None:
        _SCOREFXN = pyrosetta.get_fa_scorefxn()
    return _SCOREFXN


def load_pose(path: Path):
    pyrosetta = init_pyrosetta()
    return pyrosetta.pose_from_pdb(str(path))


def score_pose(pose, *, scorefxn_name: str = "ref2015") -> float:
    scorefxn = get_scorefxn(scorefxn_name)
    return float(scorefxn(pose))


def score_pdb(path: Path, *, scorefxn_name: str = "ref2015") -> PyRosettaScore:
    pose = load_pose(path)
    total_score = score_pose(pose, scorefxn_name=scorefxn_name)
    length = pose.total_residue()
    if length <= 0:
        raise ValueError(f"{path} has no residues after PyRosetta import")
    return PyRosettaScore(
        structure_id=path.stem,
        source=path,
        length=length,
        total_score=total_score,
    )


def _pdb_residue_label(pose, residue_index: int) -> tuple[str, str]:
    pdb_info = pose.pdb_info()
    if pdb_info is None:
        return "", str(residue_index)
    chain = str(pdb_info.chain(residue_index)).strip()
    number = str(pdb_info.number(residue_index)).strip()
    insertion_code = ""
    if hasattr(pdb_info, "icode"):
        insertion_code = str(pdb_info.icode(residue_index)).strip()
    return chain, f"{number}{insertion_code}"


def _xyz_tuple(xyz) -> tuple[float, float, float]:
    def component(name: str) -> float:
        value = getattr(xyz, name)
        return float(value() if callable(value) else value)

    return component("x"), component("y"), component("z")


def _distance(a, b) -> float:
    ax, ay, az = _xyz_tuple(a)
    bx, by, bz = _xyz_tuple(b)
    return sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)


def _rmsd_from_distances(distances: Sequence[float]) -> float:
    if not distances:
        return 0.0
    return sqrt(sum(d * d for d in distances) / len(distances))


def ca_residue_indices(pose) -> list[int]:
    return [
        idx for idx in range(1, pose.total_residue() + 1)
        if pose.residue(idx).has("CA")
    ]


def ca_displacements(before_pose, after_pose, residue_indices: Sequence[int]) -> list[float]:
    distances: list[float] = []
    for idx in residue_indices:
        before = before_pose.residue(idx)
        after = after_pose.residue(idx)
        if not before.has("CA") or not after.has("CA"):
            continue
        distances.append(_distance(before.xyz("CA"), after.xyz("CA")))
    return distances


def ca_rmsd(before_pose, after_pose, residue_indices: Sequence[int]) -> float:
    return _rmsd_from_distances(ca_displacements(before_pose, after_pose, residue_indices))


def all_atom_rmsd(before_pose, after_pose) -> float:
    distances: list[float] = []
    for idx in range(1, before_pose.total_residue() + 1):
        before = before_pose.residue(idx)
        after = after_pose.residue(idx)
        for atom_idx in range(1, before.natoms() + 1):
            atom_name = before.atom_name(atom_idx).strip()
            if not after.has(atom_name):
                continue
            distances.append(
                _distance(before.xyz(atom_idx), after.xyz(after.atom_index(atom_name)))
            )
    return _rmsd_from_distances(distances)


def motif_pose_indices(pose, motif_residues: str | None) -> list[int]:
    """Map a motif spec to pose residue indices.

    Accepted tokens:
      * pose indices: `10,15,20` or `10-15`
      * PDB residue labels: `A:160,A:206,A:237` or `A:160-165`
    """
    if motif_residues is None or not motif_residues.strip():
        return []

    label_to_index = {}
    for idx in range(1, pose.total_residue() + 1):
        chain, pdb_residue = _pdb_residue_label(pose, idx)
        label_to_index[(chain, pdb_residue)] = idx

    indices: list[int] = []
    seen: set[int] = set()
    tokens = [
        token.strip()
        for token in motif_residues.replace(";", ",").split(",")
        if token.strip()
    ]

    def add_index(idx: int, token: str) -> None:
        if idx < 1 or idx > pose.total_residue():
            raise ValueError(
                f"motif residue {token!r} maps outside pose length {pose.total_residue()}"
            )
        if idx not in seen:
            indices.append(idx)
            seen.add(idx)

    for token in tokens:
        if ":" in token:
            chain, residue_text = token.split(":", 1)
            chain = chain.strip()
            residue_text = residue_text.strip()
            if not chain or not residue_text:
                raise ValueError(f"invalid motif residue token {token!r}")
            if "-" in residue_text:
                start_text, end_text = residue_text.split("-", 1)
                start = int(start_text)
                end = int(end_text)
                step = 1 if end >= start else -1
                for resnum in range(start, end + step, step):
                    key = (chain, str(resnum))
                    if key not in label_to_index:
                        raise ValueError(f"motif residue {chain}:{resnum} not found")
                    add_index(label_to_index[key], f"{chain}:{resnum}")
            else:
                key = (chain, residue_text)
                if key not in label_to_index:
                    raise ValueError(f"motif residue {token!r} not found")
                add_index(label_to_index[key], token)
        elif "-" in token:
            start_text, end_text = token.split("-", 1)
            start = int(start_text)
            end = int(end_text)
            step = 1 if end >= start else -1
            for idx in range(start, end + step, step):
                add_index(idx, str(idx))
        else:
            add_index(int(token), token)

    return indices


def score_pdb_residues(
    path: Path,
    *,
    scorefxn_name: str = "ref2015",
) -> tuple[PyRosettaScore, list[ResidueScore]]:
    pose = load_pose(path)
    total_score = score_pose(pose, scorefxn_name=scorefxn_name)
    length = pose.total_residue()
    if length <= 0:
        raise ValueError(f"{path} has no residues after PyRosetta import")

    residues: list[ResidueScore] = []
    for idx in range(1, length + 1):
        chain, pdb_residue = _pdb_residue_label(pose, idx)
        residues.append(
            ResidueScore(
                structure_id=path.stem,
                source=path,
                residue_index=idx,
                chain=chain,
                pdb_residue=pdb_residue,
                aa=str(pose.residue(idx).name1()),
                residue_total_score=float(pose.energies().residue_total_energy(idx)),
            )
        )
    score = PyRosettaScore(
        structure_id=path.stem,
        source=path,
        length=length,
        total_score=total_score,
    )
    return score, residues


def rolling_score_windows(
    residue_scores: Sequence[ResidueScore],
    *,
    window_size: int = 9,
) -> list[WindowScore]:
    if window_size <= 0:
        raise ValueError("window_size must be positive")
    if not residue_scores:
        return []
    effective_size = min(window_size, len(residue_scores))
    windows: list[WindowScore] = []
    for start in range(0, len(residue_scores) - effective_size + 1):
        window = residue_scores[start:start + effective_size]
        windows.append(
            WindowScore(
                structure_id=window[0].structure_id,
                source=window[0].source,
                window_start=window[0].residue_index,
                window_end=window[-1].residue_index,
                window_size=effective_size,
                window_total_score=sum(r.residue_total_score for r in window),
            )
        )
    return windows


def restrained_relax_pdb(
    path: Path,
    output_path: Path,
    *,
    scorefxn_name: str = "ref2015",
    coord_constraint_weight: float = 1.0,
    max_iter: int = 200,
    motif_residues: str | None = None,
) -> RelaxResult:
    if max_iter <= 0:
        raise ValueError("max_iter must be positive")

    pyrosetta = init_pyrosetta()
    pose = load_pose(path)
    length = pose.total_residue()
    if length <= 0:
        raise ValueError(f"{path} has no residues after PyRosetta import")
    pre_relax_pose = pose.clone()
    motif_indices = motif_pose_indices(pre_relax_pose, motif_residues)

    base_scorefxn = get_scorefxn(scorefxn_name)
    pre_score = float(base_scorefxn(pose))

    relax_scorefxn = base_scorefxn.clone()
    relax_scorefxn.set_weight(
        pyrosetta.rosetta.core.scoring.coordinate_constraint,
        coord_constraint_weight,
    )

    relax = pyrosetta.rosetta.protocols.relax.FastRelax()
    relax.set_scorefxn(relax_scorefxn)
    relax.constrain_relax_to_start_coords(True)
    relax.coord_constrain_sidechains(True)
    relax.ramp_down_constraints(False)
    relax.max_iter(max_iter)
    relax.apply(pose)

    post_score = float(base_scorefxn(pose))
    ca_indices = ca_residue_indices(pre_relax_pose)
    ca_distances = ca_displacements(pre_relax_pose, pose, ca_indices)
    motif_rmsd = ca_rmsd(pre_relax_pose, pose, motif_indices) if motif_indices else None
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pose.dump_pdb(str(output_path))

    return RelaxResult(
        structure_id=path.stem,
        source=path,
        relaxed_pdb=output_path,
        length=length,
        pre_relax_score=pre_score,
        post_relax_score=post_score,
        relax_CA_RMSD=_rmsd_from_distances(ca_distances),
        relax_all_atom_RMSD=all_atom_rmsd(pre_relax_pose, pose),
        max_CA_displacement=max(ca_distances) if ca_distances else 0.0,
        relax_motif_RMSD=motif_rmsd,
        motif_residue_count=len(motif_indices),
        motif_residues=motif_residues or "",
        coord_constraint_weight=coord_constraint_weight,
        max_iter=max_iter,
    )


def pose_sequence(pose) -> str:
    return pose.sequence()
