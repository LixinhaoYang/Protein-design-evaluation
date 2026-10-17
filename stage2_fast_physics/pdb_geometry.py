"""Small PDB geometry helpers for Stage 2 fast-physics screens.

These helpers are intentionally lightweight. They are not replacements for
Rosetta, Reduce, MolProbity, packstat, or a proper SASA engine; they provide
stable, inspectable proxies for cheap triage.
"""

from __future__ import annotations

import math
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

AA3 = frozenset(AA3_TO_1)
HYDROPHOBIC_RESNAMES = frozenset({
    "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "PRO", "CYS",
})

BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})

VDW_RADII = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "H": 1.10,
}


@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    element: str
    resname: str
    chain: str
    resseq: str
    icode: str
    x: float
    y: float
    z: float
    residue_order: int

    @property
    def residue_id(self) -> tuple[str, str, str]:
        return (self.chain, self.resseq, self.icode)

    @property
    def radius(self) -> float:
        return VDW_RADII.get(self.element, 1.70)

    @property
    def is_sidechain(self) -> bool:
        return self.name not in BACKBONE_ATOMS

    @property
    def is_polar_heavy(self) -> bool:
        return self.element in {"N", "O", "S"}


def infer_element(atom_name: str, element_field: str = "") -> str:
    element = element_field.strip().upper()
    if element:
        return element
    cleaned = "".join(ch for ch in atom_name.strip().upper() if ch.isalpha())
    if not cleaned:
        return "C"
    if len(cleaned) >= 2 and cleaned[:2] in {
        "CL", "BR", "NA", "MG", "FE", "ZN", "CA", "MN", "CU", "NI", "CO",
    }:
        return cleaned[:2]
    return cleaned[0]


def parse_pdb_atoms(
    path: Path,
    *,
    standard_only: bool = True,
    heavy_only: bool = True,
) -> list[Atom]:
    atoms: list[Atom] = []
    residue_orders: dict[tuple[str, str, str], int] = {}

    with path.open() as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resname = line[17:20].strip().upper()
            if standard_only and resname not in AA3:
                continue
            name = line[12:16].strip()
            element = infer_element(name, line[76:78] if len(line) >= 78 else "")
            if heavy_only and element in {"H", "D"}:
                continue
            chain = line[21].strip()
            resseq = line[22:26].strip()
            icode = line[26].strip()
            residue_id = (chain, resseq, icode)
            if residue_id not in residue_orders:
                residue_orders[residue_id] = len(residue_orders) + 1
            atoms.append(Atom(
                serial=int(line[6:11]),
                name=name,
                element=element,
                resname=resname,
                chain=chain,
                resseq=resseq,
                icode=icode,
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
                residue_order=residue_orders[residue_id],
            ))

    if not atoms:
        raise ValueError(f"{path} has no standard protein atoms")
    return atoms


def residue_ids(atoms: Iterable[Atom]) -> list[tuple[str, str, str]]:
    seen: set[tuple[str, str, str]] = set()
    ordered: list[tuple[str, str, str]] = []
    for atom in atoms:
        if atom.residue_id not in seen:
            seen.add(atom.residue_id)
            ordered.append(atom.residue_id)
    return ordered


def residue_count(atoms: Iterable[Atom]) -> int:
    return len(residue_ids(atoms))


def sequence_from_atoms(atoms: Iterable[Atom]) -> str:
    by_residue: dict[tuple[str, str, str], str] = {}
    for atom in atoms:
        by_residue.setdefault(atom.residue_id, atom.resname)
    return "".join(AA3_TO_1[by_residue[rid]] for rid in residue_ids(atoms))


def distance(a: Atom, b: Atom) -> float:
    return math.sqrt(
        (a.x - b.x) * (a.x - b.x)
        + (a.y - b.y) * (a.y - b.y)
        + (a.z - b.z) * (a.z - b.z)
    )


def _cell(coord: tuple[float, float, float], size: float) -> tuple[int, int, int]:
    return (
        math.floor(coord[0] / size),
        math.floor(coord[1] / size),
        math.floor(coord[2] / size),
    )


def nearby_atom_pairs(
    atoms: list[Atom],
    cutoff: float,
) -> Iterable[tuple[int, int, float]]:
    if cutoff <= 0:
        return
    grid: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for idx, atom in enumerate(atoms):
        c = _cell((atom.x, atom.y, atom.z), cutoff)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for other_idx in grid.get((c[0] + dx, c[1] + dy, c[2] + dz), []):
                        other = atoms[other_idx]
                        dist = distance(atom, other)
                        if dist <= cutoff:
                            yield other_idx, idx, dist
        grid[c].append(idx)


def nonlocal_residue_pair(a: Atom, b: Atom, *, skip_adjacent: bool = True) -> bool:
    if a.residue_id == b.residue_id:
        return False
    if skip_adjacent and a.chain == b.chain:
        try:
            a_resseq = int(a.resseq)
            b_resseq = int(b.resseq)
        except ValueError:
            if abs(a.residue_order - b.residue_order) <= 1:
                return False
        else:
            if a.icode == b.icode == "" and abs(a_resseq - b_resseq) <= 1:
                return False
    return True


def sphere_points(n_points: int) -> list[tuple[float, float, float]]:
    if n_points <= 0:
        raise ValueError("n_points must be positive")
    points: list[tuple[float, float, float]] = []
    offset = 2.0 / n_points
    increment = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n_points):
        y = (i * offset) - 1.0 + (offset / 2.0)
        r = math.sqrt(max(0.0, 1.0 - y * y))
        phi = i * increment
        points.append((math.cos(phi) * r, y, math.sin(phi) * r))
    return points


def atom_sasa(
    atoms: list[Atom],
    *,
    probe_radius: float = 1.4,
    n_points: int = 64,
) -> list[float]:
    """Return Shrake-Rupley-like SASA per atom in square Angstroms."""
    points = sphere_points(n_points)
    inflated = [atom.radius + probe_radius for atom in atoms]
    max_inflated = max(inflated)
    cell_size = max_inflated

    grid: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for idx, atom in enumerate(atoms):
        grid[_cell((atom.x, atom.y, atom.z), cell_size)].append(idx)

    values: list[float] = []
    for idx, atom in enumerate(atoms):
        radius = inflated[idx]
        accessible = 0
        for px, py, pz in points:
            sx = atom.x + px * radius
            sy = atom.y + py * radius
            sz = atom.z + pz * radius
            c = _cell((sx, sy, sz), cell_size)
            occluded = False
            for dx in (-1, 0, 1):
                if occluded:
                    break
                for dy in (-1, 0, 1):
                    if occluded:
                        break
                    for dz in (-1, 0, 1):
                        for other_idx in grid.get((c[0] + dx, c[1] + dy, c[2] + dz), []):
                            if other_idx == idx:
                                continue
                            other = atoms[other_idx]
                            cutoff = inflated[other_idx]
                            if (
                                (sx - other.x) * (sx - other.x)
                                + (sy - other.y) * (sy - other.y)
                                + (sz - other.z) * (sz - other.z)
                            ) < cutoff * cutoff:
                                occluded = True
                                break
                        if occluded:
                            break
            if not occluded:
                accessible += 1
        values.append(4.0 * math.pi * radius * radius * accessible / n_points)
    return values


def residue_sasa(atoms: list[Atom], sasa: list[float]) -> dict[tuple[str, str, str], float]:
    totals: dict[tuple[str, str, str], float] = defaultdict(float)
    for atom, value in zip(atoms, sasa):
        totals[atom.residue_id] += value
    return dict(totals)


def sidechain_sasa(atoms: list[Atom], sasa: list[float]) -> dict[tuple[str, str, str], float]:
    totals: dict[tuple[str, str, str], float] = defaultdict(float)
    for atom, value in zip(atoms, sasa):
        if atom.is_sidechain:
            totals[atom.residue_id] += value
    return dict(totals)


def sidechain_capacity(
    atoms: list[Atom],
    *,
    probe_radius: float = 1.4,
) -> dict[tuple[str, str, str], float]:
    totals: dict[tuple[str, str, str], float] = defaultdict(float)
    for atom in atoms:
        if atom.is_sidechain:
            radius = atom.radius + probe_radius
            totals[atom.residue_id] += 4.0 * math.pi * radius * radius
    return dict(totals)


def residue_resnames(atoms: Iterable[Atom]) -> dict[tuple[str, str, str], str]:
    names: dict[tuple[str, str, str], str] = {}
    for atom in atoms:
        names.setdefault(atom.residue_id, atom.resname)
    return names


def hydrophobic_residue_components(
    atoms: list[Atom],
    exposed_residues: set[tuple[str, str, str]],
    *,
    contact_cutoff: float = 6.0,
) -> list[set[tuple[str, str, str]]]:
    graph: dict[tuple[str, str, str], set[tuple[str, str, str]]] = {
        rid: set() for rid in exposed_residues
    }
    sidechain_indices = [
        idx for idx, atom in enumerate(atoms)
        if atom.is_sidechain and atom.residue_id in exposed_residues
    ]
    sidechain_atoms = [atoms[idx] for idx in sidechain_indices]
    for i, j, _dist in nearby_atom_pairs(sidechain_atoms, contact_cutoff):
        a = sidechain_atoms[i]
        b = sidechain_atoms[j]
        if a.residue_id == b.residue_id:
            continue
        graph[a.residue_id].add(b.residue_id)
        graph[b.residue_id].add(a.residue_id)

    components: list[set[tuple[str, str, str]]] = []
    remaining = set(exposed_residues)
    while remaining:
        start = remaining.pop()
        component = {start}
        queue = deque([start])
        while queue:
            current = queue.popleft()
            for nxt in graph[current]:
                if nxt in remaining:
                    remaining.remove(nxt)
                    component.add(nxt)
                    queue.append(nxt)
        components.append(component)
    return components


def round_float(value: float, digits: int = 3) -> float:
    return round(float(value), digits)
