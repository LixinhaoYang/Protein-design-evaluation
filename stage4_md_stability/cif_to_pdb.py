#!/usr/bin/env python3
"""Convert a single-chain protein mmCIF (e.g. an AF3 model) to PDB for Stage 4.

Stage 4 preparation accepts PDB only. This tool performs that one conversion and
nothing else: it does not select, score, or repair structures. Anything outside a
single-model, single-conformer protein chain is rejected rather than quietly
dropped, so an unexpected input fails here instead of inside PSFGEN.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from .common import identifier
except ImportError:
    from common import identifier


# _atom_site columns this converter reads. Others may follow them; the required
# ones must be present, since mmCIF does not fix their order.
REQUIRED_COLUMNS = (
    "group_PDB", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
    "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv",
    "auth_seq_id", "auth_asym_id", "pdbx_PDB_model_num",
)

STANDARD_RESIDUES = frozenset("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())


def _atom_site_rows(text: str, source: Path) -> tuple[dict[str, int], list[list[str]]]:
    """Return the _atom_site column index map and its whitespace-delimited rows."""
    order: list[str] = []
    rows: list[list[str]] = []
    in_loop = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("_atom_site."):
            order.append(stripped.split(".", 1)[1].split()[0])
            in_loop = True
            continue
        if not in_loop:
            continue
        if stripped.startswith(("ATOM", "HETATM")):
            if any(quote in stripped for quote in ("'", '"')):
                raise ValueError(f"{source}: quoted _atom_site value is not supported")
            fields = stripped.split()
            if len(fields) != len(order):
                raise ValueError(f"{source}: _atom_site row does not match its column count")
            rows.append(fields)
        elif rows:
            break
    if not order:
        raise ValueError(f"{source}: no _atom_site loop found")
    missing = [column for column in REQUIRED_COLUMNS if column not in order]
    if missing:
        raise ValueError(f"{source}: _atom_site is missing {', '.join(missing)}")
    if not rows:
        raise ValueError(f"{source}: _atom_site loop contains no atoms")
    return {name: position for position, name in enumerate(order)}, rows


def convert(source: Path, output: Path, chain: str = "A") -> int:
    """Write one protein chain of `source` as PDB and return the atom count."""
    if len(chain) != 1:
        raise ValueError("chain ID must be one character")
    index, rows = _atom_site_rows(source.read_text(), source)
    get = lambda row, name: row[index[name]]  # noqa: E731

    chains = {get(row, "auth_asym_id") for row in rows}
    if chain not in chains:
        raise ValueError(f"{source}: chain {chain!r} not present; found {', '.join(sorted(chains))}")

    written = 0
    serial = 0
    with output.open("w") as handle:
        for row in rows:
            if get(row, "auth_asym_id") != chain:
                continue
            if get(row, "group_PDB") != "ATOM":
                raise ValueError(f"{source}: chain {chain} contains a HETATM record")
            if get(row, "pdbx_PDB_model_num") != "1":
                raise ValueError(f"{source}: more than one model; split it first")
            if get(row, "label_alt_id") not in (".", "?", "A"):
                raise ValueError(f"{source}: alternate conformers are not supported")
            residue = get(row, "label_comp_id")
            if residue not in STANDARD_RESIDUES:
                raise ValueError(f"{source}: non-standard residue {residue} in chain {chain}")
            element = get(row, "type_symbol")
            name = get(row, "label_atom_id")
            # PSFGEN rebuilds hydrogens, and its CTER patch builds OT1/OT2 itself.
            if element == "H" or name == "OXT":
                continue
            serial += 1
            written += 1
            # PDB columns: name occupies 13-16, with 1-3 character names indented
            # by one so the element letter sits in column 14.
            atom_name = f"{name:<4s}" if len(name) == 4 else f" {name:<3s}"
            handle.write(
                f"ATOM  {serial:5d} {atom_name} {residue:>3s} {chain:1s}"
                f"{int(get(row, 'auth_seq_id')):4d}    "
                f"{float(get(row, 'Cartn_x')):8.3f}"
                f"{float(get(row, 'Cartn_y')):8.3f}"
                f"{float(get(row, 'Cartn_z')):8.3f}"
                f"{float(get(row, 'occupancy')):6.2f}"
                f"{float(get(row, 'B_iso_or_equiv')):6.2f}"
                f"          {element:>2s}\n")
        if not written:
            raise ValueError(f"{source}: no protein atoms written for chain {chain}")
        handle.write("TER\nEND\n")
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cifs", nargs="+", type=Path)
    parser.add_argument("--chain", default="A", help="mmCIF auth_asym_id to extract")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="directory for the converted PDBs; named <sequence_id>.pdb")
    parser.add_argument("--sequence-id",
                        help="output name for a single mmCIF; default: filename stem without '_model'")
    args = parser.parse_args(argv)

    if args.sequence_id and len(args.cifs) != 1:
        parser.error("--sequence-id accepts exactly one mmCIF")
    try:
        args.outdir.mkdir(parents=True, exist_ok=True)
        for cif in args.cifs:
            stem = args.sequence_id or cif.stem.removesuffix("_model")
            output = args.outdir / f"{identifier(stem)}.pdb"
            if output.exists():
                raise ValueError(f"{output}: already exists; remove it or choose another --outdir")
            print(f"{output}: {convert(cif, output, args.chain)} atoms")
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
