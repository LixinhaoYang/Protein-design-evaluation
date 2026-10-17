#!/usr/bin/env python3
"""Prepare selected single-chain PDBs with VMD, explicit water/KCl, and optional HMR.

This experimental tool runs system preparation only, not molecular dynamics.
Topology files and VMD are supplied by the user; no force-field files are bundled.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

try:
    from .common import TEMPLATES, identifier, nonnegative_float, positive_float, require_file, tcl_list, tcl_word, write_json
except ImportError:
    from common import TEMPLATES, identifier, nonnegative_float, positive_float, require_file, tcl_list, tcl_word, write_json


STANDARD_RESIDUES = frozenset("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())


def normalize_pdb(source: Path, output: Path, chain: str, sequence_id: str) -> list[dict]:
    """Extract one heavy-atom protein chain and record all residue renumbering."""
    if source.suffix.lower() != ".pdb":
        raise ValueError(f"{source}: Stage 4 currently requires PDB; convert AF3 mmCIF first")
    if len(chain) != 1:
        raise ValueError("chain ID must be one character; use '_' for a blank PDB chain")
    chain = " " if chain == "_" else chain
    residues: dict[tuple[str, str], dict] = {}
    models = 0
    chain_ended = False
    for line in source.read_text().splitlines():
        if line.startswith("MODEL "):
            models += 1
            if models > 1:
                raise ValueError(f"{source}: multiple models are not supported")
        if line.startswith("SSBOND"):
            raise ValueError(f"{source}: disulfide patching is not supported by this preparer")
        if line.startswith("TER") and residues:
            chain_ended = True
        if not line.startswith(("ATOM  ", "HETATM")) or line[21:22] != chain:
            continue
        resname = line[17:20].strip()
        if resname in {"HOH", "WAT", "TIP3"}:
            continue
        if resname not in STANDARD_RESIDUES:
            raise ValueError(f"{source}: unsupported selected-chain residue {resname!r}")
        if chain_ended:
            raise ValueError(f"{source}: multiple fragments in the selected chain are not supported")
        key = (line[22:26].strip(), line[26:27].strip())
        int(key[0])
        residue = residues.setdefault(key, {"resname": resname, "atoms": {}})
        if residue["resname"] != resname:
            raise ValueError(f"{source}: alternate residue identities are not supported")
        name = line[12:16].strip()
        element = line[76:78].strip().upper() or next((c for c in name if c.isalpha()), "").upper()
        if element in {"H", "D"}:
            continue
        coords = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        if not all(math.isfinite(value) for value in coords):
            raise ValueError(f"{source}: nonfinite coordinates")
        altloc = line[16:17]
        occupancy = float(line[54:60].strip() or "0")
        priority = (altloc == " ", occupancy, altloc == "A")
        old = residue["atoms"].get(name)
        if old is None or priority > old[0]:
            residue["atoms"][name] = (priority, name, coords, element)
    if not residues:
        raise ValueError(f"{source}: no protein residues in selected chain {chain!r}")
    if len(residues) > 9999:
        raise ValueError("too many residues for normalized PDB numbering")
    lines, mapping = [], []
    serial = 0
    for position, ((number, insertion), residue) in enumerate(residues.items(), start=1):
        if "CA" not in residue["atoms"]:
            raise ValueError(f"{source}: residue {number}{insertion} has no CA atom")
        mapping.append({"sequence_id": sequence_id, "input_chain": chain.strip(),
                        "input_residue_number": number, "input_insertion_code": insertion,
                        "prepared_residue_number": position, "residue_name": residue["resname"]})
        for _, name, (x, y, z), element in residue["atoms"].values():
            serial += 1
            if serial > 99999:
                raise ValueError("too many protein atoms for PDB serial numbering")
            atom_field = f" {name:<3}" if len(name) < 4 else name
            lines.append(f"ATOM  {serial:5d} {atom_field} {residue['resname']:3s} A{position:4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n")
    output.write_text("".join(lines) + "TER\nEND\n")
    return mapping


def read_inputs(args) -> list[dict]:
    if args.case_csv:
        if args.pdbs or args.source or args.sequence_id:
            raise ValueError("use either --case-csv or PDB paths with --source/--sequence-id")
        with args.case_csv.open(newline="") as handle:
            reader = csv.DictReader(handle)
            required = {"sequence_id", "input_structure", "structure_source"}
            if not required <= set(reader.fieldnames or []):
                raise ValueError("case CSV requires sequence_id,input_structure,structure_source")
            cases = [{"sequence_id": row["sequence_id"],
                      "input_structure": (args.case_csv.parent / row["input_structure"]).resolve(),
                      "structure_source": row["structure_source"],
                      "chain": row.get("chain_id") or args.chain}
                     for row in reader]
    else:
        if not args.pdbs or not args.source:
            raise ValueError("provide PDB paths with --source design|af3, or --case-csv")
        if args.sequence_id and len(args.pdbs) != 1:
            raise ValueError("--sequence-id requires exactly one PDB")
        cases = [{"sequence_id": args.sequence_id or path.stem,
                  "input_structure": path.resolve(), "structure_source": args.source,
                  "chain": args.chain} for path in args.pdbs]
    if not cases:
        raise ValueError("no candidate systems supplied")
    seen = set()
    for case in cases:
        name = identifier(case["sequence_id"])
        if name in seen:
            raise ValueError(f"duplicate sequence_id {name!r}")
        seen.add(name)
        if case["structure_source"] not in {"design", "af3"}:
            raise ValueError("structure_source must be design or af3")
        require_file(case["input_structure"])
    return cases


def prepare(args) -> Path:
    cases = read_inputs(args)
    vmd = shutil.which(args.vmd)
    if not vmd:
        raise ValueError("VMD executable not found; set --vmd or VMD_BIN")
    topologies = [path.resolve() for path in args.topology]
    for path in topologies:
        require_file(path)
    root = args.outdir.resolve()
    if root.exists() and (not root.is_dir() or any(root.iterdir())):
        raise ValueError(f"output directory must be new or empty: {root}")
    root.mkdir(parents=True, exist_ok=True)
    systems = []
    for case in cases:
        case_dir = root / case["sequence_id"]
        case_dir.mkdir()
        for subdir in ("conf", "out", "logs"):
            (case_dir / subdir).mkdir()
        mapping = normalize_pdb(case["input_structure"], case_dir / "protein.pdb", case["chain"], case["sequence_id"])
        with (case_dir / "residue_map.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(mapping[0]))
            writer.writeheader()
            writer.writerows(mapping)
        driver = [f"set topology_files {tcl_list(topologies)}",
                  f"set water_padding {args.padding}", f"set salt_molar {args.salt_molar}",
                  f"set use_hmr {int(args.hmr)}",
                  f"source {tcl_word(TEMPLATES / 'geometry.tcl')}",
                  f"source {tcl_word(TEMPLATES / 'hmr.tcl')}",
                  f"if {{[catch {{source {tcl_word(TEMPLATES / 'prepare.tcl')}}} message]}} {{",
                  '    puts stderr "Stage 4 preparation failed: $message"', "    exit 1", "}",
                  'puts "STAGE4_PREPARED"', "exit 0"]
        (case_dir / "build.tcl").write_text("\n".join(driver) + "\n")
        log = case_dir / "logs" / "prepare.log"
        with log.open("w") as handle:
            result = subprocess.run([vmd, "-dispdev", "text", "-e", "build.tcl"],
                                    cwd=case_dir, stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode or "STAGE4_PREPARED" not in log.read_text():
            raise ValueError(f"{case['sequence_id']}: VMD preparation failed; inspect {log}")
        for name in ("system.psf", "system.pdb", "box.dat"):
            require_file(case_dir / name)
        box = [float(value) for value in (case_dir / "box.dat").read_text().split()]
        if len(box) != 6 or not all(math.isfinite(v) for v in box) or min(box[:3]) <= 0:
            raise ValueError(f"{case['sequence_id']}: invalid prepared box metadata")
        metadata = {"schema_version": 1, "state": "prepared", "sequence_id": case["sequence_id"],
                    "structure_source": case["structure_source"],
                    "input_structure": str(case["input_structure"]), "input_chain": case["chain"],
                    "normalized_chain": "A", "protein_segid": "PRO", "protein_residues": len(mapping),
                    "box_lengths": box[:3], "box_origin": box[3:], "hmr": args.hmr,
                    "hydrogen_mass": 3.024 if args.hmr else None,
                    "water_padding_angstrom": args.padding, "salt_molar": args.salt_molar,
                    "topology_files": [str(path) for path in topologies]}
        write_json(case_dir / "case.json", metadata)
        systems.append({"sequence_id": case["sequence_id"]})
        print(f"Prepared {case['sequence_id']}: {case_dir}")
    manifest = root / "systems.json"
    write_json(manifest, {"schema_version": 1, "systems": systems})
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdbs", nargs="*", type=Path)
    parser.add_argument("--case-csv", type=Path)
    parser.add_argument("--source", choices=("design", "af3"))
    parser.add_argument("--sequence-id", help="candidate ID for a single PDB; default: filename stem")
    parser.add_argument("--chain", default="A", help="PDB chain to extract; '_' selects a blank chain")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--vmd", default=os.environ.get("VMD_BIN", "vmd"))
    parser.add_argument("--topology", type=Path, action="append", required=True,
                        help="external CHARMM protein topology; repeat for additional files")
    parser.add_argument("--padding", type=positive_float, default=25.0, help="water padding on each face in Angstrom")
    parser.add_argument("--salt-molar", type=nonnegative_float, default=0.15, help="KCl concentration passed to autoionize")
    parser.add_argument("--hmr", action="store_true", help="repartition solute hydrogen masses to 3.024 Da")
    args = parser.parse_args(argv)
    try:
        print(f"Manifest: {prepare(args)}")
    except (OSError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
