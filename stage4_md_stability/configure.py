#!/usr/bin/env python3
"""Generate NAMD initial/restart segment configurations without running NAMD.

Run length is always explicit. The NPT template is experimental; segment names
such as 'equilibration' or 'production' do not establish scientific equilibration.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

try:
    from .common import TEMPLATES, identifier, load_cases, nonnegative_int, positive_float, positive_int, require_file, restart_step, tcl_list, tcl_word
except ImportError:
    from common import TEMPLATES, identifier, load_cases, nonnegative_int, positive_float, positive_int, require_file, restart_step, tcl_list, tcl_word


def configure(args) -> list[Path]:
    identifier(args.segment)
    cases = load_cases(args.manifest, args.sequence_id)
    if args.restart_prefix and len(cases) != 1:
        raise ValueError("--restart-prefix requires exactly one selected system")
    if args.restart_from:
        identifier(args.restart_from)
        if args.restart_from == args.segment:
            raise ValueError("restart input and new output segment must differ")
    restarting = bool(args.restart_from or args.restart_prefix)
    minimize = args.minimize_steps if args.minimize_steps is not None else (0 if restarting else 200)
    if restarting and minimize:
        raise ValueError("restart segments must not minimize again; use --minimize-steps 0")
    if args.run_steps % 10 or minimize % 10 or args.output_every % 10 or args.energy_every % 10:
        raise ValueError("step counts and output intervals must be multiples of the template's 10-step cycle")
    template = (TEMPLATES / "namd.tcl").read_text()
    plans = []
    for case_dir, metadata in cases:
        for name in ("system.psf", "system.pdb"):
            require_file(case_dir / name)
        if type(metadata.get("hmr")) is not bool:
            raise ValueError("prepared metadata must explicitly record hmr as true or false")
        timestep = args.timestep_fs if args.timestep_fs is not None else (4.0 if metadata["hmr"] else 2.0)
        if timestep > 4.0 or (timestep > 2.0 and not metadata["hmr"]):
            raise ValueError("timesteps above 2 fs require an HMR-prepared system; maximum template timestep is 4 fs")
        lengths, origin = metadata["box_lengths"], metadata["box_origin"]
        if (len(lengths) != 3 or len(origin) != 3 or
                not all(math.isfinite(float(v)) for v in lengths + origin) or min(lengths) <= 0):
            raise ValueError("invalid prepared box dimensions or origin")
        config = case_dir / "conf" / f"{args.segment}.conf"
        output = case_dir / "out" / args.segment
        if config.exists() or any(output.parent.glob(output.name + ".*")):
            raise ValueError(f"configuration or output already exists for {metadata['sequence_id']}/{args.segment}; choose a new segment")
        prefix = ""
        first_step = 0
        if args.restart_from:
            prefix = f"out/{args.restart_from}.restart"
            first_step = restart_step(case_dir / prefix)
        elif args.restart_prefix:
            prefix = str(args.restart_prefix.resolve())
            if Path(prefix) == output.resolve() or Path(prefix) == Path(str(output.resolve()) + ".restart"):
                raise ValueError("restart prefix would collide with new output files")
            first_step = restart_step(Path(prefix))
        settings = {
            "segment": tcl_word(args.segment), "parameter_files": tcl_list(path.resolve() for path in args.parameter),
            "box_lengths": tcl_list(lengths), "box_origin": tcl_list(origin),
            "restart_prefix": tcl_word(prefix), "restart_relative": int(bool(args.restart_from)),
            "timestep_fs": timestep, "run_steps": args.run_steps, "minimize_steps": minimize,
            "temperature_k": args.temperature, "pressure_bar": args.pressure_bar,
            "output_every": args.output_every, "energy_every": args.energy_every,
        }
        header = (f"# Generated experimental Stage 4 configuration for {metadata['sequence_id']}.\n"
                  f"# Additional MD time: {args.run_steps * timestep / 1e6:g} ns.\n"
                  f"# Trajectory/restart interval: {args.output_every * timestep / 1000:g} ps.\n"
                  f"# Restart step at configuration time: {first_step}.\n")
        content = header + "".join(f"set {key} {value}\n" for key, value in settings.items()) + "\n" + template
        plans.append((config, content))
    for config, content in plans:
        for name in ("conf", "out", "logs"):
            (config.parent.parent / name).mkdir(parents=True, exist_ok=True)
        config.write_text(content)
        print(f"Wrote {config}")
    return [path for path, _ in plans]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--sequence-id", action="append", help="select one or more IDs; default: all systems")
    parser.add_argument("--segment", required=True, help="new segment name, such as equilibration or production")
    parser.add_argument("--run-steps", type=positive_int, required=True, help="additional integration steps for this segment")
    parser.add_argument("--minimize-steps", type=nonnegative_int, help="initial default: 200; restart default: 0")
    parser.add_argument("--parameter", type=Path, action="append", required=True,
                        help="external CHARMM parameter file; repeat in loading order; must exist on the NAMD host")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--restart-from", help="previous segment in each system's out directory; reads SEGMENT.restart.{coor,vel,xsc}")
    group.add_argument("--restart-prefix", type=Path, help="explicit existing restart prefix for one selected system")
    parser.add_argument("--timestep-fs", type=positive_float, help="default: 2 fs without HMR, 4 fs with HMR")
    parser.add_argument("--temperature", type=positive_float, default=300.0, help="Kelvin")
    parser.add_argument("--pressure-bar", type=positive_float, default=1.01325)
    parser.add_argument("--output-every", type=positive_int, default=25000, help="DCD/restart/cell output interval in steps")
    parser.add_argument("--energy-every", type=positive_int, default=5000)
    args = parser.parse_args(argv)
    try:
        configure(args)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
