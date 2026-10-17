#!/usr/bin/env python3
"""Generate site-configurable Slurm scripts. Submission requires explicit --submit.

Targets a single-node multicore NAMD executable: one process with +p CPU threads.
An afterok dependency orders jobs; restart inputs must be configured separately.
"""

from __future__ import annotations

import argparse
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

try:
    from .common import DEFAULT_NAMD, identifier, nonnegative_int, positive_int, require_file
except ImportError:
    from common import DEFAULT_NAMD, identifier, nonnegative_int, positive_int, require_file


def generate(args) -> list[Path]:
    if not re.fullmatch(r"(?:\d+-)?\d{1,3}:\d{2}:\d{2}", args.time):
        raise ValueError("--time must use [days-]hours:minutes:seconds")
    for name in ("account", "partition", "qos", "gpu_type"):
        value = getattr(args, name)
        if value and not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
            raise ValueError(f"invalid --{name.replace('_', '-')} token")
    if args.gpu_type and not args.gpus:
        raise ValueError("--gpu-type requires --gpus greater than zero")
    if args.gpu_devices and not args.gpus:
        raise ValueError("--gpu-devices requires --gpus greater than zero")
    gpu_devices = args.gpu_devices or ",".join(str(index) for index in range(args.gpus))
    if args.gpus and (not re.fullmatch(r"\d+(?:,\d+)*", gpu_devices) or
                      len(gpu_devices.split(",")) != args.gpus or
                      len(set(gpu_devices.split(","))) != args.gpus):
        raise ValueError("--gpu-devices must list each allocated NAMD GPU once, such as 0 or 0,1")
    if args.dependency and not re.fullmatch(r"\d+(?::\d+)*", args.dependency):
        raise ValueError("--dependency must contain job IDs separated by colons")
    if any(c in args.namd for c in "\n\r"):
        raise ValueError("--namd must name an executable, not a multiline command")
    plans = []
    for path in args.configs:
        config = path.resolve()
        require_file(config)
        identifier(config.stem)
        script = config.with_suffix(".slurm")
        logs = config.parent.parent / "logs"
        job_name = f"{identifier(config.parent.parent.name)}.{config.stem}"[:80]
        options = [("job-name", job_name), ("nodes", "1"), ("ntasks", "1"),
                   ("cpus-per-task", str(args.cpus)), ("time", args.time),
                   ("chdir", str(config.parent)),
                   ("output", str(logs / f"{config.stem}.%j.out")),
                   ("error", str(logs / f"{config.stem}.%j.err"))]
        for name in ("account", "partition", "qos"):
            if getattr(args, name):
                options.append((name, getattr(args, name)))
        if args.gpus:
            resource = f"gpu:{args.gpu_type}:{args.gpus}" if args.gpu_type else f"gpu:{args.gpus}"
            options.append(("gres", resource))
        if args.dependency:
            options.append(("dependency", "afterok:" + args.dependency))
        body = "#!/usr/bin/env bash\n" + "".join(f"#SBATCH --{key}={shlex.quote(value)}\n" for key, value in options)
        namd_command = [args.namd, f"+p{args.cpus}"]
        if args.gpus:
            namd_command.extend(("+setcpuaffinity", "+devices", gpu_devices))
        namd_command.append(str(config))
        body += "set -euo pipefail\nexec " + shlex.join(namd_command) + "\n"
        if script.exists() and script.read_text() != body:
            raise ValueError(f"existing Slurm script has different contents: {script}")
        plans.append((script, logs, body))
    if len({path for path, _, _ in plans}) != len(plans):
        raise ValueError("duplicate configurations supplied")
    if args.submit and not shutil.which("sbatch"):
        raise ValueError("sbatch not found; omit --submit to generate scripts locally")
    for script, logs, body in plans:
        logs.mkdir(parents=True, exist_ok=True)
        if script.exists():
            print(f"Reusing {script}")
        else:
            script.write_text(body)
            print(f"Wrote {script}")
        if args.submit:
            result = subprocess.run(["sbatch", "--parsable", str(script)], capture_output=True, text=True, check=True)
            print(f"Submitted: {result.stdout.strip()}")
        else:
            print(f"Preview only. Submit explicitly with: {shlex.join(['sbatch', str(script)])}")
    return [script for script, _, _ in plans]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configs", nargs="+", type=Path)
    parser.add_argument("--namd", default=str(DEFAULT_NAMD), help="compute-node executable path or command name; default: local CUDA NAMD 3")
    parser.add_argument("--time", required=True, help="Slurm wall time in [days-]hours:minutes:seconds")
    parser.add_argument("--cpus", type=positive_int, default=6)
    parser.add_argument("--gpus", type=nonnegative_int, default=1)
    parser.add_argument("--gpu-type", default="H200", help="Slurm GPU type; default: H200")
    parser.add_argument("--gpu-devices", help="NAMD-visible GPU IDs, e.g. 0 or 0,1; default: 0 through gpus-1")
    parser.add_argument("--account", default="gts-jcgumbart-coda20")
    parser.add_argument("--partition", default="embers")
    parser.add_argument("--qos")
    parser.add_argument("--dependency", help="afterok job IDs separated by colons; does not configure restart")
    parser.add_argument("--submit", action="store_true", help="actually call sbatch; default only writes scripts")
    args = parser.parse_args(argv)
    try:
        generate(args)
    except subprocess.CalledProcessError as exc:
        print(f"error: sbatch failed: {exc.stderr.strip()}", file=sys.stderr)
        return 2
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
