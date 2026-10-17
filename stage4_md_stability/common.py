"""Shared helpers for the experimental MD stage (standard library only)."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


TEMPLATES = Path(__file__).resolve().parent / "templates"


def identifier(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", value):
        raise ValueError(f"invalid identifier {value!r}; use 1–80 letters, digits, dots, underscores or hyphens, starting with a letter or digit")
    return value


def positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return number


def nonnegative_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError("must be finite and nonnegative")
    return number


def positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def nonnegative_int(value: str) -> int:
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return number


def tcl_word(value: object) -> str:
    """Quote one Tcl word, including paths containing spaces or substitutions."""
    text = str(value)
    for old, new in (("\\", "\\\\"), ('"', '\\"'), ("$", "\\$"),
                     ("[", "\\["), ("]", "\\]"), ("\n", "\\n"),
                     ("\r", "\\r"), ("\t", "\\t")):
        text = text.replace(old, new)
    return f'"{text}"'


def tcl_list(values) -> str:
    return "[list " + " ".join(tcl_word(value) for value in values) + "]"


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty file: {path}")


def load_cases(manifest: Path, selected: list[str] | None = None) -> list[tuple[Path, dict]]:
    payload = json.loads(manifest.read_text())
    if payload.get("schema_version") != 1 or not payload.get("systems"):
        raise ValueError("expected a nonempty Stage 4 systems.json manifest (schema_version=1)")
    cases = []
    seen = set()
    for entry in payload["systems"]:
        sequence_id = identifier(entry["sequence_id"])
        if sequence_id in seen:
            raise ValueError(f"duplicate sequence_id {sequence_id!r} in manifest")
        seen.add(sequence_id)
        case_dir = manifest.resolve().parent / sequence_id
        metadata = json.loads((case_dir / "case.json").read_text())
        if metadata.get("sequence_id") != sequence_id or metadata.get("state") != "prepared":
            raise ValueError(f"{sequence_id}: case metadata is not a matching prepared system")
        if selected is None or sequence_id in selected:
            cases.append((case_dir, metadata))
    if selected and (set(selected) - seen):
        raise ValueError(f"unknown sequence IDs: {', '.join(sorted(set(selected) - seen))}")
    return cases


def restart_step(prefix: Path) -> int:
    for suffix in (".coor", ".vel", ".xsc"):
        require_file(Path(str(prefix) + suffix))
    rows = [line.split() for line in Path(str(prefix) + ".xsc").read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    if len(rows) != 1 or len(rows[0]) < 13:
        raise ValueError(f"invalid NAMD restart cell file: {prefix}.xsc")
    step = int(rows[0][0])
    if step < 0 or not all(math.isfinite(float(value)) for value in rows[0][1:13]):
        raise ValueError(f"invalid NAMD restart cell values: {prefix}.xsc")
    return step
