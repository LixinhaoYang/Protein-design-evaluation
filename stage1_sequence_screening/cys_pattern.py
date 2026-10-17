#!/usr/bin/env python3
"""Step 1 subcheck -- cysteine pattern annotation.

Cys is context-dependent. A free cysteine can cause oxidation, disulfide
scrambling, oligomerisation, purification headaches, or unwanted cross-linking;
but cysteines are also correct in disulfide-stabilised, metal-binding, redox, or
catalytic proteins. Sequence alone cannot see exposure, pairing, or function.

So this module is deliberately an ANNOTATION step, not a gate. Every flag is a
WARN, nothing fails, and the measured Cys pattern is the real output.

Measures
  * Cys count, fraction, and 1-based positions
  * odd/even parity
  * closest Cys-Cys spacing
  * maximum local Cys count in a window
  * dense Cys regions

Task metadata matters:
  * --disulfide-expected suppresses high-count / clustered-Cys warnings, while
    still warning on an odd Cys count because one Cys cannot be paired.
  * --cys-expected suppresses Cys liability flags entirely for designs where
    free/reactive/metal-binding/catalytic Cys residues are intended. The pattern
    is still measured and written to CSV/JSON.

Usage
    conda run -n base python stage1_sequence_screening/cys_pattern.py designs.fasta
    conda run -n base python stage1_sequence_screening/cys_pattern.py --seq MKVLAA...
    conda run -n base python stage1_sequence_screening/cys_pattern.py designs.fasta --disulfide-expected
    conda run -n base python stage1_sequence_screening/cys_pattern.py designs.fasta --cys-expected
    conda run -n base python stage1_sequence_screening/cys_pattern.py designs.fasta --csv cys.csv --json cys.json
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path

try:
    from .basic_validity import (
        Issue,
        Record,
        STANDARD_AA,
        _WHITESPACE,
        read_fasta,
        render_positions,
    )
except ImportError:  # pragma: no cover - direct script execution path
    from basic_validity import (
        Issue,
        Record,
        STANDARD_AA,
        _WHITESPACE,
        read_fasta,
        render_positions,
    )

# --------------------------------------------------------------------------
# Defaults.
#
# These are pragmatic review thresholds, not biochemical constants. They are
# tuned to stay silent on the current worked examples:
#
#     5OJ8 natural      4 Cys / 210 aa; max 1 per 30-aa window
#     5OJ8 design       2 Cys / 210 aa; max 1 per 30-aa window
#     1Q0S natural      2 Cys / 241 aa; max 1 per 30-aa window
#     1Q0S design       1 Cys / 247 aa; odd-count warning
#
# Hen lysozyme has 8 Cys and real disulfides; it should warn in the default
# soluble/cytosolic context, and become clean when disulfides are expected.
# --------------------------------------------------------------------------

DEFAULT_MAX_CYS_COUNT = 6
DEFAULT_MAX_CYS_FRACTION = 0.06
DEFAULT_CLUSTER_WINDOW = 30
DEFAULT_MAX_CYS_IN_WINDOW = 3
DEFAULT_CLOSE_PAIR_DISTANCE = 3


@dataclass
class CysPatternResult:
    seq_id: str
    length: int
    sequence: str
    metrics: dict
    cys_positions: list[int] = field(default_factory=list)
    dense_regions: list[tuple[int, int, int]] = field(default_factory=list)
    close_pairs: list[tuple[int, int]] = field(default_factory=list)
    issues: list[Issue] = field(default_factory=list)

    @property
    def status(self) -> str:
        return "WARN" if self.issues else "PASS"


# --------------------------------------------------------------------------
# Measurements
# --------------------------------------------------------------------------


def cys_positions(seq: str) -> list[int]:
    """1-based positions of every cysteine."""
    return [i + 1 for i, aa in enumerate(seq) if aa == "C"]


def closest_cys_spacing(positions: list[int]) -> int | None:
    """Minimum sequence distance between neighbouring cysteines, or None."""
    if len(positions) < 2:
        return None
    return min(b - a for a, b in zip(positions, positions[1:]))


def close_cys_pairs(
    positions: list[int], max_distance: int = DEFAULT_CLOSE_PAIR_DISTANCE
) -> list[tuple[int, int]]:
    """Adjacent Cys pairs separated by at most `max_distance` residues."""
    return [(a, b) for a, b in zip(positions, positions[1:]) if b - a <= max_distance]


def cys_window_counts(seq: str, window: int = DEFAULT_CLUSTER_WINDOW) -> list[int]:
    """Number of Cys residues in every window; short sequences score whole."""
    if not seq:
        return []
    if len(seq) <= window:
        return [seq.count("C")]

    values = [1 if aa == "C" else 0 for aa in seq]
    running = sum(values[:window])
    out = [running]
    for i in range(1, len(values) - window + 1):
        running += values[i + window - 1] - values[i - 1]
        out.append(running)
    return out


def max_cys_in_window(seq: str, window: int = DEFAULT_CLUSTER_WINDOW) -> tuple[int, int]:
    """Maximum Cys count in a window, and the 1-based window start."""
    counts = cys_window_counts(seq, window)
    if not counts:
        return 0, 0
    best = max(range(len(counts)), key=lambda i: counts[i])
    return counts[best], best + 1


def dense_cys_regions(
    seq: str,
    window: int = DEFAULT_CLUSTER_WINDOW,
    min_count: int = DEFAULT_MAX_CYS_IN_WINDOW + 1,
) -> list[tuple[int, int, int]]:
    """Merged spans where a window contains at least `min_count` Cys.

    Returns (start, end, max_count) with 1-based inclusive coordinates.
    """
    counts = cys_window_counts(seq, window)
    if not counts:
        return []
    if len(seq) <= window:
        return [(1, len(seq), counts[0])] if counts[0] >= min_count else []

    regions: list[list[int]] = []
    for i, count in enumerate(counts):
        if count < min_count:
            continue
        start, end = i + 1, i + window
        if regions and start <= regions[-1][1] + 1:
            regions[-1][1] = max(regions[-1][1], end)
            regions[-1][2] = max(regions[-1][2], count)
        else:
            regions.append([start, end, count])
    return [(a, b, c) for a, b, c in regions]


def measure(seq: str, window: int = DEFAULT_CLUSTER_WINDOW) -> dict:
    positions = cys_positions(seq)
    n_cys = len(positions)
    max_local, max_local_at = max_cys_in_window(seq, window)
    spacing = closest_cys_spacing(positions)
    return {
        "length": len(seq),
        "n_cys": n_cys,
        "cys_fraction": round(n_cys / len(seq), 4) if seq else 0.0,
        "cys_parity": "odd" if n_cys % 2 else "even",
        "closest_cys_spacing": spacing,
        f"max_cys_in_w{window}": max_local,
        f"max_cys_in_w{window}_at": max_local_at,
        "n_close_cys_pairs": len(close_cys_pairs(positions)),
    }


# --------------------------------------------------------------------------
# Checks -- advisory only
# --------------------------------------------------------------------------


def check_odd_cys_count(seq: str) -> list[Issue]:
    positions = cys_positions(seq)
    if len(positions) == 0 or len(positions) % 2 == 0:
        return []
    return [
        Issue(
            check="cys_odd_count",
            label=f"odd Cys count ({len(positions)})",
            reason="Odd Cys count means at least one cysteine cannot be paired if "
                   "the design is meant to use disulfides. In a soluble/cytosolic "
                   "design it is a free-thiol caution, not a defect by itself: "
                   "buried, catalytic, metal-binding or intentionally reactive Cys "
                   "can be correct.",
            positions=positions,
            severity="WARN",
        )
    ]


def check_high_cys_count(seq: str, max_count: int) -> list[Issue]:
    positions = cys_positions(seq)
    if len(positions) < max_count:
        return []
    return [
        Issue(
            check="cys_high_count",
            label=f"{len(positions)} Cys residues",
            reason=f"Cys count is at or above {max_count}. For ordinary soluble "
                   f"designs this raises redox, mispairing and aggregation risk. "
                   f"This is context-dependent: disulfide-stabilised or "
                   f"metal-binding proteins may legitimately be Cys-rich.",
            positions=positions,
            severity="WARN",
        )
    ]


def check_high_cys_fraction(seq: str, max_fraction: float) -> list[Issue]:
    positions = cys_positions(seq)
    if not seq or len(positions) / len(seq) <= max_fraction:
        return []
    return [
        Issue(
            check="cys_high_fraction",
            label=f"Cys fraction is high ({len(positions) / len(seq):.1%})",
            reason=f"Above the provisional {max_fraction:.0%} fraction. Cys-rich "
                   f"sequence can be correct in disulfide or metal-binding "
                   f"contexts, but is a liability for ordinary soluble/cytosolic "
                   f"expression.",
            positions=positions,
            severity="WARN",
        )
    ]


def check_dense_cys_cluster(seq: str, window: int, max_count: int) -> list[Issue]:
    regions = dense_cys_regions(seq, window, max_count + 1)
    if not regions:
        return []
    positions = [
        pos for pos in cys_positions(seq)
        if any(start <= pos <= end for start, end, _ in regions)
    ]
    spans = ", ".join(f"{a}-{b}" for a, b, _ in regions[:6])
    more = f" (+{len(regions) - 6} more)" if len(regions) > 6 else ""
    worst = max(regions, key=lambda r: r[2])
    return [
        Issue(
            check="cys_dense_cluster",
            label=f"dense Cys region(s): {spans}{more}",
            reason=f"Peak window contains {worst[2]} Cys in {window} aa, above "
                   f"the provisional limit of {max_count}. Dense Cys clusters "
                   f"are often disulfide, metal-binding, or redox motifs; flagging "
                   f"them is an annotation to check task intent, not a rejection.",
            positions=positions,
            severity="WARN",
        )
    ]


def check_close_cys_pairs(seq: str, max_distance: int) -> list[Issue]:
    pairs = close_cys_pairs(cys_positions(seq), max_distance)
    if not pairs:
        return []
    positions = sorted({p for pair in pairs for p in pair})
    pair_text = ", ".join(f"{a}/{b}" for a, b in pairs[:6])
    more = f" (+{len(pairs) - 6} more)" if len(pairs) > 6 else ""
    return [
        Issue(
            check="cys_close_pair",
            label=f"nearby Cys pair(s): {pair_text}{more}",
            reason=f"Cys pairs within {max_distance} residues can indicate redox, "
                   f"metal-binding or other functional motifs. They can be exactly "
                   f"right when intended, but should not appear by accident in a "
                   f"generic soluble scaffold.",
            positions=positions,
            severity="WARN",
        )
    ]


def analyse(
    record: Record,
    max_cys_count: int = DEFAULT_MAX_CYS_COUNT,
    max_cys_fraction: float = DEFAULT_MAX_CYS_FRACTION,
    window: int = DEFAULT_CLUSTER_WINDOW,
    max_cys_in_window: int = DEFAULT_MAX_CYS_IN_WINDOW,
    close_pair_distance: int = DEFAULT_CLOSE_PAIR_DISTANCE,
    expected_disulfides: bool = False,
    cys_expected: bool = False,
) -> CysPatternResult:
    """Profile one sequence.

    `cys_expected=True` is the broad task-specific exemption: cysteines are
    meant to be there, so the module reports the pattern but raises no flags.

    `expected_disulfides=True` is narrower: Cys-rich and clustered patterns are
    expected, but an odd count still warns because one Cys cannot be paired
    unless a free cysteine is also intended.
    """
    seq = record.sequence
    positions = cys_positions(seq)
    dense = dense_cys_regions(seq, window, max_cys_in_window + 1)
    pairs = close_cys_pairs(positions, close_pair_distance)

    if cys_expected:
        issues: list[Issue] = []
    elif expected_disulfides:
        issues = [*check_odd_cys_count(seq)]
    else:
        issues = [
            *check_odd_cys_count(seq),
            *check_high_cys_count(seq, max_cys_count),
            *check_high_cys_fraction(seq, max_cys_fraction),
            *check_dense_cys_cluster(seq, window, max_cys_in_window),
            *check_close_cys_pairs(seq, close_pair_distance),
        ]

    return CysPatternResult(
        seq_id=record.seq_id,
        length=len(seq),
        sequence=seq,
        metrics=measure(seq, window),
        cys_positions=positions,
        dense_regions=dense,
        close_pairs=pairs,
        issues=issues,
    )


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def _short(seq_id: str, width: int = 34) -> str:
    head = seq_id.split("|")[0]
    return head if len(head) <= width else head[: width - 1] + "..."


def print_report(
    results: list[CysPatternResult],
    window: int = DEFAULT_CLUSTER_WINDOW,
    max_cys_count: int = DEFAULT_MAX_CYS_COUNT,
    max_cys_fraction: float = DEFAULT_MAX_CYS_FRACTION,
    max_cys_in_window: int = DEFAULT_MAX_CYS_IN_WINDOW,
    close_pair_distance: int = DEFAULT_CLOSE_PAIR_DISTANCE,
    show_context: bool = True,
) -> None:
    print(f"\nStep 1 subcheck -- cysteine pattern   ({len(results)} sequence(s))")

    for result in results:
        m = result.metrics
        n = len(result.issues)
        summary = "no Cys liability signal" if not n else f"{n} flag(s)"
        print(f"\n  {_short(result.seq_id)}   {result.length} aa   "
              f"{result.status}   {summary}")

        def fmt(value):
            return "n/a" if value is None else str(value)

        print(f"\n      {'measure':<26}{'value':>10}   threshold")
        rows = [
            ("Cys count", str(m["n_cys"]), f"< {max_cys_count}"),
            ("Cys fraction", f"{m['cys_fraction']:.3f}",
             f"<= {max_cys_fraction:.2f}"),
            ("Cys parity", m["cys_parity"], "even if all Cys pair"),
            ("closest Cys spacing", fmt(m["closest_cys_spacing"]),
             f"> {close_pair_distance}"),
            (f"max Cys in w={window}", str(m[f"max_cys_in_w{window}"]),
             f"<= {max_cys_in_window}"),
        ]
        for name, value, bound in rows:
            print(f"      {name:<26}{value:>10}   {bound}")

        if result.cys_positions:
            shown = " ".join(str(p) for p in result.cys_positions[:20])
            more = f" +{len(result.cys_positions) - 20} more" if len(result.cys_positions) > 20 else ""
            print(f"\n      Cys positions: {shown}{more}")

        for issue in result.issues:
            print(f"\n    [{issue.severity}] {issue.label}")
            for line in textwrap.wrap(issue.reason, width=76,
                                      initial_indent="           ",
                                      subsequent_indent="           "):
                print(line)
            if show_context and issue.positions:
                for line in render_positions(result.sequence, issue.positions, flank=5):
                    print(line)

    n_warn = sum(1 for r in results if r.issues)
    print(f"\n  {'-' * 44}")
    print(f"  {len(results) - n_warn} clean   {n_warn} with a Cys pattern signal")
    print("  Advisory only: sequence alone cannot decide whether Cys is intended.\n")


def write_csv(results: list[CysPatternResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not results:
        return
    columns = list(results[0].metrics)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sequence_id", "status", "n_flags", *columns,
                         "cys_positions", "dense_regions", "close_pairs", "flags"])
        for r in results:
            writer.writerow([
                r.seq_id, r.status, len(r.issues),
                *["" if r.metrics[c] is None else r.metrics[c] for c in columns],
                " ".join(str(p) for p in r.cys_positions),
                "; ".join(f"{a}-{b}" for a, b, _ in r.dense_regions),
                "; ".join(f"{a}/{b}" for a, b in r.close_pairs),
                "; ".join(i.label for i in r.issues),
            ])


def write_json(results: list[CysPatternResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "caveat": "Cys pattern is sequence-only annotation. It cannot infer exposure, "
                  "pairing, metal binding, redox role or catalytic intent.",
        "thresholds_caveat": "All thresholds are provisional review hints for ordinary "
                             "soluble designs; use task metadata for disulfide or "
                             "functional-Cys cases.",
        "summary": {
            "n_sequences": len(results),
            "n_flagged": sum(1 for r in results if r.issues),
        },
        "sequences": [
            {
                "sequence_id": r.seq_id,
                "length": r.length,
                "status": r.status,
                "metrics": r.metrics,
                "cys_positions": r.cys_positions,
                "dense_regions": [
                    {"start": a, "end": b, "max_cys": c}
                    for a, b, c in r.dense_regions
                ],
                "close_pairs": [{"first": a, "second": b} for a, b in r.close_pairs],
                "flags": [
                    {"check": i.check, "severity": i.severity, "label": i.label,
                     "reason": i.reason, "positions": i.positions}
                    for i in r.issues
                ],
            }
            for r in results
        ],
    }
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("fasta", nargs="*", type=Path, help="FASTA file(s)")
    parser.add_argument("--seq", help="check a single sequence from the command line")
    parser.add_argument("--max-cys-count", type=int, default=DEFAULT_MAX_CYS_COUNT,
                        help=f"flag Cys count at or above this "
                             f"(default {DEFAULT_MAX_CYS_COUNT})")
    parser.add_argument("--max-cys-fraction", type=float,
                        default=DEFAULT_MAX_CYS_FRACTION,
                        help=f"flag Cys fraction above this "
                             f"(default {DEFAULT_MAX_CYS_FRACTION})")
    parser.add_argument("--window", type=int, default=DEFAULT_CLUSTER_WINDOW,
                        help=f"window for local Cys density "
                             f"(default {DEFAULT_CLUSTER_WINDOW})")
    parser.add_argument("--max-cys-in-window", type=int,
                        default=DEFAULT_MAX_CYS_IN_WINDOW,
                        help=f"flag windows with more than this many Cys "
                             f"(default {DEFAULT_MAX_CYS_IN_WINDOW})")
    parser.add_argument("--close-pair-distance", type=int,
                        default=DEFAULT_CLOSE_PAIR_DISTANCE,
                        help=f"flag adjacent Cys pairs within this sequence distance "
                             f"(default {DEFAULT_CLOSE_PAIR_DISTANCE})")
    parser.add_argument("--disulfide-expected", action="store_true",
                        help="this design is meant to be disulfide-stabilised or "
                             "secreted: still measure Cys-rich patterns, but only "
                             "flag odd Cys count")
    parser.add_argument("--cys-expected", action="store_true",
                        help="free/reactive/metal-binding/catalytic Cys is part of "
                             "the task: still measure, but do not flag Cys patterns")
    parser.add_argument("--no-context", action="store_true")
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 if any sequence is flagged")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    if not args.fasta and not args.seq:
        parser.error("give at least one FASTA file, or --seq")
    if args.window < 2:
        parser.error("--window must be at least 2")
    if args.max_cys_count < 1:
        parser.error("--max-cys-count must be at least 1")
    if args.max_cys_in_window < 1:
        parser.error("--max-cys-in-window must be at least 1")
    if args.close_pair_distance < 1:
        parser.error("--close-pair-distance must be at least 1")

    records: list[Record] = []
    if args.seq:
        records.append(
            Record(seq_id="<--seq>", sequence=_WHITESPACE.sub("", args.seq).upper())
        )
    for path in args.fasta:
        try:
            records.extend(read_fasta(path))
        except (OSError, ValueError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2

    for record in records:
        bad = set(record.sequence) - STANDARD_AA
        if bad:
            print(f"error: {record.seq_id} contains non-standard residue(s) "
                  f"{''.join(sorted(bad))} -- run basic_validity.py first",
                  file=sys.stderr)
            return 2

    results = [
        analyse(
            r,
            args.max_cys_count,
            args.max_cys_fraction,
            args.window,
            args.max_cys_in_window,
            args.close_pair_distance,
            args.disulfide_expected,
            args.cys_expected,
        )
        for r in records
    ]

    if not args.quiet:
        print_report(results, args.window, args.max_cys_count,
                     args.max_cys_fraction, args.max_cys_in_window,
                     args.close_pair_distance, show_context=not args.no_context)
    if args.csv:
        write_csv(results, args.csv)
    if args.json:
        write_json(results, args.json)

    if args.strict:
        return 1 if any(r.issues for r in results) else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
