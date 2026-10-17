#!/usr/bin/env python3
"""Step 1 subcheck -- hydrophobicity and hydrophobic stretches.

Global average hydrophobicity plus, more usefully, LOCAL maxima. A soluble
protein can have a perfectly ordinary GRAVY while hiding a single very
hydrophobic 20-residue stretch, and that stretch is what causes inclusion
bodies, aggregation, or an accidental membrane-spanning segment. The global
number averages it away; the sliding window does not.

Advisory only, like the other Step 1 liability subchecks: every flag is a WARN,
nothing fails, and the measured values are the real output. Use --strict to make
flags fatal for a run.

Measures
  * max windowed hydropathy at w=19 -- the PRIMARY signal, flagged as a TM-LIKE
    stretch (a hint, never a prediction), with merged regions and positions
  * min windowed hydropathy -- an extreme-hydrophilic OUTLIER flag, which has no
    equivalent literature basis and is weak evidence on its own
  * GRAVY -- descriptive only; it dilutes with length and cannot find a patch
  * longest hydrophobic run -- an auxiliary composition descriptor, explicitly
    NOT a membrane-helix detector (real TM helices reach runs of only 3-5)

Set expected_soluble=False (CLI: --membrane-expected) for designs that are meant
to contain a membrane anchor or signal peptide; the stretch is then measured and
reported but not treated as a liability.

Usage
    python stage1_sequence_screening/hydrophobicity.py designs.fasta
    python stage1_sequence_screening/hydrophobicity.py --seq MKVLAA...
    python stage1_sequence_screening/hydrophobicity.py designs.fasta --window 19 --max-hydropathy 1.6
    python stage1_sequence_screening/hydrophobicity.py designs.fasta --csv h.csv --json h.json
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
# Kyte & Doolittle (1982) J. Mol. Biol. 157:105-132, hydropathy index.
# The conventional scale for GRAVY and for hydrophobic-window analysis.
# --------------------------------------------------------------------------
KYTE_DOOLITTLE = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

# --------------------------------------------------------------------------
# Thresholds, and what each one is actually worth.
#
# These four measurements are NOT equally well founded. Ranked by how much
# interpretive weight they can carry:
#
# 1. max hydropathy over w=19  -- the primary signal, and the only one with a
#    real literature basis. Kyte & Doolittle's own reading is that a 19-residue
#    mean above roughly +1.6 is LIKELY to correspond to a membrane-spanning
#    segment -- a strong hint, not a determination. The classic rule has known
#    false positives and false negatives, so this flag is deliberately named a
#    "TM-like hydrophobic stretch" and never a TM prediction. It must not become
#    a hard cutoff.
#
# 2. min hydropathy over w=19  -- an outlier flag with NO equivalent literature
#    interpretation. There is no published "-3.0 means X". Extremely hydrophilic
#    is not the same as wrong: natural acidic regions, disordered tails and
#    charged linkers are all legitimately very hydrophilic. Hence the name
#    `extreme_hydrophilic_outlier`, and hence it should be read as one piece of
#    evidence among several -- convincing when it coincides with an extreme
#    negative fraction (composition subcheck) and low entropy (low-complexity
#    subcheck), weak on its own.
#
# 3. GRAVY -- descriptive. Reported because it is conventional and cheap, but the
#    host-length experiment in the tests shows it dilutes with length and so
#    cannot be relied on to find a local patch.
#
# 4. longest hydrophobic run -- an AUXILIARY COMPOSITION DESCRIPTOR, and
#    explicitly NOT a membrane-helix detector. Measured against four real TM
#    helices:
#
#        glycophorin A               max19 2.67   run 4
#        glycophorin A (full TM)     max19 2.67   run 4
#        bacteriorhodopsin helix A   max19 1.84   run 3
#        influenza M2 TM             max19 2.35   run 5
#
#    Every one is caught by max19; NOT ONE reaches a run of 8. Membrane helices
#    are routinely broken by Gly, Thr and Ser -- glycophorin A's is interrupted
#    by T, G, G, G, T, S -- so consecutive-run length simply does not measure
#    this. It stays because uninterrupted runs are a genuine aggregation
#    liability in their own right, but it carries no TM interpretation.
#
# Calibration context for the two judgement-based bounds (run, hydrophilic
# floor), from the four real sequences in this repo:
#
#                          GRAVY          max(w=19)   longest hydrophobic run
#   four real sequences   -1.08..-0.25    0.52-1.18            3-5
#
# Both inherit ../README.md -> PENDING -> P1/P2: provisional, pending a matched
# background.
# --------------------------------------------------------------------------

DEFAULT_WINDOW = 19             # Kyte-Doolittle membrane window
SECONDARY_WINDOW = 9            # short patches, reported for context
DEFAULT_MAX_HYDROPATHY = 1.6    # KD: a 19-aa mean above this is LIKELY membrane-spanning
DEFAULT_MAX_GRAVY = 0.0         # soluble proteins are normally negative
DEFAULT_MAX_HYDROPHOBIC_RUN = 8   # aggregation descriptor only -- NOT a TM detector
DEFAULT_MIN_HYDROPATHY = -3.0     # outlier hint only; no literature meaning

# A residue counts as hydrophobic when its hydropathy is positive:
# I V L F C M A. Spelled out rather than left implicit.
HYDROPHOBIC_AA = frozenset(aa for aa, h in KYTE_DOOLITTLE.items() if h > 0)


@dataclass
class HydrophobicityResult:
    seq_id: str
    length: int
    sequence: str
    metrics: dict
    regions: list[tuple[int, int, float]] = field(default_factory=list)
    issues: list[Issue] = field(default_factory=list)

    @property
    def status(self) -> str:
        return "WARN" if self.issues else "PASS"


# --------------------------------------------------------------------------
# Measurements
# --------------------------------------------------------------------------


def gravy(seq: str) -> float:
    """Grand average of hydropathy: (1/L) * sum(h_i)."""
    if not seq:
        return 0.0
    return sum(KYTE_DOOLITTLE[aa] for aa in seq) / len(seq)


def window_hydropathy(seq: str, window: int) -> list[float]:
    """Mean hydropathy of every window, via a rolling sum.

    Sequences shorter than the window are scored whole rather than skipped, so a
    15-aa input is still measured instead of silently passing.
    """
    if len(seq) <= window:
        return [gravy(seq)] if seq else []
    values = [KYTE_DOOLITTLE[aa] for aa in seq]
    running = sum(values[:window])
    out = [running / window]
    for i in range(1, len(values) - window + 1):
        running += values[i + window - 1] - values[i - 1]
        out.append(running / window)
    return out


def max_window_hydropathy(seq: str, window: int = DEFAULT_WINDOW) -> tuple[float, int]:
    """Highest windowed mean hydropathy and its 1-based start position."""
    values = window_hydropathy(seq, window)
    if not values:
        return 0.0, 0
    best = max(range(len(values)), key=lambda i: values[i])
    return values[best], best + 1


def min_window_hydropathy(seq: str, window: int = DEFAULT_WINDOW) -> tuple[float, int]:
    """Lowest windowed mean hydropathy -- the strongly hydrophilic extreme."""
    values = window_hydropathy(seq, window)
    if not values:
        return 0.0, 0
    best = min(range(len(values)), key=lambda i: values[i])
    return values[best], best + 1


def hydrophobic_regions(
    seq: str,
    window: int = DEFAULT_WINDOW,
    threshold: float = DEFAULT_MAX_HYDROPATHY,
) -> list[tuple[int, int, float]]:
    """Merged spans whose windows exceed `threshold`.

    Returns (start, end, max_hydropathy) with 1-based inclusive bounds.
    Overlapping windows merge so one stretch is reported once rather than once
    per offset.
    """
    values = window_hydropathy(seq, window)
    if not values:
        return []
    if len(seq) <= window:
        return [(1, len(seq), values[0])] if values[0] > threshold else []

    regions: list[list] = []
    for i, value in enumerate(values):
        if value <= threshold:
            continue
        start, end = i + 1, i + window
        if regions and start <= regions[-1][1] + 1:
            regions[-1][1] = max(regions[-1][1], end)
            regions[-1][2] = max(regions[-1][2], value)
        else:
            regions.append([start, end, value])
    return [(int(a), int(b), float(c)) for a, b, c in regions]


def longest_hydrophobic_run(seq: str) -> tuple[int, int]:
    """Longest run of consecutive hydrophobic residues: (length, 1-based start)."""
    best_len, best_start = 0, 0
    cur_len, cur_start = 0, 0
    for i, aa in enumerate(seq):
        if aa in HYDROPHOBIC_AA:
            if cur_len == 0:
                cur_start = i + 1
            cur_len += 1
            if cur_len > best_len:
                best_len, best_start = cur_len, cur_start
        else:
            cur_len = 0
    return best_len, best_start


def measure(seq: str, window: int) -> dict:
    """Every raw number for one sequence. This table is the real output."""
    max_h, max_at = max_window_hydropathy(seq, window)
    max_short, max_short_at = max_window_hydropathy(seq, SECONDARY_WINDOW)
    min_h, min_at = min_window_hydropathy(seq, window)
    run_len, run_at = longest_hydrophobic_run(seq)
    return {
        "length": len(seq),
        "gravy": round(gravy(seq), 3),
        f"max_hydropathy_w{window}": round(max_h, 3),
        f"max_hydropathy_w{window}_at": max_at,
        f"max_hydropathy_w{SECONDARY_WINDOW}": round(max_short, 3),
        f"max_hydropathy_w{SECONDARY_WINDOW}_at": max_short_at,
        f"min_hydropathy_w{window}": round(min_h, 3),
        f"min_hydropathy_w{window}_at": min_at,
        "longest_hydrophobic_run": run_len,
        "longest_hydrophobic_run_at": run_at,
        "hydrophobic_fraction": round(
            sum(1 for aa in seq if aa in HYDROPHOBIC_AA) / len(seq), 4
        ) if seq else 0.0,
    }


# --------------------------------------------------------------------------
# Checks -- advisory only
# --------------------------------------------------------------------------


def check_hydrophobic_stretch(
    seq: str, window: int, threshold: float, expected_soluble: bool = True
) -> list[Issue]:
    """The check this step exists for: a local maximum a global average hides.

    Reported as a TM-LIKE stretch, never as a TM prediction -- the underlying
    rule has known false positives and false negatives. When the design is not
    expected to be soluble (a membrane anchor, a signal peptide, an engineered
    secretion sequence), the stretch is still measured and reported but is not
    a liability, so no flag is raised.
    """
    if not expected_soluble:
        return []
    regions = hydrophobic_regions(seq, window, threshold)
    if not regions:
        return []
    worst = max(regions, key=lambda r: r[2])
    spans = ", ".join(f"{a}-{b}" for a, b, _ in regions[:6])
    more = f" (+{len(regions) - 6} more)" if len(regions) > 6 else ""
    positions = [p for a, b, _ in regions for p in range(a, b + 1)]
    return [
        Issue(
            check="tm_like_hydrophobic_stretch",
            label=f"{len(regions)} TM-like hydrophobic stretch(es): {spans}{more}",
            reason=f"Peak mean hydropathy {worst[2]:.2f} over {window} aa, above "
                   f"{threshold:.1f} -- Kyte & Doolittle's threshold for a segment "
                   f"LIKELY to be membrane-spanning. This is a hint, not a prediction: "
                   f"the rule has known false positives and false negatives. In a "
                   f"design meant to be soluble it flags aggregation and "
                   f"inclusion-body risk, and a global GRAVY averages it away.",
            positions=positions,
            severity="WARN",
        )
    ]


def check_gravy(seq: str, limit: float) -> list[Issue]:
    value = gravy(seq)
    if value <= limit:
        return []
    return [
        Issue(
            check="hydrophobicity_gravy",
            label=f"GRAVY is {value:+.2f}, above {limit:+.1f}",
            reason="Soluble globular proteins are normally negative overall. A "
                   "positive GRAVY suggests a membrane protein, a poorly soluble "
                   "design, or an exposed hydrophobic surface. Descriptive only: "
                   "GRAVY dilutes with sequence length, so it cannot be relied on "
                   "to find a local patch -- that is what the window scan is for.",
            severity="WARN",
        )
    ]


def check_hydrophobic_run(seq: str, limit: int) -> list[Issue]:
    length, start = longest_hydrophobic_run(seq)
    if length < limit:
        return []
    return [
        Issue(
            check="hydrophobic_run",
            label=f"{length} consecutive hydrophobic residues",
            reason=f"An unbroken run of {length} residues from "
                   f"{''.join(sorted(HYDROPHOBIC_AA))} starting at {start}. Real "
                   f"sequences in this repo reach 5. Uninterrupted runs are an "
                   f"aggregation liability in their own right -- but this is a "
                   f"composition descriptor, NOT a membrane-helix detector: real TM "
                   f"helices are broken by Gly/Thr/Ser and reach runs of only 3-5.",
            positions=list(range(start, start + length)),
            severity="WARN",
        )
    ]


def check_hydrophilic_outlier(seq: str, window: int, floor: float) -> list[Issue]:
    """The opposite extreme -- an OUTLIER flag, not a defect.

    Unlike the +1.6 criterion this has no published interpretation behind it.
    Natural acidic regions, disordered tails and charged linkers are all
    legitimately very hydrophilic. Treat it as one piece of evidence: strong
    when it coincides with an extreme negative fraction (composition subcheck)
    and low entropy (low-complexity subcheck), weak on its own.
    """
    value, start = min_window_hydropathy(seq, window)
    if value >= floor:
        return []
    end = min(len(seq), start + window - 1)
    return [
        Issue(
            check="extreme_hydrophilic_outlier",
            label=f"extreme hydrophilic outlier, mean hydropathy {value:.2f}",
            reason=f"Positions {start}-{end} average {value:.2f} over {window} aa, "
                   f"below {floor:.1f}. Usually a charged or low-complexity tract. "
                   f"NOT a defect on its own -- natural acidic regions, disordered "
                   f"tails and charged linkers are legitimately this hydrophilic. "
                   f"Convincing only alongside an extreme negative fraction "
                   f"(composition subcheck) or low entropy "
                   f"(low-complexity subcheck).",
            positions=list(range(start, end + 1)),
            severity="WARN",
        )
    ]


def analyse(
    record: Record,
    window: int = DEFAULT_WINDOW,
    max_hydropathy: float = DEFAULT_MAX_HYDROPATHY,
    max_gravy: float = DEFAULT_MAX_GRAVY,
    max_hydrophobic_run: int = DEFAULT_MAX_HYDROPHOBIC_RUN,
    min_hydropathy: float = DEFAULT_MIN_HYDROPATHY,
    expected_soluble: bool = True,
) -> HydrophobicityResult:
    """Profile one sequence.

    `expected_soluble` is task metadata: set it False for a design that is meant
    to contain a membrane anchor, signal peptide or secretion sequence. The
    TM-like stretch is then still measured and reported, but not flagged -- for
    such a case a hydrophobic peak may be exactly correct.
    """
    seq = record.sequence
    return HydrophobicityResult(
        seq_id=record.seq_id,
        length=len(seq),
        sequence=seq,
        metrics=measure(seq, window),
        regions=hydrophobic_regions(seq, window, max_hydropathy),
        issues=[
            *check_hydrophobic_stretch(seq, window, max_hydropathy, expected_soluble),
            *check_gravy(seq, max_gravy),
            *check_hydrophobic_run(seq, max_hydrophobic_run),
            *check_hydrophilic_outlier(seq, window, min_hydropathy),
        ],
    )


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def _short(seq_id: str, width: int = 34) -> str:
    head = seq_id.split("|")[0]
    return head if len(head) <= width else head[: width - 1] + "…"


def print_report(
    results: list[HydrophobicityResult],
    window: int = DEFAULT_WINDOW,
    max_hydropathy: float = DEFAULT_MAX_HYDROPATHY,
    show_context: bool = True,
) -> None:
    print(f"\nStep 1 subcheck -- hydrophobicity   ({len(results)} sequence(s))")

    for result in results:
        m = result.metrics
        n = len(result.issues)
        summary = "no hydrophobic stretch found" if not n else f"{n} flag(s)"
        print(f"\n  {_short(result.seq_id)}   {result.length} aa   "
              f"{result.status}   {summary}")

        print(f"\n      {'measure':<28}{'value':>8}   threshold")
        rows = [
            ("GRAVY (whole sequence)", f"{m['gravy']:+.2f}",
             f"<= {DEFAULT_MAX_GRAVY:+.1f}"),
            (f"max hydropathy (w={window})", f"{m[f'max_hydropathy_w{window}']:+.2f}",
             f"<= {max_hydropathy:.1f}  (Kyte-Doolittle TM)"),
            (f"max hydropathy (w={SECONDARY_WINDOW})",
             f"{m[f'max_hydropathy_w{SECONDARY_WINDOW}']:+.2f}", "context only"),
            (f"min hydropathy (w={window})", f"{m[f'min_hydropathy_w{window}']:+.2f}",
             f">= {DEFAULT_MIN_HYDROPATHY:.1f}"),
            ("longest hydrophobic run", str(m["longest_hydrophobic_run"]),
             f"< {DEFAULT_MAX_HYDROPHOBIC_RUN}"),
            ("hydrophobic fraction", f"{m['hydrophobic_fraction']:.3f}", ""),
        ]
        for name, value, bound in rows:
            print(f"      {name:<28}{value:>8}   {bound}")
        print(f"\n      scale: Kyte & Doolittle (1982); w={window} is their "
              f"membrane window")

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
    print(f"  {len(results) - n_warn} clean   {n_warn} with a hydrophobicity signal")
    print("  Advisory only: nothing here fails a design.\n")


def write_csv(results: list[HydrophobicityResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not results:
        return
    columns = list(results[0].metrics)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sequence_id", "status", "n_flags", *columns,
                         "n_regions", "regions", "flags"])
        for r in results:
            writer.writerow([
                r.seq_id, r.status, len(r.issues),
                *[r.metrics[c] for c in columns],
                len(r.regions),
                "; ".join(f"{a}-{b}" for a, b, _ in r.regions),
                "; ".join(i.label for i in r.issues),
            ])


def write_json(results: list[HydrophobicityResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "scale": "Kyte & Doolittle (1982) J. Mol. Biol. 157:105-132",
        "membrane_criterion": f"mean hydropathy > {DEFAULT_MAX_HYDROPATHY} over "
                              f"{DEFAULT_WINDOW} residues",
        "thresholds_caveat": "the w=19 / 1.6 criterion is Kyte & Doolittle's own; the "
                             "hydrophobic-run and hydrophilic-floor bounds are "
                             "provisional, calibrated on four sequences. "
                             "See ../README.md PENDING P1/P2.",
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
                "hydrophobic_regions": [
                    {"start": a, "end": b, "max_hydropathy": round(c, 3)}
                    for a, b, c in r.regions
                ],
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
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW,
                        help=f"hydropathy window (default {DEFAULT_WINDOW}, the "
                             f"Kyte-Doolittle membrane window)")
    parser.add_argument("--max-hydropathy", type=float, default=DEFAULT_MAX_HYDROPATHY,
                        help=f"flag windows above this (default "
                             f"{DEFAULT_MAX_HYDROPATHY}, KD's membrane criterion)")
    parser.add_argument("--max-gravy", type=float, default=DEFAULT_MAX_GRAVY,
                        help=f"flag GRAVY above this (default {DEFAULT_MAX_GRAVY})")
    parser.add_argument("--max-hydrophobic-run", type=int,
                        default=DEFAULT_MAX_HYDROPHOBIC_RUN,
                        help=f"flag runs at or above this "
                             f"(default {DEFAULT_MAX_HYDROPHOBIC_RUN})")
    parser.add_argument("--min-hydropathy", type=float, default=DEFAULT_MIN_HYDROPATHY,
                        help=f"flag windows below this (default {DEFAULT_MIN_HYDROPATHY})")
    parser.add_argument("--membrane-expected", action="store_true",
                        help="this design is meant to contain a membrane anchor, "
                             "signal peptide or secretion sequence: still measure "
                             "TM-like stretches, but do not flag them")
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

    # The scale is defined only over the 20 standard residues; step 1 owns the rest.
    for record in records:
        bad = set(record.sequence) - STANDARD_AA
        if bad:
            print(f"error: {record.seq_id} contains non-standard residue(s) "
                  f"{''.join(sorted(bad))} -- run basic_validity.py first",
                  file=sys.stderr)
            return 2

    results = [
        analyse(r, args.window, args.max_hydropathy, args.max_gravy,
                args.max_hydrophobic_run, args.min_hydropathy,
                expected_soluble=not args.membrane_expected)
        for r in records
    ]

    if not args.quiet:
        print_report(results, args.window, args.max_hydropathy,
                     show_context=not args.no_context)
    if args.csv:
        write_csv(results, args.csv)
    if args.json:
        write_json(results, args.json)

    if args.strict:
        return 1 if any(r.issues for r in results) else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
