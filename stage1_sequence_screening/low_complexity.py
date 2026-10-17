#!/usr/bin/env python3
"""Step 1 subcheck -- low-complexity detection.

Finds stretches that are repetitive or built from too few residue types:
homopolymer runs, periodic repeats like GSGSGSGS, and windows with low Shannon
entropy.

Like the other Step 1 liability subchecks this is a PROFILING step, not a gate.
Every flag is a WARN, nothing fails, and the measured values are the real output.

LOW COMPLEXITY MUST NEVER BE A HARD FAIL, and not merely because this layer is
advisory by convention -- because low complexity is often CORRECT:

  * natural proteins genuinely contain low-complexity segments, and two of the
    real proteins in this repo sit below SEG's own default trigger;
  * designed constructs legitimately contain linkers -- (GGGGS)n and (EAAAK)n
    are standard practice and are maximally periodic by design;
  * flexible loops, disordered termini and some binding motifs are low-complexity
    for functional reasons.

So the question this step answers is "does this sequence contain an unusual
region worth looking at?" -- never "is this design unacceptable?". Use --strict
only when a specific run wants flags to be fatal.

Measures
  * sliding-window Shannon entropy, H = -sum(p_i * log2 p_i), at two window sizes
  * the lowest-entropy window and the regions below a threshold
  * longest homopolymer run
  * strongest periodic repeat, with its period (2 = dipeptide, 3 = tripeptide)
  * fewest distinct residues in any 15-aa window

Usage
    conda run -n base python stage1_sequence_screening/low_complexity.py designs.fasta
    conda run -n base python stage1_sequence_screening/low_complexity.py --seq MKVLAA...
    conda run -n base python stage1_sequence_screening/low_complexity.py designs.fasta --csv lc.csv --json lc.json
    conda run -n base python stage1_sequence_screening/low_complexity.py designs.fasta --window 20 --min-entropy 2.0
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import textwrap
from collections import Counter
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
# Defaults
#
# CALIBRATION, and the reason these are not tighter. Measured on the four real
# sequences in this repo -- three well-behaved (one design, two natural parent
# proteins) and one pathological design with a 35-residue poly-E tract:
#
#                    minH(12)   minH(20)   run   repeat   distinct(15)
#   three good        2.05-2.42  2.80-3.01   2-3  .19-.29       6-7
#   poly-E design      0.00       0.00       35    1.00          1
#
# Every default below sits in the gap, with margin on the natural side. The
# previous version of this layer warned below 3.0 bits at window 20, which is
# ABOVE the natural minimum of 2.80 -- it flagged real, crystallised, expressed
# proteins. That single number is why this module exists in its current form.
#
# Note these are deliberately more conservative than SEG, whose default trigger
# is 2.2 bits at window 12. SEG's purpose is to MASK low-complexity regions for
# database search, and it finds them in most proteins -- two of our three good
# sequences sit below its trigger. Masking and rejecting are different jobs and
# want different thresholds.
#
# Like the composition bands, these remain provisional pending a matched background;
# see ../README.md -> "PENDING -- blocked on the reference database" -> P1.
# --------------------------------------------------------------------------

DEFAULT_WINDOW = 12             # primary window for entropy flagging
SECONDARY_WINDOW = 20           # also reported, for context
DEFAULT_MIN_ENTROPY = 1.5       # bits; naturals >= 2.05 at window 12
DEFAULT_MAX_HOMOPOLYMER = 6     # naturals reach 3
DEFAULT_MAX_REPEAT = 0.60       # fraction of a window matching at some period
DEFAULT_REPEAT_WINDOW = 24

# Periodicities checked: 1 through DEFAULT_MAX_PERIOD inclusive.
#   1 homopolymer   AAAAAA
#   2 dipeptide     GSGSGS
#   3 tripeptide    GPGGPGGPG, collagen-like GXY
#   4 tetrapeptide  GGGS linkers
#   5 pentapeptide  (GGGGS)n and (EAAAK)n -- the two commonest engineered linkers
#   6 hexapeptide   headroom, and (GGGGS)+1 framings
#
# Checking only short periods misses real cases: (EAAAK)x8 scores 0.43 at p<=4
# but 1.00 at p<=6, so a period-4 limit would catch it only incidentally via
# entropy and would misreport the mechanism. Raising the limit is close to free:
# across the four real sequences in this repo the strongest spurious repeat rises
# only from 0.23 to 0.26 as the limit goes from 2 to 8, against a 0.60 threshold.
DEFAULT_MAX_PERIOD = 6
DEFAULT_DISTINCT_WINDOW = 15
DEFAULT_MIN_DISTINCT = 4        # naturals have 6-7 distinct residues per 15 aa

PERIOD_NAMES = {1: "homopolymer", 2: "dipeptide", 3: "tripeptide",
                4: "tetrapeptide", 5: "pentapeptide", 6: "hexapeptide"}


@dataclass
class LowComplexityResult:
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


def sequence_entropy(seq: str) -> float:
    """Shannon entropy in bits: H = -sum(p_i * log2 p_i).

    `+ 0.0` normalises the -0.0 that a single-residue window otherwise produces,
    which would print as "-0.00".
    """
    if not seq:
        return 0.0
    counts = Counter(seq)
    n = len(seq)
    return -sum((c / n) * math.log2(c / n) for c in counts.values()) + 0.0


def window_entropies(seq: str, window: int) -> list[float]:
    """Entropy of every window. Sequences shorter than the window score whole.

    Scoring a short sequence whole rather than skipping it means a 15-residue
    input is still measured instead of silently passing.
    """
    if len(seq) <= window:
        return [sequence_entropy(seq)]
    return [sequence_entropy(seq[i : i + window]) for i in range(len(seq) - window + 1)]


def min_window_entropy(seq: str, window: int = DEFAULT_WINDOW) -> tuple[float, int]:
    """Lowest windowed entropy and its 1-based start position."""
    values = window_entropies(seq, window)
    best = min(range(len(values)), key=lambda i: values[i])
    return values[best], best + 1


def low_complexity_regions(
    seq: str, window: int = DEFAULT_WINDOW, threshold: float = DEFAULT_MIN_ENTROPY
) -> list[tuple[int, int, float]]:
    """Merged spans whose windows fall below `threshold`.

    Returns (start, end, min_entropy) with 1-based inclusive bounds. Overlapping
    low windows are merged so one tract is reported once rather than once per
    offset.
    """
    if len(seq) <= window:
        h = sequence_entropy(seq)
        return [(1, len(seq), h)] if h < threshold and seq else []

    regions: list[list] = []
    for i, h in enumerate(window_entropies(seq, window)):
        if h >= threshold:
            continue
        start, end = i + 1, i + window
        if regions and start <= regions[-1][1] + 1:
            regions[-1][1] = max(regions[-1][1], end)
            regions[-1][2] = min(regions[-1][2], h)
        else:
            regions.append([start, end, h])
    return [(int(a), int(b), float(c)) for a, b, c in regions]


def longest_homopolymer(seq: str) -> tuple[int, str, int]:
    """Longest run of one residue: (length, residue, 1-based start)."""
    if not seq:
        return 0, "", 0
    best_len, best_aa, best_start = 1, seq[0], 1
    cur_len, cur_start = 1, 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur_len += 1
        else:
            cur_len, cur_start = 1, i + 1
        if cur_len > best_len:
            best_len, best_aa, best_start = cur_len, seq[i], cur_start
    return best_len, best_aa, best_start


def max_periodic_repeat(
    seq: str,
    window: int = DEFAULT_REPEAT_WINDOW,
    max_period: int = DEFAULT_MAX_PERIOD,
) -> tuple[float, int, int]:
    """Strongest periodic repeat: (fraction, period, 1-based start).

    For each period p, mark every position matching the residue p ahead, then
    slide a rolling sum to find the densest window. This catches GSGSGSGS
    (period 2) and similar, which homopolymer and entropy checks both miss --
    GSGSGSGS has no run longer than 1 and a respectable two-symbol entropy.

    O(max_period * n) via the rolling sum, rather than O(max_period * n * window).
    """
    n = len(seq)
    if n < 4:
        return 0.0, 0, 0
    window = min(window, n)
    best = (0.0, 0, 0)

    for period in range(1, min(max_period, window - 1) + 1):
        matches = [1 if seq[i] == seq[i + period] else 0 for i in range(n - period)]
        if not matches:
            continue
        span = min(window - period, len(matches))
        if span <= 0:
            continue
        running = sum(matches[:span])
        candidates = [(running / span, 1)]
        for i in range(1, len(matches) - span + 1):
            running += matches[i + span - 1] - matches[i - 1]
            candidates.append((running / span, i + 1))
        fraction, start = max(candidates)
        if fraction > best[0]:
            best = (fraction, period, start)
    return best


def min_distinct_residues(
    seq: str, window: int = DEFAULT_DISTINCT_WINDOW
) -> tuple[int, int]:
    """Fewest distinct residue types in any window: (count, 1-based start)."""
    if not seq:
        return 0, 0
    if len(seq) <= window:
        return len(set(seq)), 1
    best, best_i = len(STANDARD_AA) + 1, 1
    for i in range(len(seq) - window + 1):
        distinct = len(set(seq[i : i + window]))
        if distinct < best:
            best, best_i = distinct, i + 1
    return best, best_i


def measure(seq: str, window: int, distinct_window: int,
            repeat_window: int, max_period: int) -> dict:
    """Every raw number for one sequence. This table is the real output."""
    h_primary, h_at = min_window_entropy(seq, window)
    h_secondary, _ = min_window_entropy(seq, SECONDARY_WINDOW)
    run_len, run_aa, run_at = longest_homopolymer(seq)
    rep_frac, rep_period, rep_at = max_periodic_repeat(seq, repeat_window, max_period)
    distinct, distinct_at = min_distinct_residues(seq, distinct_window)
    return {
        "length": len(seq),
        "entropy_whole": round(sequence_entropy(seq), 3),
        f"min_entropy_w{window}": round(h_primary, 3),
        f"min_entropy_w{window}_at": h_at,
        f"min_entropy_w{SECONDARY_WINDOW}": round(h_secondary, 3),
        "max_entropy_possible": round(math.log2(min(window, len(seq) or 1)), 3),
        "longest_homopolymer": run_len,
        "longest_homopolymer_aa": run_aa,
        "longest_homopolymer_at": run_at,
        "max_repeat_fraction": round(rep_frac, 3),
        "max_repeat_period": rep_period,
        "max_repeat_at": rep_at,
        f"min_distinct_w{distinct_window}": distinct,
        f"min_distinct_w{distinct_window}_at": distinct_at,
    }


# --------------------------------------------------------------------------
# Checks -- advisory only
# --------------------------------------------------------------------------


def check_entropy(seq: str, window: int, threshold: float) -> list[Issue]:
    regions = low_complexity_regions(seq, window, threshold)
    if not regions:
        return []
    worst = min(regions, key=lambda r: r[2])
    spans = ", ".join(f"{a}-{b}" for a, b, _ in regions[:6])
    more = f" (+{len(regions) - 6} more)" if len(regions) > 6 else ""
    positions = [p for a, b, _ in regions for p in range(a, b + 1)]
    return [
        Issue(
            check="low_complexity_entropy",
            label=f"{len(regions)} low-complexity region(s): {spans}{more}",
            reason=f"Lowest windowed entropy {worst[2]:.2f} bits over {window} aa "
                   f"(max possible {math.log2(window):.2f}), below the provisional "
                   f"{threshold:.2f}. Low-entropy tracts are poorly expressed, often "
                   f"disordered, and are a common artefact of generative models.",
            positions=positions,
            severity="WARN",
        )
    ]


def check_homopolymer(seq: str, limit: int) -> list[Issue]:
    length, aa, start = longest_homopolymer(seq)
    if length < limit:
        return []
    return [
        Issue(
            check="low_complexity_homopolymer",
            label=f"{length}-residue poly-{aa} run",
            reason=f"A run of {length} identical residues starting at {start}. Natural "
                   f"proteins in this repo reach 3. Long homopolymers impair expression "
                   f"and are frequently a generative artefact.",
            positions=list(range(start, start + length)),
            severity="WARN",
        )
    ]


def check_repeat(seq: str, window: int, max_period: int, limit: float) -> list[Issue]:
    fraction, period, start = max_periodic_repeat(seq, window, max_period)
    if fraction < limit or period == 0:
        return []
    name = PERIOD_NAMES.get(period, f"period-{period}")
    end = min(len(seq), start + window - 1)
    return [
        Issue(
            check="low_complexity_repeat",
            label=f"{name} repeat, {fraction:.0%} of a {window}-aa window",
            reason=f"Period {period} ({name}) repetition starting at {start}. Periodic "
                   f"repeats such as GSGSGSGS evade both the homopolymer and entropy "
                   f"checks -- no run exceeds one residue and two symbols give "
                   f"respectable entropy -- so they need their own test.",
            positions=list(range(start, end + 1)),
            severity="WARN",
        )
    ]


def check_distinct(seq: str, window: int, minimum: int) -> list[Issue]:
    distinct, start = min_distinct_residues(seq, window)
    if distinct > minimum:
        return []
    end = min(len(seq), start + window - 1)
    alphabet = "".join(sorted(set(seq[start - 1 : end])))
    return [
        Issue(
            check="low_complexity_alphabet",
            label=f"only {distinct} distinct residue(s) in a {window}-aa window",
            reason=f"Positions {start}-{end} use just {alphabet}. Natural proteins in "
                   f"this repo use 6-7 distinct residues per {window} aa.",
            positions=list(range(start, end + 1)),
            severity="WARN",
        )
    ]


def analyse(
    record: Record,
    window: int = DEFAULT_WINDOW,
    min_entropy: float = DEFAULT_MIN_ENTROPY,
    max_homopolymer: int = DEFAULT_MAX_HOMOPOLYMER,
    max_repeat: float = DEFAULT_MAX_REPEAT,
    repeat_window: int = DEFAULT_REPEAT_WINDOW,
    max_period: int = DEFAULT_MAX_PERIOD,
    distinct_window: int = DEFAULT_DISTINCT_WINDOW,
    min_distinct: int = DEFAULT_MIN_DISTINCT,
) -> LowComplexityResult:
    seq = record.sequence
    return LowComplexityResult(
        seq_id=record.seq_id,
        length=len(seq),
        sequence=seq,
        metrics=measure(seq, window, distinct_window, repeat_window, max_period),
        regions=low_complexity_regions(seq, window, min_entropy),
        issues=[
            *check_entropy(seq, window, min_entropy),
            *check_homopolymer(seq, max_homopolymer),
            *check_repeat(seq, repeat_window, max_period, max_repeat),
            *check_distinct(seq, distinct_window, min_distinct),
        ],
    )


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def _short(seq_id: str, width: int = 34) -> str:
    head = seq_id.split("|")[0]
    return head if len(head) <= width else head[: width - 1] + "…"


def print_report(
    results: list[LowComplexityResult],
    window: int = DEFAULT_WINDOW,
    min_entropy: float = DEFAULT_MIN_ENTROPY,
    show_context: bool = True,
) -> None:
    print(f"\nStep 1 subcheck -- low complexity   ({len(results)} sequence(s))")

    for result in results:
        m = result.metrics
        n = len(result.issues)
        summary = "no low-complexity signal" if not n else f"{n} flag(s)"
        print(f"\n  {_short(result.seq_id)}   {result.length} aa   "
              f"{result.status}   {summary}")

        print(f"\n      {'measure':<26}{'value':>9}   provisional")
        rows = [
            (f"min entropy (w={window})", f"{m[f'min_entropy_w{window}']:.2f} bits",
             f">= {min_entropy:.2f}"),
            (f"min entropy (w={SECONDARY_WINDOW})",
             f"{m[f'min_entropy_w{SECONDARY_WINDOW}']:.2f} bits", "context only"),
            ("whole-sequence entropy", f"{m['entropy_whole']:.2f} bits", ""),
            ("longest homopolymer",
             f"{m['longest_homopolymer']}x {m['longest_homopolymer_aa']}",
             f"< {DEFAULT_MAX_HOMOPOLYMER}"),
            ("strongest repeat",
             f"{m['max_repeat_fraction']:.0%} p={m['max_repeat_period']}",
             f"< {DEFAULT_MAX_REPEAT:.0%}"),
            (f"min distinct (w={DEFAULT_DISTINCT_WINDOW})",
             str(m[f'min_distinct_w{DEFAULT_DISTINCT_WINDOW}']),
             f"> {DEFAULT_MIN_DISTINCT}"),
        ]
        for name, value, bound in rows:
            print(f"      {name:<26}{value:>9}   {bound}")
        print(f"\n      thresholds are provisional outlier hints, not pass/fail")

        for issue in result.issues:
            print(f"\n    [{issue.severity}] {issue.label}")
            for line in textwrap.wrap(issue.reason, width=76,
                                      initial_indent="           ",
                                      subsequent_indent="           "):
                print(line)
            if show_context and issue.positions:
                for line in render_positions(result.sequence, issue.positions, flank=6):
                    print(line)

    n_warn = sum(1 for r in results if r.issues)
    print(f"\n  {'-' * 44}")
    print(f"  {len(results) - n_warn} clean   {n_warn} with low-complexity signal")
    print("  Advisory only: nothing here fails a design.\n")


def write_csv(results: list[LowComplexityResult], path: Path) -> None:
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


def write_json(results: list[LowComplexityResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "method": "sliding-window Shannon entropy H = -sum(p_i log2 p_i), plus "
                  "homopolymer, periodic-repeat and alphabet-size checks",
        "thresholds_caveat": "provisional outlier hints, not pass/fail; calibrated "
                             "against four sequences only. See ../README.md PENDING P1.",
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
                "low_complexity_regions": [
                    {"start": a, "end": b, "min_entropy": round(c, 3)}
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
                        help=f"entropy window (default {DEFAULT_WINDOW}; "
                             f"w={SECONDARY_WINDOW} is always reported too)")
    parser.add_argument("--min-entropy", type=float, default=DEFAULT_MIN_ENTROPY,
                        help=f"flag windows below this many bits "
                             f"(default {DEFAULT_MIN_ENTROPY})")
    parser.add_argument("--max-homopolymer", type=int, default=DEFAULT_MAX_HOMOPOLYMER,
                        help=f"flag runs at or above this (default {DEFAULT_MAX_HOMOPOLYMER})")
    parser.add_argument("--max-repeat", type=float, default=DEFAULT_MAX_REPEAT,
                        help=f"flag periodic repeats at or above this fraction "
                             f"(default {DEFAULT_MAX_REPEAT})")
    parser.add_argument("--min-distinct", type=int, default=DEFAULT_MIN_DISTINCT,
                        help=f"flag windows with this many distinct residues or fewer "
                             f"(default {DEFAULT_MIN_DISTINCT})")
    parser.add_argument("--no-context", action="store_true",
                        help="do not print sequence excerpts")
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

    results = [
        analyse(r, window=args.window, min_entropy=args.min_entropy,
                max_homopolymer=args.max_homopolymer, max_repeat=args.max_repeat,
                min_distinct=args.min_distinct)
        for r in records
    ]

    if not args.quiet:
        print_report(results, args.window, args.min_entropy,
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
