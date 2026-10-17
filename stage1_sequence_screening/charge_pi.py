#!/usr/bin/env python3
"""Step 1 subcheck -- net charge and isoelectric point.

Two tiers, deliberately:

  SIMPLE (always available, stdlib only)
      net charge   = (N_K + N_R) - (N_D + N_E)
      normalised   = that, divided by length

    Crude -- it ignores His, Cys, Tyr and the termini, and assumes full
    ionisation -- but it needs no dependency and is already enough to find
    population outliers, which is all this layer is for.

  STANDARD (requires Biopython)
      isoelectric point and charge at pH 7.0, via Bio.SeqUtils.ProtParam

    Proper Henderson-Hasselbalch with a standard pKa set. We deliberately do NOT
    hand-roll acid-base chemistry: a subtly wrong pKa table or a sign error is
    easy to write and hard to notice, and ProtParam is the conventional
    reference implementation. If Biopython is absent the simple tier still runs
    and the report says which numbers are missing -- it never silently
    substitutes a worse estimate for a better one.

Advisory only, like the other Step 1 liability subchecks: every flag is a WARN,
nothing fails.

Usage
    python stage1_sequence_screening/charge_pi.py designs.fasta
    python stage1_sequence_screening/charge_pi.py --seq MKVLAA...
    python stage1_sequence_screening/charge_pi.py designs.fasta --csv q.csv --json q.json
    python stage1_sequence_screening/charge_pi.py designs.fasta --require-biopython    # fail if absent
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import textwrap
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

try:
    from .basic_validity import Issue, Record, STANDARD_AA, _WHITESPACE, read_fasta
except ImportError:  # pragma: no cover - direct script execution path
    from basic_validity import Issue, Record, STANDARD_AA, _WHITESPACE, read_fasta

# Biopython is optional. Import once, record the outcome, and let every caller
# branch on the flag rather than re-trying the import.
try:  # pragma: no cover - trivial, and depends on the environment
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    HAS_BIOPYTHON = True
except ImportError:  # pragma: no cover
    ProteinAnalysis = None
    HAS_BIOPYTHON = False

BIOPYTHON_HINT = "pip install biopython"

POSITIVE = frozenset("KR")
NEGATIVE = frozenset("DE")

# --------------------------------------------------------------------------
# Thresholds.
#
# Calibrated on the real and control sequences in this repo:
#
#     well-behaved (4 real + lysozyme + ubiquitin)   q_norm  -0.024 .. +0.062
#     19_1Q0S design (poly-E tract)                  q_norm          -0.166
#     poly-E(40) / poly-K(40) inserts                q_norm  -0.262 / +0.226
#
# 0.10 sits in the gap with margin on the well-behaved side. Provisional like
# every other hand-set bound here: ../README.md -> PENDING -> P1/P2.
#
# The pI window is different in kind. A protein is often least soluble near its
# own pI, so a pI sitting on top of the working pH is a condition-specific buffer
# caution rather than a sequence-population outlier. It is not a design failure
# by itself: natural soluble proteins can trigger it.
# --------------------------------------------------------------------------

DEFAULT_MAX_ABS_CHARGE_DENSITY = 0.10
DEFAULT_WORKING_PH = 7.0
DEFAULT_PI_WINDOW = 0.5


@dataclass
class ChargeResult:
    seq_id: str
    length: int
    metrics: dict
    issues: list[Issue] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def status(self) -> str:
        return "WARN" if self.issues else "PASS"


# --------------------------------------------------------------------------
# Simple tier -- always available
# --------------------------------------------------------------------------


def net_charge_simple(seq: str) -> int:
    """(N_K + N_R) - (N_D + N_E).

    Counts only the four unambiguously charged residues at physiological pH.
    His is excluded on purpose: at ~10% protonated at pH 7 it is neither, and
    assigning it either way would bias the result.
    """
    counts = Counter(seq)
    return (counts["K"] + counts["R"]) - (counts["D"] + counts["E"])


def normalized_charge(seq: str) -> float:
    """Net charge per residue -- the length-independent form, for comparison."""
    return net_charge_simple(seq) / len(seq) if seq else 0.0


def charged_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    counts = Counter(seq)
    return sum(counts[aa] for aa in "DEKR") / len(seq)


# --------------------------------------------------------------------------
# Standard tier -- Biopython
# --------------------------------------------------------------------------


def biopython_properties(seq: str, ph: float = DEFAULT_WORKING_PH) -> dict | None:
    """pI and charge at `ph` via ProtParam, or None when Biopython is absent.

    Returning None rather than a fallback estimate is deliberate: a worse number
    silently standing in for a better one is how a report becomes misleading.
    """
    if not HAS_BIOPYTHON or not seq:
        return None
    analysis = ProteinAnalysis(seq)
    return {
        "isoelectric_point": round(analysis.isoelectric_point(), 3),
        f"charge_at_ph{ph:g}": round(analysis.charge_at_pH(ph), 3),
        "gravy_protparam": round(analysis.gravy(), 4),
    }


def measure(seq: str, ph: float = DEFAULT_WORKING_PH) -> tuple[dict, list[str]]:
    """Every raw number, plus notes about anything unavailable."""
    metrics = {
        "length": len(seq),
        "net_charge_simple": net_charge_simple(seq),
        "charge_density": round(normalized_charge(seq), 4),
        "abs_charge_density": round(abs(normalized_charge(seq)), 4),
        "charged_fraction": round(charged_fraction(seq), 4),
        "n_pos_KR": sum(1 for aa in seq if aa in POSITIVE),
        "n_neg_DE": sum(1 for aa in seq if aa in NEGATIVE),
        "n_his": seq.count("H"),
    }
    notes: list[str] = []
    extra = biopython_properties(seq, ph)
    if extra is None:
        metrics.update({"isoelectric_point": None,
                        f"charge_at_ph{ph:g}": None,
                        "gravy_protparam": None})
        notes.append(
            f"pI and charge at pH {ph:g} unavailable -- Biopython not installed "
            f"({BIOPYTHON_HINT}). The simple net-charge metrics above are unaffected."
        )
    else:
        metrics.update(extra)
    return metrics, notes


# --------------------------------------------------------------------------
# Checks -- advisory only
# --------------------------------------------------------------------------


def check_charge_density(seq: str, limit: float) -> list[Issue]:
    density = normalized_charge(seq)
    if abs(density) <= limit:
        return []
    direction = "acidic" if density < 0 else "basic"
    return [
        Issue(
            check="charge_density_outlier",
            label=f"strongly {direction}: net charge {net_charge_simple(seq):+d} "
                  f"over {len(seq)} aa ({density:+.3f}/residue)",
            reason=f"Beyond the provisional +/-{limit:.2f} per residue. Well-behaved "
                   f"sequences in this repo sit within -0.024..+0.062. A strong net "
                   f"charge usually comes from a low-complexity charged tract and "
                   f"often accompanies disorder; it also complicates ion-exchange "
                   f"purification.",
            severity="WARN",
        )
    ]


def check_pi_near_working_ph(
    metrics: dict, ph: float = DEFAULT_WORKING_PH, window: float = DEFAULT_PI_WINDOW
) -> list[Issue]:
    """Flag pI close to the requested working pH.

    Unlike the charge-density check this is mechanistic rather than a population
    comparison: close to pI, net charge is near zero and colloidal stability can
    be worse in that buffer. It remains advisory and context-specific; natural
    soluble proteins can trigger it. Skipped entirely when Biopython is absent,
    rather than guessed at.
    """
    pi = metrics.get("isoelectric_point")
    if pi is None:
        return []
    if abs(pi - ph) > window:
        return []
    return [
        Issue(
            check="pi_near_working_ph",
            label=f"pI {pi:.2f} is within {window:g} of pH {ph:g}",
            reason=f"Near the isoelectric point, net charge is close to zero and "
                   f"solubility or colloidal stability can be worse in this buffer. "
                   f"This is a condition-specific caution, not a design failure by "
                   f"itself; many natural soluble proteins sit near neutral pI. "
                   f"Consider checking an alternate buffer pH before redesign.",
            severity="WARN",
        )
    ]


def analyse(
    record: Record,
    max_abs_charge_density: float = DEFAULT_MAX_ABS_CHARGE_DENSITY,
    ph: float = DEFAULT_WORKING_PH,
    pi_window: float = DEFAULT_PI_WINDOW,
) -> ChargeResult:
    seq = record.sequence
    metrics, notes = measure(seq, ph)
    return ChargeResult(
        seq_id=record.seq_id,
        length=len(seq),
        metrics=metrics,
        issues=[
            *check_charge_density(seq, max_abs_charge_density),
            *check_pi_near_working_ph(metrics, ph, pi_window),
        ],
        notes=notes,
    )


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def _short(seq_id: str, width: int = 34) -> str:
    head = seq_id.split("|")[0]
    return head if len(head) <= width else head[: width - 1] + "…"


def print_report(
    results: list[ChargeResult],
    max_abs_charge_density: float = DEFAULT_MAX_ABS_CHARGE_DENSITY,
    ph: float = DEFAULT_WORKING_PH,
) -> None:
    print(f"\nStep 1 subcheck -- charge and pI   ({len(results)} sequence(s))")
    if not HAS_BIOPYTHON:
        print(f"  note: Biopython not installed -- pI and charge at pH {ph:g} are "
              f"unavailable ({BIOPYTHON_HINT}).")
        print("        the simple net-charge metrics below are unaffected.")

    for result in results:
        m = result.metrics
        n = len(result.issues)
        summary = "within reference range" if not n else f"{n} flag(s)"
        print(f"\n  {_short(result.seq_id)}   {result.length} aa   "
              f"{result.status}   {summary}")

        def fmt(key, spec="{:+.3f}"):
            value = m.get(key)
            return "n/a" if value is None else spec.format(value)

        print(f"\n      {'measure':<26}{'value':>10}   threshold")
        rows = [
            ("net charge (K+R-D-E)", f"{m['net_charge_simple']:+d}", ""),
            ("charge density", fmt("charge_density"),
             f"|.| <= {max_abs_charge_density:.2f}"),
            ("charged fraction", f"{m['charged_fraction']:.3f}", ""),
            ("K+R / D+E", f"{m['n_pos_KR']}/{m['n_neg_DE']}", ""),
            ("isoelectric point", fmt("isoelectric_point", "{:.2f}"),
             f"not within {DEFAULT_PI_WINDOW:g} of pH {ph:g}"),
            (f"charge at pH {ph:g}", fmt(f"charge_at_ph{ph:g}"), ""),
        ]
        for name, value, bound in rows:
            print(f"      {name:<26}{value:>10}   {bound}")

        for note in result.notes:
            for line in textwrap.wrap(note, width=74,
                                      initial_indent="      ~ ",
                                      subsequent_indent="        "):
                print(line)

        for issue in result.issues:
            print(f"\n    [{issue.severity}] {issue.label}")
            for line in textwrap.wrap(issue.reason, width=76,
                                      initial_indent="           ",
                                      subsequent_indent="           "):
                print(line)

    n_warn = sum(1 for r in results if r.issues)
    print(f"\n  {'-' * 44}")
    print(f"  {len(results) - n_warn} clean   {n_warn} flagged")
    print("  Advisory only: nothing here fails a design.\n")


def write_csv(results: list[ChargeResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not results:
        return
    columns = list(results[0].metrics)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sequence_id", "status", "n_flags", *columns, "flags"])
        for r in results:
            writer.writerow([
                r.seq_id, r.status, len(r.issues),
                *["" if r.metrics[c] is None else r.metrics[c] for c in columns],
                "; ".join(i.label for i in r.issues),
            ])


def write_json(results: list[ChargeResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "simple_formula": "net charge = (N_K + N_R) - (N_D + N_E); density = that / length",
        "simple_caveat": "ignores His, Cys, Tyr and the termini, and assumes full "
                         "ionisation; sufficient for population outliers",
        "biopython_available": HAS_BIOPYTHON,
        "standard_source": "Bio.SeqUtils.ProtParam (isoelectric_point, charge_at_pH)",
        "thresholds_caveat": "charge-density bound is provisional, calibrated on a "
                             "handful of sequences. See ../README.md PENDING P1/P2.",
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
                "notes": r.notes,
                "flags": [
                    {"check": i.check, "severity": i.severity,
                     "label": i.label, "reason": i.reason}
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
    parser.add_argument("--max-charge-density", type=float,
                        default=DEFAULT_MAX_ABS_CHARGE_DENSITY,
                        help=f"flag |net charge / length| above this "
                             f"(default {DEFAULT_MAX_ABS_CHARGE_DENSITY})")
    parser.add_argument("--ph", type=float, default=DEFAULT_WORKING_PH,
                        help=f"working pH for charge and the pI check "
                             f"(default {DEFAULT_WORKING_PH})")
    parser.add_argument("--pi-window", type=float, default=DEFAULT_PI_WINDOW,
                        help=f"flag a pI within this of the working pH "
                             f"(default {DEFAULT_PI_WINDOW})")
    parser.add_argument("--require-biopython", action="store_true",
                        help="exit 2 if Biopython is unavailable, instead of "
                             "reporting only the simple metrics")
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 if any sequence is flagged")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    if not args.fasta and not args.seq:
        parser.error("give at least one FASTA file, or --seq")
    if args.require_biopython and not HAS_BIOPYTHON:
        print(f"error: --require-biopython given but Biopython is not installed "
              f"({BIOPYTHON_HINT})", file=sys.stderr)
        return 2

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
        analyse(r, args.max_charge_density, args.ph, args.pi_window) for r in records
    ]

    if not args.quiet:
        print_report(results, args.max_charge_density, args.ph)
    if args.csv:
        write_csv(results, args.csv)
    if args.json:
        write_json(results, args.json)

    if args.strict:
        return 1 if any(r.issues for r in results) else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
