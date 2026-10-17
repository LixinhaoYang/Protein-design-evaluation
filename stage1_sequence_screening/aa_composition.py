#!/usr/bin/env python3
"""Step 1 subcheck -- amino-acid composition.

Profiles what a sequence is made of, and flags compositions that sit far outside
what natural proteins look like.

This is a PROFILING step, not a gate. Composition is never "broken input" the
way a stop codon is -- it is a design-quality signal. So every flag here is a
WARN, nothing fails, and the raw fractions are the real output. Use --strict if
you want warnings to be fatal for a particular run.

Reports
  * the fraction of each of the 20 standard residues
  * grouped fractions: hydrophobic, charged, positive, negative, gly_pro, cys
  * charge balance |positive - negative|
  * each value against a natural-protein reference band

Usage
    python stage1_sequence_screening/aa_composition.py designs.fasta
    python stage1_sequence_screening/aa_composition.py --seq MKVLAA...
    python stage1_sequence_screening/aa_composition.py designs.fasta --csv comp.csv --json comp.json
    python stage1_sequence_screening/aa_composition.py designs.fasta --all-residues
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import textwrap
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

try:
    from .basic_validity import Issue, Record, STANDARD_AA, _WHITESPACE, read_fasta
except ImportError:  # pragma: no cover - direct script execution path
    from basic_validity import Issue, Record, STANDARD_AA, _WHITESPACE, read_fasta

# --------------------------------------------------------------------------
# Reference composition
# --------------------------------------------------------------------------

# Amino-acid background frequencies, as percentages.
#
# Source: UniProtKB/Swiss-Prot release statistics, release 2026_03 (02-Sep-2026),
#         575,748 sequence entries, 209,017,843 amino acids.
#         https://web.expasy.org/docs/relnotes/relstat.html
#
# VERSION-LOCKED ON PURPOSE. Swiss-Prot composition drifts slightly with every
# release as the database grows, so "the Swiss-Prot average" without a release
# number is not reproducible. UniProt publishes current and previous release
# statistics, so pinning the release removes the drift problem entirely. To
# update, replace this table AND the release string together -- never one
# without the other.
SWISSPROT_RELEASE = "2026_03 (02-Sep-2026); 575,748 entries; 209,017,843 aa"

SWISSPROT_PERCENT = {
    "A": 8.25, "R": 5.52, "N": 4.06, "D": 5.46, "C": 1.38,
    "Q": 3.93, "E": 6.71, "G": 7.07, "H": 2.27, "I": 5.90,
    "L": 9.64, "K": 5.79, "M": 2.41, "F": 3.86, "P": 4.75,
    "S": 6.66, "T": 5.36, "W": 1.10, "Y": 2.92, "V": 6.85,
}

# Residue groups. Membership is a CONVENTION, not a fact -- "hydrophobic" and
# "charged" each have several defensible definitions -- so every set is spelled
# out here, named in the output, and used to DERIVE the reference percentages
# below, so a group can never drift out of step with its own reference.
#
# NOTE: these groups OVERLAP and do not partition the sequence. `positive` and
# `negative` are subsets of `charged`; `cys` is a subset of `hydrophobic`. Their
# percentages therefore do not sum to 100%, by construction.
RESIDUE_GROUPS = {
    "hydrophobic": frozenset("AVILMFWYC"),
    "charged": frozenset("DEKR"),
    "positive": frozenset("KR"),
    "negative": frozenset("DE"),
    "gly_pro": frozenset("GP"),
    "cys": frozenset("C"),
}

# His is ~10% protonated at pH 7, so it is neither cleanly charged nor cleanly
# neutral. It is excluded from the charged/positive groups above and reported in
# the per-residue table instead of being silently assigned to one side.
PARTIALLY_CHARGED = frozenset("H")


def reference_group_fractions() -> dict[str, float]:
    """Expected group fractions, derived from SWISSPROT_PERCENT.

    Computed rather than hard-coded so the reference can never disagree with the
    group definitions above.
    """
    return {
        name: sum(SWISSPROT_PERCENT[aa] for aa in members) / 100.0
        for name, members in RESIDUE_GROUPS.items()
    }


REFERENCE_GROUPS = reference_group_fractions()

# PROVISIONAL outlier bands -- NOT pass/fail thresholds.
#
# The primary output of this step is the measured fraction and its deviation
# from the reference. These bands only draw attention to extreme values, and
# they are the weakest part of the module:
#
#   * judgement-based, anchored on the reference means and sanity-checked
#     against four real sequences in this repo -- four is nowhere near enough to
#     justify a general range;
#   * Swiss-Prot pools soluble, membrane, secreted and disordered proteins
#     across all species, so "deviates from the global average" does NOT mean
#     "unreasonable sequence". A genuine membrane protein SHOULD be an extreme
#     hydrophobic outlier.
#
# The proper replacement is a matched background -- experimentally solved,
# soluble, single-chain globular proteins in our length range -- scored by
# EMPIRICAL PERCENTILE, flagging only the far tail. Not a plain z-score: these
# fractions are bounded on [0,1] and non-Gaussian, and rare-residue fractions are
# badly skewed, so (x-mu)/sigma misfits them. A robust z, (x-median)/(1.4826*MAD),
# is the fallback when one signed number is needed. Until then, read a flag as
# "look at this", never as "reject this".
#
# Concretely: this design sits at hydrophobic 0.324 against a 0.33 floor. That
# 0.006 gap means nothing biologically, and no design should flip from flagged to
# clean by crossing it. A percentile removes the artificial boundary instead of
# re-tuning it.
#
# PARKED until the reference database is assembled. Full write-up, the decisions
# to make, and the acceptance test any replacement must keep passing:
#     ../README.md  ->  "PENDING -- blocked on the reference database"  ->  P1
# (at the very end of that document).
PROVISIONAL_BANDS = {
    "hydrophobic": (0.33, 0.52),
    "charged": (0.15, 0.35),
    "positive": (0.05, 0.20),
    "negative": (0.05, 0.20),
    "gly_pro": (0.04, 0.22),
    "cys": (0.00, 0.06),
}

# |positive - negative|. A large imbalance means a strongly charged protein,
# which is a solubility and expression consideration.
DEFAULT_MAX_CHARGE_IMBALANCE = 0.10

# A single residue above this fraction is over-represented however you measure
# it: the most common residue in Swiss-Prot is Leu at 9.7%, and the highest we
# see in a real sequence here is 13.8%.
DEFAULT_MAX_SINGLE_RESIDUE = 0.15

# Fold-change vs. the Swiss-Prot mean at which a residue is called out in the
# over/under-representation summary.
DEFAULT_ENRICHMENT_FOLD = 2.0

WHY = {
    "hydrophobic": {
        "high": "Hydrophobic-rich sequences tend to aggregate or express into "
                "inclusion bodies.",
        "low": "Too few hydrophobic residues to pack a stable core; often a sign "
               "of a disordered or unfolded design.",
    },
    "charged": {
        "high": "Very highly charged sequences are often disordered rather than "
                "folded.",
        "low": "Few charged residues leaves little to keep the protein soluble.",
    },
    "positive": {
        "high": "Strongly basic; may bind nucleic acids non-specifically and can "
                "complicate purification.",
        "low": "",
    },
    "negative": {
        "high": "Strongly acidic; often accompanies low-complexity acidic tracts "
                "and disorder.",
        "low": "",
    },
    "gly_pro": {
        "high": "Gly and Pro break regular secondary structure; an excess suggests "
                "a floppy or loop-dominated design.",
        "low": "",
    },
    "cys": {
        "high": "Cys-rich sequences risk incorrect disulfide pairing and "
                "aggregation, especially in the E. coli cytoplasm.",
        "low": "",
    },
}


@dataclass
class CompositionResult:
    seq_id: str
    length: int
    counts: dict[str, int]
    fractions: dict[str, float]
    groups: dict[str, float]
    deltas: dict[str, float]          # measured - Swiss-Prot reference
    charge_imbalance: float
    issues: list[Issue]

    @property
    def status(self) -> str:
        return "WARN" if self.issues else "PASS"


# --------------------------------------------------------------------------
# Computation
# --------------------------------------------------------------------------


def residue_counts(sequence: str) -> dict[str, int]:
    """Counts for all 20 standard residues, zeros included.

    Absent residues are kept at zero so every sequence yields the same columns;
    a ragged table would be unusable downstream.
    """
    counts = Counter(sequence)
    return {aa: counts.get(aa, 0) for aa in sorted(STANDARD_AA)}


def residue_fractions(sequence: str) -> dict[str, float]:
    counts = residue_counts(sequence)
    n = len(sequence)
    if n == 0:
        return {aa: 0.0 for aa in counts}
    return {aa: c / n for aa, c in counts.items()}


def group_fractions(sequence: str) -> dict[str, float]:
    """Fraction of the sequence in each named residue group."""
    n = len(sequence)
    if n == 0:
        return {name: 0.0 for name in RESIDUE_GROUPS}
    counts = Counter(sequence)
    return {
        name: sum(counts.get(aa, 0) for aa in members) / n
        for name, members in RESIDUE_GROUPS.items()
    }


def enrichment(fractions: dict[str, float]) -> dict[str, float]:
    """Fold-change of each residue against the Swiss-Prot mean."""
    return {
        aa: (frac / (SWISSPROT_PERCENT[aa] / 100.0)) if SWISSPROT_PERCENT.get(aa) else 0.0
        for aa, frac in fractions.items()
    }


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------


def check_groups(groups: dict[str, float], bands: dict[str, tuple[float, float]]) -> list[Issue]:
    """Flag group fractions outside their reference band. WARN only."""
    issues: list[Issue] = []
    for name, value in groups.items():
        band = bands.get(name)
        if band is None:
            continue
        low, high = band
        if value > high:
            direction, reason = "high", WHY.get(name, {}).get("high", "")
        elif value < low:
            direction, reason = "low", WHY.get(name, {}).get("low", "")
        else:
            continue
        members = "".join(sorted(RESIDUE_GROUPS[name]))
        ref = REFERENCE_GROUPS[name]
        detail = (f"{value:.1%} of residues vs. a Swiss-Prot reference of {ref:.1%} "
                  f"({value - ref:+.1%}), beyond the provisional {low:.0%}-{high:.0%} band")
        issues.append(
            Issue(
                check=f"composition_{name}",
                label=f"{name} fraction is {direction} ({members})",
                reason=f"{detail}. {reason}".strip(),
                severity="WARN",
            )
        )
    return issues


def check_charge_imbalance(imbalance: float, limit: float) -> list[Issue]:
    if imbalance <= limit:
        return []
    return [
        Issue(
            check="composition_charge_imbalance",
            label=f"charge imbalance is high ({imbalance:.1%})",
            reason=f"|positive - negative| exceeds {limit:.0%}. A strongly net-charged "
                   f"protein can be harder to express and purify, and the imbalance often "
                   f"comes from a low-complexity charged tract.",
            severity="WARN",
        )
    ]


def check_single_residues(
    fractions: dict[str, float], length: int, cap: float
) -> list[Issue]:
    """Flag any single residue that dominates the sequence."""
    issues: list[Issue] = []
    for aa, frac in sorted(fractions.items(), key=lambda kv: -kv[1]):
        if frac <= cap:
            break
        fold = frac / (SWISSPROT_PERCENT[aa] / 100.0)
        issues.append(
            Issue(
                check="composition_single_residue",
                label=f"{aa} makes up {frac:.1%} of the sequence",
                reason=f"{round(frac * length)} of {length} residues, {fold:.1f}x the "
                       f"natural average of {SWISSPROT_PERCENT[aa]:.1f}%. A single "
                       f"dominant residue usually means a low-complexity tract.",
                severity="WARN",
            )
        )
    return issues


def analyse(
    record: Record,
    bands: dict[str, tuple[float, float]] | None = None,
    max_charge_imbalance: float = DEFAULT_MAX_CHARGE_IMBALANCE,
    max_single_residue: float = DEFAULT_MAX_SINGLE_RESIDUE,
) -> CompositionResult:
    """Profile one sequence and flag outlying composition."""
    bands = bands or PROVISIONAL_BANDS
    sequence = record.sequence
    fractions = residue_fractions(sequence)
    groups = group_fractions(sequence)
    imbalance = abs(groups["positive"] - groups["negative"])

    issues = [
        *check_groups(groups, bands),
        *check_charge_imbalance(imbalance, max_charge_imbalance),
        *check_single_residues(fractions, len(sequence), max_single_residue),
    ]

    return CompositionResult(
        seq_id=record.seq_id,
        length=len(sequence),
        counts=residue_counts(sequence),
        fractions=fractions,
        groups=groups,
        deltas={g: groups[g] - REFERENCE_GROUPS[g] for g in groups},
        charge_imbalance=imbalance,
        issues=issues,
    )


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def _short(seq_id: str, width: int = 34) -> str:
    head = seq_id.split("|")[0]
    return head if len(head) <= width else head[: width - 1] + "…"


def print_report(
    results: list[CompositionResult],
    bands: dict[str, tuple[float, float]] | None = None,
    show_all_residues: bool = False,
) -> None:
    bands = bands or PROVISIONAL_BANDS
    print(f"\nStep 1 subcheck -- amino-acid composition   ({len(results)} sequence(s))")

    for result in results:
        n = len(result.issues)
        summary = "within reference ranges" if not n else f"{n} flag(s)"
        print(f"\n  {_short(result.seq_id)}   {result.length} aa   "
              f"{result.status}   {summary}")

        print(f"\n      {'group':<14}{'value':>8}{'ref':>8}{'delta':>9}   {'band':<13}")
        for name in RESIDUE_GROUPS:
            value = result.groups[name]
            ref = REFERENCE_GROUPS[name]
            low, high = bands.get(name, (0.0, 1.0))
            mark = "" if low <= value <= high else "  <-- outlier"
            print(f"      {name:<14}{value:>8.3f}{ref:>8.3f}{value - ref:>+9.3f}"
                  f"   {low:.2f}-{high:.2f}{mark}")
        limit = DEFAULT_MAX_CHARGE_IMBALANCE
        mark = "" if result.charge_imbalance <= limit else "  <-- outlier"
        print(f"      {'charge_imbal':<14}{result.charge_imbalance:>8.3f}{'':>8}{'':>9}"
              f"   <={limit:.2f}{mark}")
        print(f"\n      reference: UniProtKB/Swiss-Prot {SWISSPROT_RELEASE}")
        print( "      bands are provisional outlier hints, not pass/fail thresholds")

        folds = enrichment(result.fractions)
        notable = [
            (aa, result.fractions[aa], folds[aa])
            for aa in sorted(folds, key=lambda a: -folds[a])
            if folds[aa] >= DEFAULT_ENRICHMENT_FOLD and result.fractions[aa] > 0.02
        ]
        if notable:
            parts = ", ".join(f"{aa} {f:.1%} ({x:.1f}x)" for aa, f, x in notable[:6])
            print(f"\n      over-represented vs. natural: {parts}")

        if show_all_residues:
            print()
            items = sorted(result.fractions.items(), key=lambda kv: -kv[1])
            for i in range(0, len(items), 5):
                row = "  ".join(f"{aa} {v:>6.1%}" for aa, v in items[i : i + 5])
                print(f"      {row}")

        for issue in result.issues:
            print(f"\n    [{issue.severity}] {issue.label}")
            for line in textwrap.wrap(issue.reason, width=76,
                                      initial_indent="           ",
                                      subsequent_indent="           "):
                print(line)

    n_warn = sum(1 for r in results if r.issues)
    print(f"\n  {'-' * 44}")
    print(f"  {len(results) - n_warn} near reference   {n_warn} flagged as outlying")
    print("  Advisory only: deviation from a global average is not a defect.\n")


def write_csv(results: list[CompositionResult], path: Path) -> None:
    """One row per sequence, with every residue fraction as its own column."""
    path.parent.mkdir(parents=True, exist_ok=True)
    residues = sorted(STANDARD_AA)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["sequence_id", "length", "status", "n_flags",
             *RESIDUE_GROUPS,
             *[f"delta_{g}" for g in RESIDUE_GROUPS],
             "charge_imbalance",
             *[f"frac_{aa}" for aa in residues],
             *[f"count_{aa}" for aa in residues],
             "flags"]
        )
        for r in results:
            writer.writerow([
                r.seq_id, r.length, r.status, len(r.issues),
                *[f"{r.groups[g]:.5f}" for g in RESIDUE_GROUPS],
                *[f"{r.deltas[g]:+.5f}" for g in RESIDUE_GROUPS],
                f"{r.charge_imbalance:.5f}",
                *[f"{r.fractions[aa]:.5f}" for aa in residues],
                *[r.counts[aa] for aa in residues],
                "; ".join(i.label for i in r.issues),
            ])


def write_json(results: list[CompositionResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "reference_source": "UniProtKB/Swiss-Prot release statistics",
        "reference_release": SWISSPROT_RELEASE,
        "reference_percent": SWISSPROT_PERCENT,
        "reference_groups": REFERENCE_GROUPS,
        "groups": {k: "".join(sorted(v)) for k, v in RESIDUE_GROUPS.items()},
        "groups_overlap": "positive and negative are subsets of charged; cys is a "
                          "subset of hydrophobic; fractions do not sum to 1",
        "provisional_bands": {k: list(v) for k, v in PROVISIONAL_BANDS.items()},
        "bands_caveat": "provisional outlier hints, not pass/fail thresholds; "
                        "replace with percentiles against a matched background",
        "summary": {
            "n_sequences": len(results),
            "n_warn": sum(1 for r in results if r.issues),
        },
        "sequences": [
            {
                "sequence_id": r.seq_id,
                "length": r.length,
                "status": r.status,
                "groups": r.groups,
                "deltas_vs_reference": r.deltas,
                "charge_imbalance": r.charge_imbalance,
                "counts": r.counts,
                "fractions": r.fractions,
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
    parser.add_argument("fasta", nargs="*", type=Path, help="FASTA file(s) to profile")
    parser.add_argument("--seq", help="profile a single sequence given on the command line")
    parser.add_argument("--all-residues", action="store_true",
                        help="print all 20 residue fractions, not just the outliers")
    parser.add_argument("--max-charge-imbalance", type=float,
                        default=DEFAULT_MAX_CHARGE_IMBALANCE,
                        help=f"flag |positive-negative| above this "
                             f"(default {DEFAULT_MAX_CHARGE_IMBALANCE})")
    parser.add_argument("--max-single-residue", type=float,
                        default=DEFAULT_MAX_SINGLE_RESIDUE,
                        help=f"flag any residue above this fraction "
                             f"(default {DEFAULT_MAX_SINGLE_RESIDUE})")
    parser.add_argument("--csv", type=Path, help="write the full composition table here")
    parser.add_argument("--json", type=Path, help="write the full report here")
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 if any sequence is flagged")
    parser.add_argument("--quiet", action="store_true", help="suppress the report")
    args = parser.parse_args(argv)

    if not args.fasta and not args.seq:
        parser.error("give at least one FASTA file, or --seq")

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

    # Composition is only meaningful over the standard alphabet; anything else
    # should have been stopped by step 1.
    for record in records:
        bad = set(record.sequence) - STANDARD_AA
        if bad:
            print(f"error: {record.seq_id} contains non-standard residue(s) "
                  f"{''.join(sorted(bad))} -- run basic_validity.py first",
                  file=sys.stderr)
            return 2

    results = [
        analyse(r, None, args.max_charge_imbalance, args.max_single_residue)
        for r in records
    ]

    if not args.quiet:
        print_report(results, show_all_residues=args.all_residues)
    if args.csv:
        write_csv(results, args.csv)
    if args.json:
        write_json(results, args.json)

    if args.strict:
        return 1 if any(r.issues for r in results) else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
