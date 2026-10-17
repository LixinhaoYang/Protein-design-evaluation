#!/usr/bin/env python3
"""Step 1 subcheck -- basic validity QC.

The first and cheapest gate in the triage pipeline. Its only job is to catch
illegal residues and obvious input errors, so nothing downstream ever wastes
compute on a malformed sequence.

Every problem is reported with three things: WHAT is wrong, WHY it matters, and
WHERE it is -- exact 1-based positions, plus a printed excerpt of the sequence
with the offending residues marked.

Checks
  1. length            reported always; compared against an expected length
                       only when one is available (auto-skips otherwise)
  2. alphabet          only the 20 standard amino acids
  3. special codes     X / B / Z / U / O / * and alignment gaps, classified
  4. duplicates        exact duplicate sequences, and duplicate IDs
  5. length bounds     absurdly short or long

Usage
    conda run -n base python stage1_sequence_screening/basic_validity.py designs.fasta
    conda run -n base python stage1_sequence_screening/basic_validity.py --seq MKVLAA...
    conda run -n base python stage1_sequence_screening/basic_validity.py a.fasta b.fasta --csv out.csv --json out.json
    conda run -n base python stage1_sequence_screening/basic_validity.py designs.fasta --expected-length 210
    conda run -n base python stage1_sequence_screening/basic_validity.py designs.fasta --context 30
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path

# The 20 standard amino acids. Anything else is either a special code (below)
# or garbage; both are input errors as far as this step is concerned.
STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")

# Non-standard characters that carry a specific meaning. Naming them matters:
# each implies a different upstream mistake, and "invalid character" alone does
# not tell you which. (label, why it matters, severity)
#
# Severity distinguishes two genuinely different situations:
#   FAIL -- the sequence is broken. A stop codon, an alignment gap or an
#           undetermined position means the input itself is wrong.
#   WARN -- the sequence is real but awkward. Sec and Pyl are genuine residues
#           found in genuine proteins; the problem is only that downstream tools
#           mostly cannot handle them. That is a compatibility caveat, not a
#           malformed input, so it is reported and allowed through.
SPECIAL_CODES = {
    "X": ("unknown residue",
          "Position is undetermined. The sequence was never fully resolved.", "FAIL"),
    "B": ("Asx ambiguity (Asn or Asp)",
          "Unresolved ambiguity code -- the residue was never disambiguated.", "FAIL"),
    "Z": ("Glx ambiguity (Gln or Glu)",
          "Unresolved ambiguity code -- the residue was never disambiguated.", "FAIL"),
    "U": ("selenocysteine (Sec, 21st aa)",
          "A real, genetically encoded residue, but unsupported by most design "
          "and folding tools. Check that every downstream step accepts it.", "WARN"),
    "O": ("pyrrolysine (Pyl, 22nd aa)",
          "A real, genetically encoded residue, but unsupported by most design "
          "and folding tools. Check that every downstream step accepts it.", "WARN"),
    "*": ("stop codon",
          "Translated through a stop -- indicates a frameshift or truncation.", "FAIL"),
    "-": ("alignment gap",
          "Sequence was taken from an MSA row and never cleaned.", "FAIL"),
    ".": ("alignment gap",
          "Sequence was taken from an MSA row and never cleaned.", "FAIL"),
}

# Absolute sanity bounds. NOT design constraints: they exist only to catch
# obviously broken input (a 3-residue "protein", a concatenation accident).
# Real length limits for a campaign belong in --expected-length.
DEFAULT_MIN_LENGTH = 20
DEFAULT_MAX_LENGTH = 3000

# How much sequence to show either side of a flagged position.
DEFAULT_CONTEXT = 12

# Our FASTA headers carry the intended length, e.g.
#   >1_5OJ8|candidate_state=s1|novel_unique_success|length=210
# When present it is used as the expected length; when absent the comparison is
# simply skipped, which is what makes that check self-switching.
_HEADER_LENGTH = re.compile(r"\blength=(\d+)\b", re.IGNORECASE)
_WHITESPACE = re.compile(r"\s+")


@dataclass
class Record:
    seq_id: str
    sequence: str
    header: str = ""
    source: str = ""
    expected_length: int | None = None


@dataclass
class Issue:
    """One problem found in one sequence.

    `positions` holds 1-based indices into the sequence and is empty for
    problems that are not positional (length, duplication). Keeping positions
    structured rather than baked into the message is what lets the report draw
    the sequence excerpt and the CSV/JSON carry machine-readable locations.
    """

    check: str
    label: str
    reason: str
    positions: list[int] = field(default_factory=list)
    severity: str = "FAIL"

    def where(self, limit: int = 12) -> str:
        """Human-readable position list, truncated when there are many."""
        if not self.positions:
            return ""
        shown = ", ".join(str(p) for p in self.positions[:limit])
        if len(self.positions) > limit:
            shown += f", +{len(self.positions) - limit} more"
        word = "position" if len(self.positions) == 1 else "positions"
        return f"{word} {shown}"

    def as_text(self) -> str:
        where = self.where()
        return f"{self.label} ({where})" if where else self.label


@dataclass
class Result:
    seq_id: str
    length: int
    sequence: str = ""
    issues: list[Issue] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not any(i.severity == "FAIL" for i in self.issues)

    @property
    def status(self) -> str:
        if not self.issues:
            return "PASS"
        return "FAIL" if not self.ok else "WARN"


# --------------------------------------------------------------------------
# FASTA input
# --------------------------------------------------------------------------


def read_fasta(path: Path) -> list[Record]:
    """Parse FASTA. Whitespace inside sequences is stripped; case is normalised.

    Everything else is left exactly as written -- silently repairing bad
    characters here would defeat the point of the step.
    """
    records: list[Record] = []
    header: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        if header is None:
            return
        seq = _WHITESPACE.sub("", "".join(chunks)).upper()
        seq_id = header.split()[0] if header.split() else header
        match = _HEADER_LENGTH.search(header)
        records.append(
            Record(
                seq_id=seq_id,
                sequence=seq,
                header=header,
                source=str(path),
                expected_length=int(match.group(1)) if match else None,
            )
        )

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                flush()
                header, chunks = line[1:].strip(), []
            elif line.strip():
                if header is None:
                    raise ValueError(f"{path}: sequence data before the first '>' header")
                chunks.append(line)
    flush()

    if not records:
        raise ValueError(f"{path}: no FASTA records found")
    return records


# --------------------------------------------------------------------------
# The checks
# --------------------------------------------------------------------------


def check_alphabet(sequence: str) -> list[Issue]:
    """Checks 2 and 3: standard alphabet, with every offending position recorded.

    Special codes are reported separately from unrecognised characters because
    each one points at a different upstream mistake -- a '*' means a frameshift
    or read-through, an 'X' means an undetermined position, a '-' means someone
    pasted a row out of an alignment.
    """
    issues: list[Issue] = []

    positions_by_char: dict[str, list[int]] = {}
    for index, aa in enumerate(sequence):
        if aa not in STANDARD_AA:
            positions_by_char.setdefault(aa, []).append(index + 1)

    # Special codes get one issue each: the whole point is that each character
    # implies a different upstream mistake.
    for char in sorted(positions_by_char):
        if char in SPECIAL_CODES:
            label, reason, severity = SPECIAL_CODES[char]
            issues.append(
                Issue(
                    check="special_code",
                    label=f"{char!r} -- {label}",
                    reason=reason,
                    positions=positions_by_char[char],
                    severity=severity,
                )
            )

    # Unrecognised characters share one reason, so they collapse into a single
    # issue rather than one per stray character.
    unknown = {c: p for c, p in positions_by_char.items() if c not in SPECIAL_CODES}
    if unknown:
        chars = ", ".join(repr(c) for c in sorted(unknown))
        positions = sorted(p for plist in unknown.values() for p in plist)
        noun = "not amino acid codes" if len(unknown) > 1 else "not an amino acid code"
        issues.append(
            Issue(
                check="unrecognised_character",
                label=f"{chars} -- {noun}",
                reason="Not a standard residue or a recognised special code. Check for "
                       "pasted position numbers, punctuation, or file corruption.",
                positions=positions,
            )
        )

    # Report in sequence order so the excerpts read top to bottom.
    issues.sort(key=lambda i: i.positions[0] if i.positions else 0)
    return issues


def check_length_bounds(length: int, min_length: int, max_length: int) -> list[Issue]:
    """Check 5: absurdly short or long. Not positional."""
    if length == 0:
        return [Issue("length_bounds", "sequence is empty",
                      "No residues were parsed. Check the file format.")]
    if length < min_length:
        return [Issue("length_bounds", f"length {length} is below the floor of {min_length} aa",
                      "Too short to be a folded domain -- likely truncated or the wrong input.")]
    if length > max_length:
        return [Issue("length_bounds", f"length {length} is above the ceiling of {max_length} aa",
                      "Implausibly long -- likely a concatenation accident.")]
    return []


def check_expected_length(length: int, expected: int | None) -> tuple[list[Issue], bool]:
    """Check 1: compare against the expected length, if we have one.

    Returns (issues, ran). This check switches itself off when no expected
    length was supplied or discoverable, rather than inventing one.
    """
    if expected is None:
        return [], False
    if length != expected:
        delta = length - expected
        return [Issue(
            "expected_length",
            f"length {length} does not match the expected {expected} ({delta:+d})",
            "The sequence is not the length it was declared to be. Check for "
            "truncation, an extra tag, or a mismatched header.",
        )], True
    return [], True


def check_duplicates(records: list[Record]) -> dict[str, list[Issue]]:
    """Check 4: exact duplicate sequences and duplicate IDs, across the whole set.

    Duplication is a property of the set, so it cannot be evaluated one sequence
    at a time. The first occurrence is kept clean and every later copy is
    flagged, so a representative always survives.
    """
    issues: dict[str, list[Issue]] = {}

    first_by_seq: dict[str, str] = {}
    for record in records:
        if not record.sequence:
            continue
        original = first_by_seq.get(record.sequence)
        if original is None:
            first_by_seq[record.sequence] = record.seq_id
        else:
            issues.setdefault(record.seq_id, []).append(
                Issue("duplicate_sequence", f"exact duplicate of {original}",
                      "The identical sequence already appears in the input set.")
            )

    first_by_id: dict[str, Record] = {}
    for record in records:
        original = first_by_id.get(record.seq_id)
        if original is None:
            first_by_id[record.seq_id] = record
        else:
            issues.setdefault(record.seq_id, []).append(
                Issue("duplicate_id", f"duplicate sequence ID (also in {original.source})",
                      "Two records share an ID; downstream tables would silently collide.")
            )

    return issues


def validate(
    record: Record,
    min_length: int = DEFAULT_MIN_LENGTH,
    max_length: int = DEFAULT_MAX_LENGTH,
    expected_length: int | None = None,
) -> Result:
    """Run every per-sequence check. Duplicates are handled by the caller."""
    sequence = record.sequence
    result = Result(seq_id=record.seq_id, length=len(sequence), sequence=sequence)

    result.issues.extend(check_length_bounds(len(sequence), min_length, max_length))
    result.issues.extend(check_alphabet(sequence))

    expected = expected_length if expected_length is not None else record.expected_length
    issues, ran = check_expected_length(len(sequence), expected)
    result.issues.extend(issues)
    if not ran:
        result.skipped.append("expected_length (none supplied or found in header)")

    return result


def validate_all(
    records: list[Record],
    min_length: int = DEFAULT_MIN_LENGTH,
    max_length: int = DEFAULT_MAX_LENGTH,
    expected_length: int | None = None,
) -> list[Result]:
    """Validate a whole set, including the set-level duplicate check."""
    results = [validate(r, min_length, max_length, expected_length) for r in records]
    duplicates = check_duplicates(records)
    if duplicates:
        by_id: dict[str, Result] = {}
        for result in results:
            by_id.setdefault(result.seq_id, result)
        for seq_id, issues in duplicates.items():
            target = by_id.get(seq_id)
            if target is not None:
                target.issues.extend(issues)
    return results


# --------------------------------------------------------------------------
# Position rendering
# --------------------------------------------------------------------------


def _merge_windows(positions: list[int], length: int, flank: int) -> list[tuple[int, int, list[int]]]:
    """Group positions into overlapping context windows.

    Returns (start, end, positions_in_window) with 1-based inclusive bounds, so
    two nearby problems share one excerpt instead of printing twice.
    """
    windows: list[tuple[int, int, list[int]]] = []
    for pos in sorted(positions):
        start = max(1, pos - flank)
        end = min(length, pos + flank)
        if windows and start <= windows[-1][1] + 1:
            prev_start, prev_end, prev_positions = windows[-1]
            windows[-1] = (prev_start, max(prev_end, end), prev_positions + [pos])
        else:
            windows.append((start, end, [pos]))
    return windows


def render_positions(
    sequence: str, positions: list[int], flank: int = DEFAULT_CONTEXT,
    indent: str = "      ", max_windows: int = 4,
) -> list[str]:
    """Draw sequence excerpts with carets under each flagged position.

        29  KFSVSGEG*EGDATYGKLT  47
                    ^
    """
    if not sequence or not positions:
        return []

    lines: list[str] = []
    windows = _merge_windows(positions, len(sequence), flank)
    for start, end, marked in windows[:max_windows]:
        excerpt = sequence[start - 1 : end]
        prefix = f"{indent}{start}  "
        lines.append(f"{prefix}{excerpt}  {end}")
        carets = [" "] * len(excerpt)
        for pos in marked:
            carets[pos - start] = "^"
        lines.append(" " * len(prefix) + "".join(carets).rstrip())
    if len(windows) > max_windows:
        lines.append(f"{indent}... {len(windows) - max_windows} more location(s)")
    return lines


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def _short(seq_id: str, width: int = 34) -> str:
    """FASTA ids here are long pipe-delimited headers; show the leading field."""
    head = seq_id.split("|")[0]
    return head if len(head) <= width else head[: width - 1] + "…"


def print_report(
    results: list[Result], show_skipped: bool = True,
    flank: int = DEFAULT_CONTEXT, show_context: bool = True,
) -> None:
    print(f"\nStep 1 subcheck -- basic validity   ({len(results)} sequence(s))")

    for result in results:
        n = len(result.issues)
        summary = "no problems found" if not n else f"{n} issue(s)"
        print(f"\n  {_short(result.seq_id)}   {result.length} aa   "
              f"{result.status}   {summary}")

        for issue in result.issues:
            where = issue.where()
            head = f"    [{issue.severity}] {issue.label}"
            print(f"{head}  ·  {where}" if where else head)
            for line in textwrap.wrap(issue.reason, width=76,
                                      initial_indent="           ",
                                      subsequent_indent="           "):
                print(line)
            if show_context and issue.positions:
                for line in render_positions(result.sequence, issue.positions, flank):
                    print(line)

        if show_skipped:
            for note in result.skipped:
                print(f"    [skip]  {note}")

    n_fail = sum(1 for r in results if r.status == "FAIL")
    n_warn = sum(1 for r in results if r.status == "WARN")
    n_pass = sum(1 for r in results if r.status == "PASS")
    print(f"\n  {'-' * 44}")
    print(f"  PASS {n_pass}   WARN {n_warn}   FAIL {n_fail}"
          f"   ->  {n_pass + n_warn}/{len(results)} continue")
    if n_fail:
        print("  Failed sequences are malformed input and must not continue.")
    if n_warn:
        print("  Warned sequences are valid and continue; check tool compatibility.")
    print()


def write_csv(results: list[Result], path: Path) -> None:
    """One row per issue, so positions stay machine-readable.

    A clean sequence still gets a row (with an empty check) so that every input
    appears in the output and a missing id cannot be mistaken for a pass.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["sequence_id", "length", "status", "check", "severity",
             "label", "reason", "positions", "n_positions"]
        )
        for result in results:
            if not result.issues:
                writer.writerow([result.seq_id, result.length, result.status,
                                 "", "", "", "", "", 0])
                continue
            for issue in result.issues:
                writer.writerow([
                    result.seq_id, result.length, result.status,
                    issue.check, issue.severity, issue.label, issue.reason,
                    " ".join(str(p) for p in issue.positions), len(issue.positions),
                ])


def write_json(results: list[Result], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "summary": {
            "n_sequences": len(results),
            "n_pass": sum(1 for r in results if r.status == "PASS"),
            "n_warn": sum(1 for r in results if r.status == "WARN"),
            "n_fail": sum(1 for r in results if r.status == "FAIL"),
        },
        "sequences": [
            {
                "sequence_id": r.seq_id,
                "length": r.length,
                "status": r.status,
                "skipped": r.skipped,
                "issues": [
                    {
                        "check": i.check,
                        "severity": i.severity,
                        "label": i.label,
                        "reason": i.reason,
                        "positions": i.positions,
                    }
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
    parser.add_argument("fasta", nargs="*", type=Path, help="FASTA file(s) to check")
    parser.add_argument("--seq", help="check a single sequence given on the command line")
    parser.add_argument(
        "--expected-length", type=int,
        help="expected length for every sequence; overrides any 'length=' in the header. "
             "Omit it and the comparison switches itself off.",
    )
    parser.add_argument("--min-length", type=int, default=DEFAULT_MIN_LENGTH,
                        help=f"sanity floor (default {DEFAULT_MIN_LENGTH})")
    parser.add_argument("--max-length", type=int, default=DEFAULT_MAX_LENGTH,
                        help=f"sanity ceiling (default {DEFAULT_MAX_LENGTH})")
    parser.add_argument("--context", type=int, default=DEFAULT_CONTEXT,
                        help=f"residues shown either side of a flagged position "
                             f"(default {DEFAULT_CONTEXT}; 0 disables excerpts)")
    parser.add_argument("--csv", type=Path, help="write one row per issue to this CSV")
    parser.add_argument("--json", type=Path, help="write the full report to this JSON")
    parser.add_argument("--strict", action="store_true",
                        help="treat warnings as failures (exit 1 on any issue at all)")
    parser.add_argument("--quiet", action="store_true", help="suppress the report")
    parser.add_argument("--no-skipped", action="store_true",
                        help="hide notes about checks that switched themselves off")
    args = parser.parse_args(argv)

    if not args.fasta and not args.seq:
        parser.error("give at least one FASTA file, or --seq")
    if args.min_length > args.max_length:
        parser.error("--min-length cannot exceed --max-length")
    if args.context < 0:
        parser.error("--context cannot be negative")

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

    results = validate_all(records, args.min_length, args.max_length, args.expected_length)

    if not args.quiet:
        print_report(results, show_skipped=not args.no_skipped,
                     flank=args.context or DEFAULT_CONTEXT,
                     show_context=args.context > 0)
    if args.csv:
        write_csv(results, args.csv)
    if args.json:
        write_json(results, args.json)

    if args.strict:
        return 1 if any(r.issues for r in results) else 0
    return 1 if any(not r.ok for r in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
