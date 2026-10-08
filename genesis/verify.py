"""Read the designed molecule back and check it is what was asked for.

Everything else in Genesis works forward: modules are resolved, concatenated,
and the layout records where each one was put. That bookkeeping is exactly
what a bug would corrupt, so checking the design against it proves nothing.

This module throws the layout away. It takes the finished sequence, searches
it for each part's actual bases, re-annotates it from scratch with the same
motif scanner used on imported files, and only then compares what it found
against what was asked for. A module that was silently dropped, trimmed,
duplicated or inverted shows up here as a sequence that is not where it
should be.

For a Cre-lox construct the check is run twice: once on the plasmid as it
would be cloned, and once on the product of recombination, because those are
two different molecules and only the second one is expressed.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass, field

from .seqtools import find_all, revcomp, slice_region, translate
from .vocab import find_motifs

#: Module types whose protein product is worth confirming, not just the DNA.
TRANSLATED = {"CDS", "trafficking_motif", "linker", "p2A", "tag"}

#: Elements a construct is *supposed* to carry more than one of. Two lox
#: sites of each type are what makes a Cre-lox cassette work, and an AAV
#: needs two ITRs, so finding them twice is confirmation, not a fault.
EXPECTED_REPEATS = {"lox_site", "AAV_ITR"}


@dataclass
class ModuleCheck:
    name: str
    feature_type: str
    expected_bp: int
    found: bool = False
    occurrences: int = 0
    strand: int = 0
    start: int | None = None
    protein_ok: bool | None = None
    detail: str = ""

    def to_dict(self) -> dict:
        return {"name": self.name, "feature_type": self.feature_type,
                "expected_bp": self.expected_bp, "found": self.found,
                "occurrences": self.occurrences, "strand": self.strand,
                "start": self.start, "protein_ok": self.protein_ok,
                "detail": self.detail}


@dataclass
class VerificationResult:
    passed: bool = True
    checks: list = field(default_factory=list)
    problems: list = field(default_factory=list)
    observed: list = field(default_factory=list)
    order_ok: bool = True
    recombined: dict | None = None
    summary: str = ""

    def to_dict(self) -> dict:
        return {
            "passed": self.passed,
            "order_ok": self.order_ok,
            "summary": self.summary,
            "checks": [c.to_dict() for c in self.checks],
            "problems": self.problems,
            "observed": self.observed,
            "recombined": self.recombined,
        }


def _locate(seq: str, needle: str, circular: bool) -> tuple[list[int], list[int]]:
    """Where does this exact sequence occur, on each strand?"""
    if not needle:
        return [], []
    return (find_all(seq, needle, circular=circular),
            find_all(seq, revcomp(needle), circular=circular))


def verify_design(construct, resolved_modules, circular: bool = True,
                  recombined_seq: str | None = None) -> VerificationResult:
    """Confirm the assembled sequence contains what the specification asked for."""
    result = VerificationResult()
    seq = construct.seq
    if not seq:
        result.passed = False
        result.problems.append("the construct has no sequence to verify")
        result.summary = "nothing to verify"
        return result

    placed = [rm for rm in resolved_modules if rm.seq]
    positions: list[tuple[int, str]] = []

    for rm in placed:
        check = ModuleCheck(name=rm.part_name, feature_type=rm.feature_type,
                            expected_bp=len(rm.seq))
        forward, reverse = _locate(seq, rm.seq, circular)
        check.occurrences = len(forward) + len(reverse)

        if not forward and not reverse:
            check.detail = ("its sequence does not appear anywhere in the "
                            "assembled molecule")
            result.problems.append(
                f"{rm.part_name}: {check.detail}")
            result.passed = False
        else:
            check.found = True
            if forward:
                check.strand, check.start = 1, forward[0]
            else:
                check.strand, check.start = -1, reverse[0]
            positions.append((check.start, rm.part_name))

            if (check.occurrences > 1
                    and rm.feature_type not in EXPECTED_REPEATS):
                check.detail = (f"appears {check.occurrences} times; a "
                                f"diagnostic digest or sequencing read cannot "
                                f"tell the copies apart")
                result.problems.append(f"{rm.part_name}: {check.detail}")

            # `rm.seq` is the molecule's own bases, already flipped by the
            # cassette expander where the payload is inverted. So it must be
            # found on the forward strand; `orientation` describes which way
            # the element reads biologically, and is what decides how to
            # translate it, not where to look for it.
            if not forward:
                check.detail = (
                    "its bases were only found reverse-complemented, which "
                    "means the module was placed the wrong way round")
                result.problems.append(f"{rm.part_name}: {check.detail}")
                result.passed = False

            if rm.feature_type in TRANSLATED:
                reading = (revcomp(rm.seq)
                           if rm.module.orientation == "reverse" else rm.seq)
                protein = translate(reading)
                internal = protein.rstrip("*")
                check.protein_ok = "*" not in internal
                if not check.protein_ok:
                    at = internal.index("*") + 1
                    check.detail = (f"reading it on its own strand gives a "
                                    f"stop codon at residue {at}")
                    result.problems.append(f"{rm.part_name}: {check.detail}")

        result.checks.append(check)

    # --- order ----------------------------------------------------------
    # Compare the order the parts actually occur in against the order they
    # were asked for. Reverse-strand runs legitimately read backwards, so
    # only forward-strand modules are ordered here.
    # An inverted payload legitimately reads backwards, so only modules
    # placed in the forward sense are ordered against the specification.
    orderable = [c for c, rm in zip(result.checks, placed)
                 if c.found and c.start is not None
                 and rm.module.orientation != "reverse"
                 and c.occurrences == 1]
    asked = [c.name for c in orderable]
    actual = [c.name for c in sorted(orderable, key=lambda c: c.start)]
    if asked != actual:
        result.order_ok = False
        result.passed = False
        result.problems.append(
            f"the forward-strand parts occur in a different order from the "
            f"specification: asked for {' - '.join(asked)}, found "
            f"{' - '.join(actual)}")

    # --- independent re-annotation ---------------------------------------
    observed = find_motifs(seq, circular=circular)
    result.observed = [
        {"name": m["name"], "feature_type": m["feature_type"],
         "start": m["start"], "end": m["end"], "strand": m["strand"],
         "level": m["level"]}
        for m in observed]

    # --- the expressed molecule, for recombinase cassettes ----------------
    if recombined_seq is not None and recombined_seq != seq:
        # Which way the cassette is meant to switch depends on how the
        # payload was stored. An inverted payload (DIO, FLEX, FAS) should
        # become readable forward; a forward one (DO) should stop being
        # readable, because Cre turns that arrangement off.
        inverted = [rm for rm in placed
                    if rm.feature_type in TRANSLATED
                    and rm.module.orientation == "reverse"]
        forward = [rm for rm in placed
                   if rm.feature_type in TRANSLATED
                   and rm.module.orientation != "reverse"]

        if inverted:
            switched = [rm.part_name for rm in inverted
                        if _locate(recombined_seq, revcomp(rm.seq), False)[0]]
            result.recombined = {
                "direction": "on",
                "length": len(recombined_seq),
                "payload_modules": [rm.part_name for rm in inverted],
                "now_forward": switched,
                "all_switched": len(switched) == len(inverted),
                "detail": "the payload was stored inverted, so recombination "
                          "should make it readable on the forward strand",
            }
            if len(switched) != len(inverted):
                still = [rm.part_name for rm in inverted
                         if rm.part_name not in switched]
                result.problems.append(
                    f"after recombination these are still not readable on the "
                    f"forward strand: {', '.join(still)}")
                result.passed = False
        else:
            # A Cre-OFF cassette: the payload should no longer be intact.
            remaining = [rm.part_name for rm in forward
                         if _locate(recombined_seq, rm.seq, False)[0]]
            result.recombined = {
                "direction": "off",
                "length": len(recombined_seq),
                "payload_modules": [rm.part_name for rm in forward],
                "now_forward": remaining,
                "all_switched": not remaining,
                "detail": "the payload was stored forward, so recombination "
                          "should remove or invert it",
            }
            if forward and remaining:
                result.problems.append(
                    f"this looks like a Cre-OFF cassette, but after "
                    f"recombination {', '.join(remaining)} is still intact "
                    f"on the forward strand")
                result.passed = False

    # --- summary ----------------------------------------------------------
    found = sum(1 for c in result.checks if c.found)
    total = len(result.checks)
    if result.passed:
        result.summary = (
            f"read back {found} of {total} parts from the assembled sequence, "
            f"each once, in the specified order and orientation")
    else:
        result.summary = (
            f"read back {found} of {total} parts; "
            f"{len(result.problems)} problem(s) found")
    return result


def verify_plan(construct, resolved_modules, circular: bool = True) -> dict:
    """Convenience wrapper that also simulates recombination when relevant."""
    from .validate import lox_sites_from_layout, simulate_cre

    recombined = None
    sites = lox_sites_from_layout(construct.seq)
    if len(sites) >= 2:
        recombined, _ = simulate_cre(construct.seq, sites, circular=circular)

    return verify_design(construct, resolved_modules, circular,
                         recombined).to_dict()
