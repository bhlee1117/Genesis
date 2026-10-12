"""Stage 5: validate the intended construct in silico, before any planning.

Three layers of checks, as the specification requires:

* DNA level - order, orientation, duplications, ITR integrity, site uniqueness
* Protein level - reading frame, stops, linker and 2A placement
* Recombination level - simulate Cre so that DIO/FLEX designs are checked in
  their *expressed* state, not only as deposited DNA
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

from .construct import AssembledConstruct, translated_span
from .enzymes import CATALOG, find_sites
from .seqtools import (STOP_CODONS, longest_repeat, revcomp, slice_region,
                       translate)
from .vocab import PROTECTED_TYPES

#: Payload capacity between AAV ITRs. Past this, packaging efficiency falls
#: off sharply and truncated genomes appear.
AAV_PAYLOAD_LIMIT = 4700
AAV_ITR_LENGTH = 145


@dataclass
class Issue:
    level: str        # error | warning | info
    category: str     # dna | protein | recombination | packaging
    message: str
    location: tuple | None = None

    def to_dict(self) -> dict:
        return asdict(self)


# --------------------------------------------------------------------------
# Cre / lox recombination
# --------------------------------------------------------------------------

LOX_SPACER_OFFSET = 13      # 13 bp inverted repeat, then the 8 bp spacer
LOX_LENGTH = 34


@dataclass
class LoxSite:
    name: str
    start: int
    end: int
    strand: int

    @property
    def crossover(self) -> int:
        """Midpoint of the asymmetric spacer, where strand exchange occurs."""
        return self.start + LOX_LENGTH // 2


def simulate_cre(seq: str, lox_sites: list[LoxSite], circular: bool = True,
                 max_rounds: int = 12) -> tuple[str, list[str]]:
    """Apply Cre recombination until the molecule stops changing.

    Only sites with *identical* spacers recombine with each other, which is
    exactly why DIO/FLEX cassettes use two lox variants.

    Two rules drive the outcome, and both matter:

    * A pair in **opposite** orientation inverts the DNA between them. This
      reaction is reversible - Cre keeps flipping it - so it reaches an
      equilibrium rather than a product, and a given pair is therefore
      applied only once here.
    * A pair in the **same** orientation excises the DNA between them as a
      circle, which is effectively irreversible and so drives the system.
      Excision is always preferred when one is available.

    In a FLEX cassette the first inversion flips the orientation of the
    *other* lox pair, which turns it into an excisable pair and locks the
    payload in the expressed orientation. That interplay is the whole point
    of the design, so it is simulated rather than assumed.
    """
    log: list[str] = []
    current = seq
    sites = sorted(lox_sites, key=lambda s: s.start)
    inverted_pairs: set[tuple] = set()

    def pair_key(a: LoxSite, b: LoxSite) -> tuple:
        return (a.name, min(a.crossover, b.crossover), max(a.crossover, b.crossover))

    for _ in range(max_rounds):
        by_name: dict[str, list[LoxSite]] = {}
        for s in sites:
            by_name.setdefault(s.name, []).append(s)

        excision = None
        inversion = None
        for name, group in by_name.items():
            if len(group) < 2:
                continue
            group = sorted(group, key=lambda s: s.start)
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    a, b = group[i], group[j]
                    if a.strand == b.strand:
                        if excision is None:
                            excision = (name, a, b)
                    elif inversion is None and pair_key(a, b) not in inverted_pairs:
                        inversion = (name, a, b)

        if excision:
            name, a, b = excision
            left, right = sorted((a.crossover, b.crossover))
            removed = right - left
            current = current[:left] + current[right:]
            log.append(f"Cre: excised {removed} bp between the parallel "
                       f"{name} pair (one {name} site retained)")
            new_sites = []
            for s in sites:
                if left <= s.crossover < right:
                    continue
                if s.crossover >= right:
                    new_sites.append(LoxSite(s.name, s.start - removed,
                                             s.end - removed, s.strand))
                else:
                    new_sites.append(s)
            sites = sorted(new_sites, key=lambda s: s.start)
            inverted_pairs = {p for p in inverted_pairs if p[2] < left}
            continue

        if inversion:
            name, a, b = inversion
            left, right = sorted((a.crossover, b.crossover))
            current = current[:left] + revcomp(current[left:right]) + current[right:]
            log.append(f"Cre: inverted {right - left} bp between the "
                       f"antiparallel {name} pair")
            inverted_pairs.add(pair_key(a, b))
            # Sites wholly inside the segment are mirrored and flip strand.
            # The two acting sites straddle the crossover points and are
            # reconstituted in place, so they keep their orientation.
            new_sites = []
            for s in sites:
                if s is a or s is b:
                    new_sites.append(s)
                elif left <= s.start and s.end <= right:
                    new_sites.append(LoxSite(s.name, left + (right - s.end),
                                             left + (right - s.start), -s.strand))
                else:
                    new_sites.append(s)
            sites = sorted(new_sites, key=lambda s: s.start)
            continue

        break
    return current, log


def lox_sites_from_layout(seq: str) -> list[LoxSite]:
    from .vocab import find_motifs
    out = []
    for m in find_motifs(seq, include_protein=False):
        if m["feature_type"] == "lox_site" and m["name"].lower().startswith(("lox", "frt")):
            out.append(LoxSite(m["name"], m["start"], m["end"], m["strand"]))
    return out


# --------------------------------------------------------------------------
# Main validation
# --------------------------------------------------------------------------

def validate_construct(construct: AssembledConstruct,
                       expected_order: list[str] | None = None,
                       stocked_enzymes: list[str] | None = None,
                       check_aav: bool = True,
                       wants_cassette: bool = False) -> list[Issue]:
    issues: list[Issue] = []
    seq = construct.seq
    if not seq:
        return [Issue("error", "dna", "construct has no sequence")]

    for w in construct.warnings:
        issues.append(Issue("warning", "dna", w))

    # A recombinase site in a construct that never asked for one is not a
    # tidiness problem. Cre turns up in the animal, not in the plasmid: a
    # payload left sitting between somebody else's lox sites is silent
    # until it does, and an inverted one stays inverted. Whichever recipient
    # it came from and whichever path put it there, this cannot ship.
    if not wants_cassette:
        stray = [i["name"] for i in construct.layout
                 if i["type"] == "lox_site"]
        if stray:
            issues.append(Issue(
                "error", "dna",
                f"this design does not ask for a recombinase cassette, but "
                f"the construct carries {len(stray)} lox site(s) "
                f"({', '.join(dict.fromkeys(stray))}) - it would be "
                f"Cre-dependent, which is a different experiment; start "
                f"from a recipient without a cassette"))

    # A required module that resolved to nothing means the assembled molecule
    # is not the construct that was specified. For something like a DIO
    # cassette this is the difference between a Cre-dependent construct and a
    # constitutive one, so it is an error, never a quiet omission.
    for missing in getattr(construct, "missing_required", []):
        issues.append(Issue(
            "error", "dna",
            f"required module {missing['type']}"
            f"{'/' + missing['name'] if missing['name'] else ''} could not be "
            f"resolved to any part in the inventory, so it is absent from the "
            f"assembled construct; add the part, name it exactly, or allow "
            f"synthesis for it"))

    # --- order and orientation ------------------------------------------
    # Only the modules this design places are compared. A backbone brings
    # its own promoter, ITRs and origin, and those were never part of the
    # request.
    actual = [item["type"] for item in construct.layout
              if item["type"] != "backbone" and not item.get("from_backbone")]
    if expected_order and actual != expected_order:
        issues.append(Issue("error", "dna",
                            f"module order differs from the specification: "
                            f"expected {expected_order}, assembled {actual}"))
    for item in construct.layout:
        if item["strand"] == -1 and item["type"] in ("promoter", "polyA"):
            issues.append(Issue("warning", "dna",
                                f"{item['name']} is on the reverse strand; "
                                f"check this is intentional",
                                (item["start"], item["end"])))

    # --- duplications / repeats -----------------------------------------
    repeats = longest_repeat(seq, min_len=40)
    for r in repeats[:5]:
        issues.append(Issue(
            "warning", "dna",
            f"{r['length']} bp direct repeat at {r['first']} and "
            f"{r['second']}; this can misprime and can misassemble in "
            f"Gibson/HiFi", (r["first"], r["second"])))
    if len(repeats) > 5:
        issues.append(Issue("info", "dna",
                            f"{len(repeats) - 5} further direct repeats "
                            f"of 40 bp or more were found"))

    # --- AAV-specific geometry ------------------------------------------
    itrs = [i for i in construct.layout if i["type"] == "AAV_ITR"]
    if check_aav and itrs:
        if len(itrs) < 2:
            issues.append(Issue("error", "packaging",
                                f"only {len(itrs)} AAV ITR found; a packageable "
                                f"genome needs two flanking ITRs"))
        else:
            left, right = min(i["start"] for i in itrs), max(i["end"] for i in itrs)
            payload = right - left
            if payload > AAV_PAYLOAD_LIMIT:
                issues.append(Issue(
                    "error", "packaging",
                    f"payload between ITRs is {payload} bp, over the ~"
                    f"{AAV_PAYLOAD_LIMIT} bp AAV limit; expect truncated genomes",
                    (left, right)))
            elif payload > AAV_PAYLOAD_LIMIT - 300:
                issues.append(Issue(
                    "warning", "packaging",
                    f"payload between ITRs is {payload} bp, close to the ~"
                    f"{AAV_PAYLOAD_LIMIT} bp AAV limit", (left, right)))
        for itr in itrs:
            span = itr["end"] - itr["start"]
            if span < 100:
                issues.append(Issue(
                    "warning", "dna",
                    f"ITR at {itr['start']} spans only {span} bp; full ITRs are "
                    f"~{AAV_ITR_LENGTH} bp and are easily lost during cloning",
                    (itr["start"], itr["end"])))

    # --- protein level ---------------------------------------------------
    issues.extend(_check_protein(construct))

    # --- recombination ---------------------------------------------------
    issues.extend(_check_recombination(construct))

    # --- restriction site uniqueness ------------------------------------
    if stocked_enzymes:
        unique = []
        for name in stocked_enzymes:
            enz = CATALOG.get(name)
            if not enz:
                continue
            n = len(find_sites(seq, enz, circular=construct.circular))
            if n == 1:
                unique.append(name)
        if unique:
            issues.append(Issue("info", "dna",
                                "unique sites available for diagnostic digest: "
                                + ", ".join(sorted(unique)[:12])))
        else:
            issues.append(Issue("warning", "dna",
                                "no stocked enzyme cuts this construct exactly "
                                "once; colony screening by digest will be harder"))

    return issues


def _check_protein(construct: AssembledConstruct) -> list[Issue]:
    issues: list[Issue] = []
    runs = translated_span(construct)
    if not runs:
        return issues

    for run in runs:
        members = run["members"]
        span = construct.seq[run["start"]:run["end"]]
        if run.get("strand", 1) == -1:
            # Read the fusion off the strand it is actually encoded on.
            span = revcomp(span)
        names = " + ".join(m["name"] for m in members)

        if len(span) % 3 != 0:
            issues.append(Issue(
                "error", "protein",
                f"fusion ORF ({names}) is {len(span)} bp, not a multiple of 3; "
                f"the reading frame is broken "
                f"({len(span) % 3} extra base(s))", (run["start"], run["end"])))

        prot = translate(span)
        stop_at = prot.find("*")
        if stop_at != -1 and stop_at < len(prot) - 1:
            upstream = next((m["name"] for m in members
                             if m["start"] <= run["start"] + stop_at * 3 < m["end"]),
                            "?")
            issues.append(Issue(
                "error", "protein",
                f"premature stop codon at residue {stop_at + 1} of the fusion "
                f"({names}), inside {upstream}",
                (run["start"] + stop_at * 3, run["start"] + stop_at * 3 + 3)))

        # Per-module frame bookkeeping, walking the fusion 5' to 3'.
        offset = 0
        for m in members:
            length = m["end"] - m["start"]
            if m["type"] in ("CDS", "p2A", "trafficking_motif", "tag", "linker"):
                if offset % 3 != 0:
                    issues.append(Issue(
                        "error", "protein",
                        f"{m['name']} starts out of frame (offset {offset % 3}) "
                        f"relative to the fusion start",
                        (m["start"], m["end"])))
                if m["type"] != "CDS" and length % 3 != 0:
                    issues.append(Issue(
                        "error", "protein",
                        f"{m['name']} is {length} bp, which is not a whole "
                        f"number of codons", (m["start"], m["end"])))
            offset += length

        # Tandem 2A: two self-cleaving peptides back to back produce a tiny
        # orphan peptide between them and usually mean a CDS carried a fused
        # 2A that nobody noticed. Checked on the assembled sequence, so it is
        # caught whether the duplicate came from a module or from inside one.
        from .vocab import find_motifs
        twoa = sorted(
            (m for m in find_motifs(span, include_protein=True)
             if m["feature_type"] == "p2A" and m["strand"] == 1),
            key=lambda m: m["start"])
        for a, b in zip(twoa, twoa[1:]):
            gap = b["start"] - a["end"]
            if gap < 90:
                issues.append(Issue(
                    "error", "protein",
                    f"two 2A peptides ({a['name']} then {b['name']}) sit "
                    f"{gap} bp apart in {names}; one of the coding sequences "
                    f"almost certainly already contains a fused 2A, and the "
                    f"tandem pair would release a short orphan peptide",
                    (run["start"] + a["start"], run["start"] + b["end"])))

        # 2A placement: must sit between two coding modules, not at an end.
        for i, m in enumerate(members):
            if m["type"] != "p2A":
                continue
            before = [x for x in members[:i] if x["type"] == "CDS"]
            after = [x for x in members[i + 1:] if x["type"] == "CDS"]
            if not before or not after:
                issues.append(Issue(
                    "error", "protein",
                    f"{m['name']} needs a coding sequence on both sides; it is "
                    f"currently {'first' if not before else 'last'} in the fusion",
                    (m["start"], m["end"])))
            seg = construct.seq[m["start"]:m["end"]]
            if translate(seg).endswith("*"):
                issues.append(Issue("error", "protein",
                                    f"{m['name']} ends in a stop codon",
                                    (m["start"], m["end"])))

        # Trafficking motifs are C-terminal signals in these constructs.
        for i, m in enumerate(members):
            if m["type"] == "trafficking_motif" and i + 1 < len(members):
                nxt = members[i + 1]
                if nxt["type"] == "CDS":
                    issues.append(Issue(
                        "warning", "protein",
                        f"{m['name']} is followed by {nxt['name']}; export and "
                        f"trafficking motifs normally sit at the C-terminus of "
                        f"the protein they act on", (m["start"], m["end"])))

        # Start / stop codons for the whole fusion.
        if not span.upper().startswith("ATG") and members[0]["type"] == "CDS":
            issues.append(Issue("warning", "protein",
                                f"fusion ({names}) does not begin with ATG",
                                (run["start"], run["start"] + 3)))
        last_codon = span[-3:].upper()
        if last_codon not in STOP_CODONS:
            issues.append(Issue("warning", "protein",
                                f"fusion ({names}) does not end in a stop codon",
                                (run["end"] - 3, run["end"])))
        issues.append(Issue("info", "protein",
                            f"predicted fusion protein ({names}): "
                            f"{len(prot.rstrip('*'))} aa"))
    return issues


def _check_recombination(construct: AssembledConstruct) -> list[Issue]:
    issues: list[Issue] = []
    sites = lox_sites_from_layout(construct.seq)
    if len(sites) < 2:
        return issues

    names = sorted({s.name for s in sites})
    issues.append(Issue("info", "recombination",
                        f"{len(sites)} recombinase sites found ({', '.join(names)})"))

    recombined, log = simulate_cre(construct.seq, sites,
                                   circular=construct.circular)
    for line in log:
        issues.append(Issue("info", "recombination", line))

    if recombined == construct.seq:
        issues.append(Issue("warning", "recombination",
                            "lox sites present but Cre would not act on them; "
                            "check that the pairs are homotypic and correctly "
                            "oriented"))
        return issues

    # What matters for a Cre switch is the ORF on the *forward* strand, the
    # one the promoter reads. An inverted payload still shows a long ORF on
    # the reverse strand before Cre, so comparing both strands would report
    # a working DIO cassette as getting worse.
    from .seqtools import find_orfs
    fwd = lambda s: max((o["aa_length"] for o in find_orfs(s, min_aa=50)
                         if o["strand"] == 1), default=0)
    pre, post = fwd(construct.seq), fwd(recombined)

    issues.append(Issue(
        "info", "recombination",
        f"longest forward-strand ORF before Cre: {pre} aa; after Cre: "
        f"{post} aa ({len(recombined)} bp product)"))

    # An ORF is only measurable up to a stop codon. If the designed fusion has
    # no terminal stop, the numbers above describe whatever the ribosome runs
    # into next, not the intended protein, so no conclusion is drawn from them.
    runs = translated_span(construct)
    unterminated = [r for r in runs
                    if (revcomp(construct.seq[r["start"]:r["end"]])
                        if r.get("strand", 1) == -1
                        else construct.seq[r["start"]:r["end"]]
                        )[-3:].upper() not in STOP_CODONS]
    if unterminated:
        issues.append(Issue(
            "warning", "recombination",
            "the Cre comparison above is unreliable because the fusion has no "
            "terminal stop codon, so translation runs past the end of the "
            "payload; add a stop codon (and a polyA) and re-check"))
        return issues

    if post > pre * 1.5 and post >= 100:
        issues.append(Issue(
            "info", "recombination",
            f"Cre switches expression on ({pre} aa -> {post} aa on the "
            f"forward strand), which is what a correctly built DIO/FLEX "
            f"cassette should do"))
    elif pre >= 100 and post <= pre:
        issues.append(Issue(
            "warning", "recombination",
            f"the forward-strand ORF is already {pre} aa before Cre, so this "
            f"payload looks expressible without Cre; check the orientation of "
            f"the cassette"))
    elif post < 100:
        issues.append(Issue(
            "error", "recombination",
            f"after Cre the longest forward-strand ORF is only {post} aa; the "
            f"recombined construct would not produce a full-length protein"))
    return issues


def summarize(issues: list[Issue]) -> dict:
    return {
        "errors": sum(1 for i in issues if i.level == "error"),
        "warnings": sum(1 for i in issues if i.level == "warning"),
        "info": sum(1 for i in issues if i.level == "info"),
        "issues": [i.to_dict() for i in issues],
    }
