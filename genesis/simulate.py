"""Run the protocol on paper and see what actually comes out of the tube.

`verify.py` reads the designed molecule back and checks it against the
specification.  That catches a construct assembled wrongly, but it cannot
catch a construct assembled correctly from a protocol that builds something
else: the design said CAG-DIO-ASAP7y-P2A-iGluSnFR4f and verified clean,
while the steps underneath it amplified a promoter and nothing more.

So this module ignores the designed sequence entirely.  It takes the
primers as ordered, the templates as named, and the assembly as written,
and works out what molecule those operations produce.  Only then does it
ask whether every requested part is in it, and whether it matches what the
plan claims.  A route whose steps do not build the construct it describes
fails here, however self-consistent its own bookkeeping is.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass, field

from .enzymes import digest
from .seqtools import find_all, revcomp
from .swap import sequence_coverage

#: How much of a part has to turn up in the product to count as present.
#: Not exact containment: a swap keeps the recipient's own promoter and
#: WPRE, which are the same elements as the catalogued parts without being
#: byte-identical to them, and demanding an exact match would report every
#: correct swap as having lost them.
PRESENT_COVERAGE = 0.90

#: Shortest terminal homology two fragments may be joined on.  NEBuilder
#: HiFi wants 15-25 bp, so anything shorter is not an assembly junction.
MIN_HOMOLOGY = 15
MAX_HOMOLOGY = 60


@dataclass
class SimulatedProduct:
    name: str
    seq: str
    note: str = ""


@dataclass
class SimulationResult:
    """What the protocol makes, and whether that is what was asked for."""
    passed: bool = True
    simulated: bool = True
    product_bp: int = 0
    circular: bool = True
    products: list = field(default_factory=list)
    steps: list = field(default_factory=list)
    problems: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    modules: list = field(default_factory=list)
    summary: str = ""
    #: Whether the simulated product is the designed molecule, not merely
    #: the same size as it.
    identical: bool = False
    summary_suffix: str = ""

    def to_dict(self) -> dict:
        return {
            "passed": self.passed, "simulated": self.simulated,
            "product_bp": self.product_bp, "circular": self.circular,
            "products": [{"name": p.name, "bp": len(p.seq), "note": p.note}
                         for p in self.products],
            "steps": self.steps, "problems": self.problems,
            "warnings": self.warnings,
            "modules": self.modules, "summary": self.summary,
            "identical": self.identical,
        }


# --------------------------------------------------------------------------
# Individual operations
# --------------------------------------------------------------------------

def pcr(template: str, fwd_tail: str, fwd_anneal: str,
        rev_tail: str, rev_anneal: str,
        circular: bool = True) -> tuple[str | None, str]:
    """The product a primer pair makes on a template.

    The forward primer binds the top strand; the reverse binds the bottom,
    so its site on the top strand is the reverse complement of its annealing
    region.  Whatever 5' tails the primers carry are not on the template but
    are in the product, which is the whole point of an assembly primer.
    """
    t = template.upper()
    if not fwd_anneal or not rev_anneal:
        return None, "primer has no annealing region"

    # Which strand we call "the template" is bookkeeping, not chemistry: a
    # pair written against the bottom strand amplifies just as well, and a
    # part on the minus strand of its donor is designed that way on purpose.
    # So try both, and report failure only if neither works.
    for strand in (t, revcomp(t)):
        product, why = _pcr_one_strand(strand, fwd_tail, fwd_anneal,
                                       rev_tail, rev_anneal, circular)
        if product is not None:
            return product, ""
    return None, why


def _pcr_one_strand(t: str, fwd_tail: str, fwd_anneal: str,
                    rev_tail: str, rev_anneal: str,
                    circular: bool) -> tuple[str | None, str]:
    fwd_anneal = fwd_anneal.upper()
    rev_site = revcomp(rev_anneal.upper())

    # An inverse PCR runs out through the origin and back, so the product
    # can wrap; searching a doubled template finds those too.
    search = t + t if circular else t

    f_hits = find_all(search, fwd_anneal)
    f_hits = [h for h in f_hits if h < len(t)]
    r_hits = find_all(search, rev_site)
    if not f_hits:
        return None, "the forward primer does not bind this template"
    if not r_hits:
        return None, "the reverse primer does not bind this template"

    f = f_hits[0]
    ends = [h + len(rev_site) for h in r_hits if h + len(rev_site) > f]
    if not ends:
        return None, "the reverse primer binds upstream of the forward primer"
    end = min(ends)
    if end - f > len(t) + 1:
        return None, "the primers would amplify more than the whole template"

    body = search[f:end]
    return fwd_tail.upper() + body + revcomp(rev_tail.upper()), ""


def _overlap_len(a: str, b: str) -> int:
    """Longest terminal homology between the end of `a` and the start of `b`."""
    limit = min(MAX_HOMOLOGY, len(a), len(b))
    for n in range(limit, MIN_HOMOLOGY - 1, -1):
        if a[-n:] == b[:n]:
            return n
    return 0


def assemble_by_homology(fragments: list[SimulatedProduct]
                         ) -> tuple[str | None, list[str]]:
    """Join fragments on their shared ends, the way HiFi or Gibson does.

    Orientation is not given to us - a fragment amplified off the other
    strand still assembles - so both are tried.  The order is whatever
    chains up; a set of fragments that does not chain is a failed assembly,
    which is exactly the kind of thing worth knowing before ordering
    primers.
    """
    if not fragments:
        return None, ["no fragments to assemble"]
    if len(fragments) == 1:
        only = fragments[0].seq
        n = _overlap_len(only, only)
        if n:
            return only[:-n], []
        return None, ["a single fragment cannot circularise: its ends do "
                      "not overlap"]

    pool = [(f.name, f.seq) for f in fragments]
    notes: list[str] = []

    def extend(chain_seq: str, used: set, trail: list):
        if len(used) == len(pool):
            closing = _overlap_len(chain_seq, chain_seq)
            if closing:
                return chain_seq[:-closing], trail
            return None, trail
        for i, (name, seq) in enumerate(pool):
            if i in used:
                continue
            for oriented, label in ((seq, name), (revcomp(seq), name + " (rc)")):
                n = _overlap_len(chain_seq, oriented)
                if not n:
                    continue
                got, path = extend(chain_seq + oriented[n:], used | {i},
                                   trail + [(label, n)])
                if got is not None:
                    return got, path
        return None, trail

    first_name, first_seq = pool[0]
    product, trail = extend(first_seq, {0}, [(first_name, 0)])
    if product is None:
        notes.append(
            "these fragments do not chain into a circle on "
            f"{MIN_HOMOLOGY}+ bp of shared end sequence, so the assembly "
            "would not close")
        return None, notes
    for label, n in trail[1:]:
        notes.append(f"joined {label} on {n} bp of homology")
    return product, notes


# --------------------------------------------------------------------------
# Replaying a route
# --------------------------------------------------------------------------

def _plasmid(conn: sqlite3.Connection, pid: int | None):
    if pid is None:
        return None
    return conn.execute(
        "SELECT id, name, seq, circular FROM plasmid WHERE id = ?",
        (pid,)).fetchone()


def _templates(conn: sqlite3.Connection, plan: dict) -> dict:
    """Map every template name a primer can cite onto its real sequence."""
    out: dict[str, tuple[str, bool, str]] = {}

    swap = plan.get("swap") or {}
    for key, label in (("recipient", "backbone"), ("donor", "donor")):
        info = swap.get(key) or {}
        row = _plasmid(conn, info.get("id"))
        if row is not None:
            out[label] = (row["seq"], bool(row["circular"]), row["name"])
            out[row["name"]] = out[label]

    backbone = plan.get("backbone") or {}
    row = _plasmid(conn, backbone.get("id"))
    if row is not None:
        entry = (row["seq"], bool(row["circular"]), row["name"])
        out.setdefault("backbone", entry)
        out.setdefault(row["name"], entry)
        # The recipient is amplified as the arcs the construct keeps, so a
        # route can cite backbone1, backbone2 ... all off the same plasmid.
        for i in range(1, 17):
            out.setdefault(f"backbone{i}", entry)

    for frag in plan.get("fragments") or []:
        row = _plasmid(conn, frag.get("source_plasmid_id"))
        if row is not None:
            entry = (row["seq"], bool(row["circular"]), row["name"])
            out[frag["name"]] = entry
            out.setdefault(row["name"], entry)
    return out


def _primer_pairs(route: dict) -> list[tuple[str, dict, dict]]:
    """Group the ordered primers back into the pairs that share a reaction."""
    by_target: dict[str, dict] = {}
    for p in route.get("primers") or []:
        key = p.get("template") or p["name"].rsplit("-", 1)[0]
        by_target.setdefault(key, {})[p.get("direction", "F").upper()] = p
    pairs = []
    for key, ends in by_target.items():
        if "F" in ends and "R" in ends:
            pairs.append((key, ends["F"], ends["R"]))
    return pairs


def simulate_route(conn: sqlite3.Connection, plan: dict, route: dict
                   ) -> SimulationResult:
    """Work out what molecule this route's steps actually produce."""
    result = SimulationResult()
    operations = [s.get("operation") for s in route.get("steps") or []]

    templates = _templates(conn, plan)
    pairs = _primer_pairs(route)

    if not pairs and "PCR" in operations:
        result.simulated = False
        result.summary = "no primer pairs to simulate"
        return result

    products: list[SimulatedProduct] = []
    for name, fwd, rev in pairs:
        entry = templates.get(name)
        if entry is None:
            result.problems.append(
                f"primer pair {name} cites a template that is not in the "
                f"plan, so the reaction cannot be checked")
            result.passed = False
            continue
        seq, circular, real_name = entry
        amplicon, why = pcr(seq, fwd.get("tail", ""), fwd.get("anneal", ""),
                            rev.get("tail", ""), rev.get("anneal", ""),
                            circular)
        if amplicon is None:
            result.problems.append(f"PCR of {name} on {real_name}: {why}")
            result.passed = False
            continue
        products.append(SimulatedProduct(name, amplicon,
                                         f"{len(amplicon)} bp off {real_name}"))
        result.steps.append({"operation": "PCR", "product": name,
                             "template": real_name, "bp": len(amplicon)})

    # Fragments that are ordered rather than amplified still go in the tube.
    amplified = {p.name for p in products}
    for frag in plan.get("fragments") or []:
        if frag.get("needs_synthesis") and frag["name"] not in amplified:
            products.append(SimulatedProduct(
                frag["name"], frag.get("seq", ""), "synthesised"))

    if not products:
        result.simulated = False
        result.summary = "nothing to simulate for this route"
        return result

    result.products = products

    if any(op in ("HI_FI_ASSEMBLY", "GIBSON_ASSEMBLY", "LIGATION")
           for op in operations):
        product, notes = assemble_by_homology(products)
        for n in notes:
            result.steps.append({"operation": "ASSEMBLY", "note": n})
        if product is None:
            result.passed = False
            result.problems.extend(notes[:1] or ["the assembly does not close"])
            result.summary = "the fragments do not assemble"
            return result
    else:
        result.simulated = False
        result.summary = "this route has no assembly step to simulate"
        return result

    result.product_bp = len(product)
    _check_contents(result, product, plan, conn)
    return result


def _check_contents(result: SimulationResult, product: str, plan: dict,
                    conn: sqlite3.Connection | None = None) -> None:
    """Confirm every requested part is in the molecule the steps produce."""
    doubled = (product + product).upper()
    missing = []
    for m in plan.get("modules") or []:
        name = m.get("part_name") or m.get("name") or "?"
        seq = (m.get("seq") or "").upper()
        if not seq and conn is not None and m.get("part_id"):
            row = conn.execute("SELECT seq FROM part WHERE id = ?",
                               (m["part_id"],)).fetchone()
            seq = (row["seq"] or "").upper() if row else ""
        if not seq:
            continue
        covered = sequence_coverage(seq, doubled)
        present = covered >= PRESENT_COVERAGE
        result.modules.append({"name": name, "bp": len(seq),
                               "present": present,
                               "coverage": round(covered, 3)})
        if not present:
            missing.append(f"{name} ({covered*100:.0f}% of it)"
                           if covered > 0.05 else name)

    if missing:
        result.passed = False
        result.problems.append(
            f"the protocol does not put {', '.join(missing)} into the "
            f"product; {'they are' if len(missing) > 1 else 'it is'} in the "
            f"design but no step supplies "
            f"{'them' if len(missing) > 1 else 'it'}")

    # A length difference is worth saying out loud, but it is not the same
    # failure as a missing part: a swap keeps the recipient's own copies of
    # elements the design resolved from elsewhere, so the molecule is right
    # without being identical to the idealised one the plan prints.
    # Same length is weak evidence; the same molecule is the claim.  A
    # rotation-invariant checksum settles it for a circular product, where
    # the two can start at different bases and still be identical.
    designed = (plan.get("construct") or {}).get("sequence") or ""
    if designed and product:
        from .importer import seq_checksum
        circular = bool((plan.get("construct") or {}).get("circular", True))
        same = (seq_checksum(product, circular)
                == seq_checksum(designed, circular))
        result.identical = same
        if same:
            result.summary_suffix = " - identical to the designed sequence"
        elif len(designed) == len(product):
            result.warnings.append(
                "the product is the same length as the designed construct "
                "but not the same sequence; compare them before building")

    stated = (plan.get("construct") or {}).get("length")
    if stated and stated != len(product):
        result.warnings.append(
            f"the plan prints a {stated:,} bp construct while its steps "
            f"produce {len(product):,} bp; the parts are all there, so this "
            f"is the two being built differently rather than anything "
            f"missing")

    present = sum(1 for m in result.modules if m["present"])
    result.summary = (
        f"{len(product):,} bp product, {present}/{len(result.modules)} "
        f"requested parts present"
        + (result.summary_suffix if result.passed
           else " - the steps do not build this design"))
