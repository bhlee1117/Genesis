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
    #: Where each amplicon is cut from, and where it lands in the product.
    sources: list = field(default_factory=list)
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
            "sources": self.sources,
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


def _plasmid_features(conn: sqlite3.Connection, pid: int) -> list[dict]:
    """The annotations worth drawing on a source plasmid."""
    return [{"feature_type": r["feature_type"], "name": r["name"],
             "start": r["start"], "end": r["end"], "strand": r["strand"]}
            for r in conn.execute(
                """SELECT feature_type, name, start, end, strand
                     FROM annotation WHERE plasmid_id = ?
                      AND feature_type NOT IN ('primer_binding_site', 'unknown')
                    ORDER BY start""", (pid,))]


def _templates(conn: sqlite3.Connection, plan: dict,
               route: dict | None = None) -> dict:
    """Map every template name a primer can cite onto its real sequence.

    "backbone" means different plasmids in different routes: the recipient
    being edited in a swap, and the chosen vector everywhere else. Reading
    it the same way for both simulated one route's primers against the
    other's template.
    """
    out: dict[str, tuple] = {}
    swap_route = bool(route and str(route.get("id", "")).startswith("swap"))

    swap = plan.get("swap") or {}
    for key, label in (("recipient", "backbone"), ("donor", "donor")):
        info = swap.get(key) or {}
        row = _plasmid(conn, info.get("id"))
        if row is not None:
            entry = (row["seq"], bool(row["circular"]), row["name"], row["id"])
            if label != "backbone" or swap_route:
                out[label] = entry
            out[row["name"]] = entry

    backbone = plan.get("backbone") or {}
    row = _plasmid(conn, backbone.get("id"))
    if row is not None:
        entry = (row["seq"], bool(row["circular"]), row["name"], row["id"])
        if not swap_route:
            out["backbone"] = entry
        out.setdefault("backbone", entry)
        out.setdefault(row["name"], entry)
        # The recipient is amplified as the arcs the construct keeps, so a
        # route can cite backbone1, backbone2 ... all off the same plasmid.
        for i in range(1, 17):
            out.setdefault(f"backbone{i}", entry)

    for frag in plan.get("fragments") or []:
        row = _plasmid(conn, frag.get("source_plasmid_id"))
        if row is not None:
            entry = (row["seq"], bool(row["circular"]), row["name"], row["id"])
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

    templates = _templates(conn, plan, route)
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
        seq, circular, real_name, template_id = entry
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
        result.sources.append({
            "name": name, "template": real_name, "template_id": template_id,
            "template_length": len(seq), "template_circular": circular,
            "from_start": fwd.get("target_start"),
            "from_end": fwd.get("target_end"),
            "bp": len(amplicon),
            "features": _plasmid_features(conn, template_id),
        })

    # Fragments that are ordered rather than amplified still go in the tube.
    # Only the fragments this route actually orders. The module grouping
    # belongs to the all-PCR route; a route that orders the whole insert
    # instead must not also be handed the pieces it replaced.
    ordered_here = {str(step.get("output") or "")
                    for step in route.get("steps") or []
                    if step.get("operation") == "SYNTHESIS"}
    amplified = {p.name for p in products}
    for frag in plan.get("fragments") or []:
        if not frag.get("needs_synthesis") or frag["name"] in amplified:
            continue
        if f"{frag['name']} gBlock" not in ordered_here:
            continue
        products.append(SimulatedProduct(
            frag["name"], frag.get("seq", ""), "synthesised"))

    # A route may order a fragment of its own rather than one the module
    # grouping produced - the gene-block route orders the whole insert -
    # and that is just as much a thing in the tube.
    for step in route.get("steps") or []:
        if step.get("operation") != "SYNTHESIS":
            continue
        ordered = (step.get("details") or {}).get("sequence")
        if ordered and not any(p.seq == ordered for p in products):
            products.append(SimulatedProduct(
                step.get("output") or "ordered fragment", ordered,
                f"{len(ordered)} bp ordered"))
            result.steps.append({"operation": "SYNTHESIS",
                                 "product": step.get("output") or "insert",
                                 "bp": len(ordered)})

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

    # Where each amplicon ended up, so the assembly can be drawn rather than
    # only tabulated.
    doubled = (product + product).upper()
    for src, made in zip(result.sources, products):
        body = made.seq.upper()
        if len(body) < 60:
            anchor, offset = body, 0
        else:
            offset = len(body) // 2 - 30
            anchor = body[offset:offset + 60]
        at = doubled.find(anchor)
        if at < 0:
            at = doubled.find(revcomp(anchor))
            src["flipped"] = True
        if at >= 0:
            src["at_start"] = (at - offset) % len(product)
            src["at_end"] = (src["at_start"] + len(body)) % len(product)

    _check_contents(result, product, plan, conn, templates)
    return result


#: Feature types whose presence defines what a construct *is*, as opposed
#: to what carries it. A recipient legitimately supplies the vector parts -
#: WPRE, polyA, ITRs, origin, resistance - without any of them being named
#: in the request, so only these are leftovers worth failing over.
LEFTOVER_TYPES = {"CDS"}

#: Below this a match is not worth reporting: a 12 bp Factor Xa site turns
#: up everywhere and says nothing about whether the payload was replaced.
LEFTOVER_MIN_BP = 150


def _check_leftovers(result: "SimulationResult", product: str, plan: dict,
                     conn: sqlite3.Connection | None = None,
                     templates: dict | None = None) -> None:
    """Fail when the product still carries what it was meant to replace.

    Checking that every requested part is present says nothing about what
    else came along. A swap that removes too little leaves the old payload
    in the plasmid, and every presence check still passes: the promoter is
    there, the new coding sequence is there, and so is the indicator the
    whole point was to get rid of.

    That is not a cosmetic difference. It produced a plan whose steps
    build pAAV-CAG-<new CDS>-Voltron2-Ace2N-HaloTag and reported it as
    "2 of 2 requested parts present", with the 1,601 bp of leftover
    payload explained away as the design "being built differently".
    """
    if conn is None or not product:
        return
    # Ask of every plasmid this route actually amplifies from, not of
    # whichever one the plan calls the backbone: a swap edits the
    # recipient while the generic planner names a different vector, and
    # checking the wrong one of the two finds nothing.
    parents: dict[int, str] = {}
    for entry in (templates or {}).values():
        try:
            _seq, _circ, pname, ppid = entry
        except (TypeError, ValueError):
            continue
        if ppid:
            parents[ppid] = pname
    for key in ("backbone",):
        info = (plan.get(key) or {})
        if info.get("id"):
            parents.setdefault(info["id"], info.get("name") or "")
    swap_recipient = ((plan.get("swap") or {}).get("recipient") or {})
    if swap_recipient.get("id"):
        parents.setdefault(swap_recipient["id"],
                           swap_recipient.get("name") or "")
    if not parents:
        return
    doubled = (product + product).upper()

    # What the construct says it is made of. A recipient's own copy of a
    # part that was asked for is not a leftover.
    wanted = []
    for m in plan.get("modules") or []:
        seq = (m.get("seq") or "").upper()
        if not seq and m.get("part_id"):
            got = conn.execute("SELECT seq FROM part WHERE id = ?",
                               (m["part_id"],)).fetchone()
            seq = (got["seq"] or "").upper() if got else ""
        if seq:
            wanted.append(seq)

    marks = ",".join("?" * len(LEFTOVER_TYPES))
    left: list[str] = []
    blame = ""
    for ppid, pname in parents.items():
        prow = conn.execute("SELECT seq FROM plasmid WHERE id = ?",
                            (ppid,)).fetchone()
        if prow is None or not prow["seq"]:
            continue
        host = prow["seq"].upper()
        for a in conn.execute(
                f"""SELECT name, start, end FROM annotation
                    WHERE plasmid_id = ? AND feature_type IN ({marks})
                      AND (end - start) >= ?
                    ORDER BY (end - start) DESC""",
                (ppid, *sorted(LEFTOVER_TYPES), LEFTOVER_MIN_BP)):
            seg = host[a["start"]:a["end"]]
            if len(seg) < LEFTOVER_MIN_BP:
                continue
            if sequence_coverage(seg, doubled) < PRESENT_COVERAGE:
                continue                  # removed, as intended
            if any(sequence_coverage(seg, w) >= PRESENT_COVERAGE
                   for w in wanted):
                continue                  # it is one of the parts asked for
            name = a["name"] or "an unnamed coding sequence"
            if name not in left:
                left.append(name)
                blame = blame or pname

    if not left:
        return
    result.passed = False
    listed = ", ".join(left[:4]) + (", ..." if len(left) > 4 else "")
    result.problems.append(
        f"the product still carries {listed} from {blame or 'the recipient'}, "
        f"which this design does not ask for: the region removed was too "
        f"small, so the new coding sequence goes in beside the old payload "
        f"instead of replacing it")


def _check_contents(result: SimulationResult, product: str, plan: dict,
                    conn: sqlite3.Connection | None = None,
                    templates: dict | None = None) -> None:
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

    _check_leftovers(result, product, plan, conn, templates)

    stated = (plan.get("construct") or {}).get("length")
    if stated and stated != len(product):
        result.warnings.append(
            f"the plan prints a {stated:,} bp construct while its steps "
            f"produce {len(product):,} bp")

    present = sum(1 for m in result.modules if m["present"])
    result.summary = (
        f"{len(product):,} bp product, {present}/{len(result.modules)} "
        f"requested parts present"
        + (result.summary_suffix if result.passed
           else " - the steps do not build this design"))
