"""Stages 7, 8 and 10: generate cloning routes and rank them.

Two route families are generated as first-class alternatives:

* **Restriction route** - digest backbone and donor with compatible enzymes,
  gel purify, ligate.
* **PCR / HiFi route** - amplify each fragment with overlap tails, remove
  parental template, assemble isothermally.

Both are emitted as explicit operation-coded steps, and both are costed with
the same configurable weights so they can be compared honestly.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass, field

from .config import Config
from .construct import AssembledConstruct, ResolvedModule
from .enzymes import CATALOG, digest, ends_compatible, find_sites, methylation_warnings
from .operations import Step, renumber
from .primers import (Junction, Primer, analyse_junctions, check_overlap_uniqueness,
                      design_gibson_pair, design_re_pair)
from .seqtools import gc_fraction, revcomp, slice_region
from .vocab import PROTECTED_TYPES


def _dedupe(items: list[str]) -> list[str]:
    """Drop repeated warnings while keeping their original order.

    The same template can back several fragments, so an unqualified warning
    about it would otherwise be repeated once per PCR.
    """
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


@dataclass
class Route:
    id: str
    strategy: str                 # restriction | hifi | hybrid
    title: str
    steps: list = field(default_factory=list)
    primers: list = field(default_factory=list)
    junctions: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    cost: float = 0.0
    cost_breakdown: dict = field(default_factory=dict)
    feasible: bool = True
    #: Drawing data for the two-parent scheme, when this route is a swap.
    scheme: dict = field(default_factory=dict)

    def to_dict(self, cfg: Config) -> dict:
        return {
            "id": self.id, "strategy": self.strategy, "title": self.title,
            "steps": [s.to_dict() for s in self.steps],
            "primers": [p.to_dict(cfg) for p in self.primers],
            "junctions": [j.to_dict() for j in self.junctions],
            "warnings": _dedupe(self.warnings), "cost": round(self.cost, 2),
            "cost_breakdown": self.cost_breakdown, "feasible": self.feasible,
            "n_pcr": sum(1 for s in self.steps if s.operation == "PCR"),
            "n_assembly": sum(1 for s in self.steps
                              if s.operation in ("GIBSON_ASSEMBLY", "HI_FI_ASSEMBLY",
                                                 "LIGATION")),
            "scheme": self.scheme,
            "n_steps": len(self.steps),
            # Active bench time, and total elapsed turnaround including
            # overnight growth, sequencing and any synthesis order.
            "hands_on_minutes": sum(s.to_dict()["hands_on"] for s in self.steps),
            "turnaround_minutes": sum(s.to_dict()["minutes"] for s in self.steps),
        }


# --------------------------------------------------------------------------
# Fragment planning
# --------------------------------------------------------------------------

@dataclass
class PlannedFragment:
    name: str
    seq: str
    source_plasmid_id: int | None
    source_name: str
    start: int
    end: int
    strand: int
    circular_source: bool
    modules: list = field(default_factory=list)
    needs_synthesis: bool = False
    verified: str = "unverified"
    #: Orientation this fragment occupies in the construct, which is not
    #: the orientation it has in its donor.  A Cre-lox payload is stored
    #: inverted, so a part read forward from its donor goes in backwards.
    placed_strand: int = 1

    @property
    def amplified_strand(self) -> int:
        """Which strand of the donor to amplify to get construct orientation.

        Only a fallback: the frame is chosen by which orientation of the
        donor slice actually equals the sequence wanted.
        """
        return self.strand

    @property
    def placed_seq(self) -> str:
        """This fragment as the finished molecule reads it.

        Which is `seq` itself.  A Cre-lox payload is inverted by expand_dio,
        which hands back a module whose sequence is already the other strand
        and merely labels it as reverse; flipping it again here counted the
        inversion twice and cancelled it out, so the payload was amplified
        the right way up and assembled into a construct that wanted it
        upside down.
        """
        return self.seq


def _pair_with_layout(placed: list, items: list[dict]) -> list[tuple]:
    """Match each module to its own entry in the finished layout.

    Position is not enough.  The layout holds everything the molecule ends
    up with - the lox sites a Cre-lox cassette adds, and the ITRs, origin
    and resistance that come from the backbone - while the module list holds
    only what was asked for, so the two run at different lengths and zipping
    them walks out of step.  Every fragment after the first then inherits a
    neighbour's orientation, which is how an inverted payload came to be
    amplified the right way up and assembled into junctions that share
    nothing.

    Matching on identity instead, left to right, keeps a part that appears
    twice - a WPRE in the payload and another in the backbone - paired with
    the occurrence that is actually its own.
    """
    pairs: list[tuple] = []
    cursor = 0
    for rm in placed:
        for j in range(cursor, len(items)):
            it = items[j]
            if it.get("supplied_by_recipient"):
                # The recipient's own copy, kept rather than replaced. No
                # fragment, no reaction: the arc around it is amplified
                # anyway.
                continue
            same = (it.get("part_id") is not None
                    and it.get("part_id") == rm.part_id)
            if same or (it.get("name") == rm.part_name
                        and it.get("type") == rm.feature_type):
                pairs.append((rm, it))
                cursor = j + 1
                break
        else:
            # Nothing in the layout claims this module.  Pairing it with a
            # stranger would be worse than leaving it unplaced.
            continue
    return pairs


def group_contiguous(resolved: list[ResolvedModule], construct: AssembledConstruct,
                     cfg: Config) -> list[PlannedFragment]:
    """Merge neighbouring modules that come from the same source plasmid and
    sit adjacently there, so they can be amplified as one PCR product.

    This is where most of the laboratory saving comes from: a promoter and a
    DIO cassette already next to each other in a donor plasmid should be one
    amplicon, not two.
    """
    fragments: list[PlannedFragment] = []
    current: PlannedFragment | None = None

    # assemble() drops modules that resolved to no sequence, so the layout is
    # shorter than `resolved` whenever a part could not be found. Filter the
    # same way before pairing them up, or every module after the missing one
    # would be matched to the wrong layout entry.
    placed = [rm for rm in resolved if rm.seq]
    layout_items = [i for i in construct.layout if i["type"] != "backbone"]

    for rm, item in _pair_with_layout(placed, layout_items):
        source = rm.sources[0] if rm.sources else None
        if source is None:
            fragments.append(PlannedFragment(
                rm.part_name, rm.seq, None, "(synthesis)", 0, len(rm.seq), 1,
                False, [item], needs_synthesis=True))
            current = None
            continue

        # Which way along the template consecutive modules travel. It is
        # not the source strand alone: a Cre-lox payload is placed
        # inverted, so its modules run backwards through a donor that
        # holds them the right way round, and testing only the source
        # strand left P2A and EGFP - adjacent bases in one plasmid - as
        # two reactions, one of them 66 bp long.
        placed_strand = item.get("strand", 1)
        forward = source["strand"] * placed_strand == 1
        can_extend = (
            current is not None
            and current.source_plasmid_id == source["id"]
            and current.strand == source["strand"]
            and current.placed_strand == placed_strand
            and not current.needs_synthesis
            and (
                (forward and source["start"] == current.end)
                or (not forward and source["end"] == current.start)
            )
            and (len(current.seq) + len(rm.seq)) <= cfg.max_pcr_length
        )
        if can_extend and current is not None:
            current.seq += rm.seq
            current.name += f"+{rm.part_name}"
            current.modules.append(item)
            if forward:
                current.end = source["end"]
            else:
                current.start = source["start"]
        else:
            current = PlannedFragment(
                rm.part_name, rm.seq, source["id"], source["name"],
                source["start"], source["end"], source["strand"],
                bool(source.get("seq_length")), [item],
                verified=source.get("verified", "unverified"),
                placed_strand=item.get("strand", 1))
            fragments.append(current)
    return fragments


# --------------------------------------------------------------------------
# Route 1: PCR + HiFi / Gibson  (Stage 8)
# --------------------------------------------------------------------------

def plan_hifi_route(conn: sqlite3.Connection, construct: AssembledConstruct,
                    fragments: list[PlannedFragment], backbone: dict | None,
                    cfg: Config, kit: str = "HI_FI_ASSEMBLY") -> Route:
    route = Route(id="hifi", strategy="hifi",
                  title=f"PCR + {cfg.preferred_assembly_kit} assembly")
    steps: list[Step] = []
    primers: list[Primer] = []

    # What each PCR will actually yield, which is not always the catalogued
    # part: the recorded boundaries of an occurrence can run a few bases wide
    # of the part sequence.  Overlaps have to be taken from the amplicon that
    # will exist, or a junction is designed against bases the neighbouring
    # tube does not contain.
    realised: dict[int, tuple[str, int, int, str, bool]] = {}
    for idx, f in enumerate(fragments):
        if f.needs_synthesis or f.source_plasmid_id is None:
            realised[idx] = ("", 0, 0, f.placed_seq, False)
            continue
        row = conn.execute(
            "SELECT seq, circular FROM plasmid WHERE id = ?",
            (f.source_plasmid_id,)).fetchone()
        if row is None:
            realised[idx] = ("", 0, 0, f.placed_seq, False)
            continue
        circ = bool(row["circular"])
        tpl, st, en = _amplification_frame(
            row["seq"], f.start, f.end, f.amplified_strand,
            wanted=f.placed_seq, circular=circ)
        realised[idx] = (tpl, st, en,
                         slice_region(tpl, st, en, circular=circ), circ)

    # The finished molecule in order: the arcs of the recipient it keeps,
    # and the fragments supplying everything else.  Amplifying the whole
    # recipient instead of the arcs it keeps is how a product came to carry
    # the payload it was supposed to replace.
    ordered: list[dict] = []
    segments = construct.backbone_segments if backbone else []
    for i, seg in enumerate(segments):
        ordered.append({
            "kind": "backbone",
            "name": "backbone" if len(segments) == 1 else f"backbone{i + 1}",
            "seq": seg["seq"], "at": seg["construct_start"],
            "source_start": seg["source_start"], "source_end": seg["source_end"],
        })
    if backbone and not segments:
        # Nothing recorded which arcs survive, so the whole plasmid stands in.
        ordered.append({"kind": "backbone", "name": "backbone",
                        "seq": backbone["seq"], "at": -1,
                        "source_start": 0, "source_end": len(backbone["seq"])})
    bb_id = backbone["id"] if backbone else None
    bb_len = len(backbone["seq"]) if backbone else 0
    for idx, f in enumerate(fragments):
        tpl = realised[idx][0]
        entry = {
            "kind": "fragment", "name": f.name, "seq": realised[idx][3],
            "at": f.modules[0]["start"] if f.modules else 0, "index": idx,
            "synthetic": f.needs_synthesis,
        }
        # A module the recipient itself supplies is, as far as the bench is
        # concerned, part of the recipient: same template, same strand,
        # neighbouring coordinates.
        if (bb_id is not None and f.source_plasmid_id == bb_id
                and tpl == backbone["seq"]):
            entry["source_start"] = realised[idx][1]
            entry["source_end"] = realised[idx][2]
        ordered.append(entry)
    ordered.sort(key=lambda o: o["at"])

    if backbone:
        ordered, absorbed = _merge_contiguous(ordered, bb_len)
    else:
        absorbed = set()
    absorbed = set(absorbed)

    # A kept arc can be a handful of bases - the gap between a promoter and
    # the coding sequence after it - and there is no such thing as a 25 bp
    # PCR.  Anything too short to amplify is carried on the next piece's
    # forward primer instead, which puts the same bases in the same place
    # for the cost of a longer oligo.
    carried: dict[int, str] = {}
    keep: list[dict] = []
    pending = ""
    for o in ordered:
        # Too short to be its own reaction, whichever side it came from.
        # A Cre-lox cassette's sites are 34 bp: ordering four of those as
        # separate synthetic fragments and assembling them is not how
        # anyone builds a cassette, and they cannot carry an overlap on
        # both sides while also being the lox site.
        short = len(o["seq"]) < cfg.min_pcr_length
        if short and (o["kind"] == "backbone" or o.get("synthetic")):
            pending += o["seq"]
            if o["kind"] == "fragment":
                absorbed.add(o["index"])
            continue
        if pending:
            carried[len(keep)] = pending
            pending = ""
        keep.append(o)
    if keep:
        if pending:
            # Arcs after the last amplified piece wrap onto the first.
            carried[0] = pending + carried.get(0, "")
        ordered = keep
    else:
        carried = {}

    pieces: list[tuple[str, str]] = [
        (o["name"], carried.get(i, "") + o["seq"])
        for i, o in enumerate(ordered)]
    if not pieces:
        route.feasible = False
        route.warnings.append("nothing to assemble")
        return route

    n = len(pieces)
    templates_needing_dpni: list[str] = []

    # --- backbone amplification ------------------------------------------
    for k, piece in enumerate(ordered):
        if piece["kind"] != "backbone":
            continue
        bb_seq = piece["seq"]
        up = pieces[k - 1][1]
        down = pieces[(k + 1) % n][1]
        fwd, rev = design_gibson_pair(
            0, piece["name"], backbone["seq"],
            piece["source_start"], piece["source_end"], up, down, cfg,
            circular_template=bool(backbone.get("circular")),
            carry_downstream=False, unique_in=construct.seq)
        _carry(fwd, carried.get(k, ""), route, cfg)
        _check_unique(fwd, construct.seq, route)
        primers += [fwd, rev]
        if len(bb_seq) > cfg.max_pcr_length:
            route.warnings.append(
                f"the {piece['name']} amplicon is {len(bb_seq)} bp, above the "
                f"{cfg.max_pcr_length} bp PCR limit in the configuration; "
                f"linearising the backbone with a restriction digest is "
                f"usually more reliable than amplifying it")
        steps.append(Step(
            0, "PCR",
            f"Amplify the {len(bb_seq):,} bp {piece['name']} arc "
            f"({piece['source_start']:,}-{piece['source_end']:,}) from "
            f"{backbone['name']} with primers {fwd.name} / {rev.name}",
            inputs=[backbone["name"], f"{fwd.name}+{rev.name}"],
            output=f"{piece['name']} amplicon",
            details={"length": len(bb_seq), "template": backbone["name"]}))
        if backbone["name"] not in templates_needing_dpni:
            templates_needing_dpni.append(backbone["name"])

    # --- insert amplification ---------------------------------------------
    where = {o["index"]: k for k, o in enumerate(ordered)
             if o["kind"] == "fragment"}
    for i, frag in enumerate(fragments):
        if i in absorbed:
            continue            # amplified as part of the recipient's arc
        piece_index = where.get(i, i)
        up = pieces[piece_index - 1][1]
        down = pieces[(piece_index + 1) % n][1]

        if frag.needs_synthesis:
            if not any(m.get("allow_synthesis") for m in frag.modules) and \
                    not cfg.allow_unverified_templates:
                route.warnings.append(
                    f"{frag.name} has no lab source and would need synthesis")
            steps.append(Step(
                0, "SYNTHESIS",
                f"Order {frag.name} as a synthetic fragment "
                f"({len(frag.seq)} bp) with {cfg.overlap_length} bp overlaps",
                output=f"{frag.name} gBlock",
                details={"length": len(frag.seq)}))
            continue

        row = conn.execute("SELECT seq, circular, name, verified FROM plasmid WHERE id = ?",
                           (frag.source_plasmid_id,)).fetchone()
        tpl, start, end, _body, circ = realised[i]
        if not tpl:
            tpl, start, end = frag.seq, 0, len(frag.seq)

        fwd, rev = design_gibson_pair(
            i, frag.name, tpl, start, end, up, down, cfg,
            circular_template=circ, carry_downstream=False,
            unique_in=construct.seq)
        _carry(fwd, carried.get(piece_index, ""), route, cfg)
        # Report where the part is in the donor, not where it is in the
        # strand the primers were written against.
        for pr in (fwd, rev):
            pr.target_start, pr.target_end = frag.start, frag.end
        primers += [fwd, rev]

        if len(frag.seq) > cfg.max_pcr_length:
            route.warnings.append(
                f"{frag.name} amplicon is {len(frag.seq)} bp, above the "
                f"{cfg.max_pcr_length} bp limit set in the configuration")
        if len(frag.seq) < cfg.min_pcr_length:
            # Short elements are cheaper and more reliable to carry on the
            # neighbouring primer's tail, or to order as annealed oligos,
            # than to amplify and purify as their own fragment.
            route.warnings.append(
                f"{frag.name} is only {len(frag.seq)} bp; rather than a "
                f"separate PCR, add it to an adjacent fragment's primer tail "
                f"or order it as annealed oligos")
        if row and row["verified"] == "unverified" and not cfg.allow_unverified_templates:
            route.warnings.append(
                f"template {frag.source_name} is not sequence-verified")

        steps.append(Step(
            0, "PCR",
            f"Amplify {frag.name} ({len(frag.seq)} bp) from {frag.source_name} "
            f"with primers {fwd.name} / {rev.name}",
            inputs=[frag.source_name, f"{fwd.name}+{rev.name}"],
            output=f"{frag.name} amplicon",
            details={"length": len(frag.seq), "template": frag.source_name,
                     "gc": round(gc_fraction(frag.seq), 3)}))
        templates_needing_dpni.append(frag.source_name)

    # --- template removal and cleanup -------------------------------------
    if templates_needing_dpni:
        steps.append(Step(
            0, "DPNI_DIGEST",
            "Digest the PCR reactions with DpnI to remove the parental "
            "plasmid templates (" + ", ".join(sorted(set(templates_needing_dpni))[:4]) + ")",
            output="template-free amplicons"))

    # Gel extraction only where it is actually warranted (Stage 8 rule).
    gel_needed = [f.name for f in fragments if _expects_extra_bands(f, cfg)]
    if gel_needed:
        steps.append(Step(
            0, "GEL_EXTRACTION",
            "Gel purify " + ", ".join(gel_needed) +
            " (multiple products or a size-selection requirement)",
            output="size-selected fragments"))
        remaining = [p[0] for p in pieces if p[0] not in gel_needed]
        if remaining:
            steps.append(Step(0, "DNA_CLEANUP",
                              "Column purify the remaining PCR products",
                              output="clean amplicons"))
    else:
        steps.append(Step(0, "DNA_CLEANUP",
                          "Column purify all PCR products",
                          output="clean amplicons"))

    # --- assembly ----------------------------------------------------------
    names = [p[0] for p in pieces]
    steps.append(Step(
        0, kit,
        f"Assemble {len(names)} fragments ({' + '.join(names)}) using "
        f"{cfg.preferred_assembly_kit}, {cfg.overlap_length} bp overlaps",
        inputs=names, output="assembled plasmid",
        details={"fragments": len(names), "overlap": cfg.overlap_length}))

    _append_finishing_steps(steps, cfg, construct)

    route.steps = renumber(steps)
    route.primers = primers
    route.junctions = analyse_junctions(construct.layout, construct.seq, "hifi", cfg)
    route.warnings += check_overlap_uniqueness(construct.layout, construct.seq, cfg)
    for j in route.junctions:
        route.warnings += [f"junction {j.index} ({j.left_name}|{j.right_name}): {w}"
                           for w in j.warnings]
    for p in primers:
        route.warnings += [f"primer {p.name}: {w}" for w in p.warnings]
    return route


def _expects_extra_bands(frag: PlannedFragment, cfg: Config) -> bool:
    """Decide whether a PCR genuinely needs gel purification."""
    if frag.needs_synthesis:
        return False
    if len(frag.seq) > 5000:
        return True            # long amplicons usually carry shorter by-products
    if gc_fraction(frag.seq) > 0.70:
        return True            # GC-rich reactions are the messy ones
    return False


def _append_finishing_steps(steps: list[Step], cfg: Config,
                            construct: AssembledConstruct) -> None:
    if cfg.include_transformation_steps:
        steps.append(Step(0, "TRANSFORMATION",
                          "Transform the assembly into competent cells and "
                          "plate on the appropriate antibiotic",
                          output="colonies"))
        screen = _screening_enzyme(construct, cfg)
        detail = (f"colony PCR, or diagnostic digest with {screen}"
                  if screen else "colony PCR")
        steps.append(Step(0, "COLONY_SCREEN",
                          f"Screen 4-8 colonies by {detail}",
                          output="candidate clones",
                          details={"diagnostic_enzyme": screen}))
    if cfg.include_sequencing_steps:
        steps.append(Step(0, "SEQUENCE_VERIFY",
                          "Sanger sequence across every junction, then confirm "
                          "the full insert",
                          output="verified plasmid"))


def _screening_enzyme(construct: AssembledConstruct, cfg: Config) -> str | None:
    """Pick a stocked enzyme that cuts the product a small number of times."""
    best = None
    for name in cfg.stocked_enzymes:
        enz = CATALOG.get(name)
        if not enz:
            continue
        n = len(find_sites(construct.seq, enz, circular=construct.circular))
        if 2 <= n <= 3:
            return name
        if n == 1 and best is None:
            best = name
    return best


# --------------------------------------------------------------------------
# Route 1b: order the insert, amplify the vector  (Stage 8)
# --------------------------------------------------------------------------

def plan_gblock_route(conn: sqlite3.Connection, construct: AssembledConstruct,
                      backbone: dict | None, cfg: Config,
                      kit: str = "HI_FI_ASSEMBLY") -> Route:
    """Order everything that is not the recipient as one synthetic fragment.

    Some things are cheaper to buy than to build. A Cre-lox cassette on a
    recipient that has none is the clearest case: its four lox sites are 34
    bp each, too short to amplify and too short to carry an assembly
    overlap on both sides while also being the lox site, so the PCR route
    ends up ordering them as separate oligos and joining five pieces to
    make a stretch somebody will sell as one.

    This takes the whole non-recipient region instead - cassette, payload
    and all - and orders it. One vector PCR, one fragment, one assembly.
    """
    route = Route(id="gblock", strategy="gblock",
                  title="Order the insert, amplify the vector")
    if not backbone or not construct.backbone_segments:
        route.feasible = False
        route.warnings.append(
            "there is no recipient to keep, so there is nothing to order "
            "the rest of: this construct is built from parts")
        return route

    # Keep the longest arc of the recipient and buy everything else, which
    # leaves exactly one junction at each end.
    keep = max(construct.backbone_segments, key=lambda s: len(s["seq"]))
    kept = keep["seq"]
    start = (keep["construct_start"] + len(kept)) % construct.length
    insert_len = construct.length - len(kept)
    if insert_len <= 0:
        route.feasible = False
        route.warnings.append("the recipient already is this construct")
        return route

    insert = slice_region(construct.seq, start, keep["construct_start"],
                          circular=True)
    arm = cfg.overlap_length
    ordered_seq = (slice_region(construct.seq, start - arm, start,
                                circular=True)
                   + insert
                   + slice_region(construct.seq, keep["construct_start"],
                                  keep["construct_start"] + arm,
                                  circular=True))

    if len(ordered_seq) > cfg.max_synthetic_fragment:
        route.warnings.append(
            f"the fragment to order is {len(ordered_seq):,} bp, past the "
            f"{cfg.max_synthetic_fragment:,} bp a gBlock runs to; an eBlock "
            f"reaches further, or order it in two overlapping pieces")

    steps: list[Step] = []
    primers: list[Primer] = []

    # The ordered fragment already carries homology to the vector at both
    # ends, so the vector primers carry none: a tail here would put the
    # start of the arc in twice and leave neither junction able to anneal.
    bb_f, bb_r = design_gibson_pair(
        0, "backbone", backbone["seq"], keep["source_start"],
        keep["source_end"], "", "", cfg,
        circular_template=bool(backbone.get("circular")),
        carry_downstream=False, unique_in=construct.seq)
    primers += [bb_f, bb_r]
    _check_unique(bb_f, construct.seq, route)

    steps.append(Step(
        0, "SYNTHESIS",
        f"Order the insert as a synthetic fragment: {len(ordered_seq):,} bp, "
        f"carrying {arm} bp of homology to the vector at each end",
        output="synthetic insert",
        details={"length": len(ordered_seq), "sequence": ordered_seq,
                 "homology_arm": arm}))
    steps.append(Step(
        0, "PCR",
        f"Amplify the {len(kept):,} bp vector arc "
        f"({keep['source_start']:,}-{keep['source_end']:,}) from "
        f"{backbone['name']} with primers {bb_f.name} / {bb_r.name}",
        inputs=[backbone["name"], f"{bb_f.name}+{bb_r.name}"],
        output="linearised vector",
        details={"length": len(kept), "template": backbone["name"]}))
    steps.append(Step(
        0, "DPNI_DIGEST",
        f"Digest the PCR with DpnI to remove the {backbone['name']} template",
        inputs=["vector PCR"], output="template-free vector"))
    steps.append(Step(
        0, "DNA_CLEANUP", "Column purify the vector PCR",
        output="clean vector"))
    steps.append(Step(
        0, kit,
        f"Assemble the vector with the ordered insert using "
        f"{cfg.preferred_assembly_kit} ({arm} bp overlaps)",
        inputs=["linearised vector", "synthetic insert"],
        output="assembled plasmid"))
    _append_finishing_steps(steps, cfg, construct)

    route.steps = renumber(steps)
    route.primers = primers
    route.junctions = analyse_junctions(construct.layout, construct.seq,
                                        "gblock", cfg)
    route.warnings.append(
        f"one fragment to order and one reaction to run, against "
        f"{len(construct.backbone_segments)} arc(s) of recipient kept")
    return route


# --------------------------------------------------------------------------
# Route 2: restriction digest + ligation  (Stage 7)
# --------------------------------------------------------------------------

def protected_spans(construct: AssembledConstruct) -> list[tuple[int, int, str]]:
    spans = []
    for item in construct.layout:
        if item["type"] in PROTECTED_TYPES:
            spans.append((item["start"], item["end"], item["name"]))
    return spans


def plan_restriction_route(conn: sqlite3.Connection, construct: AssembledConstruct,
                           fragments: list[PlannedFragment], backbone: dict | None,
                           cfg: Config) -> Route:
    """Find an enzyme pair that opens the backbone and excises the insert."""
    route = Route(id="restriction", strategy="restriction",
                  title="Restriction digest + ligation")

    if backbone is None:
        route.feasible = False
        route.warnings.append("a restriction route needs an explicit backbone")
        return route

    # A restriction route can only excise an insert that already exists as one
    # contiguous piece in a single donor. Fragments here have already been
    # merged by donor and adjacency, so more than one fragment means the
    # insert would have to be rebuilt from separate pieces - which for
    # restriction cloning means a multi-fragment ligation.
    real = [f for f in fragments if not f.needs_synthesis]
    donors = {f.source_plasmid_id for f in real}
    MAX_LIGATION_INSERTS = 2

    if len(fragments) > MAX_LIGATION_INSERTS:
        route.feasible = False
        route.warnings.append(
            f"the insert is made of {len(fragments)} pieces from "
            f"{len(donors)} donor plasmid(s); a ligation of that many "
            f"fragments is not practical, so this construct needs the "
            f"PCR/HiFi route (or a staged build of two pieces at a time)")
        return route
    if len(fragments) > 1:
        route.warnings.append(
            f"{len(fragments)}-piece ligation; expect low efficiency and screen "
            f"more colonies than usual")

    single_donor_contiguous = len(fragments) == 1 and bool(real)
    insert_seq = "".join(f.seq for f in fragments)
    bb_seq = backbone["seq"]
    bb_circ = bool(backbone.get("circular"))
    protected = protected_spans(construct)

    candidates = _find_enzyme_pairs(bb_seq, insert_seq, construct, protected, cfg,
                                    bb_circ)
    if not single_donor_contiguous:
        # Without a single contiguous donor the sites have to be added by PCR.
        candidates = [c for c in candidates if c["mode"] == "add_by_pcr"] or [
            dict(c, mode="add_by_pcr",
                 warnings=c["warnings"] + [
                     f"{c['enzymes'][0]} sites must be added by PCR because the "
                     f"insert is not contiguous in a single donor"])
            for c in candidates]
    if not candidates:
        route.feasible = False
        route.warnings.append(
            "no stocked enzyme pair cuts the backbone uniquely while leaving "
            "the insert and all protected features intact; use the HiFi route "
            "or add sites by PCR")
        return route

    best = candidates[0]
    e1, e2 = best["enzymes"]
    steps: list[Step] = []
    primers: list[Primer] = []

    if best["mode"] == "native":
        steps.append(Step(
            0, "RE_DIGEST",
            f"Digest {backbone['name']} with {e1} + {e2}",
            inputs=[backbone["name"], f"{e1}, {e2}"], output="linearised backbone",
            details={"enzymes": [e1, e2], "buffer": CATALOG[e1].buffer}))
        donor = fragments[0].source_name if fragments else "donor plasmid"
        steps.append(Step(
            0, "RE_DIGEST",
            f"Digest {donor} with {e1} + {e2} to release the insert",
            inputs=[donor, f"{e1}, {e2}"], output="digested donor",
            details={"enzymes": [e1, e2]}))
    else:
        # Sites must be added by PCR first.
        for frag in fragments:
            row = conn.execute("SELECT seq, circular FROM plasmid WHERE id = ?",
                               (frag.source_plasmid_id,)).fetchone()
            tmpl = row["seq"] if row else frag.seq
            fwd, rev = design_re_pair(frag.name, tmpl, frag.start, frag.end,
                                      e1, e2, cfg,
                                      circular_template=bool(row and row["circular"]))
            primers += [fwd, rev]
            steps.append(Step(
                0, "PCR",
                f"Amplify {frag.name} from {frag.source_name} with "
                f"{fwd.name} / {rev.name}, adding {e1} and {e2} sites",
                inputs=[frag.source_name, f"{fwd.name}+{rev.name}"],
                output=f"{frag.name} amplicon with sites",
                details={"length": len(frag.seq)}))
        steps.append(Step(0, "DNA_CLEANUP", "Column purify the PCR products",
                          output="clean amplicons"))
        steps.append(Step(
            0, "RE_DIGEST", f"Digest the amplicons with {e1} + {e2}",
            inputs=[f"{e1}, {e2}"], output="digested insert"))
        steps.append(Step(
            0, "RE_DIGEST", f"Digest {backbone['name']} with {e1} + {e2}",
            inputs=[backbone["name"]], output="linearised backbone"))

    if e1 == e2:
        steps.append(Step(
            0, "DEPHOSPHORYLATION",
            "Dephosphorylate the cut vector to suppress self-ligation "
            "(single-enzyme digest leaves compatible ends)",
            output="dephosphorylated vector"))

    steps.append(Step(
        0, "GEL_EXTRACTION",
        f"Gel purify the {best['backbone_band']} bp linearised backbone band",
        output="purified backbone", details={"size": best["backbone_band"]}))
    steps.append(Step(
        0, "GEL_EXTRACTION",
        f"Gel purify the {best['insert_band']} bp insert band",
        output="purified insert", details={"size": best["insert_band"]}))
    steps.append(Step(
        0, "LIGATION",
        f"Ligate the insert into the backbone "
        f"({best['left_end']} / {best['right_end']} ends), 3:1 insert:vector",
        inputs=["purified backbone", "purified insert"], output="ligated plasmid"))

    _append_finishing_steps(steps, cfg, construct)

    route.steps = renumber(steps)
    route.primers = primers
    route.junctions = analyse_junctions(construct.layout, construct.seq,
                                        "ligation", cfg)
    route.warnings += best["warnings"]
    for p in primers:
        route.warnings += [f"primer {p.name}: {w}" for w in p.warnings]
    route.cost_breakdown["enzymes"] = [e1, e2]
    return route


def _find_enzyme_pairs(bb_seq: str, insert_seq: str,
                       construct: AssembledConstruct,
                       protected: list[tuple[int, int, str]], cfg: Config,
                       bb_circular: bool) -> list[dict]:
    """Search for enzyme pairs that cleanly open the backbone.

    Implements the Stage 7 checklist: uniqueness, absence from the insert,
    no cuts inside protected features, end compatibility, and whether a
    single or double digest is possible.
    """
    results: list[dict] = []
    names = [n for n in cfg.stocked_enzymes if n in CATALOG]

    bb_counts, insert_counts, product_counts = {}, {}, {}
    for n in names:
        enz = CATALOG[n]
        bb_counts[n] = len(find_sites(bb_seq, enz, bb_circular))
        insert_counts[n] = len(find_sites(insert_seq, enz, False))
        product_counts[n] = find_sites(construct.seq, enz, construct.circular)

    def hits_protected(name: str) -> str | None:
        for site in product_counts[name]:
            for start, end, label in protected:
                if start <= site.top_cut < end:
                    return label
        return None

    usable = []
    for n in names:
        if bb_counts[n] != 1:
            continue
        if insert_counts[n] != 0:
            continue
        clash = hits_protected(n)
        if clash:
            continue
        usable.append(n)

    for i, a in enumerate(usable):
        for b in usable[i:]:
            if a == b and len(usable) > 1:
                continue
            enz_a, enz_b = CATALOG[a], CATALOG[b]
            warnings: list[str] = []
            if enz_a.end_type != enz_b.end_type and "blunt" not in (
                    enz_a.end_type, enz_b.end_type):
                pass  # different overhang classes are fine on opposite ends
            if ends_compatible("N" * max(enz_a.overhang_len, 0), enz_a.end_type,
                               "N" * max(enz_b.overhang_len, 0), enz_b.end_type) \
                    and a != b:
                warnings.append(
                    f"{a} and {b} leave compatible ends, so the vector can "
                    f"self-ligate; dephosphorylate or pick another pair")
            for name in (a, b):
                warnings += methylation_warnings(bb_seq, CATALOG[name], bb_circular)
                if CATALOG[name].star_activity:
                    warnings.append(f"{name} shows star activity; avoid long "
                                    f"digests and high enzyme excess")
            if enz_a.buffer != enz_b.buffer:
                warnings.append(f"{a} and {b} prefer different buffers; "
                                f"consider sequential digestion")

            insert_band = len(insert_seq)
            backbone_band = len(bb_seq)
            if insert_band < 100:
                warnings.append(f"insert band is only {insert_band} bp and will "
                                f"be hard to recover from a gel")
            if abs(insert_band - backbone_band) < 200:
                warnings.append(
                    f"insert ({insert_band} bp) and backbone ({backbone_band} bp) "
                    f"bands are too close in size to separate cleanly")

            results.append({
                "enzymes": (a, b), "mode": "native",
                "warnings": warnings,
                "insert_band": insert_band, "backbone_band": backbone_band,
                "left_end": enz_a.end_type, "right_end": enz_b.end_type,
                "score": len(warnings) + (2 if a == b else 0),
            })

    if not results:
        # Fall back to adding sites by PCR, which only needs the enzyme to be
        # absent from the insert and the backbone to be cuttable.
        for n in names:
            if insert_counts[n] == 0 and bb_counts[n] == 1 and not hits_protected(n):
                results.append({
                    "enzymes": (n, n), "mode": "add_by_pcr",
                    "warnings": [f"{n} sites must be added by PCR"],
                    "insert_band": len(insert_seq), "backbone_band": len(bb_seq),
                    "left_end": CATALOG[n].end_type, "right_end": CATALOG[n].end_type,
                    "score": 5,
                })

    results.sort(key=lambda r: r["score"])
    return results


# --------------------------------------------------------------------------
# Stage 10: ranking
# --------------------------------------------------------------------------

def score_route(route: Route, cfg: Config, history: dict | None = None) -> Route:
    """Cost a route with the configurable weights. Lower is better."""
    w = cfg.cost_weights
    breakdown: dict[str, float] = {}

    counts: dict[str, int] = {}
    for s in route.steps:
        counts[s.operation] = counts.get(s.operation, 0) + 1

    def add(key: str, n: float, weight_key: str):
        if n:
            breakdown[key] = round(n * w.get(weight_key, 0.0), 2)

    add("assembly reactions", counts.get("GIBSON_ASSEMBLY", 0)
        + counts.get("HI_FI_ASSEMBLY", 0) + counts.get("LIGATION", 0),
        "assembly_reaction")
    add("PCRs", counts.get("PCR", 0), "pcr")
    add("primers", len(route.primers), "primer")
    add("digests", counts.get("RE_DIGEST", 0), "restriction_digest")
    add("gel extractions", counts.get("GEL_EXTRACTION", 0), "gel_extraction")
    # Synthesis is sold by the base, and that is the whole trade here: a
    # 2.7 kb fragment is not the same purchase as a 200 bp one, and costing
    # them alike made ordering the insert look cheaper than any number of
    # PCRs.
    add("synthesis", counts.get("SYNTHESIS", 0), "synthesis")
    ordered_kb = sum(s.details.get("length", 0) for s in route.steps
                     if s.operation == "SYNTHESIS") / 1000
    if ordered_kb:
        # A configuration saved before this weight existed still has to
        # cost ordered DNA, or buying looks free to anyone who has run
        # Genesis before today.
        breakdown["synthesised DNA"] = round(
            ordered_kb * w.get("synthesis_per_kb", 90.0), 2)

    long_kb = sum(max(s.details.get("length", 0) - 2000, 0) / 1000
                  for s in route.steps if s.operation == "PCR")
    add("long PCR", long_kb, "long_pcr_per_kb")

    difficult = sum(1 for j in route.junctions if j.warnings)
    add("difficult junctions", difficult, "difficult_junction")

    # Deduplicated, so one shared template is costed once rather than once
    # per fragment it supplies.
    unverified = sum(1 for warn in _dedupe(route.warnings)
                     if "not sequence-verified" in warn)
    add("unverified templates", unverified, "unverified_template")

    non_stocked = sum(1 for e in route.cost_breakdown.get("enzymes", [])
                      if e not in cfg.stocked_enzymes)
    add("non-stocked enzymes", non_stocked, "non_stocked_enzyme")

    if history:
        fails = history.get(route.strategy, 0)
        add("past failures", fails, "historical_failure")

    route.cost_breakdown.update(breakdown)
    route.cost = sum(breakdown.values())
    if not route.feasible:
        route.cost += 1000.0
    return route


def rank_routes(routes: list[Route], cfg: Config,
                history: dict | None = None) -> list[Route]:
    for r in routes:
        score_route(r, cfg, history)
    return sorted(routes, key=lambda r: (not r.feasible, r.cost))


def failure_history(conn: sqlite3.Connection) -> dict:
    """Stage 11 feedback: how often has each strategy failed in this lab?"""
    rows = conn.execute(
        """SELECT d.plan, o.result FROM outcome o JOIN design d ON d.id = o.design_id
           WHERE o.result = 'failure'""").fetchall()
    hist: dict[str, int] = {}
    import json
    for r in rows:
        try:
            plan = json.loads(r["plan"])
        except (json.JSONDecodeError, TypeError):
            continue
        for route in plan.get("routes", []):
            hist[route.get("strategy", "?")] = hist.get(route.get("strategy", "?"), 0) + 1
    return hist


# --------------------------------------------------------------------------
# Two-parent swap routes
#
# The common case in practice: take the closest existing plasmid, cut one
# region out of it, and put a piece from a second plasmid in its place.
# Two parents, one junction pair, one assembly reaction.
# --------------------------------------------------------------------------

def _check_unique(primer: Primer, construct_seq: str, route: "Route") -> None:
    """Say so when an overlap arm cannot point at one place only.

    An arm has to say *where* a fragment goes.  When the same bases occur
    twice in the finished molecule - two copies of a WPRE, the pair of lox
    sites a Cre-lox cassette needs - the arm matches both and the fragment
    can go in either place, either way round.  That is not something the
    assembly resolves on its own, so it has to be said before anyone orders
    the oligo.
    """
    from .primers import _occurrences
    tail = primer.tail or ""
    if not tail or not construct_seq:
        return
    n = _occurrences(tail, construct_seq)
    if n > 1:
        route.warnings.append(
            f"primer {primer.name}: its {len(tail)} bp overlap occurs {n} "
            f"times in the construct, so this fragment has more than one "
            f"place it can join and the assembly may put it in either; "
            f"a repeated element here needs a different backbone, or the "
            f"junction moved")


def _merge_contiguous(ordered: list[dict], template_length: int
                      ) -> tuple[list[dict], set]:
    """Amplify neighbouring stretches of the same plasmid in one reaction.

    The recipient's kept arcs and any module the recipient itself supplies
    are the same DNA, lying next to each other on the same plasmid. Planning
    them as separate reactions asks for three PCRs and two assembly
    junctions to rebuild a stretch that was never taken apart - and it was
    only separate because one of them answers a module in the specification
    and the others do not.

    Returns the merged list, and the indices of the fragments swallowed, so
    nothing designs primers for a piece that no longer exists.
    """
    absorbed: set = set()
    if template_length <= 0 or len(ordered) < 2:
        return ordered, absorbed

    def joins(a: dict, b: dict) -> bool:
        for piece in (a, b):
            if not isinstance(piece.get("source_start"), int):
                return False
            if not isinstance(piece.get("source_end"), int):
                return False
        return a["source_end"] % template_length == \
            b["source_start"] % template_length

    pieces = list(ordered)
    changed = True
    while changed and len(pieces) > 1:
        changed = False
        for i in range(len(pieces)):
            j = (i + 1) % len(pieces)
            if i == j:
                break
            a, b = pieces[i], pieces[j]
            if not joins(a, b):
                continue
            merged = {
                "kind": "backbone",
                "name": a["name"] if a["kind"] == "backbone" else b["name"],
                "seq": a["seq"] + b["seq"],
                "at": a["at"],
                "source_start": a["source_start"],
                "source_end": b["source_end"],
            }
            for piece in (a, b):
                if piece["kind"] == "fragment":
                    absorbed.add(piece["index"])
            if j == 0:
                # Wrapped: the merged piece takes the first slot.
                pieces = [merged] + pieces[1:i]
            else:
                pieces = pieces[:i] + [merged] + pieces[j + 1:]
            changed = True
            break
    return pieces, absorbed


def _carry(fwd: Primer, extra: str, route: "Route", cfg: Config) -> None:
    """Put an un-amplifiable arc onto the primer that makes the next piece.

    The bases still have to exist in some tube.  Adding them after the
    overlap means the amplicon reads overlap-arc-body, which is the same
    molecule a separate reaction would have given, at the price of a longer
    oligo rather than a PCR nobody can run.
    """
    if not extra:
        return
    fwd.tail = fwd.tail + extra
    fwd.sequence = fwd.tail + fwd.anneal
    if len(fwd.sequence) > 100:
        route.warnings.append(
            f"primer {fwd.name} is {len(fwd.sequence)} nt because it carries "
            f"a {len(extra)} bp stretch of the recipient that is too short to "
            f"amplify; order it as an ultramer, or synthesise that junction")


def _amplification_frame(template: str, start: int, end: int,
                         strand: int, wanted: str = "",
                         circular: bool = False) -> tuple[str, int, int]:
    """Template and coordinates for amplifying a region in construct orientation.

    A part sits on whichever strand of its donor it happens to sit on, but a
    PCR product's top strand is whatever the template reads forward.  Amplify
    a minus-strand part off the plus strand and it comes out backwards, so
    its assembly overlaps get computed in one frame and attached in the
    other: the junctions then share no sequence and the assembly will not
    close.  Designing against the strand the construct uses keeps one frame
    throughout.  The primers are the same molecules either way; only which
    strand they are written against changes.

    The recorded strand is not reliable - an occurrence can be stored plus
    while the donor actually reads the part backwards - and believing it
    inverts the fragment.  So when the sequence the construct needs is
    known, take the orientation whose slice is actually equal to it and
    treat the record as a fallback.
    """
    n = len(template)
    # Mapping a region onto the other strand is a subtraction, and on a
    # region that wraps the origin - end past the last base, which is how
    # a wrapped occurrence is recorded - the plain subtraction lands
    # before position zero and the slice comes back empty or short.
    # Reduce the far end into the sequence first, then measure from it.
    low = (n - end) % n if n else 0
    flipped = (revcomp(template), low, low + (end - start))
    if wanted:
        here = slice_region(template, start, end, circular=circular)
        if here == wanted:
            return template, start, end
        if revcomp(here) == wanted:
            return flipped
    if strand != -1:
        return template, start, end
    return flipped


def _insert_sequence(conn: sqlite3.Connection, analysis) -> str:
    if analysis.donor is None or analysis.insert is None:
        return ""
    row = conn.execute("SELECT seq, circular FROM plasmid WHERE id = ?",
                       (analysis.donor.plasmid_id,)).fetchone()
    if row is None:
        return ""
    sub = slice_region(row["seq"], analysis.insert.start, analysis.insert.end,
                       circular=bool(row["circular"]))
    return revcomp(sub) if analysis.insert.strand == -1 else sub


def plan_swap_hifi_route(conn: sqlite3.Connection, analysis, cfg: Config,
                         scheme: dict) -> Route:
    """Linearise the recipient by inverse PCR, amplify the insert, assemble."""
    route = Route(id="swap_hifi", strategy="swap_hifi",
                  title="Swap by PCR + HiFi assembly (2 parents)",
                  scheme=scheme)
    if analysis is None or not analysis.feasible:
        route.feasible = False
        route.warnings.append(
            (analysis.notes[0] if analysis and analysis.notes
             else "no suitable pair of parent plasmids was found"))
        return route

    rec = conn.execute("SELECT seq, circular, name, verified FROM plasmid WHERE id = ?",
                       (analysis.recipient.plasmid_id,)).fetchone()
    don = conn.execute("SELECT seq, circular, name, verified FROM plasmid WHERE id = ?",
                       (analysis.donor.plasmid_id,)).fetchone()
    if rec is None or don is None:
        route.feasible = False
        route.warnings.append("parent plasmid sequence missing from the database")
        return route

    rec_seq, rec_circ = rec["seq"], bool(rec["circular"])
    remove = analysis.remove
    insert_seq = _insert_sequence(conn, analysis)

    # The retained backbone runs from the end of the removed region, around
    # the plasmid, back to its start.
    backbone_seq = slice_region(rec_seq, remove.end, remove.start,
                                circular=rec_circ)
    if not backbone_seq or not insert_seq:
        route.feasible = False
        route.warnings.append("could not derive the backbone or insert sequence")
        return route

    steps: list[Step] = []
    primers: list[Primer] = []

    # Inverse PCR across the recipient, skipping the region being removed.
    bb_f, bb_r = design_gibson_pair(
        0, "backbone", rec_seq, remove.end, remove.start, insert_seq,
        insert_seq, cfg, circular_template=rec_circ, carry_downstream=False)
    primers += [bb_f, bb_r]
    steps.append(Step(
        0, "PCR",
        f"Inverse PCR {rec['name']} with {bb_f.name} / {bb_r.name} to "
        f"linearise it and delete {remove.label} "
        f"({remove.length:,} bp removed)",
        inputs=[rec["name"], f"{bb_f.name}+{bb_r.name}"],
        output="linearised backbone",
        details={"length": len(backbone_seq), "template": rec["name"],
                 "gc": round(gc_fraction(backbone_seq), 3)}))

    ins_tpl, ins_start, ins_end = _amplification_frame(
        don["seq"], analysis.insert.start, analysis.insert.end,
        analysis.insert.strand, wanted=insert_seq,
        circular=bool(don["circular"]))
    ins_f, ins_r = design_gibson_pair(
        1, analysis.insert.label or "insert", ins_tpl,
        ins_start, ins_end, backbone_seq, backbone_seq,
        cfg, circular_template=bool(don["circular"]), carry_downstream=False)
    for pr in (ins_f, ins_r):
        pr.target_start, pr.target_end = (analysis.insert.start,
                                          analysis.insert.end)
    primers += [ins_f, ins_r]
    steps.append(Step(
        0, "PCR",
        f"Amplify {analysis.insert.label} ({len(insert_seq):,} bp) from "
        f"{don['name']} with {ins_f.name} / {ins_r.name}",
        inputs=[don["name"], f"{ins_f.name}+{ins_r.name}"],
        output=f"{analysis.insert.label} amplicon",
        details={"length": len(insert_seq), "template": don["name"],
                 "gc": round(gc_fraction(insert_seq), 3)}))

    if len(backbone_seq) > cfg.max_pcr_length:
        route.warnings.append(
            f"the inverse-PCR product is {len(backbone_seq):,} bp, above the "
            f"{cfg.max_pcr_length:,} bp limit; the restriction swap avoids "
            f"amplifying the backbone at all")

    dpni = dpni_decision(conn, [analysis.recipient.plasmid_id,
                                analysis.donor.plasmid_id],
                         analysis.recipient.plasmid_id, gel_purified=False)
    if dpni["needed"]:
        steps.append(Step(
            0, "DPNI_DIGEST",
            f"Digest both PCR reactions with DpnI to remove the parental "
            f"{rec['name']} and {don['name']} templates",
            output="template-free amplicons"))
    else:
        route.warnings.append("DpnI is not needed here: " + dpni["reason"])
    steps.append(Step(0, "DNA_CLEANUP",
                      "Column purify both PCR products",
                      output="clean amplicons"))
    steps.append(Step(
        0, "HI_FI_ASSEMBLY",
        f"Assemble the linearised backbone with {analysis.insert.label} "
        f"using {cfg.preferred_assembly_kit} ({cfg.overlap_length} bp overlaps)",
        inputs=["linearised backbone", f"{analysis.insert.label} amplicon"],
        output="assembled plasmid",
        details={"fragments": 2, "overlap": cfg.overlap_length}))

    for name, row in ((rec["name"], rec), (don["name"], don)):
        if row["verified"] == "unverified" and not cfg.allow_unverified_templates:
            route.warnings.append(f"parent {name} is not sequence-verified")

    _append_swap_finishing(steps, cfg, analysis, conn)
    route.steps = renumber(steps)
    route.primers = primers
    for p in primers:
        route.warnings += [f"primer {p.name}: {w}" for w in p.warnings]
    route.warnings += analysis.notes
    return route


def plan_swap_restriction_route(conn: sqlite3.Connection, analysis, cfg: Config,
                                scheme: dict) -> Route:
    """Cut the region out of the recipient with enzymes and ligate the insert."""
    route = Route(id="swap_re", strategy="swap_restriction",
                  title="Swap by restriction digest + ligation (2 parents)",
                  scheme=scheme)
    if analysis is None or not analysis.feasible:
        route.feasible = False
        route.warnings.append(
            (analysis.notes[0] if analysis and analysis.notes
             else "no suitable pair of parent plasmids was found"))
        return route

    rec = conn.execute("SELECT seq, circular, name FROM plasmid WHERE id = ?",
                       (analysis.recipient.plasmid_id,)).fetchone()
    don = conn.execute("SELECT seq, circular, name FROM plasmid WHERE id = ?",
                       (analysis.donor.plasmid_id,)).fetchone()
    rec_seq, rec_circ = rec["seq"], bool(rec["circular"])
    remove = analysis.remove
    insert_seq = _insert_sequence(conn, analysis)

    left_bound = analysis.left_anchor["end"] if analysis.left_anchor else 0
    right_bound = (analysis.right_anchor["start"] if analysis.right_anchor
                   else len(rec_seq))

    pair = _find_swap_enzymes(rec_seq, rec_circ, insert_seq, left_bound,
                              right_bound, remove, cfg)
    if not pair:
        route.feasible = False
        route.warnings.append(
            f"no stocked enzyme pair cuts {rec['name']} only within the "
            f"region being replaced while leaving {analysis.insert.label} "
            f"intact; use the PCR swap, or add sites by PCR")
        return route

    e1, e2 = pair["enzymes"]
    steps: list[Step] = []
    primers: list[Primer] = []

    steps.append(Step(
        0, "RE_DIGEST",
        f"Digest {rec['name']} with {e1} + {e2} to drop out "
        f"{remove.label} ({pair['dropout']:,} bp)",
        inputs=[rec["name"], f"{e1}, {e2}"], output="cut backbone",
        details={"enzymes": [e1, e2]}))

    donor_native = pair.get("donor_native")
    if donor_native:
        steps.append(Step(
            0, "RE_DIGEST",
            f"Digest {don['name']} with {e1} + {e2} to release "
            f"{analysis.insert.label}",
            inputs=[don["name"], f"{e1}, {e2}"], output="released insert",
            details={"enzymes": [e1, e2]}))
    else:
        f, r = design_re_pair(analysis.insert.label or "insert", don["seq"],
                              analysis.insert.start, analysis.insert.end,
                              e1, e2, cfg,
                              circular_template=bool(don["circular"]))
        primers += [f, r]
        steps.append(Step(
            0, "PCR",
            f"Amplify {analysis.insert.label} from {don['name']} with "
            f"{f.name} / {r.name}, adding {e1} and {e2} sites",
            inputs=[don["name"], f"{f.name}+{r.name}"],
            output=f"{analysis.insert.label} amplicon",
            details={"length": len(insert_seq), "template": don["name"]}))
        steps.append(Step(0, "DNA_CLEANUP", "Column purify the PCR product",
                          output="clean amplicon"))
        steps.append(Step(
            0, "RE_DIGEST", f"Digest the amplicon with {e1} + {e2}",
            output="digested insert", details={"enzymes": [e1, e2]}))

    if e1 == e2:
        steps.append(Step(
            0, "DEPHOSPHORYLATION",
            "Dephosphorylate the cut backbone to suppress self-ligation "
            "(both ends were cut with the same enzyme)",
            output="dephosphorylated backbone"))

    steps.append(Step(
        0, "GEL_EXTRACTION",
        f"Gel purify the {pair['backbone_band']:,} bp backbone band, away "
        f"from the {pair['dropout']:,} bp dropout",
        output="purified backbone", details={"size": pair["backbone_band"]}))
    steps.append(Step(
        0, "GEL_EXTRACTION",
        f"Gel purify the {len(insert_seq):,} bp {analysis.insert.label} band",
        output="purified insert", details={"size": len(insert_seq)}))
    steps.append(Step(
        0, "LIGATION",
        f"Ligate {analysis.insert.label} into the cut backbone "
        f"({pair['left_end']} / {pair['right_end']} ends), 3:1 insert:vector",
        inputs=["purified backbone", "purified insert"],
        output="ligated plasmid"))

    route.warnings += pair["warnings"]
    _append_swap_finishing(steps, cfg, analysis, conn)
    route.steps = renumber(steps)
    route.primers = primers
    for p in primers:
        route.warnings += [f"primer {p.name}: {w}" for w in p.warnings]
    route.warnings += analysis.notes
    route.cost_breakdown["enzymes"] = [e1, e2]
    return route


#: How far a cut may fall into the region being replaced. An enzyme cuts
#: a base or two inside its site, and those bases are discarded with the
#: rest of the region, so insisting the cut land exactly on the boundary
#: rules out the enzymes that sit on it.
SWAP_CUT_SLACK = 60


def _find_swap_enzymes(rec_seq: str, rec_circ: bool, insert_seq: str,
                       left_bound: int, right_bound: int, remove,
                       cfg: Config,
                       insert_must_lack_sites: bool = True) -> dict | None:
    """Find two enzymes that cut at the boundaries of the region being replaced.

    Both cuts must lie between the conserved features either side, and
    neither enzyme may cut the incoming insert, or the piece would be
    chopped up in the same reaction.

    The boundary is not a line, though. An enzyme's cut sits a base or
    two inside its recognition site, so a site that begins exactly where
    the region does cuts one base into it - and demanding that the cut
    land before the region starts rejected BamHI at the very boundary of
    a Voltron-2 cassette, which with EcoRI at the far end would have
    dropped the whole thing out in one digest. The bases between the
    boundary and the cut are being discarded anyway, so a cut is allowed
    to fall a little way into the region.
    """
    names = [n for n in cfg.stocked_enzymes if n in CATALOG]

    # A usable cut sits between a conserved feature and the region being
    # replaced. The left cut must fall before the region starts and the right
    # cut after it ends: cutting *inside* the region would leave a stub of the
    # old sequence behind, which is a silent way to build the wrong plasmid.
    left_usable: list[tuple[str, int]] = []
    right_usable: list[tuple[str, int]] = []
    for n in names:
        enz = CATALOG[n]
        sites = find_sites(rec_seq, enz, rec_circ)
        if len(sites) != 1:
            continue                     # must be unique in the recipient
        cut = sites[0].top_cut
        if insert_must_lack_sites and find_sites(insert_seq, enz, False):
            # Only matters when the insert is itself digested. If it arrives
            # as a PCR product joined by homology, the enzyme never touches
            # it and a site inside it is irrelevant.
            continue
        if left_bound <= cut <= remove.start + SWAP_CUT_SLACK:
            left_usable.append((n, cut))
        elif remove.end - SWAP_CUT_SLACK <= cut <= right_bound:
            right_usable.append((n, cut))

    if not (left_usable and right_usable):
        return None

    best = None
    for a, cut_a in left_usable:
        for b, cut_b in right_usable:
            lo, hi = sorted((cut_a, cut_b))
            dropout = hi - lo
            if dropout < remove.length - 2 * SWAP_CUT_SLACK:
                continue                 # must excise essentially the region
            if a != b and dropout < 50:
                continue                 # ends too close to resolve on a gel
            warnings: list[str] = []
            # Whatever the cuts leave behind of the old sequence, said out
            # loud: it lands between the promoter and the new insert, or
            # just after it, and the person at the bench should know it is
            # there rather than discover it in the sequencing.
            left_stub = max(0, lo - remove.start)
            right_stub = max(0, remove.end - hi)
            if left_stub:
                warnings.append(
                    f"the cut falls {left_stub} bp inside the region being "
                    f"removed, so that much of the old sequence stays "
                    f"upstream of the insert as a spacer")
            if right_stub:
                warnings.append(
                    f"{right_stub} bp of the old sequence stays downstream "
                    f"of the insert; it is after the stop codon, but check "
                    f"it carries nothing you mind keeping")
            enz_a, enz_b = CATALOG[a], CATALOG[b]
            if a == b:
                warnings.append(
                    f"only {a} is usable here, so both ends are compatible and "
                    f"the backbone can self-ligate; dephosphorylate it")
            for name in {a, b}:
                warnings += methylation_warnings(rec_seq, CATALOG[name], rec_circ)
                if CATALOG[name].star_activity:
                    warnings.append(f"{name} shows star activity; keep the "
                                    f"digest short and avoid enzyme excess")
            if enz_a.buffer != enz_b.buffer:
                warnings.append(f"{a} and {b} prefer different buffers; "
                                f"digest sequentially")
            backbone_band = len(rec_seq) - dropout
            if abs(backbone_band - dropout) < 200:
                warnings.append("backbone and dropout bands are too close in "
                                "size to separate cleanly on a gel")
            cand = {
                "enzymes": (a, b), "warnings": warnings,
                "cut_left": lo, "cut_right": hi,
                "dropout": dropout or remove.length,
                "backbone_band": backbone_band,
                "left_end": enz_a.end_type, "right_end": enz_b.end_type,
                "donor_native": False,
                "score": len(warnings) + (3 if a == b else 0),
            }
            if best is None or cand["score"] < best["score"]:
                best = cand
    return best


def _append_swap_finishing(steps: list[Step], cfg: Config, analysis,
                           conn: sqlite3.Connection) -> None:
    if cfg.include_transformation_steps:
        steps.append(Step(0, "TRANSFORMATION",
                          "Transform into competent cells and plate on the "
                          "antibiotic the recipient backbone carries",
                          output="colonies"))
        steps.append(Step(
            0, "COLONY_SCREEN",
            f"Screen 4-8 colonies by colony PCR across the new junction "
            f"(a primer in {analysis.insert.label} paired with one in the "
            f"retained backbone distinguishes insert from empty vector)",
            output="candidate clones"))
    if cfg.include_sequencing_steps:
        steps.append(Step(
            0, "SEQUENCE_VERIFY",
            f"Sanger sequence both new junctions and the whole "
            f"{analysis.insert.label} insert",
            output="verified plasmid"))


def plan_swap_hybrid_route(conn: sqlite3.Connection, analysis, cfg: Config,
                           scheme: dict) -> Route:
    """Cut the recipient open, amplify only the insert, assemble.

    The all-PCR swap amplifies the backbone as well, which on a 7 kb plasmid
    means a long inverse PCR whose every base carries polymerase error, and
    a second pair of primers to pay for. If enzymes can open the recipient
    at the right place there is no reason to copy it at all: digest it, gel
    purify the vector, and run one PCR for the piece going in.

    One PCR instead of two, one primer pair instead of two, and the vector
    arrives as original plasmid DNA rather than an amplicon.
    """
    route = Route(id="swap_hybrid", strategy="swap_hybrid",
                  title="Swap by digest + one PCR + HiFi assembly",
                  scheme=scheme)
    if analysis is None or not analysis.feasible:
        route.feasible = False
        route.warnings.append(
            (analysis.notes[0] if analysis and analysis.notes
             else "no suitable pair of parent plasmids was found"))
        return route

    rec = conn.execute(
        "SELECT seq, circular, name, verified FROM plasmid WHERE id = ?",
        (analysis.recipient.plasmid_id,)).fetchone()
    don = conn.execute(
        "SELECT seq, circular, name, verified FROM plasmid WHERE id = ?",
        (analysis.donor.plasmid_id,)).fetchone()
    if rec is None or don is None:
        route.feasible = False
        route.warnings.append("parent plasmid sequence missing")
        return route

    rec_seq, rec_circ = rec["seq"], bool(rec["circular"])
    remove = analysis.remove
    insert_seq = _insert_sequence(conn, analysis)
    left_bound = analysis.left_anchor["end"] if analysis.left_anchor else 0
    right_bound = (analysis.right_anchor["start"] if analysis.right_anchor
                   else len(rec_seq))

    pair = _find_swap_enzymes(rec_seq, rec_circ, insert_seq, left_bound,
                              right_bound, remove, cfg,
                              insert_must_lack_sites=False)
    if not pair:
        route.feasible = False
        gap_left = remove.start - left_bound
        gap_right = right_bound - remove.end
        detail = ""
        if gap_left <= 0 or gap_right <= 0:
            side = "upstream" if gap_left <= 0 else "downstream"
            detail = (f" - the part being kept {side} begins exactly where "
                      f"the region to remove ends, leaving nowhere to cut "
                      f"without cutting into it")
        route.warnings.append(
            f"no stocked enzyme opens {rec['name']} at both ends of the "
            f"region to remove{detail}. Deleting a region by assembly needs "
            f"it physically cut out at both boundaries, so this edit has to "
            f"amplify around it instead: use the all-PCR swap")
        return route

    e1, e2 = pair["enzymes"]
    steps: list[Step] = []

    steps.append(Step(
        0, "RE_DIGEST",
        f"Digest {rec['name']} with {e1} + {e2} to open it and drop out "
        f"{remove.label} ({pair['dropout']:,} bp)",
        inputs=[rec["name"], f"{e1}, {e2}"], output="linearised vector",
        details={"enzymes": [e1, e2]}))
    if e1 == e2:
        steps.append(Step(
            0, "DEPHOSPHORYLATION",
            "Dephosphorylate the cut vector; both ends came from the same "
            "enzyme and would otherwise re-close on themselves",
            output="dephosphorylated vector"))
    steps.append(Step(
        0, "GEL_EXTRACTION",
        f"Gel purify the {pair['backbone_band']:,} bp vector band, away from "
        f"the {pair['dropout']:,} bp dropout",
        output="purified vector", details={"size": pair["backbone_band"]}))

    # The single PCR: the insert, with homology to the cut vector ends.
    cut_lo = min(pair.get("cut_left", remove.start), pair.get("cut_right", remove.end))
    cut_hi = max(pair.get("cut_left", remove.start), pair.get("cut_right", remove.end))
    upstream = slice_region(rec_seq, max(cut_lo - 60, 0), cut_lo, circular=rec_circ)
    downstream = slice_region(rec_seq, cut_hi, cut_hi + 60, circular=rec_circ)

    fwd, rev = design_gibson_pair(
        0, analysis.insert.label or "insert", don["seq"],
        analysis.insert.start, analysis.insert.end, upstream, downstream,
        cfg, circular_template=bool(don["circular"]))
    steps.append(Step(
        0, "PCR",
        f"Amplify {analysis.insert.label} ({len(insert_seq):,} bp) from "
        f"{don['name']} with {fwd.name} / {rev.name}, adding homology to the "
        f"cut vector ends",
        inputs=[don["name"], f"{fwd.name}+{rev.name}"],
        output=f"{analysis.insert.label} amplicon",
        details={"length": len(insert_seq), "template": don["name"],
                 "gc": round(gc_fraction(insert_seq), 3)}))
    dpni = dpni_decision(conn, [analysis.donor.plasmid_id],
                         analysis.recipient.plasmid_id, gel_purified=False)
    if dpni["needed"]:
        steps.append(Step(
            0, "DPNI_DIGEST",
            f"Digest the PCR with DpnI to remove the parental {don['name']}",
            output="template-free amplicon"))
    else:
        route.warnings.append("DpnI is not needed here: " + dpni["reason"])
    steps.append(Step(0, "DNA_CLEANUP", "Column purify the PCR product",
                      output="clean amplicon"))
    steps.append(Step(
        0, "HI_FI_ASSEMBLY",
        f"Assemble the cut vector with {analysis.insert.label} using "
        f"{cfg.preferred_assembly_kit}",
        inputs=["purified vector", f"{analysis.insert.label} amplicon"],
        output="assembled plasmid", details={"fragments": 2}))

    route.warnings += pair["warnings"]
    _append_swap_finishing(steps, cfg, analysis, conn)
    route.steps = renumber(steps)
    route.primers = [fwd, rev]
    for p in route.primers:
        route.warnings += [f"primer {p.name}: {w}" for w in p.warnings]
    route.warnings += analysis.notes
    route.cost_breakdown["enzymes"] = [e1, e2]
    return route


# --------------------------------------------------------------------------
# When DpnI is actually needed
# --------------------------------------------------------------------------

#: Resistance markers, normalised. Two templates that select differently can
#: be told apart by the plate, which is what makes DpnI skippable.
#: Matched as substrings of the lower-cased label, longest first, so that
#: "cmr" is recognised as well as "chloramphenicol" and "neoR/KanR" resolves
#: to KanR rather than stopping at the first thing that looks close.
_MARKERS = (("chloramphenicol", "CmR"), ("ampicillin", "AmpR"),
            ("kanamycin", "KanR"), ("spectinomycin", "SpecR"),
            ("streptomycin", "SmR"), ("tetracycline", "TetR"),
            ("blasticidin", "BsdR"), ("puromycin", "PuroR"),
            ("hygromycin", "HygR"), ("gentamicin", "GentR"),
            ("zeocin", "ZeoR"), ("carbenicillin", "AmpR"),
            ("amp", "AmpR"), ("bla", "AmpR"), ("carb", "AmpR"),
            ("kan", "KanR"), ("neo", "KanR"),
            ("chlor", "CmR"), ("cmr", "CmR"), ("cam", "CmR"),
            ("spec", "SpecR"), ("strep", "SmR"), ("smr", "SmR"),
            ("tet", "TetR"), ("zeo", "ZeoR"),
            ("blast", "BsdR"), ("bsd", "BsdR"),
            ("puro", "PuroR"), ("hygro", "HygR"), ("hygr", "HygR"),
            ("gent", "GentR"), ("aacc1", "GentR"))


def resistance_marker(conn: sqlite3.Connection, plasmid_id: int) -> str | None:
    """Which antibiotic a plasmid is selected with, as far as it is annotated."""
    rows = conn.execute(
        """SELECT name FROM annotation
           WHERE plasmid_id = ? AND feature_type = 'resistance'""",
        (plasmid_id,)).fetchall()
    for row in rows:
        label = (row["name"] or "").lower()
        for needle, marker in _MARKERS:
            if needle in label:
                return marker
    return None


def dpni_decision(conn: sqlite3.Connection, template_ids: list,
                  product_id: int | None, gel_purified: bool) -> dict:
    """Whether a DpnI step earns its place.

    DpnI cuts dam-methylated DNA, so its only job is destroying the circular
    parental plasmid that came in with the template and would otherwise
    transform directly and grow on the plate. That is a question about the
    template, not about which assembly chemistry follows: Gibson and HiFi
    are no more tolerant of a surviving parental plasmid than ligation is.

    It can be dropped when the parent cannot produce a colony anyway:

    * the amplicon was gel purified, which physically separates it from the
      template;
    * the template selects on a different antibiotic from the product, so
      the plate removes it;
    * the template was not a plasmid grown in a dam+ strain, in which case
      DpnI has nothing to cut.
    """
    reasons: list[str] = []
    if gel_purified:
        return {"needed": False,
                "reason": "the amplicon is gel purified, which separates it "
                          "from the parental plasmid; DpnI would be "
                          "belt-and-braces"}

    product_marker = (resistance_marker(conn, product_id)
                      if product_id is not None else None)
    same_selection = []
    for tid in template_ids:
        marker = resistance_marker(conn, tid)
        if marker is None or product_marker is None or marker == product_marker:
            same_selection.append(tid)
        else:
            reasons.append(f"one template selects on {marker} and the product "
                           f"on {product_marker}")

    if not same_selection:
        return {"needed": False,
                "reason": "every template selects on a different antibiotic "
                          "from the product, so the plate removes them: "
                          + "; ".join(dict.fromkeys(reasons))}

    return {"needed": True,
            "reason": "a template carries the same selection marker as the "
                      "product, so parental plasmid surviving the reaction "
                      "would grow on the plate and look like a clone"}
