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

    for rm, item in zip(placed, layout_items):
        source = rm.sources[0] if rm.sources else None
        if source is None:
            fragments.append(PlannedFragment(
                rm.part_name, rm.seq, None, "(synthesis)", 0, len(rm.seq), 1,
                False, [item], needs_synthesis=True))
            current = None
            continue

        can_extend = (
            current is not None
            and current.source_plasmid_id == source["id"]
            and current.strand == source["strand"]
            and not current.needs_synthesis
            and (
                (source["strand"] == 1 and source["start"] == current.end)
                or (source["strand"] == -1 and source["end"] == current.start)
            )
            and (len(current.seq) + len(rm.seq)) <= cfg.max_pcr_length
        )
        if can_extend and current is not None:
            current.seq += rm.seq
            current.name += f"+{rm.part_name}"
            current.modules.append(item)
            if source["strand"] == 1:
                current.end = source["end"]
            else:
                current.start = source["start"]
        else:
            current = PlannedFragment(
                rm.part_name, rm.seq, source["id"], source["name"],
                source["start"], source["end"], source["strand"],
                bool(source.get("seq_length")), [item],
                verified=source.get("verified", "unverified"))
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

    pieces: list[tuple[str, str]] = []
    if backbone:
        pieces.append(("backbone", backbone["seq"]))
    for f in fragments:
        pieces.append((f.name, f.seq))
    if not pieces:
        route.feasible = False
        route.warnings.append("nothing to assemble")
        return route

    n = len(pieces)
    templates_needing_dpni: list[str] = []

    # --- backbone amplification ------------------------------------------
    if backbone:
        bb_seq = backbone["seq"]
        up = pieces[-1][1]
        down = pieces[1][1] if n > 1 else bb_seq
        fwd, rev = design_gibson_pair(
            0, "backbone", bb_seq, 0, len(bb_seq), up, down, cfg,
            circular_template=bool(backbone.get("circular")))
        primers += [fwd, rev]
        if len(bb_seq) > cfg.max_pcr_length:
            route.warnings.append(
                f"the backbone amplicon is {len(bb_seq)} bp, above the "
                f"{cfg.max_pcr_length} bp PCR limit in the configuration; "
                f"linearising the backbone with a restriction digest is "
                f"usually more reliable than amplifying it")
        steps.append(Step(
            0, "PCR",
            f"Amplify the {len(bb_seq)} bp backbone from {backbone['name']} "
            f"with primers {fwd.name} / {rev.name}",
            inputs=[backbone["name"], f"{fwd.name}+{rev.name}"],
            output="backbone amplicon",
            details={"length": len(bb_seq), "template": backbone["name"]}))
        templates_needing_dpni.append(backbone["name"])

    # --- insert amplification ---------------------------------------------
    for i, frag in enumerate(fragments):
        piece_index = i + (1 if backbone else 0)
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
        template_seq = row["seq"] if row else frag.seq
        circ = bool(row["circular"]) if row else False
        start = frag.start if frag.strand == 1 else frag.start
        end = frag.end

        fwd, rev = design_gibson_pair(
            i, frag.name, template_seq, start, end, up, down, cfg,
            circular_template=circ)
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
    add("synthesis", counts.get("SYNTHESIS", 0), "synthesis")

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
        insert_seq, cfg, circular_template=rec_circ)
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

    ins_f, ins_r = design_gibson_pair(
        1, analysis.insert.label or "insert", don["seq"],
        analysis.insert.start, analysis.insert.end, backbone_seq, backbone_seq,
        cfg, circular_template=bool(don["circular"]))
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

    steps.append(Step(
        0, "DPNI_DIGEST",
        f"Digest both PCR reactions with DpnI to remove the parental "
        f"{rec['name']} and {don['name']} templates",
        output="template-free amplicons"))
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


def _find_swap_enzymes(rec_seq: str, rec_circ: bool, insert_seq: str,
                       left_bound: int, right_bound: int, remove,
                       cfg: Config) -> dict | None:
    """Find two enzymes that cut only inside the region being replaced.

    Both cuts must lie between the conserved features either side, and
    neither enzyme may cut the incoming insert, or the piece would be
    chopped up in the same reaction.
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
        if find_sites(insert_seq, enz, False):
            continue                     # would cut the incoming piece
        if left_bound <= cut <= remove.start:
            left_usable.append((n, cut))
        elif remove.end <= cut <= right_bound:
            right_usable.append((n, cut))

    if not (left_usable and right_usable):
        return None

    best = None
    for a, cut_a in left_usable:
        for b, cut_b in right_usable:
            lo, hi = sorted((cut_a, cut_b))
            dropout = hi - lo
            if dropout < remove.length:
                continue                 # must excise the whole region
            if a != b and dropout < 50:
                continue                 # ends too close to resolve on a gel
            warnings: list[str] = []
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
