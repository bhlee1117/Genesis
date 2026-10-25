"""The cloning figure's data model: one lane per piece that goes in the tube.

A route is a list of bench operations. A figure of that route is a list of
*pieces*, each with the plasmid it comes out of, the stages it passes
through, and the ends it carries into the join. Turning the first into the
second is the only thinking the drawing needs, and both the web app and the
protocol PDF need exactly the same answer - so it lives here rather than
twice, in JavaScript and in Python, drifting apart.

``web/app.js`` holds the JavaScript port. ``tests`` runs them against the
same plans and compares, because two implementations of one idea are worth
nothing if nobody checks they agree.
"""

from __future__ import annotations

import re

#: One colour per lane, carried from the parent circle through to the arc
#: it becomes in the product, so "where did this part come from" is
#: answered by looking rather than by reading a table.
LANE_COLOURS = [
    (0.18, 0.48, 0.83),
    (0.83, 0.34, 0.18),
    (0.18, 0.62, 0.37),
    (0.56, 0.31, 0.79),
    (0.71, 0.54, 0.00),
]

#: An enzyme list written into a step's inputs rather than its details,
#: as "HindIII, EcoRV".
_ENZYME_LIST = re.compile(r"^[A-Za-z]+[IVX0-9]+(\s*,\s*[A-Za-z]+[IVX0-9]+)*$")

_AMPLICON_SUFFIX = re.compile(r"\s+(amplicon|fragment)(\s+with sites)?$", re.I)
_PREP_PREFIX = re.compile(r"^(purified|clean|digested)\s+", re.I)

_HOMOLOGY_OPS = ("HI_FI_ASSEMBLY", "GIBSON_ASSEMBLY")


def piece_label(output: str) -> str:
    """``"EGFP amplicon with sites"`` -> ``"EGFP"``.

    What a step calls its output is bookkeeping - "clean amplicon",
    "purified insert" - and the figure wants the name of the DNA.
    """
    name = str(output or "fragment")
    return _PREP_PREFIX.sub("", _AMPLICON_SUFFIX.sub("", name))


def step_enzymes(step: dict) -> list[str]:
    """The enzymes one digest step names, wherever it put them."""
    enzymes = list((step.get("details") or {}).get("enzymes") or [])
    if enzymes:
        return enzymes
    for item in step.get("inputs") or []:
        if _ENZYME_LIST.match(str(item)):
            enzymes = [x.strip() for x in str(item).split(",")]
    return enzymes


def route_enzymes(route: dict) -> list[str]:
    """Every enzyme this route mentions anywhere.

    A restriction route cuts insert and vector with the same pair but
    only spells it out on one of the two steps, so the other was drawn
    as a bare "Digest" with nothing to say which enzymes.
    """
    seen: list[str] = []
    for step in route.get("steps") or []:
        for enzyme in step_enzymes(step):
            if enzyme not in seen:
                seen.append(enzyme)
    return seen


def step_primers(step: dict, route: dict) -> list[str]:
    """Split ``"Kozak+ASAP7y-F+Kozak+ASAP7y-R"`` back into two primers.

    Splitting on the first ``+`` is wrong whenever a part's own name
    contains one, which the fused names routinely do.
    """
    names = [p.get("name") for p in route.get("primers") or []]
    pair = next((i for i in step.get("inputs") or [] if "+" in str(i)), "")
    for a in names:
        for b in names:
            if a is not None and b is not None and a != b \
                    and f"{a}+{b}" == pair:
                return [a, b]
    m = re.match(r"^(.*-F)\+(.*-R)$", str(pair))
    return [m.group(1), m.group(2)] if m else []


def primer_span(route: dict, names: list[str]) -> dict | None:
    """Where a pair sits on its template, when the plan printed it."""
    by = {p.get("name"): p for p in route.get("primers") or []}
    fwd = by.get(names[0]) if names else None
    rev = by.get(names[1]) if len(names) > 1 else None
    if not fwd or not isinstance(fwd.get("target_start"), int):
        return None
    end = fwd.get("target_end")
    if not isinstance(end, int) and rev:
        end = rev.get("target_end")
    return {"start": fwd["target_start"], "end": end}


def join_verb(route: dict) -> str:
    """The reaction that puts the pieces together."""
    ops = [s.get("operation") for s in route.get("steps") or []]
    if "HI_FI_ASSEMBLY" in ops:
        return "NEBuilder HiFi"
    if "GIBSON_ASSEMBLY" in ops:
        return "Gibson assembly"
    if "LIGATION" in ops:
        return "T4 ligase"
    return "Assemble"


def route_lanes(route: dict, sim: dict | None = None,
                plan: dict | None = None) -> list[dict]:
    """One lane per piece that enters the final reaction.

    Each lane is ``{origin, template, src, piece, primers, stages, ends,
    end_kind, colour}``. ``origin`` is ``"plasmid"`` for DNA that comes
    out of something on the shelf and ``"order"`` for DNA nobody has.
    ``stages`` is the sequence of bench operations that piece passes
    through, so a PCR followed by a digest is two stages and not one
    hyphenated verb.
    """
    sim = sim or {}
    plan = plan or {}
    pool = [s for s in (sim.get("sources") or []) if s]
    atlas = plan.get("templates") or {}

    def take_source(template):
        for i, s in enumerate(pool):
            ok = s.get("synthesised") if template is None \
                else s.get("template") == template
            if ok:
                return pool.pop(i)
        return None

    def stand_in(template, span):
        """A route the planner rejected is never simulated, but it still
        knows its plasmids and its primers - which is all a drawing
        needs. Refusing to draw it was how the restriction route came up
        blank on screen."""
        t = atlas.get(template) or {}
        start = span and span.get("start")
        end = span and span.get("end")
        bp = (end - start) if isinstance(start, int) and isinstance(end, int) \
            and end > start else 0
        return {"name": "", "template": t.get("name") or template,
                "template_id": t.get("id"), "template_length": t.get("bp"),
                "template_circular": t.get("circular", True),
                "features": t.get("features") or [],
                "from_start": start, "from_end": end,
                "bp": bp, "estimated": True}

    lanes: list[dict] = []
    last: dict | None = None
    for step in route.get("steps") or []:
        op = step.get("operation")
        details = step.get("details") or {}
        if op == "PCR":
            template = details.get("template") \
                or next(iter(step.get("inputs") or []), "")
            names = step_primers(step, route)
            src = take_source(template) \
                or stand_in(template, primer_span(route, names))
            if not src.get("bp"):
                src["bp"] = details.get("length") or src.get("bp") or 0
            last = {"origin": "plasmid", "template": template, "src": src,
                    "piece": piece_label(step.get("output")),
                    "primers": names, "colour": None,
                    "stages": [{"verb": "PCR", "sub": _bp(src.get("bp"))}]}
            lanes.append(last)
        elif op == "RE_DIGEST":
            enzymes = step_enzymes(step)
            host = next((i for i in step.get("inputs") or []
                         if i not in enzymes and i != ", ".join(enzymes)),
                        None)
            if not enzymes:
                enzymes = route_enzymes(route)
            if host:
                src = take_source(host) or stand_in(host, None)
                last = {"origin": "plasmid", "template": host, "src": src,
                        "piece": piece_label(step.get("output")),
                        "primers": [], "colour": None, "cuts_parent": enzymes,
                        "stages": [{"verb": "Digest",
                                    "sub": " + ".join(enzymes)}]}
                lanes.append(last)
            elif last is not None:
                # A bare digest cuts whatever came off the bench last -
                # almost always the amplicon, which is how a hybrid route
                # puts sticky ends on a PCR product. That is a second
                # stage, not a footnote.
                last["stages"].append({"verb": "Digest",
                                       "sub": " + ".join(enzymes)})
                last["cuts_fragment"] = enzymes
        elif op == "SYNTHESIS":
            src = take_source(None) or {"bp": details.get("length") or 0,
                                        "synthesised": True}
            last = {"origin": "order", "template": None, "src": src,
                    "piece": piece_label(step.get("output")),
                    "primers": [], "colour": None,
                    "stages": [{"verb": "Order", "sub": _bp(src.get("bp"))}]}
            lanes.append(last)

    ops = [s.get("operation") for s in route.get("steps") or []]
    homology = any(op in _HOMOLOGY_OPS for op in ops)
    for i, lane in enumerate(lanes):
        lane["colour"] = LANE_COLOURS[i % len(LANE_COLOURS)]
        enzymes = lane.get("cuts_fragment") or lane.get("cuts_parent")
        if enzymes:
            lane["ends"] = [enzymes[0],
                            enzymes[1] if len(enzymes) > 1 else enzymes[0]]
            lane["end_kind"] = "sticky"
        elif homology:
            lane["ends"] = ["overlap", "overlap"]
            lane["end_kind"] = "homology"
        else:
            lane["ends"] = []
            lane["end_kind"] = ""
        if not lane.get("src"):
            lane["src"] = {"bp": 0, "template": lane["template"]}
    return lanes


def _bp(n) -> str:
    return f"{n:,} bp" if isinstance(n, int) and n else "? bp"


#: The feature types worth naming on a parent in a cloning figure. A full
#: plasmid map draws every annotation the file happens to carry; the
#: question here is narrower.
FIGURE_TYPES = ("CDS", "promoter", "resistance", "origin")


def outermost(items: list, start="start", end="end") -> list:
    """Drop anything that sits entirely inside something else kept.

    Real inventory files carry nested annotations - the chicken beta-actin
    promoter inside CAG - and auto-annotation leaves junk inside them
    too: YD030 has a 609 bp reverse-strand "Voltron2" buried in its CAG
    promoter that has nothing to do with the real Voltron-2 CDS further
    round the plasmid. At figure scale a feature inside another feature
    is not a part of the construct, it is a detail of one; whether it is
    genuine or junk, promoting it over its container misdescribes the
    molecule. The full map in the protocol still shows everything.
    """
    kept = []
    for a in items:
        s0, e0 = a.get(start), a.get(end)
        if not isinstance(s0, int) or not isinstance(e0, int) or e0 <= s0:
            kept.append(a)
            continue
        inside = any(
            b is not a and isinstance(b.get(start), int)
            and isinstance(b.get(end), int) and b[end] > b[start]
            and b[start] <= s0 and e0 <= b[end]
            and (b[end] - b[start]) > (e0 - s0)
            for b in items)
        if not inside:
            kept.append(a)
    return kept


def figure_features(src: dict, limit: int = 3) -> list[dict]:
    """The handful of annotations that place a parent, biggest first."""
    span = src.get("template_length") or 1
    keep = [a for a in (src.get("features") or [])
            if isinstance(a.get("start"), int)
            and a.get("end", 0) > a["start"]
            and (a["end"] - a["start"]) > span * 0.05
            and a.get("feature_type") in FIGURE_TYPES]
    keep = outermost(keep)
    keep.sort(key=lambda a: a["end"] - a["start"], reverse=True)
    return keep[:limit]


#: What the product circle names. A construct has one promoter and one
#: payload worth labelling at this size; the full map in the protocol
#: carries the rest.
PRODUCT_TYPES = ("promoter", "CDS")


def product_features(layout: list, span: int, limit: int = 4) -> list:
    """The features worth naming on the assembled plasmid."""
    span = span or 1
    keep = [a for a in layout or []
            if isinstance(a.get("start"), int)
            and a.get("type") in PRODUCT_TYPES
            and (a.get("end", 0) - a["start"]) > span * 0.03]
    keep = outermost(keep)
    keep.sort(key=lambda a: a["end"] - a["start"], reverse=True)
    return keep[:limit]


#: One colour per junction, so the two halves of a join match by eye.
#: Deliberately away from LANE_COLOURS: a lane colour says "this piece",
#: a junction colour says "these two ends meet".
JUNCTION_COLOURS = [
    (0.85, 0.65, 0.05),
    (0.74, 0.20, 0.60),
    (0.05, 0.62, 0.66),
    (0.45, 0.45, 0.48),
]
JUNCTION_LETTERS = "ABCDEF"


def junctions(lanes: list) -> list:
    """Which end meets which, as (lane index, side, junction id) triples.

    Sticky ends pair by enzyme - a HindIII end only goes into a HindIII
    end - so the enzyme name is the junction's identity. Homology arms
    have no name of their own, so they are numbered by the order the
    pieces sit in the product: piece k's right end meets piece k+1's
    left end, all the way round the circle.
    """
    out = []
    sticky = [l for l in lanes if l.get("end_kind") == "sticky"]
    if sticky and len(sticky) == len(lanes):
        for i, lane in enumerate(lanes):
            for side, name in zip(("left", "right"), lane.get("ends") or []):
                out.append((i, side, name))
        return out

    placed = [(i, l) for i, l in enumerate(lanes)
              if isinstance((l.get("src") or {}).get("at_start"), int)]
    if len(placed) != len(lanes) or len(lanes) < 2:
        return out
    placed.sort(key=lambda t: t[1]["src"]["at_start"])
    n = len(placed)
    for k in range(n):
        here = placed[k][0]
        nxt = placed[(k + 1) % n][0]
        label = JUNCTION_LETTERS[k % len(JUNCTION_LETTERS)]
        out.append((here, "right", label))
        out.append((nxt, "left", label))
    return out


def junction_colours(lanes: list) -> dict:
    """{(lane index, side): (label, colour)} for every end that has a partner."""
    ids = junctions(lanes)
    order: list = []
    for _i, _side, name in ids:
        if name not in order:
            order.append(name)
    return {(i, side): (name,
                        JUNCTION_COLOURS[order.index(name)
                                         % len(JUNCTION_COLOURS)])
            for i, side, name in ids}
