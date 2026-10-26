"""Bench-ready protocol generation: concise steps up front, detail at the end.

The protocol is built once as a structured object and then rendered to PDF
or Markdown, so the two never drift apart.

Layout, in the order a person actually uses it:

1. What you are building (one block of facts)
2. Reagents and consumables to have out
3. Primers to order
4. The numbered protocol, with reaction set-ups inline
5. Everything else - validation, junctions, sources, alternatives, sequence

Reaction set-ups are the standard supplier conditions for the kits named in
the configuration. They are a starting point to check against the current
package insert, not a substitute for it.
"""

from __future__ import annotations

import datetime as _dt
import re

from .config import Config
from .pdf import (PDF, A4, circular_map, cloning_figure,
                  figure_height, plasmid_track)
from .seqtools import gc_fraction, tm

# --------------------------------------------------------------------------
# Construct naming
# --------------------------------------------------------------------------

#: Module types that belong in a construct name, in the order people read them.
NAMED_TYPES = ("promoter", "enhancer", "Cre-lox", "DIO", "CDS",
               "trafficking_motif", "p2A", "IRES", "tag")

#: Types whose own name is the architecture, not a part label.
ARCHITECTURAL_TYPES = {"Cre-lox", "DIO"}

_DROP_TOKENS = {"cds", "sequence", "seq", "fragment", "insert", "gene",
                "new", "dna", "part"}

#: Canonical spellings for elements that lab files capitalise inconsistently.
#: A plasmid name outlives the file it came from, so it is worth getting the
#: casing right rather than inheriting whatever the annotation happened to use.
_CANONICAL = {
    "hsyn": "hSyn", "syn": "Syn", "camkii": "CaMKII", "camk2a": "CaMKIIa",
    "ef1a": "EF1a", "ef1alpha": "EF1a", "cag": "CAG", "cmv": "CMV",
    "pgk": "PGK", "tre": "TRE", "u6": "U6", "gfap": "GFAP", "mecp2": "MeCP2",
    "egfp": "EGFP", "eyfp": "EYFP", "ecfp": "ECFP", "gfp": "GFP",
    "mcherry": "mCherry", "morange2": "mOrange2", "mscarlet": "mScarlet",
    "tdtomato": "tdTomato", "citrine": "Citrine", "halotag": "HaloTag",
    "quasar2": "QuasAr2", "quasar3": "QuasAr3", "cheriff": "CheRiff",
    "chr2": "ChR2", "archon1": "Archon1", "archon": "Archon",
    "voltron": "Voltron", "voltron2": "Voltron2", "voltron-2": "Voltron2",
    "er2": "ER2", "kgc": "KGC", "ts": "TS", "p2a": "P2A", "t2a": "T2A",
    "e2a": "E2A", "f2a": "F2A", "ires": "IRES", "wpre": "WPRE",
    "dio": "DIO", "flex": "FLEX", "cre": "Cre", "flp": "Flp",
    "gcamp6f": "GCaMP6f", "gcamp7f": "GCaMP7f", "gcamp8f": "GCaMP8f",
    "gcamp8m": "GCaMP8m", "gcamp8s": "GCaMP8s", "jrgeco1a": "jRGECO1a",
    "iglusnfr3": "iGluSnFR3", "iglusnfr": "iGluSnFR",
}


def _canonical(token: str) -> str:
    return _CANONICAL.get(token.lower().replace("-", ""), token)


#: Descriptive words that belong in a part's label but not in a plasmid name.
_NOISE = (r"motif|signal|peptide|promoter|element|tag|export|sequence|site|"
          r"box|linker|region|domain|protein|terminal")


def _clean_token(name: str) -> str:
    """Turn a part label into something that reads well inside a name.

    Lab labels carry descriptions that a plasmid name does not want:
    "ER2 export motif" is just ER2, and "KGC / Kir2.1 TS" is just KGC.
    """
    token = name.split("/")[0]                   # "KGC / Kir2.1 TS" -> "KGC"
    token = re.sub(r"\(.*?\)", " ", token)       # drop parentheticals
    token = re.sub(rf"\b({_NOISE})\b", " ", token, flags=re.I)
    # Collapse what is left; join any remaining words without spaces.
    words = [w for w in re.split(r"[\s_]+", token) if w]
    token = re.sub(r"[^A-Za-z0-9.+-]+", "", "".join(words))
    if token.lower() in _DROP_TOKENS:
        return ""
    return _canonical(token.strip("-."))


def suggest_name(plan: dict) -> str:
    """Propose a construct name from the resolved architecture.

    Follows the convention these plasmids are already named by: a vector
    prefix, then the expressed elements in order, joined with hyphens.
    """
    construct = plan.get("construct", {})
    layout = construct.get("layout", [])
    types = {item["type"] for item in layout}

    # Name the architecture the user asked for. After a DIO expansion the
    # assembled layout holds the payload inverted and in reverse order, which
    # is correct DNA but the wrong thing to name a plasmid after.
    source = construct.get("spec_order") or layout

    # Vector family, from the ITRs if present and otherwise from the backbone
    # name. No word boundaries here: these markers are normally glued to a
    # prefix, as in "pAAV_backbone" or "pLJM1".
    # The vector family decides the prefix. When one was asked for it wins
    # outright: if the user said this is an AAV, the name says pAAV.
    prefix = ""
    family = (plan.get("vector") or {}).get("family") or ""
    backbone = (plan.get("backbone") or {}).get("name", "") or ""
    if family == "AAV":
        prefix = "pAAV"
    elif family == "lentiviral":
        prefix = "pLenti"
    elif "AAV_ITR" in types or re.search(r"AAV", backbone, re.I):
        prefix = "pAAV"
    elif re.search(r"lenti|fugw|pLJM|pHR|FCK|FUW", backbone, re.I):
        prefix = "pLenti"

    tokens: list[str] = []
    for item in source:
        if item["type"] not in NAMED_TYPES:
            continue
        if item["type"] in ARCHITECTURAL_TYPES:
            # Name the arrangement people recognise - DIO, FLEX, FAS -
            # rather than the family it belongs to.
            tokens.append(item.get("mode") or item.get("name") or "DIO")
            continue
        token = _clean_token(item.get("name") or "")
        if token and (not tokens or tokens[-1].lower() != token.lower()):
            tokens.append(token)

    if not tokens:
        return construct.get("name", "") or "untitled"
    parts = ([prefix] if prefix else []) + tokens
    return "-".join(parts)[:120]


# --------------------------------------------------------------------------
# Reaction set-ups
# --------------------------------------------------------------------------

def _pcr_recipe(cfg: Config) -> dict:
    return {
        "title": "PCR reaction (50 uL, high-fidelity polymerase)",
        "items": [
            ("5x high-fidelity reaction buffer", "10 uL"),
            ("10 mM dNTPs", "1 uL"),
            ("10 uM forward primer", "2.5 uL"),
            ("10 uM reverse primer", "2.5 uL"),
            ("template plasmid (1-10 ng)", "1 uL"),
            ("high-fidelity DNA polymerase", "0.5 uL"),
            ("nuclease-free water", "to 50 uL"),
        ],
    }


RECIPES = {
    "RE_DIGEST": {
        "title": "Restriction digest (50 uL)",
        "items": [("DNA", "1 ug"), ("10x digestion buffer", "5 uL"),
                  ("each enzyme", "1 uL (10-20 U)"),
                  ("nuclease-free water", "to 50 uL")],
        "conditions": "37 C for 1 h. Verify complete linearisation on a gel "
                      "before proceeding.",
    },
    "DPNI_DIGEST": {
        "title": "DpnI treatment",
        "items": [("completed PCR reaction", "all"), ("DpnI", "1 uL")],
        "conditions": "37 C for 30-60 min, added directly to the PCR. DpnI "
                      "reads methylation rather than a rare sequence: it "
                      "cuts GATC only where the A is dam-methylated, which "
                      "is every GATC in a plasmid grown in a standard "
                      "(dam+) strain and none in a PCR product, built from "
                      "unmethylated dNTPs. Cutting the template at many "
                      "sites is how the step works.",
    },
    "DNA_CLEANUP": {
        "title": "Column or bead clean-up",
        "items": [("PCR or digest reaction", "all"),
                  ("spin column or SPRI beads", "per kit")],
        "conditions": "A silica column or SPRI beads, not a gel: it "
                      "removes enzymes, buffer, dNTPs and primers but "
                      "keeps every DNA species in the tube whatever its "
                      "size. Separating one band from another is the gel "
                      "extraction step, listed separately when the plan "
                      "needs it. Elute in 20-30 uL; quantify before "
                      "assembly.",
    },
    "GEL_EXTRACTION": {
        "title": "Gel extraction",
        "items": [("agarose gel", "0.8-1 % (1.5 % below 500 bp)"),
                  ("DNA stain", "per kit"), ("gel extraction kit", "per kit")],
        "conditions": "Run far enough to resolve the bands, excise the "
                      "correct band on a blue/long-wave transilluminator, "
                      "elute in 20-30 uL.",
    },
    "HI_FI_ASSEMBLY": {
        "title": "HiFi / Gibson assembly (20 uL)",
        "items": [("2x assembly master mix", "10 uL"),
                  ("linearised vector", "50-100 ng"),
                  ("each insert", "2:1 molar ratio over vector"),
                  ("nuclease-free water", "to 20 uL")],
        "conditions": "50 C for 15 min (2-3 fragments) or 60 min (4 or more). "
                      "Keep total fragment input near 0.03-0.2 pmol.",
    },
    "GIBSON_ASSEMBLY": {
        "title": "Gibson assembly (20 uL)",
        "items": [("2x Gibson master mix", "10 uL"),
                  ("linearised vector", "50-100 ng"),
                  ("each insert", "2:1 molar ratio over vector"),
                  ("nuclease-free water", "to 20 uL")],
        "conditions": "50 C for 15-60 min depending on fragment count.",
    },
    "LIGATION": {
        "title": "T4 ligation (20 uL)",
        "items": [("10x T4 ligase buffer (fresh ATP)", "2 uL"),
                  ("digested vector", "50 ng"),
                  ("insert", "3:1 molar ratio over vector"),
                  ("T4 DNA ligase", "1 uL"),
                  ("nuclease-free water", "to 20 uL")],
        "conditions": "Room temperature for 10 min (sticky ends) or 16 C "
                      "overnight. Include a vector-only no-insert control.",
    },
    "DEPHOSPHORYLATION": {
        "title": "Vector dephosphorylation",
        "items": [("digested vector", "all"),
                  ("phosphatase (rSAP or CIP)", "1 uL")],
        "conditions": "37 C for 30 min, then heat inactivate per the supplier. "
                      "Suppresses vector self-ligation.",
    },
    "TRANSFORMATION": {
        "title": "Transformation",
        "items": [("competent cells", "50 uL"),
                  ("assembly or ligation product", "2 uL"),
                  ("SOC medium", "950 uL")],
        "conditions": "Timings below are for the strain this lab stocks. "
                      "A construct with no repeats can go into NEB "
                      "5-alpha instead, recovered at 37 C in SOC.",
    },
    "COLONY_SCREEN": {
        "title": "Colony screen",
        "items": [("colonies", "4-8"), ("screening primer pair", "per reaction"),
                  ("Taq master mix", "per reaction")],
        "conditions": "Patch each colony before screening. Confirm candidates "
                      "by diagnostic digest and set up overnight cultures.",
    },
    "SEQUENCE_VERIFY": {
        "title": "Sequence verification",
        "items": [("miniprep DNA", "per facility spec"),
                  ("sequencing primers", "one per junction")],
        "conditions": "Sequence across every junction first, then the full "
                      "insert. Check the translated fusion, not only the DNA.",
    },
    "SYNTHESIS": {
        "title": "Synthetic fragment order",
        "items": [("fragment sequence with overlaps", "see detail section")],
        "conditions": "Check the vendor's complexity screen for repeats and "
                      "GC-rich regions before ordering.",
    },
}


def _pcr_conditions(step: dict, primers: list[dict], cfg: Config) -> str:
    """Annealing temperature and extension time for one PCR step."""
    length = step.get("details", {}).get("length", 0)
    name_hint = step.get("output", "").replace(" amplicon", "").strip()
    pair = [p for p in primers if p["template"] == name_hint] or \
           [p for p in primers if p["name"].startswith(name_hint)]
    tms = [p["anneal_tm"] for p in pair if p["anneal_tm"] == p["anneal_tm"]]

    bits = []
    if tms:
        ta = min(tms)
        bits.append(f"annealing {ta:.0f} C")
        if len(tms) == 2 and abs(tms[0] - tms[1]) > 5:
            bits.append(f"primer Tm differ by {abs(tms[0] - tms[1]):.0f} C, "
                        f"consider a gradient")
    if length:
        from .reagents import for_operation
        product = for_operation("PCR")
        rate = (product.extension_s_per_kb if product
                and product.extension_s_per_kb else 30)
        seconds = max(int(round(length / 1000 * rate)), 15)
        bits.append(f"extension {seconds // 60}:{seconds % 60:02d} for "
                    f"{length:,} bp at {rate} s/kb")
    gc = step.get("details", {}).get("gc")
    if gc and gc > 0.68:
        bits.append(f"{gc:.0%} GC, add the high-GC enhancer")
    bits.append("98 C 30 s initial, 30 cycles, 72 C 2 min final")
    return "; ".join(bits) + "."


# --------------------------------------------------------------------------
# Protocol assembly
# --------------------------------------------------------------------------

def build_protocol(plan: dict, cfg: Config | None = None,
                   route_id: str | None = None) -> dict:
    """Collect everything the renderers need, once."""
    cfg = cfg or Config()
    routes = plan.get("routes", [])
    if not routes:
        return {"empty": True, "construct": plan.get("construct", {})}
    route = next((r for r in routes if r["id"] == route_id), routes[0])
    # If nothing is feasible, prefer the best-ranked route anyway so its
    # reason is the one reported.
    usable = [r for r in routes if r["feasible"]]
    if route_id is None and usable:
        route = usable[0]

    construct = plan["construct"]
    primers = route["primers"]

    # --- reagents actually required by this route ----------------------
    used_ops = []
    for step in route["steps"]:
        if step["operation"] not in used_ops:
            used_ops.append(step["operation"])

    reagents: list[tuple[str, str]] = []
    seen: set[str] = set()

    def add_reagent(item: str, amount: str):
        if item.lower() not in seen:
            seen.add(item.lower())
            reagents.append((item, amount))

    # What to order is what the steps actually call for. The protocol
    # runs on master mixes, so listing "5x reaction buffer, 10 mM dNTPs,
    # high-fidelity DNA polymerase" as three line items is a shopping
    # list for a reaction nobody is setting up - those are inside the
    # Q5 tube. Where a product is named, its own components are the
    # reagents; the generic recipe is the fallback for the rest.
    from .reagents import for_operation as _product_for
    for op in used_ops:
        product = _product_for(op)
        if product is not None and product.setup:
            for amount, item in product.setup:
                add_reagent(item, amount)
            continue
        if op == "PCR":
            for item, amount in _pcr_recipe(cfg)["items"]:
                add_reagent(item, amount)
            continue
        for item, amount in RECIPES.get(op, {}).get("items", []):
            add_reagent(item, amount)

    # --- steps, with set-ups attached ----------------------------------
    steps = []
    for step in route["steps"]:
        op = step["operation"]
        detail = ""
        recipe = None
        if op == "PCR":
            detail = _pcr_conditions(step, primers, cfg)
            recipe = _pcr_recipe(cfg)
        else:
            r = RECIPES.get(op)
            if r:
                detail = r.get("conditions", "")
                recipe = r if r.get("items") else None
        # A step may carry a fact about this plan rather than about the
        # operation in general - how many sites DpnI has to cut here, for
        # instance - and that belongs under the step, not in a table.
        own = (step.get("details") or {}).get("note")
        if own:
            detail = (detail + " " + own).strip()
        # A generic recipe is right for every kit and exact for none.
        # Where the purchasing record names the product this bench runs,
        # that product's own published setup replaces the generic one.
        from .reagents import for_operation
        product = for_operation(op)
        if product is not None and product.setup:
            recipe = {"title": product.title,
                      "items": [(item, amount)
                                for amount, item in product.setup]}
        # A double digest is one tube at 37 C until the pair does not
        # share a temperature, and then it is two stages. Which it is
        # depends on the enzymes this route chose, not on the operation.
        digest_stages, red = [], []
        if op == "RE_DIGEST":
            from .enzymes import digest_plan
            from .figure import step_enzymes
            # Routes put the pair in details on some steps and in the
            # inputs on others; figure.step_enzymes already knows both.
            names = step_enzymes(step) or route_enzyme_names(route)
            plan_d = digest_plan(list(names))
            for i, stage in enumerate(plan_d["stages"], start=1):
                label = (f"stage {i}: {' + '.join(stage['enzymes'])}"
                         if plan_d["sequential"]
                         else " + ".join(stage["enzymes"]))
                digest_stages.append(
                    (label, f"{stage['temperature']} C, {stage['minutes']}"))
                if stage["heat_kill"]:
                    digest_stages.append(("stop", stage["heat_kill"]))
            red = list(plan_d["warnings"])
        steps.append({
            "index": step["index"], "operation": op,
            "description": step["description"], "detail": detail,
            "recipe": recipe, "notes": step.get("notes", ""),
            "product": product, "digest_stages": digest_stages,
            "red_notes": red,
            # The step's own findings - the checkpoint, the DpnI site
            # count - travel with it; the renderers read them from here.
            "details": step.get("details") or {},
            "optional": bool(step.get("optional")),
        })

    return {
        "empty": False,
        "blocked": not route.get("feasible", True),
        "construct": construct,
        "suggested_name": suggest_name(plan),
        "backbone": plan.get("backbone"),
        "route": route,
        "routes": routes,
        "steps": steps,
        "primers": primers,
        "reagents": reagents,
        "validation": plan.get("validation", {}),
        "modules": plan.get("modules", []),
        "fragments": plan.get("fragments", []),
        "generated": _dt.date.today().isoformat(),
    }


# --------------------------------------------------------------------------
# PDF rendering
# --------------------------------------------------------------------------

def _anneal_site(pr: dict) -> str:
    """Where a primer sits on its template, as a position and a strand.

    A primer sequence on its own does not say where it lands, and the one
    question somebody checks by hand before ordering is whether it sits
    inside the region being cut out.
    """
    start, end = pr.get("target_start"), pr.get("target_end")
    anneal = len(pr.get("anneal") or "")
    if start is None or end is None or start == end or not anneal:
        return "-"
    if pr.get("direction") == "F":
        lo, hi, arrow = start, start + anneal, "->"
    else:
        lo, hi, arrow = end - anneal, end, "<-"
    return f"{lo:,}-{hi:,} {arrow}"


def _primer_marks(route: dict, template_name: str) -> list[dict]:
    """Where this route's primers anneal on one named template.

    Drawn on the parent, not on the product, because that is the plasmid
    somebody is holding when they want to know whether a primer sits
    inside the region they are cutting out.
    """
    # Which fragment names are amplified off this plasmid. A primer only
    # belongs on the map of the plasmid it anneals to; drawing every pair
    # on every parent says the backbone primers bind the donor, which
    # would send somebody to the wrong tube.
    mine = set()
    for st in route.get("steps") or []:
        if st.get("operation") != "PCR":
            continue
        if (st.get("details") or {}).get("template") != template_name:
            continue
        for item in st.get("inputs") or []:
            if "+" in item:
                f, _, r = item.partition("+")
                mine.update((f.strip(), r.strip()))

    marks = []
    for pr in route.get("primers") or []:
        if mine and pr.get("name") not in mine:
            continue
        start, end = pr.get("target_start"), pr.get("target_end")
        if start is None or end is None or start == end:
            continue
        anneal = len(pr.get("anneal") or "")
        if not anneal:
            continue
        if pr.get("direction") == "F":
            lo, hi, strand = start, start + anneal, 1
        else:
            lo, hi, strand = end - anneal, end, -1
        marks.append({"name": pr.get("name", ""), "start": lo, "end": hi,
                      "strand": strand})
    return marks


#: Phrases that cost a week of cloning when they are skimmed past.
_RED = ("not 37 C", "do not heat inactivate", "Do not heat inactivate",
        "not to substitute SOC", "two stages", "column purify between",
        "star activity", "exactly 30 s", "cannot be heat inactivated")


def route_enzyme_names(route: dict) -> list:
    """Every enzyme the route mentions, for a step that names none."""
    from .figure import route_enzymes
    return route_enzymes(route)


def _vendor_box(doc, recipe: dict, product, step: dict) -> None:
    """One kit's own setup, ruled off and set in a serif.

    Reagent names in a monospace column had to be squeezed to fit; a
    serif at the same point size fits them and reads as what it is - a
    quotation of somebody else's protocol inside a document otherwise
    set in the interface's own face.
    """
    from .pdf import WARN_RED

    def body():
        doc.paragraph(recipe["title"], 9, font="F6", indent=4)
        if product is not None:
            doc.paragraph(product.label, 8, font="F7", indent=4, gray=0.4)
        width = max((len(a) for a, _ in
                     [(x[1], x[0]) for x in recipe["items"]]), default=0)
        for item, amount in recipe["items"]:
            doc.paragraph(f"{amount}   {item}", 8.8, indent=14, font="F5",
                          gray=0.1)
        stages = (step.get("digest_stages") or
                  (list(product.program) if product is not None else []))
        if stages:
            doc.spacer(2)
            for stage, condition in stages:
                doc.paragraph(f"{stage}   {condition}", 8.8, indent=14,
                              font="F5", gray=0.1)
        notes = list(step.get("red_notes") or ())
        rest = [n for n in (product.notes if product is not None else ())]
        for note in notes:
            doc.bullet(note, 8.4, indent=10, marker="!", font="F5",
                       rgb=WARN_RED)
        for note in rest:
            red = any(k in note for k in _RED)
            doc.bullet(note, 8.4, indent=10, font="F5",
                       rgb=WARN_RED if red else None,
                       gray=None if red else 0.35)
        if product is not None and product.url:
            doc.link(product.url, product.url, 7.6, indent=4, font="F7")

    doc.boxed(body, pad=7, fill=(0.985, 0.985, 0.99))


def _top_level_layout(construct: dict) -> list:
    """The construct's parts, without the ones nested inside others."""
    from .figure import outermost
    return outermost(list(construct.get("layout") or []))


def _checkpoint_box(doc, step: dict) -> None:
    """What a good result looks like, and what each bad one means.

    A protocol that says "run a gel" has given advice. One that says
    "one band at 5,791 bp, and a band at 6,539 bp is uncut template"
    has given a test - which is the difference between finding a failed
    PCR now and finding it on a plate in two days.
    """
    from .pdf import WARN_RED
    check = (step.get("details") or {}).get("checkpoint")
    if not check:
        return

    def body():
        doc.paragraph(check["title"], 9.4, font="F2")
        for line in check.get("how") or ():
            doc.paragraph(line, 8.8, indent=6, gray=0.25)
        doc.spacer(1)
        doc.paragraph("Go on if you see", 8.6, font="F2", indent=6,
                      gray=0.2)
        for line in check.get("expect") or ():
            doc.bullet(line, 8.6, indent=10, gray=0.25)
        if check.get("troubles"):
            doc.spacer(1)
            doc.paragraph("If not", 8.6, font="F2", indent=6, rgb=WARN_RED)
            doc.table(["What you see", "What it means", "What to do"],
                      [list(t) for t in check["troubles"]],
                      widths=[0.24, 0.34, 0.42], size=8.0)

    doc.boxed(body, pad=7, rgb=(0.52, 0.68, 0.52),
              fill=(0.975, 0.99, 0.975))


def _already_made_banner(doc, plan: dict) -> None:
    """Say it before the protocol, not after it.

    Two weeks of cloning to arrive at a plasmid that is already in the
    freezer is the most expensive mistake this tool can fail to
    prevent, and the inventory is the only thing that knows.
    """
    from .pdf import WARN_RED
    hits = plan.get("already_made") or []
    made = [h for h in hits if h["exact"] or h["coverage"] >= 0.98]
    close = [h for h in hits if h not in made]
    if not made and not close:
        return

    def body():
        if made:
            doc.paragraph("The lab already has this construct.", 11,
                          font="F2", rgb=WARN_RED)
            for h in made:
                doc.bullet(
                    f"{h['name']} - {h['length']:,} bp, {h['owner']}"
                    + (f" ({h['verified'].replace('_', ' ')})"
                       if h.get("verified") else "")
                    + f". {h['why']}.", 9, rgb=WARN_RED)
            doc.paragraph(
                "Check one of these by sequencing before cloning it "
                "again; the protocol below is only worth running if it "
                "turns out not to be what the name says.", 8.8, gray=0.3)
        if close:
            doc.paragraph(
                "Close relatives already in the inventory - different "
                "constructs, but worth a look first:", 8.8,
                font="F2" if not made else "F1", gray=0.25)
            for h in close:
                doc.bullet(f"{h['name']} - {h['length']:,} bp, {h['owner']}; "
                           f"{h['coverage'] * 100:.0f}% shared", 8.6,
                           gray=0.35)

    doc.boxed(body, pad=8,
              rgb=WARN_RED if made else (0.80, 0.82, 0.85),
              fill=(1.0, 0.97, 0.97) if made else (0.985, 0.985, 0.99))


def _kit_protocols(doc, steps: list) -> None:
    """What to buy, and where each set of conditions came from.

    The conditions above are reaction parameters - volumes, temperatures,
    times. The vendors' own protocol pages carry their figures and their
    wording, which are theirs; this points at them rather than copying
    them, and names the catalogue number so the page that opens is the
    right one.
    """
    from .reagents import used_by, alternatives
    ops = [s["operation"] for s in steps]
    products = used_by(ops)
    if not products:
        return
    doc.heading("Kits and catalogue numbers", 2)
    doc.paragraph(
        "The products this lab buys, with the page each set of conditions "
        "above was taken from. Open the vendor page for the "
        "manufacturer's own figures and for the current revision.", 8.6,
        gray=0.4)
    doc.table(["Catalogue", "Product", "Protocol"],
              [[f"{x.vendor} #{x.catalogue}", x.name, x.url]
               for x in products],
              widths=[0.17, 0.35, 0.48], size=7.4)
    extra = []
    for op in dict.fromkeys(ops):
        for alt in alternatives(op):
            if alt not in extra:
                extra.append(alt)
    if extra:
        doc.paragraph("Also on the shelf for these steps:", 8.4, gray=0.4)
        for alt in extra:
            doc.bullet(alt, 8.2, gray=0.4)
    doc.spacer(4)


def protocol_pdf(plan: dict, cfg: Config | None = None,
                 route_id: str | None = None) -> bytes:
    cfg = cfg or Config()
    p = build_protocol(plan, cfg, route_id)
    name = p["construct"].get("name") or "construct"

    doc = PDF(title=f"{name} - cloning protocol", page_size=A4,
              footer=f"{name}  -  generated by Genesis {p.get('generated','')}")

    if p["empty"]:
        doc.heading("No route could be generated", 1)
        doc.paragraph("The design produced no feasible cloning route. See the "
                      "validation messages in Genesis.")
        return doc.to_bytes()

    route = p["route"]
    construct = p["construct"]
    v = p["validation"]

    # --- 1. header --------------------------------------------------
    doc.heading(name, 1)
    doc.paragraph(route["title"], 10.5, gray=0.35)
    doc.spacer(4)
    _already_made_banner(doc, plan)

    product = (route.get("scheme") or {}).get("product")
    if product:
        facts = [
            ("Final construct", f"{product['length']:,} bp, "
                                f"{'circular' if product['circular'] else 'linear'}"),
            ("Built from", f"{product['basis']}, with {product['removed']:,} bp "
                           f"replaced by {product['added']:,} bp of "
                           f"{product['insert_label']}"),
        ]
    else:
        facts = [
            ("Final construct", f"{construct['length']:,} bp, "
                                f"{'circular' if construct['circular'] else 'linear'}"),
        ]
        if p["backbone"]:
            facts.append(("Backbone", f"{p['backbone']['name']} "
                                      f"({p['backbone']['length']:,} bp, "
                                      f"{p['backbone']['owner']})"))
    facts += [
        # The same nesting rule the map uses, so the one-line
        # architecture and the drawing of it agree: a feature inside
        # another feature is a detail of it, junk annotation included.
        ("Architecture", " - ".join(
            i["name"] for i in _top_level_layout(construct))),
        ("Workload", f"{route['n_steps'] - route.get('n_optional', 0)} "
                      f"steps, {route['n_pcr']} PCR(s), "
                     f"{route['n_assembly']} assembly reaction(s)"
                     + (f", plus {route['n_optional']} optional check(s)"
                        if route.get("n_optional") else "")),
        ("Time", f"{route['hands_on_minutes'] / 60:.1f} h hands-on, "
                 f"{route.get('turnaround_minutes', 0) / 1440:.1f} days elapsed"
                 + (f" (+{route['optional_minutes'] / 60:.1f} h if you run "
                    f"the optional checks)"
                    if route.get("optional_minutes") else "")),
    ]
    if p["suggested_name"] and p["suggested_name"] != name:
        facts.append(("Suggested name", p["suggested_name"]))
    doc.key_values(facts)

    if v.get("errors"):
        doc.spacer(6)
        doc.heading(f"{v['errors']} validation error(s) - read the detail "
                    f"section before starting", 3)
        for issue in v.get("issues", []):
            if issue["level"] == "error":
                doc.bullet(f"{issue['category']}: {issue['message']}", 9)

    # --- 2. the designed construct ----------------------------------
    layout = construct.get("layout") or []
    marks = construct.get("codon_markers") or []
    if layout:
        doc.heading("Designed construct", 2)
        feats = [{"feature_type": i["type"], "name": i["name"],
                  "start": i["start"], "end": i["end"],
                  "strand": i.get("strand", 1)} for i in layout
                 if i["type"] not in ("primer_binding_site", "unknown")]
        if construct.get("circular"):
            # A plasmid is a circle. A linear track implies the ends are
            # ends, and a cassette that wraps the origin looks broken on one.
            circular_map(doc, name, f"{construct['length']:,} bp",
                         construct["length"], feats, codon_marks=marks)
            doc.paragraph(
                "Outer ring reads forward, inner ring reverse. Position 1 is "
                "at the top.", 8.2, gray=0.45)
        else:
            plasmid_track(doc, name, f"{construct['length']:,} bp",
                          construct["length"], feats, codon_marks=marks)
        # The map above is the architecture as specified. A swap edits a
        # real plasmid and keeps whatever spacing that plasmid has around
        # the region it replaces, so its product is a slightly different
        # molecule - and printing both sizes with no word about it is how
        # this document came to state two lengths for one construct.
        sim_bp = (route.get("simulation") or {}).get("product_bp")
        if sim_bp and sim_bp != construct["length"]:
            doc.paragraph(
                f"The route below does not build this molecule exactly: its "
                f"steps produce {sim_bp:,} bp, not {construct['length']:,} bp. "
                f"The difference is the spacing the recipient already has "
                f"around the region being replaced, which an edit keeps and "
                f"an assembly from parts does not. The simulation further "
                f"down is run on what the steps actually produce.",
                8.4, gray=0.3)

        problems = [m for m in marks
                    if m["kind"] in ("premature_stop", "missing_start",
                                     "missing_stop")]
        if problems:
            for m in problems:
                doc.bullet(f"{m['label']} at {m['start']:,} in {m['module']}",
                           8.4)
        else:
            doc.paragraph("Start and stop codons are where they should be.",
                          8.4, gray=0.4)
        doc.spacer(4)

    # --- 3. how the route builds it ---------------------------------
    # The scheme below this draws two parents and the regions they give
    # up, which only exists for a swap; and the simulation section drew
    # each source as a separate track with a table saying where it
    # landed. Neither showed the reaction, so "two PCRs and a HiFi" and
    # "a digest then a ligation" looked identical on the page.
    from .figure import route_lanes, join_verb
    lanes = route_lanes(route, route.get("simulation") or {}, plan)
    if lanes:
        # Keep the heading and its caption with the drawing: the figure
        # protects its own height, which left the title stranded at the
        # foot of the page before.
        doc._ensure(figure_height(lanes) + 56)
        doc.heading("How this route builds it", 2)
        doc.paragraph(
            "One row per piece that goes in the tube: the plasmid it comes "
            "out of, the reaction that takes it out, and the piece itself. "
            "The bracket is the single reaction that closes the circle.",
            8.8, gray=0.4)
        sim_fig = route.get("simulation") or {}
        cloning_figure(doc, lanes, join_verb(route),
                       plan.get("construct") or {},
                       sim_fig.get("product_bp") or 0)
        if not sim_fig.get("simulated"):
            doc.paragraph(
                "Drawn from the planned steps. This route was not "
                "simulated, so the product shown is the design rather "
                "than a checked molecule.", 8.2, gray=0.45)

    # The "Cloning scheme" drawing used to sit here: two parent maps
    # with the swapped regions marked. "How this route builds it" above
    # draws the same parents and more besides - which piece, by what
    # reaction, joined where - and the old one was built for the swap
    # routes only, so on a restriction route it still marked primers
    # that route never designs. Two figures of one thing, one of them
    # stale.

    # --- 5. reagents ------------------------------------------------
    doc.heading("Reagents and consumables", 2)
    doc.paragraph("Quantities are per reaction. Check against the current "
                  "package insert for the kits in use.", 8.6, gray=0.4)
    doc.table(["Item", "Amount"], [[i, a] for i, a in p["reagents"]],
              widths=[0.68, 0.32], size=8.6)

    # --- 6. primers -------------------------------------------------
    if p["primers"]:
        doc.heading("Primers to order", 2)
        doc.paragraph("Lower case is the overlap or adapter tail; upper case "
                      "binds the template. The quoted Tm is for the "
                      "template-binding region only, which is what sets the "
                      "annealing temperature.", 8.6, gray=0.4)
        rows = []
        for pr in p["primers"]:
            seq = (pr["tail"].lower() + pr["anneal"].upper()) if pr["tail"] \
                else pr["sequence"].upper()
            tm_s = "n/a" if pr["anneal_tm"] != pr["anneal_tm"] else f"{pr['anneal_tm']:.1f}"
            rows.append([pr["name"], pr["direction"], f"{pr['length']}",
                         tm_s, f"{pr['gc']:.0%}", _anneal_site(pr), seq])
        doc.table(["Name", "Dir", "nt", "Tm C", "GC", "Anneals at",
                   "Sequence (5' to 3')"],
                  rows, widths=[0.15, 0.04, 0.04, 0.06, 0.05, 0.14, 0.52],
                  size=7.8, mono_cols={6})
        doc.paragraph(
            "\"Anneals at\" is the position on the template named in the "
            "step that uses it, in that plasmid's own coordinates; the "
            "parent maps above draw the same sites in purple.",
            8.2, gray=0.45)

    # --- 7. the protocol --------------------------------------------
    doc.heading("Protocol", 2)
    # A column down the right saying what each step is doing to the
    # molecule. The figure at the front says what the route does; this
    # says which part of it you are standing in.
    from .pdf import step_thumbnail
    from .figure import route_lanes
    lanes_for_steps = route_lanes(route, route.get("simulation") or {}, plan)
    doc.gutter = 96 if lanes_for_steps else 0
    for step in p["steps"]:
        top = doc.y
        doc.step(f"{step['index']}.", step["description"],
                 badge=step["operation"],
                 optional=bool(step.get("optional")))
        if lanes_for_steps:
            # Which lane this step belongs to, not merely which kind of
            # step it is: two PCRs in a row are two different pieces.
            here = next(
                (l for l in lanes_for_steps
                 if l.get("template")
                 and str(l["template"]) in (step["description"] or "")),
                None)
            step_thumbnail(doc, lanes_for_steps, step["operation"],
                           doc.width - doc.mr - 86, top, lane=here)
        if step["detail"]:
            doc.paragraph(step["detail"], 8.8, gray=0.3, indent=22)
        doc.gutter = 0
        recipe = step["recipe"]
        product = step.get("product")
        if recipe and step["operation"] in ("PCR", "HI_FI_ASSEMBLY",
                                            "GIBSON_ASSEMBLY", "LIGATION",
                                            "RE_DIGEST", "DPNI_DIGEST",
                                            "TRANSFORMATION"):
            _vendor_box(doc, recipe, product, step)
        _checkpoint_box(doc, step)
        doc.gutter = 96 if lanes_for_steps else 0

    doc.gutter = 0
    _kit_protocols(doc, p["steps"])

    # --- 8. everything else -----------------------------------------
    doc.page_break()
    doc.heading("Detail", 1)
    doc.paragraph("Everything behind the protocol above: how the design was "
                  "validated, where each fragment comes from, what the other "
                  "routes would cost, and the final sequence.", 9, gray=0.35)

    doc.heading("Validation", 2)
    doc.paragraph(f"{v.get('errors', 0)} error(s), {v.get('warnings', 0)} "
                  f"warning(s), {v.get('info', 0)} informational.", 9)
    for level, label in (("error", "Errors"), ("warning", "Warnings"),
                         ("info", "Checks")):
        items = [i for i in v.get("issues", []) if i["level"] == level]
        if not items:
            continue
        doc.heading(label, 3)
        for issue in items:
            doc.bullet(f"{issue['category']}: {issue['message']}", 8.8,
                       gray=0.45 if level == "info" else None)

    for st in (route["steps"] or []):
        d = st.get("details") or {}
        if st.get("operation") != "SYNTHESIS" or not d.get("sequence"):
            continue
        doc.heading("Fragment to order", 2)
        doc.paragraph(
            f"{d['length']:,} bp, with {d.get('homology_arm', 0)} bp of "
            f"homology to the vector at each end. Send this as written; the "
            f"arms are part of it.", 8.8)
        ordered_seq = d["sequence"]
        for i in range(0, len(ordered_seq), 60):
            chunk = ordered_seq[i:i + 60].lower()
            doc.mono(f"{i + 1:>7}  " + " ".join(
                chunk[j:j + 10] for j in range(0, len(chunk), 10)), 7.4,
                indent=4)

    cand = plan.get("candidates") or {}
    if cand.get("recipients") or cand.get("donors"):
        doc.heading("Candidate parents", 2)
        doc.paragraph("The request read as an ordered architecture: "
                 + " - ".join(cand.get("modules", []))
                 + ". The plasmid to start from is whichever already carries "
                   "the longest leading run of it; the rest is covered by as "
                   "few donors as possible.", 8.6)
        if cand.get("recipients"):
            doc.heading("Start from", 3)
            doc.table(["Plasmid", "Size", "Owner", "Already has"],
                      [[r["name"], f"{r['length']:,} bp", r["owner"],
                        ", ".join(r["carries"])]
                       for r in cand["recipients"]],
                      widths=[0.33, 0.11, 0.16, 0.40], size=8.0)
        n = cand.get("n_donors", 0)
        doc.heading(f"Then take the rest from {n} donor"
                    f"{'s' if n != 1 else ''}", 3)
        for d in cand.get("donors", []):
            doc.bullet(" - ".join(d["covers"])
                       + (f": {d['note']}" if d.get("note") else ""), 8.8)
            for o in d.get("options", [])[:4]:
                doc.bullet(f"    {o['name']} ({o['length']:,} bp, {o['owner']})"
                           f" - has {', '.join(o['carries'])}", 8.2,
                           gray=0.45)

    sim = route.get("simulation") or {}
    if sim.get("simulated"):
        doc.heading("Cloning simulation", 2)
        doc.paragraph("The primers below were amplified off the named templates "
                 "and the products assembled, to check the steps build this "
                 "construct rather than merely describe it.", 8.6)
        doc.paragraph(
            "It replays every step that changes DNA - PCR, digestion, "
            "synthesis, DpnI and the join - against the real inventory "
            "sequences. Clean-up, transformation and screening are not "
            "replayed, so this says the intended molecule can form and "
            "that nothing else in the tube can, not what fraction of "
            "colonies will carry it.", 8.2, gray=0.45)
        doc.paragraph(f"Result: {sim.get('summary', '')}", 9.0)
        if sim.get("modules"):
            doc.table(["Part", "Length", "In the product"],
                      [[m["name"], f"{m['bp']:,} bp",
                        ("yes" if m["present"] else "NO")
                        + f" ({m['coverage']*100:.0f}%)"]
                       for m in sim["modules"]],
                      widths=[0.46, 0.18, 0.36], size=8.2)
        dpni = [x for x in sim.get("steps") or []
                if x.get("operation") == "DPNI_DIGEST" and x.get("note")]
        for x in dpni:
            doc.bullet(x["note"], 8.6)
        for w in sim.get("warnings", []):
            doc.bullet(w, 8.6)
        for pr in sim.get("problems", []):
            doc.bullet(pr, 8.6)

        sources = sim.get("sources") or []
        if sources:
            doc.heading("Where each fragment is cut from", 3)
            doc.paragraph(
                "The same pieces as the figure above, each drawn against "
                "its whole parent so the coordinates can be checked "
                "against the plasmid's own map.", 8.6)
            for src in sources:
                if src.get("synthesised"):
                    continue          # ordered DNA has no parent to draw
                hi = None
                if (isinstance(src.get("from_start"), int)
                        and isinstance(src.get("from_end"), int)):
                    hi = {"start": src["from_start"], "end": src["from_end"],
                          "label": src["name"], "length": src.get("bp")}
                sub = (f"{src.get('template_length', 0):,} bp - take "
                       f"{src.get('bp', 0):,} bp as {src['name']}")
                plasmid_track(doc, src["template"], sub,
                              src.get("template_length") or 0,
                              src.get("features") or [], hi, "take")

            placed = [x for x in sources if isinstance(x.get("at_start"), int)]
            if placed and sim.get("product_bp"):
                doc.heading("How they go together", 3)
                doc.paragraph(
                    "The simulated product, with each stretch named by the "
                    "reaction it came out of.", 8.6)
                blocks = []
                for x in placed:
                    end = x["at_end"]
                    if end < x["at_start"]:
                        end = sim["product_bp"]
                    blocks.append({"feature_type": "CDS", "name": x["name"],
                                   "start": x["at_start"], "end": end,
                                   "strand": -1 if x.get("flipped") else 1})
                plasmid_track(doc, "assembled product",
                              f"{sim['product_bp']:,} bp",
                              sim["product_bp"], blocks)
                doc.table(["Fragment", "Cut from", "Length", "Lands at"],
                          [[x["name"],
                            f"{x['template']}  {x.get('from_start')}-"
                            f"{x.get('from_end')}",
                            f"{x.get('bp', 0):,} bp",
                            f"{x['at_start']:,}-{x['at_end']:,}"
                            + (" (reversed)" if x.get("flipped") else "")]
                           for x in placed],
                          widths=[0.20, 0.42, 0.14, 0.24], size=8.0)

    ver = plan.get("verification") or {}
    if ver:
        doc.heading("Read-back", 2)
        doc.paragraph("The assembled sequence re-annotated from scratch and "
                 "compared with what was asked for, independently of the "
                 "bookkeeping that built it.", 8.6)
        doc.paragraph(ver.get("summary", ""), 9.0)
        for pr in ver.get("problems", []):
            doc.bullet(pr, 8.6)

    if route["warnings"]:
        doc.heading("Route warnings", 2)
        for w in route["warnings"]:
            doc.bullet(w, 8.8)

    if p["fragments"]:
        doc.heading("Fragments and sources", 2)
        doc.table(["Fragment", "Length", "Source template", "Synthesis"],
                  [[f["name"], f"{f['length']:,} bp", f["source"],
                    "yes" if f["needs_synthesis"] else "no"]
                   for f in p["fragments"]],
                  widths=[0.30, 0.12, 0.46, 0.12], size=8.2)

    if route["junctions"]:
        doc.heading("Junctions", 2)
        doc.table(["#", "Left", "Right", "Frame shift", "Notes"],
                  [[str(j["index"]), j["left_name"], j["right_name"],
                    str(j["frame_shift"]), "; ".join(j["warnings"]) or "ok"]
                   for j in route["junctions"]],
                  widths=[0.05, 0.19, 0.19, 0.10, 0.47], size=8.0)

    doc.heading("Route comparison", 2)
    doc.table(["Route", "Steps", "PCR", "Assembly", "Hands-on", "Elapsed",
               "Cost", "Usable"],
              [[r["title"] + (" (this one)" if r["id"] == route["id"] else ""),
                str(r["n_steps"]), str(r["n_pcr"]), str(r["n_assembly"]),
                f"{r['hands_on_minutes'] / 60:.1f} h",
                f"{r.get('turnaround_minutes', 0) / 1440:.1f} d",
                f"{r['cost']:.0f}", "yes" if r["feasible"] else "no"]
               for r in p["routes"]],
              widths=[0.30, 0.08, 0.07, 0.11, 0.11, 0.09, 0.08, 0.09], size=8.0)
    for r in p["routes"]:
        if not r["feasible"] and r["warnings"]:
            doc.paragraph(f"{r['title']} was not usable: {r['warnings'][0]}",
                          8.6, gray=0.35)

    seq = construct.get("sequence", "")
    if seq:
        doc.heading("Final sequence", 2)
        doc.paragraph(f"{len(seq):,} bp, "
                      f"{gc_fraction(seq):.1%} GC, "
                      f"{'circular' if construct['circular'] else 'linear'}.",
                      8.8, gray=0.35)
        for i in range(0, len(seq), 60):
            chunk = seq[i:i + 60].lower()
            blocks = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
            doc.mono(f"{i + 1:>7}  {blocks}", 7.4, indent=4)

    return doc.to_bytes()


# --------------------------------------------------------------------------
# Markdown rendering of the same structure
# --------------------------------------------------------------------------

def protocol_markdown(plan: dict, cfg: Config | None = None,
                      route_id: str | None = None) -> str:
    cfg = cfg or Config()
    p = build_protocol(plan, cfg, route_id)
    if p["empty"]:
        return "# No route could be generated\n"

    route, construct, v = p["route"], p["construct"], p["validation"]
    name = construct.get("name") or "construct"
    out: list[str] = [f"# {name}", "",
                      f"**{route['title']}** - generated {p['generated']} by Genesis", ""]

    out.append(f"- **Final construct:** {construct['length']:,} bp "
               f"({'circular' if construct['circular'] else 'linear'})")
    if p["backbone"]:
        out.append(f"- **Backbone:** {p['backbone']['name']} "
                   f"({p['backbone']['length']:,} bp, {p['backbone']['owner']})")
    out.append("- **Architecture:** " +
               " -> ".join(i["name"] for i in construct["layout"]))
    out.append(f"- **Workload:** {route['n_steps']} steps, {route['n_pcr']} PCR(s), "
               f"{route['n_assembly']} assembly reaction(s)")
    out.append(f"- **Time:** {route['hands_on_minutes'] / 60:.1f} h hands-on, "
               f"{route.get('turnaround_minutes', 0) / 1440:.1f} days elapsed")
    if p["suggested_name"] and p["suggested_name"] != name:
        out.append(f"- **Suggested name:** `{p['suggested_name']}`")
    out.append("")

    if v.get("errors"):
        out.append(f"> **{v['errors']} validation error(s).** Read the detail "
                   f"section before starting.")
        out.append("")
        for i in v.get("issues", []):
            if i["level"] == "error":
                out.append(f"> - {i['category']}: {i['message']}")
        out.append("")

    # The same walk the figure draws, in the form a lab notebook takes.
    from .figure import route_lanes, join_verb
    lanes = route_lanes(route, route.get("simulation") or {}, plan)
    if lanes:
        out += ["## How this route builds it", ""]
        for lane in lanes:
            src = lane.get("src") or {}
            where = ("ordered as a gBlock" if lane["origin"] == "order"
                     else f"from `{src.get('template') or lane['template']}`"
                          + (f" ({src['template_length']:,} bp)"
                             if src.get("template_length") else ""))
            # The last stage's size is the piece's size, printed at the
            # end of the line; saying it twice reads as two numbers.
            final = f"{src.get('bp', 0):,} bp"
            how = " then ".join(
                st["verb"] + (f" ({st['sub']})"
                              if st.get("sub") and st["sub"] != final else "")
                for st in lane["stages"])
            ends = (f", ends {' / '.join(lane['ends'])}"
                    if lane.get("ends") else "")
            out.append(f"- **{lane['piece']}** {where} - {how}"
                       f" -> {src.get('bp', 0):,} bp{ends}")
        sim_bp = (route.get("simulation") or {}).get("product_bp") or 0
        out += ["", f"Then {join_verb(route)}: "
                    f"{len(lanes)} piece{'' if len(lanes) == 1 else 's'}"
                    + (f" -> {sim_bp:,} bp" if sim_bp else ""), ""]

    scheme = route.get("scheme") or {}
    if scheme.get("recipient"):
        rec = scheme["recipient"]
        out += ["## Cloning scheme", "",
                "Two parent plasmids: cut a region out of the recipient and "
                "replace it with a region taken from the donor.", ""]
        out.append(f"**Recipient** - `{rec['name']}` ({rec['length']:,} bp, "
                   f"{rec['owner']})")
        if rec.get("highlight"):
            h = rec["highlight"]
            out.append(f"- Cut out **{h['start']:,}-{h['end']:,}** "
                       f"({h['length']:,} bp){': ' + h['label'] if h['label'] else ''}")
        don = scheme.get("donor")
        if don:
            out.append("")
            out.append(f"**Donor** - `{don['name']}` ({don['length']:,} bp, "
                       f"{don['owner']})")
            if don.get("highlight"):
                h = don["highlight"]
                out.append(f"- Take **{h['start']:,}-{h['end']:,}** "
                           f"({h['length']:,} bp){': ' + h['label'] if h['label'] else ''}")
        bits = []
        if scheme.get("left_anchor"):
            bits.append(f"{scheme['left_anchor']['name']} (ends "
                        f"{scheme['left_anchor']['end']:,})")
        if scheme.get("right_anchor"):
            bits.append(f"{scheme['right_anchor']['name']} (starts "
                        f"{scheme['right_anchor']['start']:,})")
        if bits:
            verb = "which is kept" if len(bits) == 1 else "which are kept"
            out += ["", "The cut is bounded by " + " and ".join(bits) +
                    f", {verb}."]
        out.append("")

    marks = construct.get("codon_markers") or []
    if marks:
        out += ["## Translation landmarks", "",
                "| Position | What | Module |", "|---|---|---|"]
        for m in marks:
            out.append(f"| {m['start']:,} | {m['label']} | {m['module']} |")
        out.append("")

    out += ["## Reagents and consumables", "",
            "| Item | Amount |", "|---|---|"]
    out += [f"| {item} | {amount} |" for item, amount in p["reagents"]]
    out.append("")

    if p["primers"]:
        out += ["## Primers to order", "",
                "| Name | Dir | nt | Tm (C) | GC | Anneals at | "
                "Sequence (5'->3') |",
                "|---|---|---|---|---|---|---|"]
        for pr in p["primers"]:
            seq = (pr["tail"].lower() + pr["anneal"].upper()) if pr["tail"] \
                else pr["sequence"].upper()
            tm_s = "n/a" if pr["anneal_tm"] != pr["anneal_tm"] else f"{pr['anneal_tm']:.1f}"
            out.append(f"| {pr['name']} | {pr['direction']} | {pr['length']} | "
                       f"{tm_s} | {pr['gc']:.0%} | {_anneal_site(pr)} | "
                       f"`{seq}` |")
        out += ["", "Lower case is the overlap or adapter tail; upper case binds "
                "the template. Tm is for the template-binding region only.", ""]

    out += ["## Protocol", ""]
    for step in p["steps"]:
        out.append(f"{step['index']}. **`{step['operation']}`**"
                   + (" *(optional)*" if step.get("optional") else "")
                   + f" {step['description']}")
        if step["detail"]:
            out.append(f"   - {step['detail']}")
        recipe = step["recipe"]
        if recipe and step["operation"] in ("PCR", "HI_FI_ASSEMBLY",
                                            "GIBSON_ASSEMBLY", "LIGATION",
                                            "RE_DIGEST", "TRANSFORMATION"):
            out.append(f"   - *{recipe['title']}*: " +
                       "; ".join(f"{a} {i}" for i, a in recipe["items"]))
        check = (step.get("details") or {}).get("checkpoint")
        if check:
            for line in check["how"]:
                out.append(f"   - {line}")
            out.append("   - **Go on if you see:** "
                       + "; ".join(check["expect"]))
            for seen, means, do in check["troubles"]:
                out.append(f"   - *{seen}* - {means}. {do}.")
    out.append("")

    out += ["---", "", "## Detail", "",
            "Everything behind the protocol above.", ""]

    out += ["### Validation", "",
            f"{v.get('errors', 0)} error(s), {v.get('warnings', 0)} warning(s), "
            f"{v.get('info', 0)} informational.", ""]
    for level, label in (("error", "Errors"), ("warning", "Warnings"),
                         ("info", "Checks")):
        items = [i for i in v.get("issues", []) if i["level"] == level]
        if items:
            out.append(f"**{label}**")
            out.append("")
            out += [f"- {i['category']}: {i['message']}" for i in items]
            out.append("")

    ordered = [st for st in route["steps"]
               if st["operation"] == "SYNTHESIS"
               and (st.get("details") or {}).get("sequence")]
    for st in ordered:
        d = st["details"]
        out += ["### Fragment to order", "",
                f"{d['length']:,} bp, with {d.get('homology_arm', 0)} bp of "
                f"homology to the vector at each end. Send this as written; "
                f"the arms are part of it.", "", "```"]
        seq = d["sequence"]
        for i in range(0, len(seq), 60):
            out.append(f"{i + 1:>7}  " + " ".join(
                seq[j:j + 10] for j in range(i, min(i + 60, len(seq)), 10)))
        out += ["```", ""]

    cand = plan.get("candidates") or {}
    if cand.get("recipients") or cand.get("donors"):
        out += ["### Candidate parents", "",
                "The request read as an ordered architecture: "
                + " - ".join(cand.get("modules", [])) + ".", "",
                "**Start from** (whichever already carries the longest "
                "leading run):", "",
                "| Plasmid | Size | Owner | Already has |",
                "|---|---|---|---|"]
        for r in cand["recipients"]:
            out.append(f"| {r['name']} | {r['length']:,} bp | {r['owner']} | "
                       f"{', '.join(r['carries'])} |")
        out.append("")
        n = cand.get("n_donors", 0)
        out.append(f"**Then take the rest from {n} donor"
                   f"{'s' if n != 1 else ''}:**")
        out.append("")
        for d in cand["donors"]:
            out.append(f"- *{' - '.join(d['covers'])}*"
                       + (f" - {d['note']}" if d.get("note") else ":"))
            for o in d.get("options", [])[:4]:
                out.append(f"    - {o['name']} ({o['length']:,} bp, "
                           f"{o['owner']}) - has {', '.join(o['carries'])}")
        out.append("")

    sim = route.get("simulation") or {}
    if sim.get("simulated"):
        out += ["### Cloning simulation", "",
                "The primers below were amplified off the named templates and "
                "the products assembled, to check the steps build this "
                "construct rather than merely describe it.", "",
                f"- **Result:** {sim.get('summary', '')}", ""]
        if sim.get("modules"):
            out += ["| Part | Length | In the product |", "|---|---|---|"]
            for m in sim["modules"]:
                mark = "yes" if m["present"] else "**NO**"
                out.append(f"| {m['name']} | {m['bp']:,} bp | {mark} "
                           f"({m['coverage']*100:.0f}%) |")
            out.append("")
        for w in sim.get("warnings", []):
            out += [f"- {w}", ""]
        for pr in sim.get("problems", []):
            out += [f"- **{pr}**", ""]

    if route["warnings"]:
        out += ["### Route warnings", ""]
        out += [f"- {w}" for w in route["warnings"]]
        out.append("")

    if p["fragments"]:
        out += ["### Fragments and sources", "",
                "| Fragment | Length | Source template | Synthesis |",
                "|---|---|---|---|"]
        out += [f"| {f['name']} | {f['length']:,} bp | {f['source']} | "
                f"{'yes' if f['needs_synthesis'] else 'no'} |"
                for f in p["fragments"]]
        out.append("")

    if route["junctions"]:
        out += ["### Junctions", "",
                "| # | Left | Right | Frame shift | Notes |", "|---|---|---|---|---|"]
        out += [f"| {j['index']} | {j['left_name']} | {j['right_name']} | "
                f"{j['frame_shift']} | {'; '.join(j['warnings']) or 'ok'} |"
                for j in route["junctions"]]
        out.append("")

    out += ["### Route comparison", "",
            "| Route | Steps | PCR | Assembly | Hands-on | Elapsed | Cost | Usable |",
            "|---|---|---|---|---|---|---|---|"]
    for r in p["routes"]:
        mark = " (this one)" if r["id"] == route["id"] else ""
        out.append(f"| {r['title']}{mark} | {r['n_steps']} | {r['n_pcr']} | "
                   f"{r['n_assembly']} | {r['hands_on_minutes'] / 60:.1f} h | "
                   f"{r.get('turnaround_minutes', 0) / 1440:.1f} d | "
                   f"{r['cost']:.0f} | {'yes' if r['feasible'] else 'no'} |")
    out.append("")
    for r in p["routes"]:
        if not r["feasible"] and r["warnings"]:
            out.append(f"*{r['title']} was not usable: {r['warnings'][0]}*")
            out.append("")

    seq = construct.get("sequence", "")
    if seq:
        out += ["### Final sequence", "",
                f"{len(seq):,} bp, {gc_fraction(seq):.1%} GC.", "", "```"]
        for i in range(0, len(seq), 60):
            chunk = seq[i:i + 60].lower()
            out.append(f"{i + 1:>7}  " +
                       " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10)))
        out += ["```", ""]

    return "\n".join(out)
