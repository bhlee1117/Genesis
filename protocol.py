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
from .pdf import PDF, A4, plasmid_track
from .seqtools import gc_fraction, tm

# --------------------------------------------------------------------------
# Construct naming
# --------------------------------------------------------------------------

#: Module types that belong in a construct name, in the order people read them.
NAMED_TYPES = ("promoter", "enhancer", "DIO", "CDS", "trafficking_motif",
               "p2A", "IRES", "tag")

#: Types whose own name is the architecture, not a part label.
ARCHITECTURAL_TYPES = {"DIO"}

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
    prefix = ""
    backbone = (plan.get("backbone") or {}).get("name", "") or ""
    if "AAV_ITR" in types or re.search(r"AAV", backbone, re.I):
        prefix = "pAAV"
    elif re.search(r"lenti|fugw|pLJM|pHR|FCK|FUW", backbone, re.I):
        prefix = "pLenti"

    tokens: list[str] = []
    for item in source:
        if item["type"] not in NAMED_TYPES:
            continue
        if item["type"] in ARCHITECTURAL_TYPES:
            tokens.append(item["type"])
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
        "conditions": "37 C for 30-60 min, added directly to the PCR. Only "
                      "cuts methylated (dam+) parental plasmid, not the "
                      "amplicon.",
    },
    "DNA_CLEANUP": {
        "title": "Column or bead clean-up",
        "items": [("PCR or digest reaction", "all"),
                  ("spin column or SPRI beads", "per kit")],
        "conditions": "Elute in 20-30 uL. Quantify before assembly.",
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
        "conditions": "Ice 30 min, 42 C for 30 s, ice 5 min, add SOC, recover "
                      "37 C for 1 h with shaking, plate on selective agar, "
                      "incubate overnight at 37 C.",
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
        seconds = max(int(round(length / 1000 * 30)), 15)
        bits.append(f"extension {seconds // 60}:{seconds % 60:02d} for "
                    f"{length:,} bp at 30 s/kb")
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

    if "PCR" in used_ops:
        for item, amount in _pcr_recipe(cfg)["items"]:
            add_reagent(item, amount)
    for op in used_ops:
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
        steps.append({
            "index": step["index"], "operation": op,
            "description": step["description"], "detail": detail,
            "recipe": recipe, "notes": step.get("notes", ""),
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

    # --- 1. header -----------------------------------------------------
    doc.heading(name, 1)
    doc.paragraph(route["title"], 10.5, gray=0.35)
    doc.spacer(4)

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
        ("Architecture", " - ".join(i["name"] for i in construct["layout"])),
        ("Workload", f"{route['n_steps']} steps, {route['n_pcr']} PCR(s), "
                     f"{route['n_assembly']} assembly reaction(s)"),
        ("Time", f"{route['hands_on_minutes'] / 60:.1f} h hands-on, "
                 f"{route.get('turnaround_minutes', 0) / 1440:.1f} days elapsed"),
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

    # --- 2. the designed construct --------------------------------------
    layout = construct.get("layout") or []
    marks = construct.get("codon_markers") or []
    if layout:
        doc.heading("Designed construct", 2)
        feats = [{"feature_type": i["type"], "name": i["name"],
                  "start": i["start"], "end": i["end"],
                  "strand": i.get("strand", 1)} for i in layout]
        plasmid_track(doc, name, f"{construct['length']:,} bp", 
                      construct["length"], feats, codon_marks=marks)
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

    # --- 3. the scheme --------------------------------------------------
    scheme = route.get("scheme") or {}
    if scheme.get("recipient"):
        doc.heading("Cloning scheme", 2)
        doc.paragraph("Two parent plasmids. Cut the marked region out of the "
                      "recipient, take the marked region from the donor, and "
                      "join them.", 8.8, gray=0.4)

        rec = scheme["recipient"]
        plasmid_track(
            doc, "RECIPIENT  " + rec["name"],
            f"{rec['length']:,} bp, {rec['owner']} - "
            f"{'cut here, keep the rest' if rec.get('highlight') else 'backbone'}",
            rec["length"], rec["features"], rec.get("highlight"), "cut")

        don = scheme.get("donor")
        if don:
            plasmid_track(
                doc, "DONOR  " + don["name"],
                f"{don['length']:,} bp, {don['owner']} - take the marked region",
                don["length"], don["features"], don.get("highlight"), "take")

        anchors = []
        if scheme.get("left_anchor"):
            anchors.append(f"{scheme['left_anchor']['name']} (ends at "
                           f"{scheme['left_anchor']['end']:,})")
        if scheme.get("right_anchor"):
            anchors.append(f"{scheme['right_anchor']['name']} (starts at "
                           f"{scheme['right_anchor']['start']:,})")
        if anchors:
            verb = "which is kept and is" if len(anchors) == 1 else \
                   "which are kept and are"
            doc.paragraph("The cut is bounded by " + " and ".join(anchors) +
                          f", {verb} what the primers anneal to or the "
                          f"enzymes cut between.", 8.6, gray=0.35)
        doc.spacer(4)

    # --- 4. reagents ---------------------------------------------------
    doc.heading("Reagents and consumables", 2)
    doc.paragraph("Quantities are per reaction. Check against the current "
                  "package insert for the kits in use.", 8.6, gray=0.4)
    doc.table(["Item", "Amount"], [[i, a] for i, a in p["reagents"]],
              widths=[0.68, 0.32], size=8.6)

    # --- 3. primers ----------------------------------------------------
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
                         tm_s, f"{pr['gc']:.0%}", seq])
        doc.table(["Name", "Dir", "nt", "Tm C", "GC", "Sequence (5' to 3')"],
                  rows, widths=[0.17, 0.05, 0.05, 0.07, 0.06, 0.60],
                  size=7.8, mono_cols={5})

    # --- 4. the protocol -----------------------------------------------
    doc.heading("Protocol", 2)
    for step in p["steps"]:
        doc.step(f"{step['index']}.", step["description"],
                 badge=step["operation"])
        if step["detail"]:
            doc.paragraph(step["detail"], 8.8, gray=0.3, indent=22)
        recipe = step["recipe"]
        if recipe and step["operation"] in ("PCR", "HI_FI_ASSEMBLY",
                                            "GIBSON_ASSEMBLY", "LIGATION",
                                            "RE_DIGEST", "TRANSFORMATION"):
            doc.paragraph(recipe["title"], 8.4, font="F2", indent=22)
            for item, amount in recipe["items"]:
                doc.paragraph(f"{amount:<22} {item}", 8.2, indent=30,
                              font="F3", gray=0.2)
            doc.spacer(2)

    # --- 5. everything else --------------------------------------------
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
                "| Name | Dir | nt | Tm (C) | GC | Sequence (5'->3') |",
                "|---|---|---|---|---|---|"]
        for pr in p["primers"]:
            seq = (pr["tail"].lower() + pr["anneal"].upper()) if pr["tail"] \
                else pr["sequence"].upper()
            tm_s = "n/a" if pr["anneal_tm"] != pr["anneal_tm"] else f"{pr['anneal_tm']:.1f}"
            out.append(f"| {pr['name']} | {pr['direction']} | {pr['length']} | "
                       f"{tm_s} | {pr['gc']:.0%} | `{seq}` |")
        out += ["", "Lower case is the overlap or adapter tail; upper case binds "
                "the template. Tm is for the template-binding region only.", ""]

    out += ["## Protocol", ""]
    for step in p["steps"]:
        out.append(f"{step['index']}. **`{step['operation']}`** "
                   f"{step['description']}")
        if step["detail"]:
            out.append(f"   - {step['detail']}")
        recipe = step["recipe"]
        if recipe and step["operation"] in ("PCR", "HI_FI_ASSEMBLY",
                                            "GIBSON_ASSEMBLY", "LIGATION",
                                            "RE_DIGEST", "TRANSFORMATION"):
            out.append(f"   - *{recipe['title']}*: " +
                       "; ".join(f"{a} {i}" for i, a in recipe["items"]))
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
