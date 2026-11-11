"""Stage 12: export a bench-ready protocol, primer table and GenBank file."""

from __future__ import annotations

import datetime as _dt
import textwrap

from .config import Config
from .seqtools import gc_fraction, tm


def _fmt_tm(value: float) -> str:
    return "n/a" if value != value else f"{value:.1f}"


# The protocol itself is built once in protocol.py and rendered to Markdown
# or PDF from the same structure, so the two cannot drift apart.
from .protocol import protocol_markdown, protocol_pdf, suggest_name  # noqa: E402,F401


def primer_table_csv(plan: dict, route_id: str | None = None) -> str:
    routes = plan.get("routes", [])
    if not routes:
        return ""
    route = next((r for r in routes if r["id"] == route_id), routes[0])
    lines = ["name,direction,length,anneal_tm_C,gc_fraction,tail,anneal,full_sequence,warnings"]
    for p in route["primers"]:
        warn = "; ".join(p["warnings"]).replace(",", ";")
        lines.append(
            f"{p['name']},{p['direction']},{p['length']},{_fmt_tm(p['anneal_tm'])},"
            f"{p['gc']},{p['tail']},{p['anneal']},{p['sequence']},\"{warn}\"")
    return "\n".join(lines)


#: IDT will not synthesise a 25 nmole oligo past this length; the next
#: scale up takes it to 90. Beyond that it is a gBlock, not a primer.
IDT_MAX_25NM = 60
IDT_MAX_100NM = 90


def idt_scale(length: int) -> str:
    """The smallest IDT scale that will make an oligo this long."""
    if length <= IDT_MAX_25NM:
        return "25nm"
    if length <= IDT_MAX_100NM:
        return "100nm"
    return "250nm"


def idt_purification(length: int) -> str:
    """Standard desalting, until the length makes it a bad idea.

    Every synthesis cycle leaves some chains unextended, so the longer
    the oligo the larger the share of the tube that is truncated
    product. Past about sixty bases that fraction stops being
    negligible for an assembly primer, where a short chain still
    anneals and then seeds a deletion.
    """
    return "STD" if length <= IDT_MAX_25NM else "PAGE"


def primer_order_csv(plan: dict, route_id: str | None = None) -> str:
    """The primers in the four columns IDT's bulk order page expects.

    Name, Sequence, Scale, Purification - nothing else, because
    anything else makes the upload fail. The analysis columns live in
    primer_table_csv; this one is for pasting into an order.
    """
    routes = plan.get("routes", [])
    if not routes:
        return ""
    route = next((r for r in routes if r["id"] == route_id), routes[0])
    lines = ["Name,Sequence,Scale,Purification"]
    for p in route.get("primers") or []:
        seq = (p.get("sequence") or "").upper()
        if not seq:
            continue
        # A comma in a name would shift every column after it.
        name = str(p.get("name") or "primer").replace(",", "-")
        lines.append(f"{name},{seq},{idt_scale(len(seq))},"
                     f"{idt_purification(len(seq))}")
    return "\n".join(lines) + "\n"


def to_genbank(plan: dict) -> str:
    """Write the designed construct as an annotated GenBank file."""
    c = plan["construct"]
    seq = c["sequence"]
    today = _dt.date.today().strftime("%d-%b-%Y").upper()
    topology = "circular" if c["circular"] else "linear"
    name = (c["name"] or "construct").replace(" ", "_")[:16]

    out = [f"LOCUS       {name:<16} {len(seq)} bp ds-DNA     {topology} SYN {today}",
           f"DEFINITION  {c['name']} (designed by Genesis).",
           "ACCESSION   .",
           "VERSION     .",
           "KEYWORDS    .",
           "SOURCE      synthetic DNA construct",
           "  ORGANISM  synthetic DNA construct",
           "FEATURES             Location/Qualifiers"]

    for item in c["layout"]:
        start, end = item["start"] + 1, item["end"]
        loc = f"{start}..{end}"
        if item["strand"] == -1:
            loc = f"complement({loc})"
        key = {"CDS": "CDS", "promoter": "promoter", "polyA": "polyA_signal",
               "lox_site": "protein_bind", "AAV_ITR": "repeat_region",
               "origin": "rep_origin"}.get(item["type"], "misc_feature")
        out.append(f"     {key:<16}{loc}")
        out.append(f"                     /label=\"{item['name']}\"")
        out.append(f"                     /note=\"{item['type']}\"")

    out.append("ORIGIN")
    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60].lower()
        blocks = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        out.append(f"{i + 1:>9} {blocks}")
    out.append("//")
    return "\n".join(out)
