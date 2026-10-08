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
