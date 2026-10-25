"""Grouping coding sequences by the protein they encode.

The part picker searched names, and a name is the worst thing to search in
this inventory. A query for "CheRiff" returned 103 rows, of which only 29
are the CheRiff part: six of them identical molecules shown six times,
two of them 624 aa soma-targeted versions wearing the same bare label, and
the rest fusions that merely contain it. Meanwhile 56 rows that do encode
CheRiff never appeared at all, because their annotation says "ChR64" or
"YQ0801L" or "Voltron2".

So the picker groups by protein instead. One row per family, the plain
protein first, and each variant named for what is fused to it - a Kv2.1
soma-targeting tag, a Lucy-Rho leader, a fluorescent protein - rather than
for whatever the file happened to call the whole thing.

What this module will not do is put two generations of a sensor in one
family. GCaMP6s and GCaMP6m differ by three residues, QuasAr2 and QuasAr3
by one, and jGCaMP7f's protein is an exact substring of GCaMP6f's. Those
distances are indistinguishable from an annotation boundary, so the
version in the name decides and sequence never overrides it.
"""

from __future__ import annotations

import sqlite3

from .swap import protein_of, same_generation
from .vocab import find_embedded_elements

#: A protein shorter than this is a fragment, not a family core.
MIN_CORE_AA = 60

#: How many residues two annotations of one protein may differ by and still
#: be the same protein. An absolute number, not a fraction: a fifth of a
#: 320 aa protein is 64 residues, which is a whole domain - the Kv2.1 soma
#: tag is 65 and the Lucy-Rho leader 41 - and absorbing those as "a
#: different boundary" is exactly the distinction the picker exists to
#: draw. The real trims measured in this catalogue span 0-7 residues.
TRIM_AA = 15

#: A floor as well, so a short peptide cannot be swallowed by a long one.
TRIM_RATIO = 0.8

#: A decoration smaller than this is not worth naming on its own.
MIN_DECORATION_AA = 5

#: Translations are expensive enough to be worth keeping (336 ms for the
#: whole table) and cheap enough to rebuild when the catalogue changes.
_CACHE: dict = {}


def _translated(conn: sqlite3.Connection, feature_type: str = "CDS") -> list:
    """Every part of this type, with its protein. Cached per catalogue state."""
    stamp = conn.execute(
        "SELECT COUNT(*), COALESCE(MAX(id), 0) FROM part WHERE feature_type = ?",
        (feature_type,)).fetchone()
    key = (feature_type, stamp[0], stamp[1])
    if _CACHE.get("key") == key:
        return _CACHE["rows"]
    rows = []
    for r in conn.execute(
            """SELECT p.id, p.name, p.seq, p.seq_length, p.contains,
                      COUNT(o.id) AS n_sources
                 FROM part p LEFT JOIN part_occurrence o ON o.part_id = p.id
                WHERE p.feature_type = ? AND p.seq IS NOT NULL
                GROUP BY p.id""", (feature_type,)):
        protein = protein_of(r["seq"])
        if len(protein) < 20:
            continue
        rows.append({
            "id": r["id"], "name": (r["name"] or "").strip(),
            "bp": r["seq_length"], "n_sources": r["n_sources"],
            "protein": protein, "seq": r["seq"],
        })
    _CACHE.clear()
    _CACHE.update(key=key, rows=rows)
    return rows


def _trim_m(protein: str) -> str:
    """Without the initiator methionine an annotation may or may not include.

    Half the duplication in this catalogue is one leading residue: a 316 aa
    "CheRiff" starting with M and a 320 aa one that does not are the same
    protein, and a plain substring test says they are unrelated.
    """
    return protein[1:] if protein[:1] == "M" else protein


def _same_protein(a: str, b: str) -> bool:
    """One protein at two annotation boundaries."""
    x, y = _trim_m(a), _trim_m(b)
    if x == y:
        return True
    short, long_ = (x, y) if len(x) <= len(y) else (y, x)
    if len(long_) - len(short) > TRIM_AA:
        return False
    return short in long_ and len(short) >= TRIM_RATIO * len(long_)


def _name_for(flank: str, rows: list) -> str:
    """What a fused stretch of protein is called, if the lab has it catalogued.

    Saying "a fused protein" of the 239 aa on the end of CheRiff_eGFP is
    technically true and useless; it is eGFP, and eGFP is in the catalogue.
    """
    target = _trim_m(flank)
    best, best_used, best_len = "", -1, 0
    for row in rows:
        other = _trim_m(row["protein"])
        if len(other) < MIN_CORE_AA or len(other) > len(target) + 12:
            continue
        # A fused partner usually arrives with a linker on it, so the flank
        # is a little longer than the protein it is. Containment with a
        # generous floor finds eGFP inside a 256 aa flank; an exact length
        # match does not.
        if other not in target:
            continue
        if len(other) < 0.6 * len(target):
            continue
        rank = (len(other), row["n_sources"])
        if rank > (best_len, best_used):
            best, best_used, best_len = row["name"], row["n_sources"], len(other)
    return best


def _describe_flank(flank: str, conn: sqlite3.Connection, where: str,
                    rows: list | None = None) -> dict:
    """Name whatever is fused on one side of a core protein."""
    if len(flank) < MIN_DECORATION_AA:
        return {}
    # Short, known peptides first: these are the motifs people name.
    hits = find_embedded_elements(_back_translate_probe(flank), "CDS")
    for hit in hits:
        return {"where": where, "name": hit["name"],
                "type": hit.get("feature_type", "trafficking_motif"),
                "aa": len(flank)}
    # Anything substantial is most likely another protein fused on - a
    # fluorescent partner, usually - and the catalogue probably knows it.
    if rows and len(flank) >= MIN_CORE_AA:
        known = _name_for(flank, rows)
        if known:
            return {"where": where, "name": known, "type": "CDS",
                    "aa": len(flank)}
    return {"where": where,
            "name": (f"a {len(flank)} aa fused protein"
                     if len(flank) >= MIN_CORE_AA
                     else f"a {len(flank)} aa peptide"),
            "type": "CDS" if len(flank) >= MIN_CORE_AA else "unknown",
            "aa": len(flank)}


def _back_translate_probe(protein: str) -> str:
    """A nucleotide string that translates to `protein`.

    find_embedded_elements works on DNA because that is what a part stores.
    The flanks here are already protein, so they are turned back into one
    arbitrary coding sequence purely so the same motif library can be used
    on them; the codons chosen are irrelevant to a protein-level match.
    """
    table = {"A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT",
             "Q": "CAA", "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT",
             "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTT", "P": "CCT",
             "S": "TCT", "T": "ACT", "W": "TGG", "Y": "TAT", "V": "GTT"}
    return "".join(table.get(aa, "NNN") for aa in protein)


def _variant(conn: sqlite3.Connection, row: dict, core: str,
             rows: list | None = None) -> dict:
    """One member of a family, described by what is attached to the core."""
    protein = row["protein"]
    at = protein.find(core)
    if at < 0:
        at = _trim_m(protein).find(_trim_m(core))
        protein = _trim_m(protein)
        core = _trim_m(core)
    decorations = []
    if at >= 0:
        n_flank, c_flank = protein[:at], protein[at + len(core):]
        for flank, where in ((n_flank, "N"), (c_flank, "C")):
            got = _describe_flank(flank, conn, where, rows)
            if got:
                decorations.append(got)
    if decorations:
        label = " + ".join(
            f"{d['name']} ({d['where']})" for d in decorations)
    elif _same_protein(protein, core):
        label = "plain"
    else:
        label = f"{len(protein)} aa form"
    return {
        "part_id": row["id"], "name": row["name"], "bp": row["bp"],
        "aa": len(protein), "n_sources": row["n_sources"],
        "decorations": decorations, "label": label,
        "is_core": len(protein) == len(core),
    }


def protein_families(conn: sqlite3.Connection, term: str,
                     feature_type: str = "CDS", limit: int = 12,
                     deep: bool = True) -> list[dict]:
    """Parts matching `term`, grouped into one entry per protein.

    `deep` also sweeps the catalogue for parts that encode the same
    protein under an unrelated name, which is where most of the missing
    recall lives - 56 of the rows encoding CheRiff are named something
    else entirely.
    """
    term = (term or "").strip()
    rows = _translated(conn, feature_type)
    if term:
        named = [r for r in rows if term.lower() in r["name"].lower()]
    else:
        named = sorted(rows, key=lambda r: -r["n_sources"])[:200]
    if not named:
        return []

    named.sort(key=lambda r: (-r["n_sources"], len(r["protein"])))

    families: list[dict] = []
    taken: set = set()

    for seed in named:
        if seed["id"] in taken:
            continue
        core = seed["protein"]
        if len(core) < MIN_CORE_AA and len(families) >= limit:
            continue
        members = [seed]
        taken.add(seed["id"])

        # Everything that carries this protein. A longer member is a
        # decorated version of the core; a slightly shorter one is the same
        # protein annotated at a different boundary.
        pool = rows if deep else named
        for other in pool:
            if other["id"] in taken:
                continue
            protein = other["protein"]
            contains = _trim_m(core) in _trim_m(protein)
            if not (contains or _same_protein(protein, core)):
                continue
            if len(protein) < TRIM_RATIO * len(core) and not contains:
                continue
            # The one thing sequence may not decide.
            if not same_generation(seed["name"], other["name"]):
                continue
            members.append(other)
            taken.add(other["id"])

        members.sort(key=lambda r: (-r["n_sources"], len(r["protein"])))

        # Which plasmids each of these parts is actually in. Four families
        # can all be called "QuasAr2" and differ by a residue or two, and
        # then the name tells you nothing - the plasmid it came out of
        # does. Fetched once for the whole family rather than per row.
        ids = [m["id"] for m in members]
        by_part: dict = {}
        if ids:
            marks = ",".join("?" * len(ids))
            for r in conn.execute(
                    f"""SELECT o.part_id AS pid, pl.name AS nm,
                               pl.seq_length AS bp
                          FROM part_occurrence o
                          JOIN plasmid pl ON pl.id = o.plasmid_id
                         WHERE o.part_id IN ({marks})
                         ORDER BY pl.seq_length""", ids):
                by_part.setdefault(r["pid"], []).append(r["nm"])
        hosts = []
        for m in members[:12]:
            for nm in by_part.get(m["id"], []):
                if nm not in hosts:
                    hosts.append(nm)
            if len(hosts) >= 6:
                break
        # One entry per distinct protein. Four rows that are the same
        # molecule trimmed four ways are one thing to choose, not four,
        # and showing them separately is the whole complaint.
        collapsed: list = []
        for m in members:
            for group in collapsed:
                if _same_protein(group["rep"]["protein"], m["protein"]):
                    group["rows"].append(m)
                    break
            else:
                collapsed.append({"rep": m, "rows": [m]})
        variants = []
        for group in collapsed:
            v = _variant(conn, group["rep"], core, rows)
            v["n_sources"] = sum(x["n_sources"] for x in group["rows"])
            v["n_rows"] = len(group["rows"])
            v["also"] = sorted({x["name"] for x in group["rows"]
                                if x["name"] != v["name"]})[:4]
            v["part_ids"] = [x["id"] for x in group["rows"]][:8]
            # The plasmids this exact protein sits in - the answer to
            # "where would I get it from", which is the reason anybody is
            # reading this list.
            seen_h: list = []
            for x in group["rows"]:
                for nm in by_part.get(x["id"], []):
                    if nm not in seen_h:
                        seen_h.append(nm)
            v["hosts"] = seen_h[:12]
            v["n_hosts"] = len(seen_h)
            variants.append(v)
        # Plain first, then the decorated ones by how widely they are used.
        variants.sort(key=lambda v: (not v["is_core"], -v["n_sources"]))
        families.append({
            "core": {"part_id": seed["id"], "name": seed["name"],
                     "aa": len(core), "bp": seed["bp"],
                     "n_sources": seed["n_sources"],
                     "protein": core},
            "hosts": hosts,
            "n_variants": len(variants),
            "n_rows": len(members),
            "sources": sum(m["n_sources"] for m in members),
            "variants": variants[:24],
        })
        if len(families) >= limit:
            break

    # A family whose name matches exactly is what was asked for.
    low = term.lower()
    families.sort(key=lambda f: (
        f["core"]["name"].lower() != low,
        not f["core"]["name"].lower().startswith(low),
        -f["sources"],
    ))

    # Where several families carry the same name - and they do; QuasAr2 is
    # four distinct proteins in this inventory - say how each differs from
    # the most-used one. Identical sequences are already one family, so a
    # second entry under one name always means the protein really differs.
    first_by_name: dict = {}
    for fam in families:
        key = fam["core"]["name"].lower()
        ref = first_by_name.setdefault(key, fam)
        if ref is fam:
            fam["differs"] = ""
            continue
        a, b = _trim_m(ref["core"]["protein"]), _trim_m(fam["core"]["protein"])
        n = min(len(a), len(b))
        subs = sum(1 for i in range(n) if a[i] != b[i])
        gap = abs(len(a) - len(b))
        bits = []
        if subs:
            bits.append(f"{subs} residue{'s' if subs != 1 else ''} different")
        if gap:
            bits.append(f"{gap} aa longer" if len(b) > len(a)
                        else f"{gap} aa shorter")
        fam["differs"] = " and ".join(bits) or "a different protein"

    for fam in families:
        fam["core"].pop("protein", None)
    return families
