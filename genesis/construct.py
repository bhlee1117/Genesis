"""Stage 4: the construct specification and in-silico assembly.

A construct is an ordered list of modules. Each module names either an exact
part, a part family, or just a functional category, and the resolver picks
candidate source parts out of the inventory.
"""

from __future__ import annotations

import json
import re
import sqlite3
from dataclasses import dataclass, field, asdict

from .seqtools import clean_seq, revcomp, translate, STOP_CODONS
from .vocab import min_part_length, find_embedded_elements, trim_embedded


@dataclass
class Module:
    type: str                       # controlled feature_type
    name: str = ""                  # exact part name or family
    required: bool = True
    orientation: str = "forward"    # forward | reverse
    translated: bool = False
    preserve_frame: bool = False
    needs_linker: bool = False
    allow_stop: bool = True
    allow_synthesis: bool = False
    preferred_sources: list = field(default_factory=list)
    part_id: int | None = None      # pinned part, if the user chose one
    seq: str = ""                   # literal sequence, if supplied directly

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ConstructSpec:
    name: str
    backbone: str = ""
    backbone_plasmid_id: int | None = None
    modules: list = field(default_factory=list)
    owner: str = ""
    circular: bool = True
    notes: str = ""
    #: "AAV", "lentiviral", "plasmid", or "" to infer from the backbone.
    vector: str = ""

    def to_dict(self) -> dict:
        d = asdict(self)
        d["modules"] = [m.to_dict() if isinstance(m, Module) else m
                        for m in self.modules]
        return d

    @classmethod
    def from_dict(cls, data: dict) -> "ConstructSpec":
        mods = []
        for m in data.get("modules", []):
            known = {f for f in Module.__dataclass_fields__}
            mods.append(Module(**{k: v for k, v in m.items() if k in known}))
        return cls(
            name=data.get("name", "untitled"),
            backbone=data.get("backbone", ""),
            backbone_plasmid_id=data.get("backbone_plasmid_id"),
            modules=mods,
            owner=data.get("owner", ""),
            circular=data.get("circular", True),
            notes=data.get("notes", ""),
            vector=_normalise_vector(data.get("vector", "")),
        )

    @classmethod
    def from_yaml(cls, text: str) -> "ConstructSpec":
        """Parse the flat YAML dialect used in the specification document.

        Supports exactly the shape shown there: top-level scalars plus a
        ``modules:`` list of ``- key: value`` blocks. Written by hand so the
        tool has no third-party dependencies.
        """
        data: dict = {"modules": []}
        current: dict | None = None
        in_modules = False
        for raw in text.splitlines():
            line = raw.split("#", 1)[0].rstrip()
            if not line.strip():
                continue
            indent = len(line) - len(line.lstrip())
            stripped = line.strip()

            if stripped == "modules:":
                in_modules = True
                continue
            if not in_modules:
                if ":" in stripped:
                    k, v = stripped.split(":", 1)
                    data[k.strip()] = _coerce(v.strip())
                continue

            if stripped.startswith("- "):
                current = {}
                data["modules"].append(current)
                stripped = stripped[2:].strip()
                if not stripped:
                    continue
            if current is not None and ":" in stripped and indent > 0:
                k, v = stripped.split(":", 1)
                current[k.strip()] = _coerce(v.strip())
            elif ":" in stripped and indent == 0:
                in_modules = False
                k, v = stripped.split(":", 1)
                data[k.strip()] = _coerce(v.strip())
        return cls.from_dict(data)


#: Accept the shorthand people actually type.
_VECTOR_ALIASES = {
    "aav": "AAV", "paav": "AAV", "aav2": "AAV",
    "lenti": "lentiviral", "lentivirus": "lentiviral",
    "lentiviral": "lentiviral", "lv": "lentiviral", "plenti": "lentiviral",
    "plasmid": "plasmid", "none": "", "auto": "", "": "",
}


def _normalise_vector(value) -> str:
    key = str(value or "").strip().lower()
    if key in _VECTOR_ALIASES:
        return _VECTOR_ALIASES[key]
    raise ValueError(
        f"unknown vector {value!r}; use AAV, lentiviral, plasmid, or leave it "
        f"blank to infer from the backbone")


def _coerce(v: str):
    v = v.strip().strip('"').strip("'")
    low = v.lower()
    if low in ("true", "yes"):
        return True
    if low in ("false", "no"):
        return False
    if re.fullmatch(r"-?\d+", v):
        return int(v)
    return v


# --------------------------------------------------------------------------
# Resolving modules against the inventory
# --------------------------------------------------------------------------

@dataclass
class ResolvedModule:
    module: Module
    part_id: int | None
    part_name: str
    seq: str
    feature_type: str
    sources: list = field(default_factory=list)   # candidate plasmid rows
    warnings: list = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.seq)


def candidate_parts(conn: sqlite3.Connection, module: Module, limit: int = 12) -> list[dict]:
    """Find inventory parts that could satisfy a module."""
    sql = """SELECT p.*, COUNT(o.id) AS n_sources
             FROM part p LEFT JOIN part_occurrence o ON o.part_id = p.id
             WHERE p.feature_type = ? AND p.seq_length >= ?"""
    args: list = [module.type, min_part_length(module.type)]
    if module.name:
        sql += " AND (p.name LIKE ? COLLATE NOCASE)"
        args.append(f"%{module.name}%")
    # Relevance before popularity. A search for "LR-Voltron2" must not land
    # on "LR-Voltron2-LR-WiChR-eGFP" just because that row has more sources:
    # exact name wins, then a name that starts with the query, then the
    # shortest name containing it, which is the one closest to what was asked.
    if module.name:
        sql += """ GROUP BY p.id ORDER BY
                     (p.name = ? COLLATE NOCASE) DESC,
                     (p.name LIKE ? COLLATE NOCASE) DESC,
                     LENGTH(p.name) ASC,
                     n_sources DESC
                   LIMIT ?"""
        args += [module.name, f"{module.name}%", limit]
    else:
        sql += " GROUP BY p.id ORDER BY n_sources DESC, p.seq_length DESC LIMIT ?"
        args.append(limit)
    return [dict(r) for r in conn.execute(sql, args).fetchall()]


def part_sources(conn: sqlite3.Connection, part_id: int) -> list[dict]:
    """Candidate source plasmids for a part, best template first.

    A PCR-product record and a maintained plasmid can both contain the part,
    but only one of them is something you can pull out of the freezer and
    amplify. Prefer circular, sequence-verified, annotated plasmids; among
    equals prefer the smaller one, which amplifies more cleanly.
    """
    rows = conn.execute(
        """SELECT pl.id, pl.name, pl.owner, pl.verified, pl.available_template,
                  pl.circular, pl.seq_length, o.start, o.end, o.strand
           FROM part_occurrence o JOIN plasmid pl ON pl.id = o.plasmid_id
           WHERE o.part_id = ? AND pl.available_template = 1
           ORDER BY
             pl.circular DESC,
             (pl.verified != 'unverified') DESC,
             (SELECT COUNT(*) FROM annotation a
              WHERE a.plasmid_id = pl.id
                AND a.feature_type IN ('origin', 'resistance')) DESC,
             pl.seq_length ASC""",
        (part_id,)).fetchall()
    return [dict(r) for r in rows]


#: Module types that are translated as part of a fusion protein.
TRANSLATED_MODULE_TYPES = {"CDS", "trafficking_motif", "linker", "p2A", "tag"}


def _module_types(spec: ConstructSpec) -> list[str]:
    out = []
    for m in spec.modules:
        out.append(m["type"] if isinstance(m, dict) else m.type)
    return out


def resolve(conn: sqlite3.Connection, spec: ConstructSpec) -> list[ResolvedModule]:
    """Turn each module into a concrete sequence plus its candidate sources."""
    resolved: list[ResolvedModule] = []
    all_types = _module_types(spec)
    for position, mod in enumerate(spec.modules):
        if isinstance(mod, dict):
            mod = Module(**mod)
        warnings: list[str] = []

        # DIO denotes an architecture that wraps the coding sequence after
        # it, not a part to look up. It is expanded into lox sites and an
        # inverted payload later, so it resolves to a marker here.
        if mod.type == "DIO" and not mod.seq and not mod.part_id:
            resolved.append(ResolvedModule(mod, None, mod.name or "DIO", "",
                                           "DIO", [], []))
            continue

        if mod.seq:
            seq = clean_seq(mod.seq)
            resolved.append(ResolvedModule(mod, None, mod.name or "(literal)",
                                           seq, mod.type, [], warnings))
            continue

        part = None
        if mod.part_id:
            row = conn.execute("SELECT * FROM part WHERE id = ?",
                               (mod.part_id,)).fetchone()
            part = dict(row) if row else None
            if part is None:
                warnings.append(f"pinned part id {mod.part_id} not found")

        if part is None:
            cands = candidate_parts(conn, mod)
            if not cands:
                warnings.append(
                    f"no inventory part matches type={mod.type!r} name={mod.name!r}"
                    + ("; synthesis allowed" if mod.allow_synthesis else
                       "; synthesis not allowed"))
                resolved.append(ResolvedModule(mod, None, mod.name or mod.type,
                                               "", mod.type, [], warnings))
                continue
            part = cands[0]

            # The last translated module has to terminate translation. Lab
            # annotations are split on this - the same CDS appears with and
            # without its stop codon - so when the best-ranked candidate has
            # no stop, prefer an equally good one that does.
            next_t = all_types[position + 1] if position + 1 < len(all_types) else None
            is_terminal = next_t not in TRANSLATED_MODULE_TYPES
            if (is_terminal and part["feature_type"] in TRANSLATED_MODULE_TYPES
                    and part["seq"][-3:].upper() not in STOP_CODONS):
                better = next(
                    (c for c in cands
                     if c["seq"][-3:].upper() in STOP_CODONS
                     and abs(c["seq_length"] - part["seq_length"]) <= 6), None)
                if better is not None:
                    warnings.append(
                        f"the best-ranked {part['name']!r} annotation "
                        f"(part {part['id']}, {part['seq_length']} bp) carries "
                        f"no stop codon and nothing translated follows it; "
                        f"used part {better['id']} "
                        f"({better['seq_length']} bp) instead, which ends in one")
                    part = better

            if len(cands) > 1:
                warnings.append(
                    f"{len(cands)} candidate parts matched; using {part['name']!r}")

        seq = part["seq"]

        # A CDS annotation often already carries a fused element such as a
        # trailing P2A. If the specification also asks for that element as
        # its own module, using the part as-is would build two of them in
        # tandem, so the fused copy is trimmed off here instead.
        embedded = find_embedded_elements(seq, part["feature_type"])
        if embedded:
            names = ", ".join(sorted({e["name"] for e in embedded}))
            neighbours = set(all_types[max(position - 1, 0):position + 2])
            clashing = {e["feature_type"] for e in embedded} & neighbours
            if clashing:
                for ftype in clashing:
                    seq, removed = trim_embedded(seq, embedded, ftype)
                    if removed:
                        warnings.append(
                            f"{part['name']!r} already contains "
                            f"{', '.join(r['name'] for r in removed)} fused to "
                            f"it, and the specification also asks for a "
                            f"{ftype} module; the fused copy was trimmed off "
                            f"so the construct does not end up with two")
                    # Anything of the same type buried inside the coding
                    # sequence cannot be removed without changing the
                    # protein, so say so rather than trimming silently.
                    stuck = [e for e in embedded
                             if e["feature_type"] == ftype
                             and e not in removed
                             and not (e["at_end"] or e["at_start"])]
                    if stuck:
                        warnings.append(
                            f"{part['name']!r} also contains "
                            f"{', '.join(sorted({e['name'] for e in stuck}))} "
                            f"internally, which cannot be trimmed without "
                            f"changing the protein; the construct will carry "
                            f"it in addition to the {ftype} module you asked "
                            f"for")
            else:
                warnings.append(
                    f"{part['name']!r} has {names} fused into it; the "
                    f"construct will contain that element even though it is "
                    f"not a separate module")

        # A coding sequence that still carries its own stop codon terminates
        # translation before anything downstream. When the next module is also
        # translated, the stop has to go or the fusion simply does not exist.
        next_type = all_types[position + 1] if position + 1 < len(all_types) else None
        if (part["feature_type"] in TRANSLATED_MODULE_TYPES
                and next_type in TRANSLATED_MODULE_TYPES
                and len(seq) >= 3 and len(seq) % 3 == 0
                and seq[-3:].upper() in STOP_CODONS):
            seq = seq[:-3]
            warnings.append(
                f"{part['name']!r} ends in a stop codon, which would stop "
                f"translation before the {next_type} that follows it; the "
                f"stop was removed so the fusion reads through")

        if mod.orientation == "reverse":
            seq = revcomp(seq)
        sources = part_sources(conn, part["id"])
        if not sources:
            warnings.append("part has no source plasmid; would require synthesis")
        resolved.append(ResolvedModule(mod, part["id"], part["name"], seq,
                                       part["feature_type"], sources, warnings))
    return resolved


# --------------------------------------------------------------------------
# In-silico assembly
# --------------------------------------------------------------------------

@dataclass
class AssembledConstruct:
    seq: str
    circular: bool
    layout: list = field(default_factory=list)   # [{name,type,start,end,strand}]
    warnings: list = field(default_factory=list)
    #: Required modules that resolved to nothing. A construct missing one of
    #: these is not the construct that was asked for, so no route may be
    #: offered for it.
    missing_required: list = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.seq)


def assemble(resolved: list[ResolvedModule], backbone_seq: str = "",
             circular: bool = True) -> AssembledConstruct:
    """Concatenate modules in order into the intended final molecule.

    The backbone is placed first so that module coordinates are reported in
    the same frame a user would read off a plasmid map.
    """
    parts: list[str] = []
    layout: list[dict] = []
    warnings: list[str] = []
    cursor = 0

    if backbone_seq:
        parts.append(backbone_seq)
        layout.append({"name": "backbone", "type": "backbone", "start": 0,
                       "end": len(backbone_seq), "strand": 1})
        cursor = len(backbone_seq)

    missing_required: list[dict] = []
    for rm in resolved:
        if not rm.seq:
            if rm.module.required:
                missing_required.append({"type": rm.module.type,
                                         "name": rm.module.name})
            continue
        parts.append(rm.seq)
        layout.append({
            "name": rm.part_name, "type": rm.feature_type,
            "start": cursor, "end": cursor + len(rm.seq),
            "strand": -1 if rm.module.orientation == "reverse" else 1,
            "part_id": rm.part_id,
        })
        cursor += len(rm.seq)
        warnings.extend(f"{rm.part_name}: {w}" for w in rm.warnings)

    return AssembledConstruct("".join(parts), circular, layout, warnings,
                              missing_required)


def translated_span(construct: AssembledConstruct) -> list[dict]:
    """Identify each contiguous run of translated modules (a fusion ORF).

    A run is broken by a non-translated element *and* by a change of strand:
    an inverted DIO payload is a real fusion, but it reads off the other
    strand, so lumping it in with forward modules would translate nonsense.
    For a reverse-strand run the members are reported in the order the
    ribosome sees them, which is the reverse of their layout order.
    """
    runs: list[dict] = []
    current: list[dict] = []
    translated_types = {"CDS", "trafficking_motif", "linker", "p2A", "tag"}

    def flush():
        nonlocal current
        if not current:
            return
        strand = current[0]["strand"]
        members = list(current) if strand == 1 else list(reversed(current))
        runs.append({"start": current[0]["start"], "end": current[-1]["end"],
                     "strand": strand, "members": members})
        current = []

    for item in construct.layout:
        if item["type"] not in translated_types:
            flush()
            continue
        if current and item["strand"] != current[-1]["strand"]:
            flush()
        current.append(item)
    flush()
    return runs


# --------------------------------------------------------------------------
# DIO / FLEX expansion
# --------------------------------------------------------------------------

LOXP_SEQ = "ATAACTTCGTATAATGTATGCTATACGAAGTTAT"
LOX2272_SEQ = "ATAACTTCGTATAAAGTATCCTATACGAAGTTAT"

#: Module types a DIO cassette encloses. The cassette stops at the elements
#: that stay outside it, such as the polyA and the WPRE.
DIO_PAYLOAD_TYPES = {"CDS", "trafficking_motif", "linker", "p2A", "tag", "IRES"}


def expand_dio(resolved: list[ResolvedModule]) -> tuple[list[ResolvedModule], list[str]]:
    """Turn a ``DIO`` module into the cassette it actually denotes.

    DIO is an architecture, not a part: it encloses the coding sequence that
    follows it in two heterotypic, antiparallel lox pairs, with the payload
    inverted so nothing is expressed until Cre flips it:

        lox2272 - loxP - [payload reversed] - lox2272' - loxP'

    The pairs are interleaved rather than nested, which is what makes the
    switch one-way: Cre first inverts the payload through one pair, and that
    flips the other pair into the same orientation so it excises, locking the
    payload in the expressed orientation. After recombination one lox2272 and
    one loxP remain, flanking the now-forward gene.
    """
    out: list[ResolvedModule] = []
    notes: list[str] = []
    i = 0
    while i < len(resolved):
        rm = resolved[i]
        if rm.feature_type != "DIO" or rm.seq:
            out.append(rm)
            i += 1
            continue

        payload = []
        j = i + 1
        while j < len(resolved) and resolved[j].feature_type in DIO_PAYLOAD_TYPES:
            if resolved[j].seq:
                payload.append(resolved[j])
            j += 1

        if not payload:
            notes.append(
                "a DIO module was given with no coding sequence after it; a "
                "DIO cassette has to enclose something, so it was dropped")
            i += 1
            continue

        def site(name: str, seq: str, strand: int) -> ResolvedModule:
            mod = Module(type="lox_site", name=name,
                         orientation="forward" if strand == 1 else "reverse")
            return ResolvedModule(mod, None, name,
                                  seq if strand == 1 else revcomp(seq),
                                  "lox_site", [], [])

        out.append(site("lox2272", LOX2272_SEQ, 1))
        out.append(site("loxP", LOXP_SEQ, 1))

        # The payload is inverted as one block, so the modules come back in
        # reverse order, each on the opposite strand.
        for rm_payload in reversed(payload):
            flipped = ResolvedModule(
                rm_payload.module, rm_payload.part_id, rm_payload.part_name,
                revcomp(rm_payload.seq), rm_payload.feature_type,
                rm_payload.sources, list(rm_payload.warnings))
            flipped.module = Module(
                **{**rm_payload.module.to_dict(), "orientation": "reverse"})
            out.append(flipped)

        out.append(site("lox2272", LOX2272_SEQ, -1))
        out.append(site("loxP", LOXP_SEQ, -1))

        notes.append(
            "DIO expanded to lox2272 - loxP - [" +
            " + ".join(p.part_name for p in payload) +
            " inverted] - lox2272' - loxP'; nothing is expressed until Cre "
            "inverts the payload and excises one pair")
        i = j
    return out, notes


# --------------------------------------------------------------------------
# Start and stop codons
# --------------------------------------------------------------------------

def codon_markers(construct: AssembledConstruct) -> list[dict]:
    """Locate the start and stop codons of every translated run.

    Returns markers in construct coordinates so a map can draw them. Four
    things are worth seeing on a map, and three of them are problems:

    * ``start``           an ATG where the fusion begins
    * ``stop``            the stop codon that ends it
    * ``premature_stop``  an in-frame stop *inside* the fusion
    * ``missing_start`` / ``missing_stop``  the codon that should be there

    Reverse-strand runs (an inverted DIO payload) are read off their own
    strand, and their coordinates are mapped back onto the construct.
    """
    markers: list[dict] = []
    for run in translated_span(construct):
        strand = run.get("strand", 1)
        start, end = run["start"], run["end"]
        span = construct.seq[start:end]
        if strand == -1:
            span = revcomp(span)
        members = run["members"]
        names = " + ".join(m["name"] for m in members)
        n_codons = len(span) // 3

        def coords(i: int) -> tuple[int, int]:
            """Construct coordinates of codon `i` of this run."""
            if strand == 1:
                return start + 3 * i, start + 3 * i + 3
            return end - 3 * (i + 1), end - 3 * i

        def owner(i: int) -> str:
            lo, _ = coords(i)
            for m in members:
                if m["start"] <= lo < m["end"]:
                    return m["name"]
            return names

        if n_codons == 0:
            continue

        lo, hi = coords(0)
        first = span[0:3].upper()
        markers.append({
            "kind": "start" if first == "ATG" else "missing_start",
            "start": lo, "end": hi, "strand": strand,
            "codon": first, "module": owner(0), "run": names,
            "label": "ATG" if first == "ATG" else f"no ATG ({first})",
        })

        for i in range(n_codons):
            codon = span[3 * i:3 * i + 3].upper()
            if codon not in STOP_CODONS:
                continue
            lo, hi = coords(i)
            terminal = (i == n_codons - 1)
            markers.append({
                "kind": "stop" if terminal else "premature_stop",
                "start": lo, "end": hi, "strand": strand,
                "codon": codon, "module": owner(i), "run": names,
                "label": codon if terminal else f"{codon} (premature)",
            })

        if span[-3:].upper() not in STOP_CODONS:
            lo, hi = coords(n_codons - 1)
            # "Just past the last codon" follows the direction of translation,
            # which runs towards *lower* construct coordinates on the reverse
            # strand. Taking `hi` there would mark the start of the last codon
            # instead, three bases out.
            past_end = hi if strand == 1 else lo
            markers.append({
                "kind": "missing_stop", "start": past_end, "end": past_end,
                "strand": strand, "codon": "", "module": owner(n_codons - 1),
                "run": names, "label": "no stop codon",
            })

    return sorted(markers, key=lambda m: m["start"])
