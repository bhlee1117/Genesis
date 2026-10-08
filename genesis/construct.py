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
from .vocab import (canonical_feature_type, find_embedded_elements,
                    min_part_length, trim_embedded)


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
    mode: str = ""                  # for Cre-lox: DIO, FLEX, FAS or DO

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
            fields = {k: v for k, v in m.items() if k in known}
            raw_type = fields.get("type", "")
            fields["type"] = canonical_feature_type(raw_type)
            # "type: DIO" names both the family and the arrangement.
            if (fields["type"] == "Cre-lox" and not fields.get("mode")
                    and raw_type != "Cre-lox"):
                fields["mode"] = raw_type
            mods.append(Module(**fields))
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
        out.append(canonical_feature_type(
            m["type"] if isinstance(m, dict) else m.type))
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
        if mod.type in ("Cre-lox", "DIO") and not mod.seq and not mod.part_id:
            resolved.append(ResolvedModule(
                mod, None, mod.name or mod.mode or "Cre-lox", "", "Cre-lox",
                [], []))
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
LOXN_SEQ = "ATAACTTCGTATAATGTATACTATACGAAGTTAT"
#: loxFAS has mutated inverted repeats, so it pairs only with itself.
LOXFAS_SEQ = "ACAACTTCGTATATACCTTTCTATACGAAGTTGT"

#: Module types a recombinase cassette encloses. The cassette stops at the
#: elements that stay outside it, such as the polyA and the WPRE.
CASSETTE_PAYLOAD_TYPES = {"CDS", "trafficking_motif", "linker", "p2A", "tag",
                          "IRES"}


@dataclass(frozen=True)
class RecombinaseMode:
    """One arrangement of recombinase sites around a payload.

    The two site types must be *heterotypic* - they must not recombine with
    each other - or the cassette would collapse the moment the recombinase
    arrived. Interleaving them rather than nesting them is what makes the
    switch one-way: the first pair inverts the payload, which flips the
    second pair into the same orientation so it excises, and the payload can
    no longer flip back.
    """
    key: str
    label: str
    outer: tuple          # (name, sequence)
    inner: tuple
    inverted: bool        # is the payload stored reverse-complemented?
    description: str


RECOMBINASE_MODES = {m.key: m for m in (
    RecombinaseMode(
        "DIO", "DIO (double-floxed inverted orientation)",
        ("lox2272", LOX2272_SEQ), ("loxP", LOXP_SEQ), True,
        "payload stored inverted; Cre turns expression ON"),
    RecombinaseMode(
        "FLEX", "FLEX (flip-excision switch)",
        ("lox2272", LOX2272_SEQ), ("loxP", LOXP_SEQ), True,
        "payload stored inverted; Cre turns expression ON"),
    RecombinaseMode(
        "FAS", "FAS (loxP / loxFAS switch)",
        ("loxFAS", LOXFAS_SEQ), ("loxP", LOXP_SEQ), True,
        "payload stored inverted; Cre turns expression ON"),
    RecombinaseMode(
        "DO", "DO (double-floxed, Cre-OFF)",
        ("lox2272", LOX2272_SEQ), ("loxP", LOXP_SEQ), False,
        "payload stored forward and expressed until Cre inverts it away"),
)}
DEFAULT_MODE = "DIO"


def resolve_mode(value: str) -> RecombinaseMode:
    key = (value or DEFAULT_MODE).strip().upper().replace("-", "").replace("_", "")
    aliases = {"CREON": "DIO", "CREOFF": "DO", "DOUBLEFLOXED": "DIO",
               "DIO": "DIO", "FLEX": "FLEX", "FAS": "FAS", "DO": "DO",
               "CRELOX": DEFAULT_MODE, "": DEFAULT_MODE}
    key = aliases.get(key, key)
    if key not in RECOMBINASE_MODES:
        raise ValueError(
            f"unknown Cre-lox mode {value!r}; use one of "
            f"{', '.join(RECOMBINASE_MODES)}")
    return RECOMBINASE_MODES[key]


def expand_dio(resolved: list[ResolvedModule]) -> tuple[list[ResolvedModule], list[str]]:
    """Turn a Cre-lox module into the cassette it actually denotes.

    Cre-lox is an architecture, not a part. It encloses the coding sequence
    that follows it in two heterotypic, antiparallel site pairs:

        outer - inner - [payload] - outer' - inner'

    In the DIO, FLEX and FAS arrangements the payload is stored inverted, so
    nothing is expressed until Cre flips it. In the DO arrangement it is
    stored forward and Cre switches it off instead. After recombination one
    site of each type remains, flanking the payload.
    """
    out: list[ResolvedModule] = []
    notes: list[str] = []
    i = 0
    while i < len(resolved):
        rm = resolved[i]
        # "DIO" is accepted as well as "Cre-lox": parts and annotations
        # stored before the rename still carry the old type.
        if rm.feature_type not in ("Cre-lox", "DIO") or rm.seq:
            out.append(rm)
            i += 1
            continue

        try:
            mode = resolve_mode(getattr(rm.module, "mode", "") or rm.module.name)
        except ValueError as exc:
            notes.append(str(exc))
            i += 1
            continue

        payload = []
        j = i + 1
        while j < len(resolved) and resolved[j].feature_type in CASSETTE_PAYLOAD_TYPES:
            if resolved[j].seq:
                payload.append(resolved[j])
            j += 1

        if not payload:
            notes.append(
                f"a {mode.label} module was given with no coding sequence "
                f"after it; the cassette has to enclose something, so it was "
                f"dropped")
            i += 1
            continue

        def site(name: str, seq: str, strand: int) -> ResolvedModule:
            mod = Module(type="lox_site", name=name,
                         orientation="forward" if strand == 1 else "reverse")
            return ResolvedModule(mod, None, name,
                                  seq if strand == 1 else revcomp(seq),
                                  "lox_site", [], [])

        out.append(site(*mode.outer, 1))
        out.append(site(*mode.inner, 1))

        # An inverted payload comes back in reverse order, each module on the
        # opposite strand, because the whole block is flipped as one.
        ordered = reversed(payload) if mode.inverted else payload
        for rm_payload in ordered:
            seq = revcomp(rm_payload.seq) if mode.inverted else rm_payload.seq
            flipped = ResolvedModule(
                rm_payload.module, rm_payload.part_id, rm_payload.part_name,
                seq, rm_payload.feature_type, rm_payload.sources,
                list(rm_payload.warnings))
            if mode.inverted:
                flipped.module = Module(
                    **{**rm_payload.module.to_dict(), "orientation": "reverse"})
            out.append(flipped)

        out.append(site(*mode.outer, -1))
        out.append(site(*mode.inner, -1))

        notes.append(
            f"{mode.label} expanded to {mode.outer[0]} - {mode.inner[0]} - [" +
            " + ".join(p.part_name for p in payload) +
            (" inverted" if mode.inverted else " forward") +
            f"] - {mode.outer[0]}' - {mode.inner[0]}'; {mode.description}")
        i = j
    return out, notes


#: Kept under the old name so existing imports continue to work.
expand_recombinase = expand_dio


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


# --------------------------------------------------------------------------
# Building the payload into a real backbone
# --------------------------------------------------------------------------

#: Element types that make up an expression cassette, and are therefore the
#: part of a recipient plasmid a new payload replaces.
CASSETTE_TYPES = {"promoter", "enhancer", "CDS", "Cre-lox", "DIO", "lox_site",
                  "p2A", "trafficking_motif", "tag", "linker", "IRES"}

#: Element types that belong to the vector itself and are never replaced:
#: the ITRs that make an AAV packageable, and the origin and resistance that
#: let the plasmid be grown at all.
VECTOR_TYPES = {"AAV_ITR", "origin", "resistance"}


def vector_arc(annotations: list[dict], length: int) -> tuple | None:
    """The arc of a circular backbone that must be kept, as ``(start, end)``.

    A transfer plasmid has two halves: the vector half, carrying the origin
    and the resistance marker that let it be grown, and the payload half,
    carrying the expression cassette. In an AAV the ITRs sit at the boundary
    between them and belong to the vector half, because the payload has to
    end up *between* the ITRs once the molecule is circular.

    The cassette frequently wraps the origin of the coordinate system - in
    one of this lab's pAAV backbones the promoter is at 9441 and its own CDS
    at 0 - so the payload cannot be described as a single ``start..end``
    slice. Describing the *kept* arc instead works in both cases: whatever
    is not kept is what the new payload replaces.
    """
    anchors = [a for a in annotations if a["feature_type"] in VECTOR_TYPES]
    if not anchors:
        return None
    start = min(a["start"] for a in anchors)
    end = max(a["end"] for a in anchors)
    if end <= start or (end - start) >= length:
        return None
    return start, end


def cassette_span(annotations: list[dict], replace_types: set) -> tuple | None:
    """Fallback for a backbone with no vector anchors: the cassette's extent."""
    candidates = [a for a in annotations
                  if a["feature_type"] in replace_types
                  and a["feature_type"] not in VECTOR_TYPES]
    if not candidates:
        return None
    start = min(a["start"] for a in candidates)
    end = max(a["end"] for a in candidates)
    return (start, end) if end > start else None


def assemble_into_backbone(resolved: list[ResolvedModule], backbone_seq: str,
                           backbone_annotations: list[dict],
                           spec_types: set, circular: bool = True
                           ) -> AssembledConstruct:
    """Put the payload where it belongs inside the recipient plasmid.

    Concatenating the payload onto the end of a backbone gives a sequence of
    the right length that is completely wrong: in an AAV it leaves the whole
    cassette outside the ITRs, where it cannot be packaged.

    What a transfer plasmid actually needs kept is everything except the
    expression cassette itself - the ITRs, the WPRE and polyA that follow
    the coding sequence, the origin and the resistance marker. Only the
    promoter-to-CDS stretch is replaced, and only up to the first WPRE or
    polyA after it, unless the specification supplies its own.

    The backbone is rotated so the payload region is contiguous before any
    of this is worked out, because in a real pAAV the cassette frequently
    wraps the origin of the coordinate system - in one of this lab's
    backbones the promoter sits at 9441 and its own CDS at 0.
    """
    warnings: list[str] = []
    missing_required: list[dict] = []

    payload_parts: list[str] = []
    payload_layout: list[dict] = []
    cursor = 0
    for rm in resolved:
        if not rm.seq:
            if rm.module.required:
                missing_required.append({"type": rm.module.type,
                                         "name": rm.module.name})
            continue
        payload_parts.append(rm.seq)
        payload_layout.append({
            "name": rm.part_name, "type": rm.feature_type,
            "start": cursor, "end": cursor + len(rm.seq),
            "strand": -1 if rm.module.orientation == "reverse" else 1,
            "part_id": rm.part_id,
        })
        cursor += len(rm.seq)
        warnings.extend(f"{rm.part_name}: {w}" for w in rm.warnings)
    payload = "".join(payload_parts)

    total = len(backbone_seq)
    arc = vector_arc(backbone_annotations, total)

    if arc is None:
        rotation = 0
        arc_end = total
        warnings.append(
            "the backbone has no ITR, origin or resistance annotated, so the "
            "insertion point was taken from its annotated cassette alone; "
            "check it before building this")
    else:
        keep_start, keep_end = arc
        # Rotate so the payload region starts at 0 and the vector arc
        # (ITRs, origin, resistance) sits at the end.
        rotation = keep_end
        arc_end = total - (keep_end - keep_start)

    def rot(pos: int) -> int:
        return (pos - rotation) % total

    rotated_seq = backbone_seq[rotation:] + backbone_seq[:rotation]
    rotated = []
    for a in backbone_annotations:
        start, end = rot(a["start"]), rot(a["end"])
        # A feature ending exactly at the rotation point maps to 0, which
        # would look like a wrap. It is the last thing in the molecule, not
        # the first - and in a pAAV it is usually an ITR, which must not be
        # lost.
        if end == 0 and start > 0:
            end = total
        if end <= start:
            # straddles the rotation point; it cannot be carried over cleanly
            warnings.append(
                f"the backbone feature {a['name']!r} spans the insertion "
                f"boundary and was not carried over")
            continue
        rotated.append({**a, "start": start, "end": end})

    # The stretch to replace: from the first cassette feature in the payload
    # region, up to the first terminator element after it.
    keeps_own_terminators = not (spec_types & {"polyA", "WPRE"})
    cassette = [a for a in rotated
                if a["feature_type"] in CASSETTE_TYPES
                and a["feature_type"] not in VECTOR_TYPES
                and a["start"] < arc_end]

    if cassette:
        replace_start = min(a["start"] for a in cassette)
        replace_end = arc_end
        if keeps_own_terminators:
            after = [a["start"] for a in rotated
                     if a["feature_type"] in ("polyA", "WPRE")
                     and a["start"] >= replace_start and a["start"] < arc_end]
            if after:
                replace_end = min(after)
    else:
        replace_start = replace_end = 0
        warnings.append(
            "the backbone has no annotated expression cassette, so the "
            "payload was inserted at the start of the payload region rather "
            "than replacing anything")

    seq = rotated_seq[:replace_start] + payload + rotated_seq[replace_end:]
    shift = len(payload) - (replace_end - replace_start)

    layout: list[dict] = []
    for a in rotated:
        if a["start"] >= replace_start and a["end"] <= replace_end:
            continue                        # replaced along with the cassette
        if a["start"] >= replace_end:
            layout.append({"name": a["name"], "type": a["feature_type"],
                           "start": a["start"] + shift, "end": a["end"] + shift,
                           "strand": a.get("strand", 1), "from_backbone": True})
        elif a["end"] <= replace_start:
            layout.append({"name": a["name"], "type": a["feature_type"],
                           "start": a["start"], "end": a["end"],
                           "strand": a.get("strand", 1), "from_backbone": True})
        else:
            warnings.append(
                f"the backbone feature {a['name']!r} straddles the insertion "
                f"point and was disrupted")

    for item in payload_layout:
        layout.append({**item, "start": item["start"] + replace_start,
                       "end": item["end"] + replace_start,
                       "from_backbone": False})

    layout.sort(key=lambda i: (i["start"], i["end"]))
    return AssembledConstruct(seq, circular, layout, warnings, missing_required)
