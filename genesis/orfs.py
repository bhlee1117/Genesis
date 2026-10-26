"""Inferring the coding sequence nobody annotated.

A third of the AAV and lentiviral plasmids on the share carry a payload
whose CDS was never annotated: the promoter is marked, the WPRE is marked,
and the 1.8 kb between them is blank.  The protein is in the file, but
Genesis only catalogues what an annotation points at, so those payloads
never become parts and never turn up in a search.

What makes them recoverable is that the lab has usually catalogued the same
protein, or a close relative, somewhere else.  The sensor in
YD043-pAAV-syn-iGluSnFR4f is unannotated, but it is 91% identical to an
iGluSnFr3 already in the catalogue, which is enough to say what it is
without reading the file name at all.

Everything inferred here is marked ``evidence='orf'`` and carries lower
confidence, so it stays visibly distinct from an annotation somebody in the
lab actually made, and `genesis infer-orfs --refresh` recomputes the lot.
"""

from __future__ import annotations

import json
import re
import sqlite3

from .seqtools import find_orfs, translate
from .vocab import classify_label

#: Shortest inferred CDS worth cataloguing.  Below this the hit rate on
#: chance reading frames climbs fast, and a real payload is never this small.
MIN_ORF_AA = 150

#: An ORF already this well covered by coding annotations is spoken for.
COVERED = 0.30

#: Least confidence that earns a place in the parts catalogue.  Below it an
#: ORF is still annotated on its plasmid, just not offered as a part.
CATALOGUE_CONFIDENCE = 0.6

#: Above this the ORF is that protein; below, it is a relative of it.
IDENTICAL = 0.95
RELATED = 0.60

#: Catalogue names that name nothing.  Inheriting one of these just moves
#: the problem, so the plasmid name gets its turn instead.
_GENERIC_NAME = re.compile(
    r"^(cds|orf|orf frame \d+|unnamed|unknown|predicted|hypothetical"
    r"|insert|product|protein|gene|sequence|translation|target|feature"
    r"|region|misc|reversed|expected|complement|template)\b", re.I)


def _useful(name: str) -> bool:
    """Whether a catalogue name is worth inheriting."""
    n = (name or "").strip()
    # Segment parts are stored as "Voltron2 [2 of 3]"; the suffix is
    # bookkeeping about the split, not part of the protein's name.
    n = re.sub(r"\s*\[\d+ of \d+\]\s*$", "", n)
    # "P (fragment)" names a fragment of something called P, which is not a
    # name anyone can search for.
    bare = re.sub(r"\s*\(fragment\)\s*$", "", n, flags=re.I)
    return len(bare) >= 3 and not _GENERIC_NAME.match(n)


def _clean_name(name: str) -> str:
    return re.sub(r"\s*\[\d+ of \d+\]\s*$", "", (name or "").strip())

#: Types that mean "something already claims this stretch".
_CODING = ("CDS", "resistance", "trafficking_motif", "tag", "linker", "p2A")

#: Elements an ORF has no business overlapping.  A reading frame that runs
#: through the origin or an ITR is an artefact of a circular sequence, not
#: a payload.
#: ``intron`` belongs here for the same reason as ``promoter``: the CAG
#: cassette's chimeric intron is a kilobase of GC-rich DNA, it reads as a
#: 202 aa Ala/Arg/Pro protein in one frame, and a spliced-out intron is
#: the one part of a transcript that is certainly not translated.
_STRUCTURAL = ("origin", "AAV_ITR", "lox_site", "WPRE", "polyA", "promoter",
               "intron")

#: Name fragments that describe the vector rather than the insert.
_NAME_NOISE = {
    "paav", "aav", "plenti", "lenti", "lentivirus", "pcdna", "pet", "pgex",
    "syn", "hsyn", "ef1a", "ef1", "cag", "cmv", "camkii", "u6", "tre", "pgk",
    "wpre", "ires", "polya", "pa", "itr", "dio", "flex", "fas", "do", "cre",
    "lox", "loxp", "flpdio", "flp", "frt", "p2a", "t2a", "e2a", "f2a", "2a",
    "ampr", "kanr", "puro", "neo", "hygro", "bla", "ori", "backbone",
    "plasmid", "vector", "construct", "final", "new", "old", "copy", "seq",
    "sequenced", "clone", "clones", "addgene", "ngr", "kgc", "er2", "ts",
    "cut", "cutvector", "predicted", "corrected", "bbcorrected", "bb",
    "insert", "fragment", "frag", "gblock", "part", "test", "temp", "tmp",
    "version", "ver", "draft", "assembly", "assembled", "map", "full",
    "sequence", "sequencing", "tbd", "fdio", "fdo", "con", "coff", "goal",
    "contig", "colony", "stats", "result", "results", "unknown", "none",
    # Status words a lab puts on the end of a file name. Without these the
    # insert name of "paav-ef1a-flex-asap6c-wpre_estimated" comes out as
    # "estimated", because it is simply the longest token left standing.
    "estimated", "confirmed", "verified", "unverified", "checked",
    "putative", "likely", "approx", "approximate", "partial", "complete",
    "cohen", "lab", "ape", "snapgene", "benchling", "expected", "reversed",
    "complete", "completed", "cropseq", "backup", "orig", "original",
    "multi", "sample", "otherwise", "fck", "fugw", "fuw", "sin", "min",
    "builder", "deisseroth", "boyden", "schnitzer", "synp", "fubi", "tra",
    "pcrprdt", "pcrprdt2", "prdt", "product", "amplicon", "digest",
    "ha", "flag", "myc", "his", "v5",
    # A named gap between two elements, which several files on the share
    # carry between the WPRE and the polyA.
    "spacer", "stuffer", "filler", "linker2", "pad",
}

#: A lab stock code: initials plus a number, like YD043 or LZF1605.
_STOCK_CODE = re.compile(r"^[A-Za-z]{1,4}[-_]?\d{1,5}[a-z]?$")


#: Vector names carry a serial ("pcDNA3a", "pAAV2"), so matching the whole
#: token is not enough.
_VECTOR_PREFIX = ("pcdna", "paav", "plenti", "pet", "pgex", "puc", "pbr",
                  "pcag", "pcmv", "pmax", "psp", "ptre")


def _is_vector_name(token: str) -> bool:
    t = token.lower()
    return any(t.startswith(v) and t[len(v):].lstrip("0123456789") in
               ("", "a", "b", "c", "d", "n", "x") for v in _VECTOR_PREFIX)


#: "pPZ021", "pP133" - a plasmid number, not a protein.
_PLASMID_CODE = re.compile(r"^p[A-Z]{1,3}\d{1,5}$")


def _is_stock_code(token: str, position: int) -> bool:
    """Whether a token is a stock number rather than a protein name.

    The pattern alone is not enough: "ASAP7y" is shaped exactly like a stock
    code.  A real one leads the plasmid name, or is written in the upper
    case the lab numbers its stocks in.
    """
    if _PLASMID_CODE.match(token):
        return True
    if not _STOCK_CODE.match(token):
        return False
    return position == 0 or token.upper() == token


def _names_an_element(token: str) -> bool:
    """Whether the vocabulary already knows this token as a non-CDS part.

    `_NAME_NOISE` is hand-kept, so it is always one word behind the share.
    "pcag-egfp-wpre-spacer-sv40dts" offered "sv40dts" as the name of an
    EGFP for exactly that reason: nobody had thought to write "sv40dts"
    down, and it was simply the longest token left standing.

    The vocabulary had thought about it, though - typing labels is what it
    does - so ask it rather than growing the list by one word per bug
    report.  Anything it reads as a promoter, a polyA, an enhancer or a
    resistance marker is a piece of the plasmid that is not the insert.
    "unknown" is the useful answer here: a protein name is precisely the
    thing the vocabulary has no rule for.
    """
    return classify_label(token) not in ("", "unknown", "CDS")


def name_from_plasmid(plasmid_name: str) -> str | None:
    """Pull the likely insert name out of a plasmid name.

    "YD043-pAAV-syn-iGluSnFR4f-NGR-WPRE" is mostly vector: a stock code, a
    serotype, a promoter, a trafficking motif and a regulatory element.  The
    one token that is none of those is the thing that was cloned.
    """
    tokens = [t for t in re.split(r"[-_\s.()\[\]{}#,;:+/\\]+",
                                 plasmid_name or "") if t]
    keep = [t for i, t in enumerate(tokens)
            if t.lower() not in _NAME_NOISE
            and not _is_vector_name(t)
            and not _is_stock_code(t, i)
            and not _names_an_element(t)
            and len(t) >= 3
            and not t.isdigit()]
    if not keep:
        return None
    # The insert name is normally the longest thing left, and ties go to
    # whichever came first.
    return max(keep, key=lambda t: (len(t), -tokens.index(t)))


class ProteinIndex:
    """One shared k-mer index over every catalogued protein.

    `protein_identity` builds an index per comparison, which is fine for a
    handful and hopeless for a thousand plasmids against fifteen hundred
    parts.  Inverting it once lets a query be scored against everything at
    the same time.
    """

    K = 8

    def __init__(self, conn: sqlite3.Connection):
        self.kmers: dict[str, set[int]] = {}
        self.parts: dict[int, tuple[str, int]] = {}
        rows = conn.execute(
            """SELECT id, name, seq FROM part
                WHERE feature_type = 'CDS' AND seq_length >= 150""")
        for r in rows:
            aa = translate(r["seq"]).rstrip("*")
            if len(aa) < self.K or "*" in aa[:-1]:
                continue
            self.parts[r["id"]] = (r["name"], len(aa))
            for i in range(len(aa) - self.K + 1):
                self.kmers.setdefault(aa[i:i + self.K], set()).add(r["id"])

    def best(self, aa: str) -> tuple[float, int, str] | None:
        """Closest catalogued protein as (identity, part id, name)."""
        if len(aa) < self.K:
            return None
        windows = [aa[i:i + self.K] for i in range(len(aa) - self.K + 1)]
        hits: dict[int, int] = {}
        for w in windows:
            for pid in self.kmers.get(w, ()):
                hits[pid] = hits.get(pid, 0) + 1
        if not hits:
            return None
        best_id, best_score = None, 0.0
        for pid, n in hits.items():
            _name, plen = self.parts[pid]
            # Scored against the LONGER of the two, unlike protein_identity.
            # That question is "are these the same part", and a cpGFP buried
            # inside a sensor answers it yes.  This question is "what is
            # this protein", where being a component of something is not
            # being it: iGluSnFR4f is 97% of a cEGFP and only 85% of the
            # iGluSnFr3 it actually is.
            denom = max(len(windows), max(plen - self.K + 1, 1))
            score = n / denom
            if score > best_score:
                best_id, best_score = pid, score
        if best_id is None:
            return None
        return min(best_score, 1.0), best_id, self.parts[best_id][0]


def _overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def candidate_orfs(seq: str, circular: bool, annotations: list[dict]) -> list[dict]:
    """ORFs in this sequence that nothing already accounts for."""
    n = len(seq)
    coding = [a for a in annotations if a["feature_type"] in _CODING]
    structural = [a for a in annotations if a["feature_type"] in _STRUCTURAL]

    kept: list[dict] = []
    for o in find_orfs(seq, min_aa=MIN_ORF_AA, circular=circular):
        # find_orfs reports wrapped ORFs in coordinates that run off the end
        # of the molecule.  Those are left alone rather than guessed at.
        if o["start"] < 0 or o["end"] > n:
            continue
        span = o["end"] - o["start"]
        if span <= 0:
            continue
        if sum(_overlap(o["start"], o["end"], a["start"], a["end"])
               for a in coding) / span >= COVERED:
            continue
        # Either way round.  Swallowing a promoter means the reading
        # frame ran straight through it; sitting *inside* one means the
        # frame is the promoter, read in a frame nobody transcribes.
        # Only the first was tested, which is why 108 plasmids carry a
        # 202 aa "CDS" that is the CAG promoter, and why one of them got
        # into the catalogue as "Voltron" and then named the same
        # stretch in thirty-four others.
        if any(ov > 0.5 * (a["end"] - a["start"]) or ov > 0.5 * span
               for a in structural
               for ov in (_overlap(o["start"], o["end"],
                                   a["start"], a["end"]),)):
            continue
        kept.append(o)

    # Nested ORFs share a stop codon; only the longest reading of each is a
    # candidate, and overlapping candidates are alternative readings of the
    # same stretch of DNA rather than two genes.
    chosen: list[dict] = []
    for o in sorted(kept, key=lambda d: -d["aa_length"]):
        if any(_overlap(o["start"], o["end"], c["start"], c["end"])
               > 0.5 * min(o["end"] - o["start"], c["end"] - c["start"])
               for c in chosen):
            continue
        chosen.append(o)
    return sorted(chosen, key=lambda d: d["start"])


def payload_span(annotations: list[dict]) -> tuple[int, int] | None:
    """The stretch an expression plasmid's insert sits in.

    Between the promoter that drives it and whatever terminates it.  The
    resistance promoter is not a candidate: it drives the backbone.
    """
    proms = [a for a in annotations if a["feature_type"] == "promoter"
             and "ampr" not in (a.get("name") or "").lower()]
    ends = [a for a in annotations
            if a["feature_type"] in ("WPRE", "polyA")]
    if not proms or not ends:
        return None
    start = min(a["end"] for a in proms)
    after = [a["start"] for a in ends if a["start"] > start]
    if not after:
        return None
    stop = min(after)
    return (start, stop) if stop - start >= 300 else None


def _names_something_else(token: str, orf: dict, anns: list) -> bool:
    """Is this token already the name of a different feature here?

    "pcag-egfp-wpre-spacer-sv40dts" has five things in its name and the
    longest leftover token is sv40dts - which is a real element of the
    plasmid, 144 bp of SV40 nuclear targeting sequence at 3983-4127,
    and not the 993 bp reading frame at 1533. Naming the ORF after it
    produced an annotation that claimed to be something the plasmid
    already had somewhere else, at seven times the length.

    If the file already says where that name lives, and it is not
    here, then it is not this ORF's name.
    """
    want = re.sub(r"[^a-z0-9]", "", (token or "").lower())
    if not want:
        return False
    for a in anns or []:
        other = re.sub(r"[^a-z0-9]", "", (a.get("name") or "").lower())
        if other != want:
            continue
        if _overlap(orf["start"], orf["end"],
                    a.get("start", 0), a.get("end", 0)) > 0:
            return False       # it is this one, under another record
        return True
    return False


def describe(orf: dict, index: ProteinIndex, plasmid_name: str,
             is_payload: bool = True, anns: list | None = None) -> dict:
    """Name an inferred ORF, and say why it got that name.

    `is_payload` says whether this ORF is the thing the plasmid was built to
    express.  Only then does the plasmid name describe it - "YD043-pAAV-syn-
    iGluSnFR4f" says nothing about the beta-lactamase in its backbone, and
    naming a backbone ORF after the payload would be worse than leaving it
    unnamed.
    """
    aa = (orf.get("protein") or "").rstrip("*")
    match = index.best(aa) if aa else None
    from_name = name_from_plasmid(plasmid_name) if is_payload else None
    # A plasmid called "...target..." names its insert no better than a
    # catalogue entry called "target" does.
    if from_name and _GENERIC_NAME.match(from_name):
        from_name = None
    # ...and a token that already names a different feature of this
    # plasmid names that feature, not this reading frame.
    if from_name and _names_something_else(from_name, orf, anns or []):
        from_name = None

    if match and match[0] >= IDENTICAL:
        ident, pid, pname = match
        if _useful(pname):
            return {"name": _clean_name(pname), "confidence": 0.8,
                    "note": f"ORF inferred; {ident*100:.0f}% identical to "
                            f"{_clean_name(pname)} (part #{pid})"}
        if from_name:
            # The match is certain but tells us nothing, so the name comes
            # from the plasmid alone and is only as good as that.
            return {"name": from_name, "confidence": 0.5,
                    "note": f"ORF inferred; named from the plasmid "
                            f"(its {ident*100:.0f}% match, part #{pid}, is "
                            f"itself unnamed)"}
    if match and match[0] >= RELATED:
        ident, pid, pname = match
        # A close but imperfect match is a variant.  The plasmid name is the
        # better source for which variant, when it offers one.
        name = from_name or (f"{_clean_name(pname)}-like"
                             if _useful(pname) else None)
        if name:
            return {"name": name, "confidence": 0.65,
                    "note": f"ORF inferred; {ident*100:.0f}% identical to "
                            f"{_clean_name(pname)} (part #{pid})"}
    if from_name:
        return {"name": from_name, "confidence": 0.5,
                "note": "ORF inferred; named from the plasmid, no close "
                        "protein in the catalogue"}
    if match:
        ident, pid, pname = match
        return {"name": f"unnamed CDS ({len(aa)} aa)", "confidence": 0.4,
                "note": f"ORF inferred; closest catalogue protein is "
                        f"{_clean_name(pname)} (part #{pid}) at "
                        f"{ident*100:.0f}%, which names nothing"}
    return {"name": f"unnamed CDS ({len(aa)} aa)", "confidence": 0.4,
            "note": "ORF inferred; no close protein in the catalogue and "
                    "nothing usable in the plasmid name"}


def infer_all(conn: sqlite3.Connection, limit: int | None = None,
              refresh: bool = False, only: list[int] | None = None,
              progress=None) -> dict:
    """Annotate unannotated payload ORFs across the inventory."""
    from .importer import extract_parts, _StoredSeq

    if refresh:
        # The annotations go, but the parts they minted used to stay: an
        # occurrence is keyed on coordinates, not on the annotation that
        # argued for it.  So a wrong guess survived the refresh meant to
        # undo it, kept its place in the protein index, and went on
        # lending its name to the next run.  Retire any occurrence that
        # nothing but an inferred ORF vouches for, then drop the parts
        # that leaves with nothing pointing at them.
        conn.execute("""CREATE TEMP TABLE _orf_orphans AS
                        SELECT DISTINCT part_id FROM part_occurrence""")
        conn.execute(
            """DELETE FROM part_occurrence
                WHERE EXISTS (SELECT 1 FROM annotation a
                               WHERE a.plasmid_id = part_occurrence.plasmid_id
                                 AND a.start = part_occurrence.start
                                 AND a.end = part_occurrence.end
                                 AND a.strand = part_occurrence.strand
                                 AND a.evidence = 'orf')
                  AND NOT EXISTS (SELECT 1 FROM annotation a
                               WHERE a.plasmid_id = part_occurrence.plasmid_id
                                 AND a.start = part_occurrence.start
                                 AND a.end = part_occurrence.end
                                 AND a.strand = part_occurrence.strand
                                 AND a.evidence <> 'orf')""")
        # Parts already unreferenced before this ran are somebody else's
        # business; only the ones this pass stranded are removed.
        conn.execute(
            """DELETE FROM part
                WHERE id IN (SELECT id FROM part WHERE id NOT IN
                                 (SELECT part_id FROM part_occurrence))
                  AND id IN (SELECT part_id FROM _orf_orphans)""")
        conn.execute("DROP TABLE _orf_orphans")
        conn.execute("DELETE FROM annotation WHERE evidence = 'orf'")

    index = ProteinIndex(conn)

    if only is not None:
        ids = list(only)
    else:
        ids = [r[0] for r in conn.execute(
            """SELECT id FROM plasmid
                WHERE seq_length BETWEEN 1000 AND 60000
                ORDER BY id""")]
    if limit is not None:
        ids = ids[:limit]

    examined = annotated = parts_made = skipped = 0
    named_by = {"catalogue": 0, "variant": 0, "plasmid_name": 0, "none": 0}

    for i, pid in enumerate(ids):
        if progress and i and i % 200 == 0:
            progress(i, len(ids), annotated, skipped)
        p = conn.execute(
            "SELECT id, name, seq, circular FROM plasmid WHERE id = ?",
            (pid,)).fetchone()
        if not p or not p["seq"]:
            continue
        anns = [dict(r) for r in conn.execute(
            """SELECT feature_type, name, start, end, strand, evidence
                 FROM annotation WHERE plasmid_id = ?""", (pid,))]
        if not refresh and any(a["evidence"] == "orf" for a in anns):
            # Already done on an earlier run; --refresh redoes them.
            skipped += 1
            continue
        examined += 1

        cands = candidate_orfs(p["seq"], bool(p["circular"]), anns)
        span = payload_span(anns)
        # A file that already annotates a coding sequence has said what
        # this plasmid expresses.  Whatever reading frame is left over is
        # not the payload, and must not be given the payload's name -
        # that is how "CAG-Voltron-ST" came to call a frame inside its
        # own chimeric intron "Voltron", and then, through the catalogue,
        # called the same stretch of CAG "Voltron" in thirteen more.
        spoken_for = any(a["feature_type"] == "CDS"
                         and a["end"] - a["start"] >= MIN_ORF_AA * 3
                         for a in anns)
        if span is None and cands and not spoken_for:
            # No promoter/terminator pair to go on, so the longest reading
            # frame is the best guess at what this plasmid is for.
            longest = max(cands, key=lambda o: o["aa_length"])
            payload_ids = {id(longest)}
        else:
            inside = [o for o in cands if span
                      and _overlap(o["start"], o["end"], *span)
                      > 0.6 * (o["end"] - o["start"])]
            # Several reading frames can fall inside the payload span, but
            # the plasmid is named after one insert.  The longest is it;
            # the rest are vector the span happened to include.
            payload_ids = ({id(max(inside, key=lambda o: o["aa_length"]))}
                           if inside else set())

        for orf in cands:
            is_payload = id(orf) in payload_ids
            d = describe(orf, index, p["name"], is_payload, anns)
            if "identical to" in d["note"]:
                named_by["catalogue" if d["confidence"] >= 0.8 else "variant"] += 1
            elif d["confidence"] >= 0.5:
                named_by["plasmid_name"] += 1
            else:
                named_by["none"] += 1
            conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        original_label, original_key, start, end, strand,
                        evidence, confidence, qualifiers)
                   VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
                (pid, "CDS", d["name"], "", "orf:inferred",
                 orf["start"], orf["end"], orf["strand"], "orf",
                 d["confidence"], json.dumps({"note": d["note"]})))
            annotated += 1
            # An ORF reaches the parts catalogue only when a protein in it
            # corroborates the name.  A guess read off the plasmid alone is
            # right often enough to record against that plasmid and wrong
            # often enough - "sequence", "TBD" - to keep out of a catalogue
            # everyone searches.
            if is_payload and d["confidence"] >= CATALOGUE_CONFIDENCE:
                parts_made += extract_parts(
                    conn, pid, _StoredSeq(p["seq"], bool(p["circular"])),
                    [{"feature_type": "CDS", "name": d["name"],
                      "start": orf["start"], "end": orf["end"],
                      "strand": orf["strand"], "confidence": d["confidence"]}])

        if examined % 200 == 0:
            conn.commit()

    # A part is shared by every plasmid carrying that sequence, and the
    # first plasmid to reach it is not necessarily the one that names it
    # best: the same ORF can arrive as "iGluSnFr3-like" from a plasmid whose
    # name says nothing, and as "iGluSnFR4f" from one whose name says
    # exactly what it is.
    upgraded = 0
    rows = conn.execute(
        """SELECT o.part_id, p.name AS part_name, a.name AS ann_name,
                  a.confidence
             FROM part_occurrence o
             JOIN part p ON p.id = o.part_id
             JOIN annotation a ON a.plasmid_id = o.plasmid_id
                              AND a.start = o.start AND a.end = o.end
                              AND a.strand = o.strand
            WHERE a.evidence = 'orf'""").fetchall()
    best: dict[int, tuple[float, str, str]] = {}
    for r in rows:
        key = (r["confidence"], not r["ann_name"].endswith("-like"))
        rank = r["confidence"] + (0.01 if not r["ann_name"].endswith("-like") else 0)
        cur = best.get(r["part_id"])
        if cur is None or rank > cur[0]:
            best[r["part_id"]] = (rank, r["ann_name"], r["part_name"])
    for part_id, (_rank, name, current) in best.items():
        if name != current:
            conn.execute("UPDATE part SET name = ? WHERE id = ?",
                         (name, part_id))
            upgraded += 1

    conn.commit()
    return {
        "plasmids_examined": examined,
        "plasmids_already_done": skipped,
        "parts_renamed": upgraded,
        "cds_inferred": annotated,
        "parts_created": parts_made,
        "named_from": named_by,
    }
