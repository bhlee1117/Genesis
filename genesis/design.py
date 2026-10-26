"""Orchestration: specification in, ranked laboratory plan out.

This is the single entry point that walks Stages 4 -> 10 in order and
returns one JSON-serialisable plan.
"""

from __future__ import annotations

import json
import sqlite3

from .config import Config
from .seqtools import find_all, revcomp, slice_region
from .construct import (VECTOR_TYPES, ConstructSpec, assemble,
                        assemble_into_backbone, codon_markers,
                        expand_dio, resolve)
from .planner import (group_contiguous, plan_gblock_route,
                      plan_hifi_route, plan_restriction_route,
                      plan_swap_hifi_route, plan_swap_hybrid_route,
                      plan_swap_restriction_route,
                      rank_routes, failure_history)
from .swap import (analyse_swap, architecture_candidates, build_scheme,
                   protein_of, rank_parents)
from .vectors import classify_vector, describe, missing_elements
from .vocab import canonical_feature_type
from .verify import verify_plan
from .validate import Issue, summarize, validate_construct


#: Shortest sequence that could plausibly act as a cloning vector. The
#: inventory contains many short fragments whose names also match a backbone
#: search term, and silently picking a 23 bp fragment as the vector produces
#: a nonsense plan.
MIN_BACKBONE_LENGTH = 1500


def _article(phrase: str) -> str:
    """'an AAV transfer plasmid' / 'a lentiviral transfer plasmid'."""
    return ("an " if phrase[:1].upper() in "AEIOU" else "a ") + phrase


def plasmid_family(conn: sqlite3.Connection, plasmid_id: int,
                   name: str = "") -> str:
    ann = [dict(r) for r in conn.execute(
        "SELECT feature_type, name FROM annotation WHERE plasmid_id = ?",
        (plasmid_id,))]
    if not name:
        row = conn.execute("SELECT name FROM plasmid WHERE id = ?",
                           (plasmid_id,)).fetchone()
        name = row["name"] if row else ""
    return classify_vector(ann, name)["family"]


#: Mutual 40-mer coverage at which two plasmids are the same construct.
#: Both directions have to clear it: one-way coverage is high whenever
#: the design is a subset of something bigger.
SAME_CONSTRUCT = 0.98
#: Below that but still close enough to be worth a look before cloning.
NEARLY_THE_SAME = 0.93


def already_made(conn: sqlite3.Connection, construct, resolved) -> list:
    """Plasmids the lab already has that are this construct.

    Cloning something that is already in the freezer is a fortnight
    spent arriving where you started, and the inventory is the only
    thing that knows. The test is mutual coverage in both directions:
    the design is inside the candidate and the candidate is inside the
    design. One direction on its own is not enough - a design is almost
    entirely contained in any larger plasmid built around the same
    backbone, which is how a Voltron2-mGluR2CT construct looks like a
    plain Voltron2 one from the inside.

    What comes back is ranked rather than binary, because the gap
    between the same construct and a close relative is a real gap and
    not always a wide one. Anything at or above SAME_CONSTRUCT is
    reported as already made; the band below it is reported as worth a
    look, which is what it is.
    """
    from .importer import seq_checksum
    from .swap import sequence_coverage

    designed = (construct.seq or "").upper()
    if len(designed) < 500:
        return []
    circular = bool(getattr(construct, "circular", True))
    want_ck = seq_checksum(designed, circular)
    doubled_design = designed + designed

    hits = []
    for row in conn.execute(
            "SELECT id, name, owner, seq, seq_length, circular, verified "
            "FROM plasmid WHERE seq_length BETWEEN ? AND ? AND seq IS NOT NULL",
            (int(len(designed) * 0.93), int(len(designed) * 1.07))):
        host = (row["seq"] or "").upper()
        if len(host) < 500:
            continue
        if seq_checksum(host, bool(row["circular"])) == want_ck:
            hits.append({"id": row["id"], "name": row["name"],
                         "owner": row["owner"], "length": row["seq_length"],
                         "verified": row["verified"], "exact": True,
                         "coverage": 1.0,
                         "why": "the same molecule, base for base"})
            continue
        forward = sequence_coverage(designed, host + host)
        if forward < NEARLY_THE_SAME:
            continue
        back = sequence_coverage(host, doubled_design)
        score = min(forward, back)
        if score < NEARLY_THE_SAME:
            continue
        hits.append({
            "id": row["id"], "name": row["name"], "owner": row["owner"],
            "length": row["seq_length"], "verified": row["verified"],
            "exact": False, "coverage": round(score, 3),
            "why": (f"{score * 100:.0f}% of this design is in it and "
                    f"{score * 100:.0f}% of it is in this design"),
        })
    hits.sort(key=lambda h: (not h["exact"], -h["coverage"], h["name"]))
    return hits[:6]


def _add_checkpoints(conn: sqlite3.Connection, plan: dict) -> None:
    """Put a QC step after every step that makes a molecule."""
    from .checkpoints import for_route
    from .operations import Step, renumber
    from .simulate import _templates

    for route in plan.get("routes") or []:
        steps = route.get("steps") or []
        if not steps or any(s.get("operation") == "QC_GEL" for s in steps):
            continue
        try:
            resolved = _templates(conn, plan, route)
        except Exception:
            resolved = {}
        lengths = {name: len(entry[0]) for name, entry in resolved.items()}
        # What the gel should show after a digest: the simulation has
        # already cut these molecules, so the sizes are known rather
        # than guessed.
        bands = {}
        from .enzymes import CATALOG
        from .figure import step_enzymes
        from .simulate import digest as cut
        for i, step in enumerate(steps):
            if step.get("operation") != "RE_DIGEST":
                continue
            enzymes = [e for e in step_enzymes(step) if e in CATALOG]
            host = next((x for x in step.get("inputs") or []
                         if x in resolved), None)
            if not (enzymes and host):
                continue
            seq, circular, _name, _tid = resolved[host]
            try:
                pieces = cut(seq, circular, enzymes)
            except Exception:
                continue
            if len(pieces) >= 2:
                bands[i] = sorted((len(x.seq) for x in pieces), reverse=True)

        marks = for_route(route, lengths, bands)
        if not marks:
            continue
        out, by_index = [], {i: c for i, c in marks}
        for i, step in enumerate(steps):
            out.append(step)
            check = by_index.get(i)
            if check is None:
                continue
            # A run of PCRs already shares one checkpoint, attached to
            # the last of them, so there is nothing to merge here.
            out.append(Step(0, "QC_GEL", check.title,
                            output=check.operation == "PCR"
                            and "a checked amplicon" or "a checked reaction",
                            details={"checkpoint": check.to_dict()}))
        route["steps"] = [s if isinstance(s, dict) else s.to_dict()
                          for s in _renumber_mixed(out)]
        # The route's counters were computed when it was serialised,
        # which was before these steps existed; recount rather than
        # print a step list of ten under a heading that says eight.
        _recount(route)


def _recount(route: dict) -> None:
    """Re-derive a route's tallies after its steps have changed."""
    steps = route.get("steps") or []
    required = [s for s in steps if not s.get("optional")]
    optional = [s for s in steps if s.get("optional")]
    route["n_steps"] = len(steps)
    route["n_optional"] = len(optional)
    route["n_pcr"] = sum(1 for s in steps if s.get("operation") == "PCR")
    route["n_assembly"] = sum(
        1 for s in steps if s.get("operation") in
        ("GIBSON_ASSEMBLY", "HI_FI_ASSEMBLY", "LIGATION"))
    route["hands_on_minutes"] = sum(s.get("hands_on", 0) for s in required)
    route["turnaround_minutes"] = sum(s.get("minutes", 0) for s in required)
    route["optional_minutes"] = sum(s.get("hands_on", 0) for s in optional)


def _renumber_mixed(steps: list) -> list:
    """Renumber a list that is part Step objects and part dicts."""
    for i, s in enumerate(steps, start=1):
        if isinstance(s, dict):
            s["index"] = i
        else:
            s.index = i
    return steps


def required_architecture(spec: ConstructSpec) -> set:
    """Structural features a recipient must already have for this design.

    A Cre-lox construct has to start from a plasmid that already carries lox
    sites, and an AAV from one that already carries ITRs. Sharing a few
    coding sequences is not enough: a plasmid with the right genes and the
    wrong architecture cannot be edited into this construct by swapping one
    region.
    """
    needs: set = set()
    types = {canonical_feature_type(
        m["type"] if isinstance(m, dict) else m.type) for m in spec.modules}
    if "Cre-lox" in types:
        needs.add("lox_site")
    if spec.vector == "AAV":
        needs.add("AAV_ITR")
    return needs


def _has_features(conn: sqlite3.Connection, plasmid_id: int,
                  needs: set) -> bool:
    if not needs:
        return True
    present = {r["feature_type"] for r in conn.execute(
        "SELECT DISTINCT feature_type FROM annotation WHERE plasmid_id = ?",
        (plasmid_id,))}
    return needs <= present


def find_backbone_for(conn: sqlite3.Connection, family: str, needs: set,
                      limit: int = 60, cfg: Config | None = None) -> dict | None:
    """Pick a backbone of the right family when none was named.

    Preferred: a plasmid of the requested vector family that already has the
    architecture this design needs, is sequence-verified, and is no larger
    than it has to be - a smaller backbone is easier to amplify and leaves
    more room under the AAV packaging limit.
    """
    # A named default first: the lab's standard backbone for this family is
    # a better starting point than whatever scores highest, because it is the
    # one people already have, have sequenced, and know works.
    preferred = ""
    if cfg is not None:
        preferred = (cfg.default_aav_backbone if family == "AAV"
                     else cfg.default_lenti_backbone if family == "lentiviral"
                     else "")
    if preferred:
        row = conn.execute(
            """SELECT * FROM plasmid WHERE name LIKE ? COLLATE NOCASE
               AND circular = 1 ORDER BY seq_length ASC LIMIT 1""",
            (f"{preferred}%",)).fetchone()
        if row is not None and _has_features(conn, row["id"], needs):
            return dict(row)

    rows = conn.execute(
        """SELECT * FROM plasmid
           WHERE seq_length >= ? AND circular = 1 AND available_template = 1
             AND EXISTS (SELECT 1 FROM annotation a
                         WHERE a.plasmid_id = plasmid.id
                           AND a.feature_type = 'origin')
             AND EXISTS (SELECT 1 FROM annotation a
                         WHERE a.plasmid_id = plasmid.id
                           AND a.feature_type = 'resistance')
           ORDER BY (verified != 'unverified') DESC, seq_length ASC
           LIMIT 400""", (MIN_BACKBONE_LENGTH,)).fetchall()
    fallback = None
    for row in rows:
        if family and plasmid_family(conn, row["id"], row["name"]) != family:
            continue
        if _has_features(conn, row["id"], needs):
            return dict(row)
        fallback = fallback or row
    return dict(fallback) if fallback else None


def _load_backbone(conn: sqlite3.Connection, spec: ConstructSpec,
                   cfg: Config | None = None,
                   preferred_id: int | None = None) -> dict | None:
    row = None
    if spec.backbone_plasmid_id:
        # An explicit id is the user's choice; honour it as given.
        row = conn.execute("SELECT * FROM plasmid WHERE id = ?",
                           (spec.backbone_plasmid_id,)).fetchone()
    elif spec.recipient_plasmid_id:
        # Choosing a recipient is choosing what to start from, and the
        # plasmid being edited is the vector the construct ends up in.
        # Letting the two differ gave one document two molecules: a swap
        # that edited the chosen plasmid, printed beside an idealised
        # construct assembled into a different backbone, with two
        # different sizes for the same product.
        row = conn.execute("SELECT * FROM plasmid WHERE id = ?",
                           (spec.recipient_plasmid_id,)).fetchone()
    elif spec.backbone:
        # Name search: a vector must be big enough to be one, and a circular,
        # verified plasmid carrying an origin beats a bare fragment.
        rows = conn.execute(
            """SELECT * FROM plasmid
               WHERE name LIKE ? COLLATE NOCASE AND seq_length >= ?
               ORDER BY
                 circular DESC,
                 (verified != 'unverified') DESC,
                 (SELECT COUNT(*) FROM annotation a
                  WHERE a.plasmid_id = plasmid.id
                    AND a.feature_type IN ('origin', 'resistance')) DESC,
                 seq_length DESC
               LIMIT 25""",
            (f"%{spec.backbone}%", MIN_BACKBONE_LENGTH)).fetchall()
        # When a vector family was asked for, prefer a backbone of that
        # family. If the name search finds none, fall back to the best match
        # anyway: "you asked for AAV but this backbone is lentiviral" is a
        # far more useful answer than silently having no backbone and a
        # cascade of missing-element errors that never names the cause.
        row = None
        if spec.vector:
            for candidate in rows:
                if plasmid_family(conn, candidate["id"],
                                  candidate["name"]) == spec.vector:
                    row = candidate
                    break
        if row is None:
            row = rows[0] if rows else None
    elif preferred_id is not None:
        # Nobody named a backbone, so start from whatever already carries the
        # longest leading run of the architecture asked for.  Picking by
        # family alone is what put a plain hSyn payload inside a DIO
        # plasmid's lox cassette.
        row = conn.execute("SELECT * FROM plasmid WHERE id = ?",
                           (preferred_id,)).fetchone()
        if row is None and spec.vector:
            return find_backbone_for(conn, spec.vector,
                                     required_architecture(spec), cfg=cfg)
    elif spec.vector:
        # No backbone named, and nothing in the inventory matches the
        # architecture: choose by family rather than building a payload with
        # no vector around it.
        return find_backbone_for(conn, spec.vector,
                                 required_architecture(spec), cfg=cfg)
    return dict(row) if row else None


#: Module types that can sit in one open reading frame together, which is
#: what a 2A peptide is for: one promoter, one transcript, several proteins.
FUSED_TYPES = {"CDS", "p2A", "linker", "tag", "trafficking_motif", "IRES"}


def _prefer_fused_source(conn: sqlite3.Connection,
                         resolved: list) -> list[str]:
    """Take a 2A fusion from a plasmid that already carries it whole.

    Asking for X - P2A - Y names three modules, and the planner will
    happily amplify three fragments and join them with two junctions it has
    to design, check and screen. But these constructs are built from each
    other: the lab has plenty of plasmids where X-P2A-Y is already sitting
    in one piece, because somebody made it that way, and a fusion taken
    whole is one reaction with no internal junction to get wrong.

    So where a run of modules that share a reading frame exists intact in
    one plasmid, every module in it is sourced from there - which is all it
    takes, because fragments from the same plasmid at adjoining
    coordinates are already merged into one amplification.
    """
    notes: list[str] = []
    runs: list[list] = []
    current: list = []
    for rm in resolved:
        if rm.seq and rm.feature_type in FUSED_TYPES:
            current.append(rm)
        else:
            if len(current) > 1:
                runs.append(current)
            current = []
    if len(current) > 1:
        runs.append(current)

    for run in runs:
        if not any(rm.feature_type == "p2A" for rm in run):
            continue
        # The longest intact stretch is worth most, so start from the whole
        # run and give up a module at a time from the end.
        for size in range(len(run), 1, -1):
            for start in range(0, len(run) - size + 1):
                piece = run[start:start + size]
                if not any(rm.feature_type == "p2A" for rm in piece):
                    continue
                found = _locate_fusion(conn, piece)
                if found is None:
                    continue
                plasmid, places = found
                for rm, where in zip(piece, places):
                    rm.sources = [{**where, "id": plasmid["id"],
                                   "name": plasmid["name"],
                                   "owner": plasmid["owner"] or "",
                                   "verified": plasmid["verified"],
                                   "available_template": 1,
                                   "circular": plasmid["circular"],
                                   "seq_length": plasmid["seq_length"]}] + [
                        x for x in rm.sources if x["id"] != plasmid["id"]]
                    rm.pinned_source = True
                notes.append(
                    " + ".join(rm.part_name for rm in piece)
                    + f" is already one piece in {plasmid['name']}, so it is "
                      f"amplified as one rather than joined from "
                      f"{len(piece)}")
                break
            else:
                continue
            break
    return notes


def _locate_fusion(conn: sqlite3.Connection, piece: list):
    """Find a plasmid carrying this run of modules end to end.

    An exact match on the whole run is the easy case and the rare one. The
    usual case is that the plasmid joins the same two proteins through its
    own 2A peptide, which is not the one a generic "P2A" resolves to -
    P2A, T2A and E2A all do the job and the lab's plasmids use all three.
    So the pair is located and whatever sits between them is read off that
    plasmid, after checking that it really is a 2A: in frame, the right
    length, and no stop codon in the middle of it.
    """
    from .seqtools import translate

    joined = "".join(rm.seq for rm in piece).upper()
    if len(joined) < 60:
        return None

    # Only a pair joined through a single 2A is handled; a longer run is
    # tried as its pairs by the caller.
    pair = (len(piece) == 3 and piece[1].feature_type == "p2A"
            and piece[0].seq and piece[2].seq)

    # A plasmid can only hold the run whole if it holds every module of
    # it, so those come first. Ranking only by the first module's sources
    # pushed the plasmids that actually carry the fusion off the end of
    # the list.
    shared = set.intersection(*[{s["id"] for s in rm.sources}
                                for rm in piece]) if piece else set()
    order = ([s for rm in piece for s in rm.sources if s["id"] in shared]
             + piece[0].sources[:40])

    seen: set = set()
    for src in order:
        if src["id"] in seen:
            continue
        seen.add(src["id"])
        if len(seen) > 60:
            break
        row = conn.execute(
            "SELECT id, name, owner, verified, circular, seq_length, seq "
            "FROM plasmid WHERE id = ?", (src["id"],)).fetchone()
        if row is None or not row["seq"]:
            continue
        host = row["seq"].upper()

        for strand, text in ((1, host), (-1, revcomp(host))):
            at = text.find(joined)
            if at >= 0:
                return dict(row), _places(piece, at, len(text), strand)

            # The same run, spelled differently. Take this plasmid's own
            # bases for each module, so what the primers amplify and what
            # the construct claims to contain are the same DNA.
            spans = _fused_in_frame(text, piece, translate)
            if spans is not None:
                for rm, (lo, hi) in zip(piece, spans):
                    bases = text[lo:hi]
                    if bases.upper() != (rm.seq or "").upper():
                        rm.warnings.append(
                            f"taken as the {len(bases)} bp {rm.part_name} "
                            f"already in {row['name']}, which encodes the "
                            f"same protein in different codons or without "
                            f"its start codon, so the run comes out of one "
                            f"reaction")
                    rm.seq = bases if strand == 1 else revcomp(bases)
                total = spans[-1][1] - spans[0][0]
                return dict(row), _places(piece, spans[0][0], len(text),
                                          strand, total)
            if not pair:
                continue
            left, right = piece[0].seq.upper(), piece[2].seq.upper()
            lo = text.find(left)
            if lo < 0:
                continue
            after = lo + len(left)
            # The window bounds where the next protein may *start*; find
            # needs room for the whole of it, or a 762 bp segment can never
            # be found inside a 120 bp window.
            gap_end = text.find(right, after, after + MAX_2A + len(right))
            if gap_end < 0:
                continue
            spacer = text[after:gap_end]
            if not _is_2a(spacer, translate):
                continue
            # Take this plasmid's own 2A, so the three pieces really are
            # one stretch of its DNA and come out of one reaction.
            piece[1].seq = spacer if strand == 1 else revcomp(spacer)
            kind = ("2A peptide" if "NPG" in translate(spacer).upper()
                    else "linker")
            piece[1].warnings.append(
                f"taken as the {len(spacer)} bp {kind} already joining "
                f"{piece[0].part_name} to {piece[2].part_name} in "
                f"{row['name']}, so the three come out of one reaction")
            total = len(left) + len(spacer) + len(right)
            return dict(row), _places(piece, lo, len(text), strand, total)
    return None


#: How far past one protein the next may start and still be a 2A fusion.
MAX_2A = 120


def _fused_in_frame(text: str, piece: list, translate) -> list | None:
    """Find this run of modules end to end in one reading frame of `text`.

    Asked at the protein level rather than on the bases, because the same
    gene is in the inventory in more than one spelling: EGFP with and
    without its ATG, and codon-optimised copies that share almost no DNA
    with the row a request resolves to. PP078 carries P2A fused straight
    to EGFP, and matching its bases against the resolved EGFP fails, which
    is what sent the planner off to amplify a 66 bp P2A on its own.

    The modules must be flush against each other, so the stretch of the
    plasmid between the first and the last is exactly the construct asked
    for and nothing is quietly left out of it.
    """
    wanted = []
    for rm in piece:
        protein = protein_of(rm.seq or "")
        if len(protein) < 6:
            return None
        wanted.append(protein[1:] if protein[:1] == "M" else protein)

    for frame in range(3):
        protein = translate(text, frame=frame)
        at = protein.find(wanted[0])
        while at >= 0:
            spans, cursor, ok = [(at, at + len(wanted[0]))], at + len(wanted[0]), True
            for want in wanted[1:]:
                if protein[cursor:cursor + len(want)] != want:
                    ok = False
                    break
                spans.append((cursor, cursor + len(want)))
                cursor += len(want)
            if ok:
                return [(frame + 3 * lo, frame + 3 * hi) for lo, hi in spans]
            at = protein.find(wanted[0], at + 1)
    return None


def _is_2a(spacer: str, translate) -> bool:
    """Whether this stretch could be the 2A peptide joining two proteins."""
    if not 45 <= len(spacer) <= MAX_2A or len(spacer) % 3:
        return False
    protein = translate(spacer)
    return "*" not in protein


def _places(piece: list, at: int, text_length: int, strand: int,
            total: int | None = None) -> list[dict]:
    """Where each module of a located run sits, in the plasmid's own frame."""
    total = total if total is not None else sum(len(rm.seq) for rm in piece)
    places, offset = [], 0
    for rm in piece:
        if strand == 1:
            places.append({"start": at + offset,
                           "end": at + offset + len(rm.seq), "strand": 1})
        else:
            # The search ran on the other strand, so map back: what starts
            # `offset` into the match there ends the same distance from its
            # far end here.
            end = text_length - (at + offset)
            places.append({"start": end - len(rm.seq), "end": end,
                           "strand": -1})
        offset += len(rm.seq)
    return places


def _locate_in(host: str, want: str, circular: bool) -> tuple[int, int]:
    """Where `want` sits in `host`, on either strand, wrap included.

    A plasmid is a circle and its sequence file has to start somewhere, so
    a part can straddle the join between the last base and the first. A
    plain substring search reads that file as a line and misses it, which
    would have Genesis tell somebody the plasmid they named does not
    contain the insert they can see in it.

    The end is returned unreduced - past the end of the sequence when the
    match wraps - which is the convention `slice_region` and the
    occurrence table already use.
    """
    hits = find_all(host, want, circular=circular)
    if hits:
        return hits[0], 1
    hits = find_all(host, revcomp(want), circular=circular)
    if hits:
        return hits[0], -1
    return -1, 0


def _prefer_backbone_as_source(resolved: list, backbone: dict) -> list[str]:
    """Take a part from the recipient when the recipient already has it.

    Modules are resolved against the whole catalogue before any backbone is
    chosen, so "hSyn" binds to whichever part the search ranks first and
    carries that part's home plasmid with it.  Nothing afterwards asks
    whether the recipient has the same bases.  It usually does - a
    pAAV-hSyn vector is being used precisely because it has hSyn - and
    amplifying them from an unrelated, often unverified plasmid instead
    adds a template, a reaction and two junctions, for a molecule that is
    identical either way.

    Only an exact match counts.  A promoter that merely resembles the
    recipient's is a different promoter, and swapping it is the whole point
    of naming it.
    """
    host = (backbone.get("seq") or "").upper()
    if not host:
        return []
    moved: list[str] = []
    for rm in resolved:
        if not rm.seq or len(rm.seq) < 40:
            continue
        if getattr(rm, "pinned_source", False):
            continue          # somebody said where this comes from
        if rm.sources and rm.sources[0]["id"] == backbone["id"]:
            continue
        wanted = rm.seq.upper()
        at, strand = _locate_in(host, wanted,
                                bool(backbone.get("circular", 1)))
        if at < 0:
            continue
        rm.sources = [{
            "id": backbone["id"], "name": backbone["name"],
            "owner": backbone.get("owner", ""),
            "verified": backbone.get("verified", "unverified"),
            "available_template": backbone.get("available_template", 1),
            "circular": backbone.get("circular", 1),
            "seq_length": backbone.get("seq_length") or len(host),
            "start": at, "end": at + len(wanted), "strand": strand,
        }] + [x for x in rm.sources if x["id"] != backbone["id"]]
        moved.append(rm.part_name)
    return moved


def rides_on_a_primer_tail(resolved: list, cfg) -> list[str]:
    """Parts the planner will write into a primer tail, not amplify.

    P2A is 66 bp. Nobody runs a PCR whose primers are longer than its
    product, so the planner splits a short element between the two
    primers that flank it - which means there is no donor plasmid for
    it, and offering a list of plasmids to take it from is offering a
    choice that changes nothing.

    The same rule the planner uses: short, with something long enough
    to amplify on each side to carry the halves.
    """
    lengths = [len(rm.seq or "") for rm in resolved]
    names = [rm.part_name or rm.module.name or rm.feature_type
             for rm in resolved]
    floor = getattr(cfg, "min_pcr_length", 100)
    out = []
    for i, length in enumerate(lengths):
        if not length or length >= floor:
            continue
        before = any(x >= floor for x in lengths[:i])
        after = any(x >= floor for x in lengths[i + 1:])
        if before and after:
            out.append(names[i])
    return out


def _drop_folded_donors(candidates: dict, resolved: list, cfg) -> None:
    """Take the tail-borne parts out of the donor picker.

    They stay named - "where did P2A go?" is a fair question - but as a
    statement rather than as a list of plasmids to choose between.
    """
    folded = set(rides_on_a_primer_tail(resolved, cfg))
    if not folded:
        return
    kept = []
    for donor in candidates.get("donors") or []:
        covers = donor.get("covers") or []
        if covers and all(c in folded for c in covers):
            continue           # nothing to choose: it rides on the tails
        kept.append(donor)
    candidates["donors"] = kept
    candidates["n_donors"] = len([d for d in kept if d.get("options")])
    candidates["on_primer_tails"] = sorted(folded)


def _offer_to_order(candidates: dict, cfg) -> None:
    """Say what could be ordered for a part no lab plasmid carries.

    Only for the parts the scan has already given up on, and only when a
    token is configured, so an ordinary scan does nothing over the network
    and works offline exactly as before. What comes back is attached
    beside the plan and never inside it: a plasmid on Addgene is weeks and
    a purchase order from being a template, and letting it compete on cost
    with something already in a freezer would make every estimate a lie.
    """
    from .addgene import token as addgene_token

    absent = [d["covers"][0] for d in candidates.get("donors", [])
              if not d.get("options") and d.get("covers")]
    if not absent:
        return
    if not addgene_token(cfg):
        # Worth saying once that this could be answered, without turning
        # a missing optional token into a warning on every scan.
        candidates["to_order_hint"] = (
            "set GENESIS_ADDGENE_TOKEN to have Genesis look these up on "
            "Addgene; free token at https://developers.addgene.org/")
        return
    try:
        from .addgene import suggest_for_missing
        candidates["to_order"] = suggest_for_missing(absent, cfg=cfg)
    except Exception as exc:                      # never fail a scan on this
        candidates["to_order_hint"] = f"Addgene could not be searched: {exc}"


def scan_sources(conn: sqlite3.Connection, spec: ConstructSpec,
                 cfg: Config | None = None) -> dict:
    """Resolve the specification and say what it could be built from.

    Everything design() does before it commits to anything: what each
    module resolves to, which plasmid carries the longest leading run of
    the architecture, and where the rest would have to come from. Designing
    primers is the expensive half and depends on those choices, so this
    stops short of it - the parents are worth settling before a protocol
    is written around them.
    """
    cfg = cfg or Config()
    resolved = resolve(conn, spec)

    architecture = [{
        "type": r.feature_type, "part_id": r.part_id,
        "part_name": ("Cre-lox" if canonical_feature_type(r.feature_type)
                      == "Cre-lox" else r.part_name),
        "seq": r.seq,
    } for r in resolved
        if r.seq or canonical_feature_type(r.feature_type) == "Cre-lox"]

    wants_cassette = any(
        canonical_feature_type(m["type"] if isinstance(m, dict) else m.type)
        == "Cre-lox" for m in spec.modules)
    candidates = architecture_candidates(
        conn, architecture, family=spec.vector or "",
        prefer_cassette=wants_cassette)
    ranked = candidates.get("recipients", [])
    backbone = _load_backbone(
        conn, spec, cfg,
        preferred_id=ranked[0]["plasmid_id"] if ranked else None)

    _drop_folded_donors(candidates, resolved, cfg)
    _offer_to_order(candidates, cfg)

    return {
        "construct": {"name": spec.name},
        "backbone": ({"id": backbone["id"], "name": backbone["name"],
                      "length": backbone["seq_length"],
                      "owner": backbone["owner"]} if backbone else None),
        "modules": [{
            "type": r.feature_type, "requested": r.module.name,
            "part_id": r.part_id, "part_name": r.part_name,
            "length": r.length,
            "sources": [{"plasmid_id": x["id"], "name": x["name"],
                         "owner": x["owner"], "verified": x["verified"]}
                        for x in r.sources[:5]],
            "warnings": r.warnings,
        } for r in resolved],
        "candidates": candidates,
    }


#: The six bases Kozak's consensus puts before an initiating ATG. The
#: position that matters most is -3: a purine there, with a G at +4, is
#: most of the effect.
KOZAK = "GCCACC"


def _looks_like_kozak(upstream: str) -> bool:
    """Whether these six bases already initiate translation well.

    Judged on the -3 purine rather than the whole consensus, because that
    is the position with the measurable effect and most lab plasmids
    carry some variant of the rest.
    """
    return len(upstream) >= 3 and upstream[-3] in "AG"


def _extend_to_kozak(conn: sqlite3.Connection, resolved: list) -> list[str]:
    """Carry the donor's Kozak sequence into the initiating coding sequence.

    A coding sequence amplified from its ATG arrives with no Kozak, and
    the recipient's own one is cut away with the payload it belonged to,
    so the new construct puts a strong promoter straight onto a bare ATG.
    That still translates; it translates less well.

    The plasmid the gene is taken from usually has one sitting right in
    front of it - paav-ef1a-dio-asap7y has GCCACC - and six more bases on
    an amplicon that is already being made costs nothing at all. Where
    the donor has no Kozak this says so rather than inventing one, since
    adding bases no template carries is a different and more expensive
    kind of change.

    Only the module that *starts* translation gets one. A coding sequence
    after a 2A peptide is reached by the ribosome skipping, not by
    initiating again, so a Kozak there would be meaningless.
    """
    from .construct import TRANSLATED_MODULE_TYPES

    notes: list[str] = []
    inserted: list = []
    types = [rm.feature_type for rm in resolved]
    for i, rm in enumerate(resolved):
        if rm.feature_type not in TRANSLATED_MODULE_TYPES or not rm.seq:
            continue
        if rm.feature_type == "Kozak":
            continue
        prev = types[i - 1] if i else None
        if prev in TRANSLATED_MODULE_TYPES or prev == "Kozak":
            continue                      # reached by skipping, not initiating
        if not rm.seq[:3].upper() == "ATG":
            continue                      # nothing here initiates anyway
        src = rm.sources[0] if rm.sources else None
        if not src or src.get("start") is None or src.get("end") is None:
            continue
        row = conn.execute("SELECT seq, circular FROM plasmid WHERE id = ?",
                           (src["id"],)).fetchone()
        if row is None or not row["seq"]:
            continue
        host = row["seq"].upper()
        circular = bool(row["circular"])
        strand = src.get("strand", 1)
        if strand == 1:
            upstream = slice_region(host, src["start"] - len(KOZAK),
                                    src["start"], circular)
        else:
            upstream = revcomp(slice_region(host, src["end"],
                                            src["end"] + len(KOZAK), circular))
        upstream = upstream.upper()
        if len(upstream) < len(KOZAK):
            continue
        if not _looks_like_kozak(upstream):
            notes.append(
                f"{rm.part_name} is amplified from its ATG and "
                f"{src.get('name', 'the donor')} has no Kozak sequence in "
                f"front of it ({upstream}); the promoter will drive a bare "
                f"start codon, which works but initiates less efficiently - "
                f"add {KOZAK} to the 5' primer tail if that matters")
            continue
        # Its own element, not six bases glued onto the front of the gene.
        # Folding it in left the coding sequence no longer starting with
        # ATG, so the map reported a missing start codon for a construct
        # that had just been given a better one.
        from .construct import Module, ResolvedModule

        if strand == 1:
            k_src = {**src, "start": src["start"] - len(KOZAK),
                     "end": src["start"]}
        else:
            k_src = {**src, "start": src["end"],
                     "end": src["end"] + len(KOZAK)}
        inserted.append((i, ResolvedModule(
            Module(type="Kozak", name="Kozak"), None, "Kozak", upstream,
            "Kozak", [k_src], [])))
        notes.append(
            f"the amplicon was extended by {len(KOZAK)} bp upstream to take "
            f"the Kozak sequence ({upstream}) that already sits in front of "
            f"{rm.part_name} in {src.get('name', 'the donor')}")

    for offset, (at, module) in enumerate(inserted):
        resolved.insert(at + offset, module)
    return notes


def _extend_to_stop(conn: sqlite3.Connection, resolved: list) -> list[str]:
    """Carry the donor's own stop codon into a terminal coding sequence.

    Some genes are catalogued without their stop - every ASAP7y row in
    this inventory is - and when nothing translated follows, the ribosome
    runs off the end of the insert and into whatever the vector puts
    there. For pAAV-CAG-ASAP7y that was the WPRE, and the promoter drove
    a 660 residue protein instead of a 427 residue one.

    The plasmid the part is taken from usually has the stop sitting right
    after it. Three more bases on the amplicon is free, so the amplified
    region is extended to include it and the construct gets a sequence
    that terminates.
    """
    from .construct import TRANSLATED_MODULE_TYPES, STOP_CODONS

    notes: list[str] = []
    types = [rm.feature_type for rm in resolved]
    for i, rm in enumerate(resolved):
        if rm.feature_type not in TRANSLATED_MODULE_TYPES or not rm.seq:
            continue
        nxt = types[i + 1] if i + 1 < len(types) else None
        if nxt in TRANSLATED_MODULE_TYPES:
            continue                       # a fusion continues; no stop here
        if rm.seq[-3:].upper() in STOP_CODONS:
            continue                       # already terminates
        src = rm.sources[0] if rm.sources else None
        if not src or src.get("start") is None or src.get("end") is None:
            continue
        row = conn.execute("SELECT seq, circular FROM plasmid WHERE id = ?",
                           (src["id"],)).fetchone()
        if row is None or not row["seq"]:
            continue
        host = row["seq"].upper()
        circular = bool(row["circular"])
        strand = src.get("strand", 1)
        if strand == 1:
            nxt_codon = slice_region(host, src["end"], src["end"] + 3,
                                     circular)
        else:
            nxt_codon = revcomp(slice_region(host, src["start"] - 3,
                                             src["start"], circular))
        if nxt_codon.upper() not in STOP_CODONS:
            continue
        rm.seq = rm.seq + nxt_codon.upper()
        if strand == 1:
            src["end"] = src["end"] + 3
        else:
            src["start"] = src["start"] - 3
        notes.append(
            f"{rm.part_name} is catalogued without a stop codon, and nothing "
            f"translated follows it; the amplicon was extended by 3 bp to "
            f"take the {nxt_codon.upper()} that follows it in "
            f"{src.get('name', 'the donor')}, so the protein terminates")
    return notes


def design(conn: sqlite3.Connection, spec: ConstructSpec,
           cfg: Config | None = None) -> dict:
    cfg = cfg or Config()

    # Stage 4: resolve modules against the inventory
    resolved = resolve(conn, spec)

    # The architecture in the order the user wrote it, kept before the
    # cassette expansion inverts the payload.
    logical_modules = [{
        "type": r.feature_type, "part_id": r.part_id,
        "part_name": r.part_name, "seq": r.seq,
    } for r in resolved if r.seq]

    # The architecture the candidate search reasons about keeps the
    # recombinase cassette in place: "hSyn - DIO - ASAP7y" is a different
    # request from "hSyn - ASAP7y", and a recipient that already has the
    # cassette is the one to start from.
    architecture = [{
        "type": r.feature_type, "part_id": r.part_id,
        "part_name": ("Cre-lox" if canonical_feature_type(r.feature_type)
                      == "Cre-lox" else r.part_name),
        "seq": r.seq,
    } for r in resolved
        if r.seq or canonical_feature_type(r.feature_type) == "Cre-lox"]

    # Read the request as an ordered architecture first, because what to
    # start from follows from it.
    wants_cassette = any(
        canonical_feature_type(m["type"] if isinstance(m, dict) else m.type)
        == "Cre-lox" for m in spec.modules)
    candidates = architecture_candidates(
        conn, architecture, family=spec.vector or "",
        prefer_cassette=wants_cassette)
    _drop_folded_donors(candidates, resolved, cfg)
    ranked = candidates.get("recipients", [])
    backbone = _load_backbone(
        conn, spec, cfg,
        preferred_id=ranked[0]["plasmid_id"] if ranked else None)

    # An explicit choice of where to take a part from outranks the search
    # that ranked the sources, and outranks the recipient too: somebody
    # saying "that tube" has a reason the inventory cannot see.
    for rm in resolved:
        wanted = spec.source_plasmid_ids.get(rm.part_id or -1)
        if not wanted or not rm.seq:
            continue
        chosen = [x for x in rm.sources if x["id"] == wanted]
        if chosen:
            rm.sources = chosen + [x for x in rm.sources if x["id"] != wanted]
            rm.pinned_source = True
            continue
        # The catalogue does not link every plasmid to every part it
        # contains - a sequence can be in there without an occurrence row -
        # and a choice that is silently ignored is worse than one that is
        # refused. Ask the DNA before giving up.
        row = conn.execute(
            "SELECT id, name, owner, verified, available_template, circular, "
            "seq_length, seq FROM plasmid WHERE id = ?", (wanted,)).fetchone()
        if row is None or not row["seq"]:
            rm.warnings.append(
                f"the plasmid chosen to take {rm.part_name} from is not in "
                f"the inventory; the usual source was used instead")
            continue
        host = row["seq"].upper()
        want = rm.seq.upper()
        at, strand = _locate_in(host, want, bool(row["circular"]))
        if at < 0:
            rm.warnings.append(
                f"{row['name']} does not contain {rm.part_name}, so it was "
                f"taken from its usual source instead")
            continue
        rm.sources = [{
            "id": row["id"], "name": row["name"], "owner": row["owner"] or "",
            "verified": row["verified"],
            "available_template": row["available_template"],
            "circular": row["circular"], "seq_length": row["seq_length"],
            "start": at, "end": at + len(want), "strand": strand,
        }] + [x for x in rm.sources if x["id"] != wanted]
        rm.pinned_source = True

    stop_notes = _extend_to_stop(conn, resolved)
    stop_notes += _extend_to_kozak(conn, resolved)
    fused_notes = _prefer_fused_source(conn, resolved)

    from_recipient = (_prefer_backbone_as_source(resolved, backbone)
                      if backbone else [])

    # Expand any DIO module into its cassette before assembling, so the
    # molecule that gets validated is the one that would actually be built.

    # When the recipient already carries lox sites, the cassette is the
    # vector's and the payload goes inside it.
    from .construct import backbone_cassette
    backbone_rows = []
    if backbone:
        backbone_rows = [
            {"feature_type": r["feature_type"], "name": r["name"],
             "start": r["start"], "end": r["end"], "strand": r["strand"]}
            for r in conn.execute(
                """SELECT feature_type, name, start, end, strand
                   FROM annotation WHERE plasmid_id = ? ORDER BY start""",
                (backbone["id"],))]
    reuse_sites = backbone_cassette(backbone_rows) is not None
    resolved, dio_notes = expand_dio(resolved, sites_from_backbone=reuse_sites)

    # Stage 4: build the intended final molecule. With a backbone the
    # payload replaces that plasmid's existing cassette, so the result is a
    # real plasmid carrying its ITRs, WPRE, polyA, origin and resistance -
    # not the payload bolted onto the end of one.
    backbone_annotations = backbone_rows
    if backbone:
        spec_types = {canonical_feature_type(
            m["type"] if isinstance(m, dict) else m.type) for m in spec.modules}
        construct = assemble_into_backbone(
            resolved, backbone["seq"], backbone_annotations, spec_types,
            circular=spec.circular)
    else:
        construct = assemble(resolved, "", circular=spec.circular)
    construct.warnings.extend(dio_notes)
    # These were being computed and dropped, so a plan could silently take
    # a fusion whole, move a source onto the recipient or extend an
    # amplicon to a stop codon without ever saying it had.
    construct.warnings.extend(stop_notes)
    construct.warnings.extend(fused_notes)
    if from_recipient:
        construct.warnings.append(
            "taken from the recipient itself rather than another plasmid: "
            + ", ".join(from_recipient))

    # A specification that never mentioned a vector takes the configured
    # default; one that explicitly asked for auto keeps inferring.
    if spec.vector is None:
        spec.vector = cfg.default_vector

    # --- vector family ---------------------------------------------------
    # Either the family was asked for, or it is inferred. Either way the
    # finished molecule is checked against what that family actually needs,
    # because a payload on its own is not a vector.
    # The layout now holds everything in the molecule, backbone included,
    # so it is the whole evidence for what family this construct belongs to.
    construct_annotations = [{"feature_type": i["type"], "name": i["name"]}
                             for i in construct.layout]

    # Only a real plasmid's name is evidence. The name someone typed for the
    # construct they intend to build is a wish, not an observation, and
    # letting it vote makes "LENTI-TEST" detect as lentiviral on its own.
    inferred = classify_vector(construct_annotations,
                               backbone["name"] if backbone else "")
    family = spec.vector or inferred["family"]
    completeness = missing_elements(
        construct_annotations, family,
        backbone["name"] if backbone else "")

    vector_info = {
        "requested": spec.vector or None,
        "detected": inferred["family"],
        "detected_evidence": inferred["evidence"],
        "family": family,
        "spec": describe(family),
        "completeness": completeness,
    }
    if inferred.get("warning"):
        vector_info["warning"] = inferred["warning"]

    # Stage 5: validate before planning anything
    # Compare against the expanded architecture: after a DIO expansion the
    # assembled order legitimately differs from the order that was typed.
    # Only the modules this design places are compared against the
    # specification; the backbone's own features are not part of the request.
    expected = [rm.feature_type for rm in resolved if rm.seq] or None
    placed_types = [i["type"] for i in construct.layout
                    if not i.get("from_backbone")]
    issues = validate_construct(construct, wants_cassette=wants_cassette,
                                expected_order=expected,
                                stocked_enzymes=list(cfg.stocked_enzymes))

    # Every architectural element that was asked for has to be in the
    # molecule. A Cre-lox module whose payload scan ended early was dropped
    # with a note nobody surfaced, and the result was a pAAV-CAG-DIO-ASAP7y
    # byte-identical to pAAV-CAG-ASAP7y: constitutive, no lox sites, no
    # error, and named for a cassette it did not contain.
    placed_set = set(placed_types)
    for m in spec.modules:
        want = canonical_feature_type(m["type"] if isinstance(m, dict)
                                      else m.type)
        if want == "Cre-lox":
            if not any(t == "lox_site" for t in placed_set):
                issues.append(Issue(
                    "error", "architecture",
                    "you asked for a Cre-lox cassette but the finished "
                    "construct has no lox sites, so it would express "
                    "constitutively rather than only in Cre-positive cells; "
                    "this is not the construct the name describes"))
        elif want and want not in placed_set and want not in VECTOR_TYPES:
            issues.append(Issue(
                "error", "architecture",
                f"you asked for a {want} but the finished construct does not "
                f"contain one, so the molecule does not match the "
                f"architecture it is named for"))

    if spec.vector and backbone is None:
        issues.append(Issue(
            "error", "packaging",
            f"you asked for {describe(spec.vector)['label']} but no backbone "
            f"was found, so the construct is only the payload; name a "
            f"backbone, or pin one with backbone_plasmid_id"))

    if spec.vector and inferred["family"] not in (spec.vector, "plasmid"):
        issues.append(Issue(
            "error", "packaging",
            f"you asked for {describe(spec.vector)['label']} but the backbone "
            f"and payload look like {inferred['label']}"
            + (f" ({'; '.join(inferred['evidence'])})"
               if inferred["evidence"] else "")))

    if family in ("AAV", "lentiviral"):
        backbone_types = {a["feature_type"] for a in backbone_annotations}
        for element in completeness["missing_required"]:
            from_backbone = backbone is not None and element not in backbone_types
            where = (f"; the backbone {backbone['name']} has none annotated "
                     f"either, so check whether it is really absent or just "
                     f"unlabelled" if from_backbone else "")
            issues.append(Issue(
                "error", "packaging",
                f"{_article(completeness['label'])} needs {element}, which "
                f"this construct does not have{where}"))
        for element in completeness["missing_named"]:
            issues.append(Issue(
                "error", "packaging",
                f"{_article(completeness['label'])} needs its {element}, "
                f"which is not annotated anywhere in this construct"))
        for element in completeness["missing_expected"]:
            issues.append(Issue(
                "warning", "packaging",
                f"no {element} found; {_article(completeness['label'])} "
                f"normally carries one and expression is usually lower "
                f"without it"))
        if completeness["complete"]:
            issues.append(Issue(
                "info", "packaging",
                f"carries the full {completeness['label']} element set"))

    # Read the finished molecule back independently of the layout that built
    # it, and check it is the construct that was asked for.
    verification = verify_plan(construct, resolved, spec.circular)
    for problem in verification["problems"]:
        issues.append(Issue("error", "verification", problem))
    if verification["passed"]:
        issues.append(Issue("info", "verification", verification["summary"]))

    validation = summarize(issues)

    # Build the module list first; the swap planner works from it.
    module_rows = [{
        "type": r.feature_type, "requested": r.module.name,
        "part_id": r.part_id, "part_name": r.part_name, "length": r.length,
        "seq": r.seq,
    } for r in resolved]

    # The usual real-world route: start from the closest existing plasmid and
    # swap one region, rather than rebuilding the construct from parts.
    analysis = analyse_swap(conn, module_rows, cfg,
                            require_types=required_architecture(spec),
                            require_family=spec.vector,
                            prefer_plasmid_id=backbone["id"] if backbone else None,
                            pin_recipient_id=spec.recipient_plasmid_id,
                            pin_donor_id=spec.donor_plasmid_id,
                            tier_modules=logical_modules)
    # Only draw a scheme for a swap that is actually on offer; a rejected
    # analysis still carries a recipient, and drawing it would illustrate
    # a plan nobody should follow.
    scheme = (build_scheme(conn, analysis)
              if analysis and analysis.feasible else {})

    # Stages 7-8: generate routes
    fragments = group_contiguous(resolved, construct, cfg)
    routes = [
        plan_hifi_route(conn, construct, fragments, backbone, cfg),
        plan_gblock_route(conn, construct, backbone, cfg),
        plan_restriction_route(conn, construct, fragments, backbone, cfg),
    ]
    if analysis is not None:
        routes.insert(0, plan_swap_hybrid_route(conn, analysis, cfg, scheme))
        routes.insert(1, plan_swap_hifi_route(conn, analysis, cfg, scheme))
        routes.insert(2, plan_swap_restriction_route(conn, analysis, cfg, scheme))

    # An unresolved required module invalidates every route: the plan would
    # build a different construct from the one that was asked for.
    if construct.missing_required:
        missing = ", ".join(
            f"{m['type']}" + (f"/{m['name']}" if m["name"] else "")
            for m in construct.missing_required)
        for route in routes:
            route.feasible = False
            route.warnings.insert(
                0, f"no route is offered because the required module(s) "
                   f"{missing} could not be resolved; the construct would be "
                   f"built without them")

    # Stage 10: rank
    routes = rank_routes(routes, cfg, failure_history(conn))

    # The order the user specified, before DIO expansion inverts the payload.
    # A construct name should describe this, not the assembled DNA order.
    spec_order = [{"type": (m["type"] if isinstance(m, dict) else m.type),
                   "name": (m["name"] if isinstance(m, dict) else m.name),
                   "mode": (m.get("mode") if isinstance(m, dict) else m.mode)}
                  for m in spec.modules]

    plan = {
        "construct": {
            "name": spec.name,
            "spec_order": spec_order,
            "length": construct.length,
            "circular": construct.circular,
            "layout": construct.layout,
            "codon_markers": codon_markers(construct),
            "sequence": construct.seq,
            "warnings": list(construct.warnings),
        },
        "backbone": {
            "id": backbone["id"], "name": backbone["name"],
            "length": backbone["seq_length"], "owner": backbone["owner"],
        } if backbone else None,
        "modules": [{
            "type": r.feature_type, "requested": r.module.name,
            "part_id": r.part_id, "part_name": r.part_name,
            "length": r.length,
            # What this module actually resolved to, which is not always the
            # catalogued part: one the recipient already carries is rebound
            # to the recipient's own copy, and anything checking the finished
            # molecule has to look for those bases rather than the part's.
            "seq": r.seq,
            "sources": [{"plasmid_id": s["id"], "name": s["name"],
                         "owner": s["owner"], "verified": s["verified"],
                         "start": s["start"], "end": s["end"],
                         "strand": s["strand"]} for s in r.sources[:5]],
            "warnings": r.warnings,
        } for r in resolved],
        "fragments": [{
            "name": f.name, "length": len(f.seq), "source": f.source_name,
            "source_plasmid_id": f.source_plasmid_id,
            "needs_synthesis": f.needs_synthesis,
            # A fragment nobody amplifies still goes in the tube, so the
            # simulation needs its bases: without them it was handed an
            # empty product and the assembly could not close.
            "seq": f.placed_seq if f.needs_synthesis else "",
            "modules": [m["name"] for m in f.modules],
        } for f in fragments],
        "validation": validation,
        "vector": vector_info,
        "verification": verification,
        "swap": analysis.to_dict() if analysis else None,
        "scheme": scheme,
        # Before anything else: is this already on the shelf?
        "already_made": already_made(conn, construct, resolved),
        "closest_parents": [p.to_dict() for p in rank_parents(
            conn, [m["part_id"] for m in module_rows if m["part_id"]], cfg,
            limit=6, require_types=required_architecture(spec),
            require_family=spec.vector)],
        "routes": [r.to_dict(cfg) for r in routes],
        # Read as an ordered architecture: which plasmids already carry the
        # longest leading run of it, and what the rest has to come from.
        "candidates": candidates,
    }

    # A name proposed from the resolved architecture, so the user does not
    # have to invent one that matches how these plasmids are already named.
    from .protocol import suggest_name
    plan["construct"]["suggested_name"] = suggest_name(plan)

    # Stage 11b: run the protocol on paper.  Verification reads the designed
    # molecule back, which cannot catch a design that verifies perfectly
    # while the steps underneath it build something else; this amplifies the
    # primers as ordered off the templates as named and sees what comes out.
    from .simulate import simulate_route, template_atlas
    # What each named template is, for every route - a route the planner
    # rejected still gets drawn, and a drawing needs the plasmid's size.
    try:
        plan["templates"] = template_atlas(conn, plan)
    except Exception:                                 # drawing data only
        plan["templates"] = {}
    for route in plan["routes"]:
        if not route.get("feasible"):
            continue
        try:
            sim = simulate_route(conn, plan, route).to_dict()
        except Exception as exc:                      # never fail a plan on this
            route["simulation"] = {"simulated": False,
                                   "summary": f"could not be simulated: {exc}"}
            continue
        route["simulation"] = sim
        if sim.get("simulated") and not sim.get("passed"):
            route["feasible"] = False
            for problem in sim.get("problems", []):
                route["warnings"].insert(0, f"simulation: {problem}")
        else:
            for w in sim.get("warnings", []):
                route["warnings"].append(f"simulation: {w}")

    # Stage 11c: a look before the next step hides the answer.
    # Spliced in here rather than in each route planner: every route
    # that makes a molecule wants the same checkpoints, and six copies
    # of that decision would drift apart.
    _add_checkpoints(conn, plan)

    # Past a few reactions, ordering the insert stops being the expensive
    # option and starts being the quick one. Cost says otherwise - a gBlock
    # is real money and a week of waiting - so this is said rather than
    # ranked, and the choice left where it belongs.
    gblock = next((r for r in plan["routes"]
                   if r["id"] == "gblock" and r.get("feasible")), None)
    if gblock:
        ordered = next((st for st in gblock["steps"]
                        if st["operation"] == "SYNTHESIS"), None)
        size = (ordered or {}).get("details", {}).get("length", 0)
        for route in plan["routes"]:
            if (route["id"] != "gblock" and route.get("feasible")
                    and route["n_pcr"] >= cfg.synthesis_beats_pcr_at):
                route["warnings"].append(
                    f"this takes {route['n_pcr']} PCRs and "
                    f"{len(route['junctions'])} junctions to get right; the "
                    f"same construct is one {size:,} bp fragment to order "
                    f"plus one vector PCR - see \"Order the insert, amplify "
                    f"the vector\" for the sequence to send")

    # Ranking happened before any of this, so a route the simulation has
    # just ruled out would otherwise keep the position it was ranked into -
    # and lead the list while being unusable. Feasible ones first; the
    # sort is stable, so the ranking survives within each group.
    def _rank(r):
        """Feasible first, then proven, then cheap.

        A route that could not be simulated is not evidence of anything,
        and it used to sort level with one whose product was checked
        against the design. That let an unverified route take the top of
        the page - and the top of the page is what gets built.
        """
        sim = r.get("simulation") or {}
        simulated = bool(sim.get("simulated"))
        return (
            not r.get("feasible"),                       # usable routes first
            not (simulated and sim.get("identical")),    # proven identical
            not simulated,                               # at least checked
            r.get("cost", 1e9),
        )

    plan["routes"].sort(key=_rank)
    return plan


def save_design(conn: sqlite3.Connection, spec: ConstructSpec, plan: dict,
                owner: str = "", mode: str = "version",
                design_id: int | None = None) -> dict:
    """Store a design (Stage 11).

    Three ways to save, because they mean different things to the person
    doing it:

    * ``version``   keep the previous one and add a new version. The default,
      and the only one that never loses work.
    * ``overwrite`` replace an existing version in place. Used when fixing a
      design that was just saved, where a new version would be noise.
    * ``save_as``   store under a new name, leaving the original untouched.

    ``save_as`` requires that ``spec.name`` is free, so it cannot silently
    merge into an existing construct's history.
    """
    if mode not in ("version", "overwrite", "save_as"):
        raise ValueError(f"unknown save mode {mode!r}")

    existing = conn.execute("SELECT id FROM construct WHERE name = ?",
                            (spec.name,)).fetchone()

    # Keep the stored plan's own label in step with the construct it is
    # saved under, or a renamed design shows its old name in every export.
    if plan.get("construct", {}).get("name") != spec.name:
        plan = {**plan, "construct": {**plan["construct"], "name": spec.name}}

    if mode == "save_as" and existing:
        raise ValueError(
            f"a construct called {spec.name!r} already exists; choose another "
            f"name, or save a new version of that one instead")

    if mode == "overwrite":
        if design_id is None:
            raise ValueError("overwrite needs the design to replace")
        row = conn.execute(
            "SELECT construct_id, version FROM design WHERE id = ?",
            (design_id,)).fetchone()
        if row is None:
            raise ValueError(f"no design with id {design_id}")
        construct_id, version = row["construct_id"], row["version"]
        conn.execute("UPDATE construct SET name = ?, spec = ? WHERE id = ?",
                     (spec.name, json.dumps(spec.to_dict()), construct_id))
        conn.execute(
            """UPDATE design SET plan = ?, final_seq = ?, validation = ?,
                      owner = ?, created_at = CURRENT_TIMESTAMP
               WHERE id = ?""",
            (json.dumps(plan), plan["construct"]["sequence"],
             json.dumps(plan["validation"]), owner or spec.owner, design_id))
        conn.commit()
        return {"design_id": design_id, "construct_id": construct_id,
                "version": version, "mode": "overwrite"}

    if existing:
        construct_id = existing["id"]
        conn.execute("UPDATE construct SET spec = ? WHERE id = ?",
                     (json.dumps(spec.to_dict()), construct_id))
    else:
        cur = conn.execute(
            "INSERT INTO construct(name, owner, spec) VALUES (?,?,?)",
            (spec.name, owner or spec.owner, json.dumps(spec.to_dict())))
        construct_id = cur.lastrowid

    version = conn.execute(
        "SELECT COALESCE(MAX(version), 0) + 1 AS v FROM design WHERE construct_id = ?",
        (construct_id,)).fetchone()["v"]

    cur = conn.execute(
        """INSERT INTO design(construct_id, version, owner, plan, final_seq,
                              validation)
           VALUES (?,?,?,?,?,?)""",
        (construct_id, version, owner or spec.owner, json.dumps(plan),
         plan["construct"]["sequence"], json.dumps(plan["validation"])))
    conn.commit()
    return {"design_id": cur.lastrowid, "construct_id": construct_id,
            "version": version, "mode": mode}


def delete_design(conn: sqlite3.Connection, design_id: int) -> dict:
    """Delete one version. The construct goes too once nothing is left."""
    row = conn.execute(
        """SELECT d.construct_id, d.version, c.name
           FROM design d JOIN construct c ON c.id = d.construct_id
           WHERE d.id = ?""", (design_id,)).fetchone()
    if row is None:
        raise ValueError(f"no design with id {design_id}")

    conn.execute("DELETE FROM outcome WHERE design_id = ?", (design_id,))
    conn.execute("DELETE FROM design WHERE id = ?", (design_id,))
    remaining = conn.execute(
        "SELECT COUNT(*) c FROM design WHERE construct_id = ?",
        (row["construct_id"],)).fetchone()["c"]
    if remaining == 0:
        conn.execute("DELETE FROM construct WHERE id = ?", (row["construct_id"],))
    conn.commit()
    return {"deleted": design_id, "name": row["name"],
            "version": row["version"], "versions_left": remaining}


def delete_construct(conn: sqlite3.Connection, construct_id: int) -> dict:
    """Delete a construct and every version of it."""
    row = conn.execute("SELECT name FROM construct WHERE id = ?",
                       (construct_id,)).fetchone()
    if row is None:
        raise ValueError(f"no construct with id {construct_id}")
    ids = [r["id"] for r in conn.execute(
        "SELECT id FROM design WHERE construct_id = ?", (construct_id,))]
    for did in ids:
        conn.execute("DELETE FROM outcome WHERE design_id = ?", (did,))
    conn.execute("DELETE FROM design WHERE construct_id = ?", (construct_id,))
    conn.execute("DELETE FROM construct WHERE id = ?", (construct_id,))
    conn.commit()
    return {"deleted_construct": construct_id, "name": row["name"],
            "versions_deleted": len(ids)}
