"""Orchestration: specification in, ranked laboratory plan out.

This is the single entry point that walks Stages 4 -> 10 in order and
returns one JSON-serialisable plan.
"""

from __future__ import annotations

import json
import sqlite3

from .config import Config
from .seqtools import revcomp
from .construct import (ConstructSpec, assemble,
                        assemble_into_backbone, codon_markers,
                        expand_dio, resolve)
from .planner import (group_contiguous, plan_hifi_route, plan_restriction_route,
                      plan_swap_hifi_route, plan_swap_hybrid_route,
                      plan_swap_restriction_route,
                      rank_routes, failure_history)
from .swap import analyse_swap, architecture_candidates, build_scheme, rank_parents
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

    seen: set = set()
    for src in piece[0].sources[:40]:
        if src["id"] in seen:
            continue
        seen.add(src["id"])
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
        at, strand = host.find(wanted), 1
        if at < 0:
            at, strand = host.find(revcomp(wanted)), -1
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
    ranked = candidates.get("recipients", [])
    backbone = _load_backbone(
        conn, spec, cfg,
        preferred_id=ranked[0]["plasmid_id"] if ranked else None)

    # An explicit choice of where to take a part from outranks the search
    # that ranked the sources, and outranks the recipient too: somebody
    # saying "that tube" has a reason the inventory cannot see.
    for rm in resolved:
        wanted = spec.source_plasmid_ids.get(rm.part_id or -1)
        if wanted and rm.sources:
            chosen = [x for x in rm.sources if x["id"] == wanted]
            if chosen:
                rm.sources = chosen + [x for x in rm.sources
                                       if x["id"] != wanted]
                rm.pinned_source = True

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
    issues = validate_construct(construct, expected_order=expected,
                                stocked_enzymes=list(cfg.stocked_enzymes))

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
    from .simulate import simulate_route
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

    # Ranking happened before any of this, so a route the simulation has
    # just ruled out would otherwise keep the position it was ranked into -
    # and lead the list while being unusable. Feasible ones first; the
    # sort is stable, so the ranking survives within each group.
    plan["routes"].sort(key=lambda r: not r.get("feasible"))
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
