"""Orchestration: specification in, ranked laboratory plan out.

This is the single entry point that walks Stages 4 -> 10 in order and
returns one JSON-serialisable plan.
"""

from __future__ import annotations

import json
import sqlite3

from .config import Config
from .construct import (ConstructSpec, assemble,
                        assemble_into_backbone, codon_markers,
                        expand_dio, resolve)
from .planner import (group_contiguous, plan_hifi_route, plan_restriction_route,
                      plan_swap_hifi_route, plan_swap_restriction_route,
                      rank_routes, failure_history)
from .swap import analyse_swap, build_scheme, rank_parents
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


def _load_backbone(conn: sqlite3.Connection, spec: ConstructSpec) -> dict | None:
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
    return dict(row) if row else None


def design(conn: sqlite3.Connection, spec: ConstructSpec,
           cfg: Config | None = None) -> dict:
    cfg = cfg or Config()

    # Stage 4: resolve modules against the inventory
    resolved = resolve(conn, spec)
    backbone = _load_backbone(conn, spec)

    # Expand any DIO module into its cassette before assembling, so the
    # molecule that gets validated is the one that would actually be built.
    resolved, dio_notes = expand_dio(resolved)

    # Stage 4: build the intended final molecule. With a backbone the
    # payload replaces that plasmid's existing cassette, so the result is a
    # real plasmid carrying its ITRs, WPRE, polyA, origin and resistance -
    # not the payload bolted onto the end of one.
    backbone_annotations = []
    if backbone:
        backbone_annotations = [
            {"feature_type": r["feature_type"], "name": r["name"],
             "start": r["start"], "end": r["end"], "strand": r["strand"]}
            for r in conn.execute(
                """SELECT feature_type, name, start, end, strand
                   FROM annotation WHERE plasmid_id = ? ORDER BY start""",
                (backbone["id"],))]
        spec_types = {canonical_feature_type(
            m["type"] if isinstance(m, dict) else m.type) for m in spec.modules}
        construct = assemble_into_backbone(
            resolved, backbone["seq"], backbone_annotations, spec_types,
            circular=spec.circular)
    else:
        construct = assemble(resolved, "", circular=spec.circular)
    construct.warnings.extend(dio_notes)

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
    expected = [rm.feature_type for rm in resolved if rm.seq] or None
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
        for element in completeness["missing_required"]:
            issues.append(Issue(
                "error", "packaging",
                f"{_article(completeness['label'])} needs {element}, which "
                f"this construct does not have"))
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
    } for r in resolved]

    # The usual real-world route: start from the closest existing plasmid and
    # swap one region, rather than rebuilding the construct from parts.
    analysis = analyse_swap(conn, module_rows, cfg)
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
        routes.insert(0, plan_swap_hifi_route(conn, analysis, cfg, scheme))
        routes.insert(1, plan_swap_restriction_route(conn, analysis, cfg, scheme))

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
            "modules": [m["name"] for m in f.modules],
        } for f in fragments],
        "validation": validation,
        "vector": vector_info,
        "verification": verification,
        "swap": analysis.to_dict() if analysis else None,
        "scheme": scheme,
        "closest_parents": [p.to_dict() for p in rank_parents(
            conn, [m["part_id"] for m in module_rows if m["part_id"]], cfg,
            limit=6)],
        "routes": [r.to_dict(cfg) for r in routes],
    }

    # A name proposed from the resolved architecture, so the user does not
    # have to invent one that matches how these plasmids are already named.
    from .protocol import suggest_name
    plan["construct"]["suggested_name"] = suggest_name(plan)
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
