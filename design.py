"""Orchestration: specification in, ranked laboratory plan out.

This is the single entry point that walks Stages 4 -> 10 in order and
returns one JSON-serialisable plan.
"""

from __future__ import annotations

import json
import sqlite3

from .config import Config
from .construct import (ConstructSpec, assemble, codon_markers,
                        expand_dio, resolve)
from .planner import (group_contiguous, plan_hifi_route, plan_restriction_route,
                      plan_swap_hifi_route, plan_swap_restriction_route,
                      rank_routes, failure_history)
from .swap import analyse_swap, build_scheme, rank_parents
from .validate import summarize, validate_construct


#: Shortest sequence that could plausibly act as a cloning vector. The
#: inventory contains many short fragments whose names also match a backbone
#: search term, and silently picking a 23 bp fragment as the vector produces
#: a nonsense plan.
MIN_BACKBONE_LENGTH = 1500


def _load_backbone(conn: sqlite3.Connection, spec: ConstructSpec) -> dict | None:
    row = None
    if spec.backbone_plasmid_id:
        # An explicit id is the user's choice; honour it as given.
        row = conn.execute("SELECT * FROM plasmid WHERE id = ?",
                           (spec.backbone_plasmid_id,)).fetchone()
    elif spec.backbone:
        # Name search: a vector must be big enough to be one, and a circular,
        # verified plasmid carrying an origin beats a bare fragment.
        row = conn.execute(
            """SELECT * FROM plasmid
               WHERE name LIKE ? COLLATE NOCASE AND seq_length >= ?
               ORDER BY
                 circular DESC,
                 (verified != 'unverified') DESC,
                 (SELECT COUNT(*) FROM annotation a
                  WHERE a.plasmid_id = plasmid.id
                    AND a.feature_type IN ('origin', 'resistance')) DESC,
                 seq_length DESC
               LIMIT 1""",
            (f"%{spec.backbone}%", MIN_BACKBONE_LENGTH)).fetchone()
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

    # Stage 4: build the intended final molecule
    construct = assemble(resolved, backbone["seq"] if backbone else "",
                         circular=spec.circular)
    construct.warnings.extend(dio_notes)

    # Stage 5: validate before planning anything
    # Compare against the expanded architecture: after a DIO expansion the
    # assembled order legitimately differs from the order that was typed.
    expected = [rm.feature_type for rm in resolved if rm.seq] or None
    issues = validate_construct(construct, expected_order=expected,
                                stocked_enzymes=list(cfg.stocked_enzymes))
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
                   "name": (m["name"] if isinstance(m, dict) else m.name)}
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
                owner: str = "") -> int:
    """Store a versioned design (Stage 11)."""
    row = conn.execute("SELECT id FROM construct WHERE name = ?",
                       (spec.name,)).fetchone()
    if row:
        construct_id = row["id"]
        conn.execute("UPDATE construct SET spec = ? WHERE id = ?",
                     (json.dumps(spec.to_dict()), construct_id))
    else:
        cur = conn.execute(
            "INSERT INTO construct(name, owner, spec) VALUES (?,?,?)",
            (spec.name, owner or spec.owner, json.dumps(spec.to_dict())))
        construct_id = cur.lastrowid

    ver = conn.execute(
        "SELECT COALESCE(MAX(version), 0) + 1 AS v FROM design WHERE construct_id = ?",
        (construct_id,)).fetchone()["v"]

    cur = conn.execute(
        """INSERT INTO design(construct_id, version, owner, plan, final_seq,
                              validation)
           VALUES (?,?,?,?,?,?)""",
        (construct_id, ver, owner or spec.owner, json.dumps(plan),
         plan["construct"]["sequence"], json.dumps(plan["validation"])))
    conn.commit()
    return cur.lastrowid
