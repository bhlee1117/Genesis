"""Autosaved working plans, and their expiry.

A design costs a minute of computation and an afternoon of thought, and
until now it lived only in the browser tab it was computed in. Closing
the tab lost it; so did the server, which is a foreground process on one
Mac and stopped when that Mac rebooted, taking a week of unsaved work
with it.

So every plan the designer computes is written down as it is produced.
These are drafts, not designs: `save_design` is the deliberate act, with
a name, a version and a history, and nothing here is allowed anywhere
near it. A draft is the opposite - nobody asked for it, so nobody will
miss it, which is the only reason it is safe to delete one on a timer.

One spec keeps one draft. Re-running the same specification replaces
its plan and pushes the expiry out, so a construct somebody comes back
to on Thursday is still there, and the week is measured from the last
time it was touched rather than the first.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3

#: How long an untouched draft survives.
TTL_DAYS = 7

#: Most drafts one person accumulates before the oldest start dropping.
#: Expiry alone is enough for a person; this is the guard against a
#: script calling /api/design in a loop.
MAX_PER_OWNER = 50


def spec_key(spec: dict) -> str:
    """A stable fingerprint of a specification.

    Keyed on the specification rather than the construct name, because
    the name is the last thing a designer settles: the same plan gets
    recomputed a dozen times while it is still called "untitled".
    """
    return hashlib.sha256(
        json.dumps(spec or {}, sort_keys=True, default=str).encode()
    ).hexdigest()


def _chosen_route(plan: dict) -> str:
    """The route a draft would reopen on: the first one that works."""
    for r in (plan or {}).get("routes") or []:
        if r.get("feasible"):
            return r.get("id") or ""
    return ""


def autosave(conn: sqlite3.Connection, spec: dict, plan: dict,
             owner: str = "", route_id: str = "") -> dict:
    """Write down a plan nobody asked to keep, and expire the stale ones.

    Returns the draft's id and when it will go, so the interface can say
    so rather than deleting the designer's work silently a week later.
    """
    purged = purge(conn)
    key = spec_key(spec)
    name = (plan.get("construct") or {}).get("name") or spec.get("name") or ""

    conn.execute(
        """INSERT INTO plan_draft(owner, name, spec_key, spec, plan, route_id,
                                  expires_at)
           VALUES (?,?,?,?,?,?, datetime('now', ?))
           ON CONFLICT(owner, spec_key) DO UPDATE SET
               name       = excluded.name,
               spec       = excluded.spec,
               plan       = excluded.plan,
               route_id   = excluded.route_id,
               saved_at   = CURRENT_TIMESTAMP,
               expires_at = excluded.expires_at""",
        (owner, name, key, json.dumps(spec, default=str),
         json.dumps(plan, default=str), route_id or _chosen_route(plan),
         f"+{TTL_DAYS} days"))

    # Oldest-first, so the one being worked on now is never the casualty.
    # `saved_at` is only second-resolution, and the case this cap exists
    # for - a caller in a loop - writes every row inside one second, so
    # ordering on it alone leaves the tie to the engine, which breaks it
    # by rowid and keeps the fifty *oldest*. The id settles it.
    conn.execute(
        """DELETE FROM plan_draft
            WHERE owner = ? AND id NOT IN (
                SELECT id FROM plan_draft WHERE owner = ?
                 ORDER BY saved_at DESC, id DESC LIMIT ?)""",
        (owner, owner, MAX_PER_OWNER))
    conn.commit()

    row = conn.execute(
        "SELECT id, expires_at FROM plan_draft WHERE owner = ? AND spec_key = ?",
        (owner, key)).fetchone()
    return {"draft_id": row["id"] if row else None,
            "expires_at": row["expires_at"] if row else None,
            "expired_removed": purged}


def purge(conn: sqlite3.Connection) -> int:
    """Delete the drafts whose week is up.

    Scoped to `plan_draft` and to nothing else, on purpose: this is the
    one place in Genesis that deletes a person's work without being
    asked, and the only thing that makes it acceptable is that it cannot
    reach anything that was deliberately saved.
    """
    cur = conn.execute(
        "DELETE FROM plan_draft WHERE expires_at <= datetime('now')")
    conn.commit()
    return cur.rowcount


def recent(conn: sqlite3.Connection, owner: str = "", limit: int = 25) -> list:
    """Drafts still in date, newest first, without their plans.

    The plan is tens of kilobytes of JSON and the list is drawn on every
    page load, so it is fetched one at a time by `load`.
    """
    purge(conn)
    rows = conn.execute(
        """SELECT id, name, route_id, created_at, saved_at, expires_at,
                  -- rounded, not truncated: a draft saved a moment ago
                  -- has 6.999 days left and should say a week, not six.
                  CAST(ROUND(julianday(expires_at) - julianday('now'))
                       AS INT) AS days_left
             FROM plan_draft WHERE owner = ?
            ORDER BY saved_at DESC, id DESC LIMIT ?""",
        (owner, limit)).fetchall()
    return [dict(r) for r in rows]


def load(conn: sqlite3.Connection, draft_id: int, owner: str = "") -> dict | None:
    """One draft, with its spec and plan, if it is still here and yours."""
    purge(conn)
    row = conn.execute(
        "SELECT * FROM plan_draft WHERE id = ? AND owner = ?",
        (draft_id, owner)).fetchone()
    if row is None:
        return None
    out = dict(row)
    for field in ("spec", "plan"):
        try:
            out[field] = json.loads(out[field])
        except (json.JSONDecodeError, TypeError):
            out[field] = None
    return out


def discard(conn: sqlite3.Connection, draft_id: int, owner: str = "") -> dict:
    """Throw one away early."""
    cur = conn.execute("DELETE FROM plan_draft WHERE id = ? AND owner = ?",
                       (draft_id, owner))
    conn.commit()
    return {"deleted": cur.rowcount}
