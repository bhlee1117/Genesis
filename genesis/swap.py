"""Two-parent swap planning: the way constructs are actually built.

Most new plasmids are not assembled from parts. Someone finds the existing
plasmid that is closest to what they want, cuts one region out of it, and
drops in one coding sequence taken from a second plasmid. Two parents, one
junction pair, one assembly.

This module finds those two parents in the inventory and works out:

* **recipient** - the closest existing plasmid, and the region to remove
* **donor**     - where the new piece comes from, and the region to take
* **boundaries** - the conserved parts either side of the cut, which is what
  the primers anneal to or the enzymes cut between

It is deliberately separate from the generic multi-part planner: that one
answers "how would I build this from scratch", this one answers "what is the
fastest edit to something we already have", which is usually the real
question.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass, field

from .config import Config
from .seqtools import revcomp, slice_region

#: Module types that define the identity of a construct. Structural elements
#: (origin, resistance, ITRs) are carried by the recipient and never swapped.
IDENTITY_TYPES = {"promoter", "enhancer", "CDS", "trafficking_motif", "p2A",
                  "IRES", "DIO", "tag", "linker", "polyA", "WPRE"}


@dataclass
class ParentMatch:
    plasmid_id: int
    name: str
    owner: str
    length: int
    circular: bool
    verified: str
    matched_part_ids: set = field(default_factory=set)
    matched_bp: int = 0
    score: float = 0.0

    def to_dict(self) -> dict:
        return {"id": self.plasmid_id, "name": self.name, "owner": self.owner,
                "length": self.length, "circular": self.circular,
                "verified": self.verified, "matched": len(self.matched_part_ids),
                "matched_bp": self.matched_bp, "score": round(self.score, 2)}


@dataclass
class SwapRegion:
    """A span on a parent plasmid, with the features it covers."""
    start: int
    end: int
    strand: int = 1
    label: str = ""
    features: list = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start

    def to_dict(self) -> dict:
        return {"start": self.start, "end": self.end, "strand": self.strand,
                "label": self.label, "length": self.length,
                "features": self.features}


@dataclass
class SwapAnalysis:
    recipient: ParentMatch
    donor: ParentMatch | None
    remove: SwapRegion                 # region cut out of the recipient
    insert: SwapRegion | None          # region taken from the donor
    left_anchor: dict | None           # conserved feature before the cut
    right_anchor: dict | None          # conserved feature after the cut
    kept: list = field(default_factory=list)
    missing: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    feasible: bool = True

    def to_dict(self) -> dict:
        return {
            "recipient": self.recipient.to_dict() if self.recipient else None,
            "donor": self.donor.to_dict() if self.donor else None,
            "remove": self.remove.to_dict() if self.remove else None,
            "insert": self.insert.to_dict() if self.insert else None,
            "left_anchor": self.left_anchor, "right_anchor": self.right_anchor,
            "kept": self.kept, "missing": self.missing,
            "notes": self.notes, "feasible": self.feasible,
        }


# --------------------------------------------------------------------------
# Finding the closest parent
# --------------------------------------------------------------------------

def rank_parents(conn: sqlite3.Connection, part_ids: list[int],
                 cfg: Config, limit: int = 8,
                 min_length: int = 1500) -> list[ParentMatch]:
    """Rank inventory plasmids by how much of the target they already contain.

    Scoring favours a plasmid that carries many of the wanted parts and few
    surprises: matched parts count most, matched base pairs break ties, and
    a verified circular plasmid is preferred over an unverified fragment.
    """
    if not part_ids:
        return []
    placeholders = ",".join("?" * len(part_ids))
    rows = conn.execute(
        f"""SELECT pl.id, pl.name, pl.owner, pl.seq_length, pl.circular,
                   pl.verified,
                   GROUP_CONCAT(DISTINCT o.part_id) AS parts,
                   SUM(o.end - o.start) AS bp
            FROM part_occurrence o
            JOIN plasmid pl ON pl.id = o.plasmid_id
            WHERE o.part_id IN ({placeholders})
              AND pl.seq_length >= ?
              AND pl.available_template = 1
            GROUP BY pl.id
            ORDER BY COUNT(DISTINCT o.part_id) DESC, bp DESC
            LIMIT 400""",
        (*part_ids, min_length)).fetchall()

    matches: list[ParentMatch] = []
    for r in rows:
        ids = {int(x) for x in (r["parts"] or "").split(",") if x}
        m = ParentMatch(r["id"], r["name"], r["owner"], r["seq_length"],
                        bool(r["circular"]), r["verified"], ids, r["bp"] or 0)
        m.score = (len(ids) * 100
                   + (m.matched_bp / 1000.0)
                   + (20 if m.circular else 0)
                   + (15 if m.verified != "unverified" else 0))
        matches.append(m)

    matches.sort(key=lambda m: -m.score)
    return matches[:limit]


def _annotations(conn: sqlite3.Connection, plasmid_id: int) -> list[dict]:
    return [dict(r) for r in conn.execute(
        """SELECT feature_type, name, start, end, strand FROM annotation
           WHERE plasmid_id = ? ORDER BY start""", (plasmid_id,)).fetchall()]


def _occurrence(conn: sqlite3.Connection, plasmid_id: int,
                part_id: int) -> dict | None:
    row = conn.execute(
        """SELECT start, end, strand FROM part_occurrence
           WHERE plasmid_id = ? AND part_id = ? ORDER BY start LIMIT 1""",
        (plasmid_id, part_id)).fetchone()
    return dict(row) if row else None


# --------------------------------------------------------------------------
# Working out what to cut and what to insert
# --------------------------------------------------------------------------

def analyse_swap(conn: sqlite3.Connection, resolved_modules: list[dict],
                 cfg: Config) -> SwapAnalysis | None:
    """Plan the shortest edit: one region out of a recipient, one piece in.

    ``resolved_modules`` is the plan's module list, in target order, each
    with ``part_id``, ``part_name`` and ``type``.
    """
    wanted = [m for m in resolved_modules
              if m.get("part_id") and m.get("type") in IDENTITY_TYPES]
    if not wanted:
        return None
    part_ids = [m["part_id"] for m in wanted]

    candidates = rank_parents(conn, part_ids, cfg)
    if not candidates:
        return None
    recipient = candidates[0]

    kept = [m for m in wanted if m["part_id"] in recipient.matched_part_ids]
    missing = [m for m in wanted if m["part_id"] not in recipient.matched_part_ids]
    notes: list[str] = []

    if not kept:
        return None
    if not missing:
        notes.append(f"{recipient.name} already contains every requested part; "
                     f"check whether this construct already exists before "
                     f"cloning anything")
        return SwapAnalysis(recipient, None, SwapRegion(0, 0), None, None, None,
                            [m["part_name"] for m in kept], [], notes,
                            feasible=False)

    new_module = missing[0]

    if len(missing) > 1:
        notes.append(
            f"{len(missing)} parts are absent from {recipient.name} "
            f"({', '.join(m['part_name'] for m in missing)}); a single swap "
            f"only replaces one region, so the extra pieces still need their "
            f"own fragments")

    # --- locate the region to remove ------------------------------------
    # The new piece goes where the target order says it goes: between the
    # kept part immediately before it and the kept part immediately after.
    target_order = wanted
    idx = target_order.index(missing[0])
    before = next((m for m in reversed(target_order[:idx])
                   if m["part_id"] in recipient.matched_part_ids), None)
    after = next((m for m in target_order[idx + 1:]
                  if m["part_id"] in recipient.matched_part_ids), None)

    annots = _annotations(conn, recipient.plasmid_id)
    left_occ = _occurrence(conn, recipient.plasmid_id, before["part_id"]) if before else None
    right_occ = _occurrence(conn, recipient.plasmid_id, after["part_id"]) if after else None

    if left_occ is None and right_occ is None:
        notes.append("could not locate the flanking parts in the recipient")
        return None

    if left_occ and right_occ and left_occ["end"] <= right_occ["start"]:
        start, end = left_occ["end"], right_occ["start"]
    elif left_occ:
        # No conserved part downstream: replace everything up to the next
        # annotated feature after the anchor.
        start = left_occ["end"]
        following = [a for a in annots if a["start"] >= start]
        end = following[0]["end"] if following else min(start + 1500,
                                                        recipient.length)
    else:
        end = right_occ["start"]
        preceding = [a for a in annots if a["end"] <= end]
        start = preceding[-1]["start"] if preceding else max(end - 1500, 0)

    covered = [a for a in annots
               if a["start"] >= start and a["end"] <= end
               and a["feature_type"] in IDENTITY_TYPES]
    remove = SwapRegion(
        start, end, 1,
        label=(", ".join(a["name"] for a in covered[:3]) or "intervening region"),
        features=covered)

    # The region being removed must actually be the kind of thing being
    # replaced. Cutting out a 21 bp export motif to put a promoter in its
    # place is arithmetically possible and biologically meaningless, so a
    # mismatch of role means there is no sensible swap here.
    if covered:
        removed_types = {a["feature_type"] for a in covered}
        if new_module["type"] not in removed_types:
            notes.append(
                f"the closest parent {recipient.name} does not hold a "
                f"{new_module['type']} in the position where "
                f"{new_module['part_name']} belongs - the region there is "
                f"{', '.join(sorted(removed_types))} - so there is no "
                f"single-region swap that builds this construct")
            return SwapAnalysis(recipient, None, remove, None, None, None,
                                [m["part_name"] for m in kept],
                                [m["part_name"] for m in missing], notes,
                                feasible=False)

    if remove.length <= 0:
        notes.append("the target order places the new part where the recipient "
                     "has no removable region; this is an insertion rather "
                     "than a swap")
        remove = SwapRegion(start, start, 1, label="insertion point")

    # --- locate the piece to take from a donor ---------------------------
    donors = rank_parents(conn, [new_module["part_id"]], cfg, limit=5)
    donors = [d for d in donors if d.plasmid_id != recipient.plasmid_id]
    donor = donors[0] if donors else None
    insert = None
    if donor:
        occ = _occurrence(conn, donor.plasmid_id, new_module["part_id"])
        if occ:
            insert = SwapRegion(occ["start"], occ["end"], occ["strand"],
                                label=new_module["part_name"])
    if insert is None:
        notes.append(f"no lab plasmid carries {new_module['part_name']}; it "
                     f"would have to be synthesised")

    anchor_left = ({"name": before["part_name"], "part_id": before["part_id"],
                    **left_occ} if (before and left_occ) else None)
    anchor_right = ({"name": after["part_name"], "part_id": after["part_id"],
                     **right_occ} if (after and right_occ) else None)

    return SwapAnalysis(
        recipient=recipient, donor=donor, remove=remove, insert=insert,
        left_anchor=anchor_left, right_anchor=anchor_right,
        kept=[m["part_name"] for m in kept],
        missing=[m["part_name"] for m in missing],
        notes=notes, feasible=insert is not None)


# --------------------------------------------------------------------------
# Drawing data for the scheme
# --------------------------------------------------------------------------

def build_scheme(conn: sqlite3.Connection, analysis: SwapAnalysis) -> dict:
    """Everything a renderer needs to draw the two parents and the product."""
    if analysis is None:
        return {}

    def track(plasmid_id: int, highlight: SwapRegion | None, kind: str) -> dict:
        row = conn.execute(
            """SELECT id, name, owner, seq_length, circular FROM plasmid
               WHERE id = ?""", (plasmid_id,)).fetchone()
        if row is None:
            return {}
        feats = [a for a in _annotations(conn, plasmid_id)
                 if a["feature_type"] in IDENTITY_TYPES
                 or a["feature_type"] in ("origin", "resistance", "AAV_ITR")]
        # Keep the map readable: largest features first, capped.
        feats.sort(key=lambda a: -(a["end"] - a["start"]))
        feats = sorted(feats[:14], key=lambda a: a["start"])
        return {
            "role": kind,
            "id": row["id"], "name": row["name"], "owner": row["owner"],
            "length": row["seq_length"], "circular": bool(row["circular"]),
            "features": feats,
            "highlight": highlight.to_dict() if highlight else None,
        }

    scheme = {
        "recipient": track(analysis.recipient.plasmid_id, analysis.remove,
                           "recipient"),
        "donor": (track(analysis.donor.plasmid_id, analysis.insert, "donor")
                  if analysis.donor else None),
        "left_anchor": analysis.left_anchor,
        "right_anchor": analysis.right_anchor,
        "remove": analysis.remove.to_dict() if analysis.remove else None,
        "insert": analysis.insert.to_dict() if analysis.insert else None,
    }
    if scheme["recipient"]:
        scheme["recipient"]["action"] = "cut out and discard this region"
    if scheme["donor"]:
        scheme["donor"]["action"] = "take this region"

    # The product of a swap is the recipient with one region exchanged, which
    # is a different molecule from the part-by-part assembly the generic
    # planner builds. Report its real size so the two are not confused.
    rec = scheme.get("recipient") or {}
    if rec and analysis.remove and analysis.insert:
        scheme["product"] = {
            "length": rec["length"] - analysis.remove.length
                      + analysis.insert.length,
            "circular": rec.get("circular", True),
            "basis": rec["name"],
            "removed": analysis.remove.length,
            "added": analysis.insert.length,
            "insert_label": analysis.insert.label,
        }
    return scheme
