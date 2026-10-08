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


#: How much of a module's sequence must already be in a plasmid for that
#: plasmid to count as having it. Well short of 1.0 on purpose: the same
#: protein is catalogued at several boundaries, so PP072's "Voltron-2" and a
#: spec's "LR-Voltron2" differ by a leader and are still the same thing to
#: anyone deciding what to clone.
VARIANT_COVERAGE = 0.6
_COVERAGE_K = 40


def sequence_coverage(query: str, host: str) -> float:
    """What fraction of `query` already appears in `host`, either strand.

    Measured over 40-mers rather than by exact containment, because the
    question being asked is "does this plasmid already have this protein",
    not "does it have this exact annotation's boundaries".
    """
    from .seqtools import revcomp

    query = (query or "").upper()
    if len(query) < _COVERAGE_K:
        return 1.0 if query and query in host else 0.0
    index = set()
    both = host + "\x00" + revcomp(host)
    for i in range(len(both) - _COVERAGE_K + 1):
        index.add(both[i:i + _COVERAGE_K])
    total = len(query) - _COVERAGE_K + 1
    hits = sum(1 for i in range(total) if query[i:i + _COVERAGE_K] in index)
    return hits / total if total else 0.0


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

    #: Short, plain reasons this plasmid was offered.
    reasons: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return {"id": self.plasmid_id, "name": self.name, "owner": self.owner,
                "length": self.length, "circular": self.circular,
                "verified": self.verified, "matched": len(self.matched_part_ids),
                "matched_bp": self.matched_bp, "score": round(self.score, 2),
                "reasons": self.reasons}


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
    donor_options: list = field(default_factory=list)
    recipient_options: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "donor_options": self.donor_options,
            "recipient_options": self.recipient_options,
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
                 min_length: int = 1500,
                 require_types: set | None = None,
                 require_family: str = "") -> list[ParentMatch]:
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
              AND pl.circular = 1
              -- A recipient has to be a plasmid that can be grown and
              -- selected. Without an origin and a resistance marker it is a
              -- fragment or an extraction, not something you can transform.
              AND EXISTS (SELECT 1 FROM annotation a
                          WHERE a.plasmid_id = pl.id
                            AND a.feature_type = 'origin')
              AND EXISTS (SELECT 1 FROM annotation a
                          WHERE a.plasmid_id = pl.id
                            AND a.feature_type = 'resistance')
            GROUP BY pl.id
            ORDER BY COUNT(DISTINCT o.part_id) DESC, bp DESC
            LIMIT 400""",
        (*part_ids, min_length)).fetchall()

    require_types = require_types or set()
    matches: list[ParentMatch] = []
    for r in rows:
        # Architecture first. A plasmid that shares three coding sequences
        # but has no lox sites cannot become a Cre-lox construct by swapping
        # one region, however well it scores on shared parts - which is
        # exactly the kind of nonsense recipient this filter exists to stop.
        if require_types:
            present = {a["feature_type"] for a in conn.execute(
                "SELECT DISTINCT feature_type FROM annotation WHERE plasmid_id = ?",
                (r["id"],))}
            if not require_types <= present:
                continue
        if require_family:
            from .vectors import classify_vector
            ann = [dict(a) for a in conn.execute(
                "SELECT feature_type, name FROM annotation WHERE plasmid_id = ?",
                (r["id"],))]
            if classify_vector(ann, r["name"])["family"] != require_family:
                continue
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


def find_donors(conn: sqlite3.Connection, part_ids: list[int],
                require_family: str = "", exclude_id: int | None = None,
                limit: int = 8) -> list[ParentMatch]:
    """Plasmids worth amplifying a part out of.

    Deliberately not ranked the way a recipient is. A recipient is scored on
    how much of the target it already holds, but a donor only has to hold
    one piece - so holding six catalogued variants of GFP does not make a
    plasmid a better place to get GFP from. What matters is how little
    trouble it is: the same vector family as the construct, sequence
    verified, and small.
    """
    if not part_ids:
        return []
    placeholders = ",".join("?" * len(part_ids))
    rows = conn.execute(
        f"""SELECT DISTINCT pl.id, pl.name, pl.owner, pl.seq_length,
                   pl.circular, pl.verified
            FROM part_occurrence o JOIN plasmid pl ON pl.id = o.plasmid_id
            WHERE o.part_id IN ({placeholders})
              AND pl.circular = 1 AND pl.available_template = 1
              AND EXISTS (SELECT 1 FROM annotation a
                          WHERE a.plasmid_id = pl.id
                            AND a.feature_type = 'origin')
              AND EXISTS (SELECT 1 FROM annotation a
                          WHERE a.plasmid_id = pl.id
                            AND a.feature_type = 'resistance')
            ORDER BY pl.seq_length ASC
            LIMIT 80""", part_ids).fetchall()

    from .vectors import classify_vector
    scored = []
    for r in rows:
        if exclude_id is not None and r["id"] == exclude_id:
            continue
        family_miss = 0
        if require_family:
            ann = [dict(a) for a in conn.execute(
                "SELECT feature_type, name FROM annotation WHERE plasmid_id = ?",
                (r["id"],))]
            family_miss = 0 if classify_vector(
                ann, r["name"])["family"] == require_family else 1
        match = ParentMatch(r["id"], r["name"], r["owner"], r["seq_length"],
                            bool(r["circular"]), r["verified"])
        match.reasons = []
        if require_family and not family_miss:
            match.reasons.append(f"same vector family ({require_family})")
        elif require_family:
            match.reasons.append("different vector family")
        if r["verified"] != "unverified":
            match.reasons.append("sequence verified")
        match.reasons.append(f"{r['seq_length']:,} bp")
        scored.append((
            (family_miss, 0 if r["verified"] != "unverified" else 1,
             r["seq_length"]), match))
    scored.sort(key=lambda item: item[0])
    return [match for _, match in scored[:limit]]


def variant_part_ids(conn: sqlite3.Connection, seq: str,
                     part_id: int | None = None) -> list[int]:
    """Every catalogued part that is substantially the same sequence.

    The same protein is catalogued several times over under different names
    and boundaries. Treating those as one thing is what lets a donor search
    find a plasmid carrying "EGFP" when the specification said "GFP".
    """
    ids = [part_id] if part_id else []
    seq = (seq or "").upper()
    if len(seq) < _COVERAGE_K:
        return ids

    # Index both strands. Inside a Cre-lox cassette the payload is stored
    # reverse-complemented, so a module's sequence arrives here flipped and
    # a one-strand index would match nothing at all.
    from .seqtools import revcomp

    index = set()
    for strand_seq in (seq, revcomp(seq)):
        for i in range(len(strand_seq) - _COVERAGE_K + 1):
            index.add(strand_seq[i:i + _COVERAGE_K])

    low, high = int(len(seq) * 0.7), int(len(seq) * 1.4)
    rows = conn.execute(
        """SELECT id, seq FROM part
           WHERE seq_length BETWEEN ? AND ? AND id != COALESCE(?, -1)""",
        (low, high, part_id)).fetchall()
    for row in rows:
        other = row["seq"].upper()
        total = len(other) - _COVERAGE_K + 1
        if total <= 0:
            continue
        hits = sum(1 for i in range(0, total, 7)
                   if other[i:i + _COVERAGE_K] in index)
        sampled = len(range(0, total, 7))
        if sampled and hits / sampled >= VARIANT_COVERAGE:
            ids.append(row["id"])
    return ids


def locate_by_coverage(query: str, host: str) -> dict | None:
    """Find the region of `host` that a variant of `query` occupies.

    Needed because a plasmid can hold a protein without holding the exact
    part that was catalogued for it: PP072 carries Voltron-2, the spec asked
    for LR-Voltron2, and they differ by a leader. The span is taken from the
    first and last 40-mer of the query that appear in the host, which locates
    the region without requiring the boundaries to agree.
    """
    from .seqtools import revcomp

    query = (query or "").upper()
    host = (host or "").upper()
    if len(query) < _COVERAGE_K or not host:
        return None

    for strand, target in ((1, host), (-1, revcomp(host))):
        positions = []
        for i in range(0, len(query) - _COVERAGE_K + 1, 10):
            at = target.find(query[i:i + _COVERAGE_K])
            if at >= 0:
                positions.append(at)
        if len(positions) < 3:
            continue
        lo, hi = min(positions), max(positions) + _COVERAGE_K
        if strand == -1:
            lo, hi = len(host) - hi, len(host) - lo
        return {"start": lo, "end": hi, "strand": strand}
    return None


def _occurrence(conn: sqlite3.Connection, plasmid_id: int, part_id: int,
                seq: str = "") -> dict | None:
    row = conn.execute(
        """SELECT start, end, strand FROM part_occurrence
           WHERE plasmid_id = ? AND part_id = ? ORDER BY start LIMIT 1""",
        (plasmid_id, part_id)).fetchone()
    if row:
        return dict(row)
    # No catalogued occurrence: the plasmid may still hold a variant of this
    # sequence, which is how it came to be counted as having the part.
    if seq:
        host = conn.execute("SELECT seq FROM plasmid WHERE id = ?",
                            (plasmid_id,)).fetchone()
        if host is not None:
            return locate_by_coverage(seq, host["seq"])
    return None


# --------------------------------------------------------------------------
# Working out what to cut and what to insert
# --------------------------------------------------------------------------

def _as_parent(conn: sqlite3.Connection, plasmid_id: int,
               part_ids: list[int]) -> ParentMatch | None:
    row = conn.execute(
        """SELECT id, name, owner, seq_length, circular, verified
           FROM plasmid WHERE id = ?""", (plasmid_id,)).fetchone()
    if row is None:
        return None
    placeholders = ",".join("?" * len(part_ids)) if part_ids else "NULL"
    held = {r["part_id"] for r in conn.execute(
        f"""SELECT DISTINCT part_id FROM part_occurrence
            WHERE plasmid_id = ? AND part_id IN ({placeholders})""",
        (plasmid_id, *part_ids))} if part_ids else set()
    bp = conn.execute(
        f"""SELECT COALESCE(SUM(end - start), 0) AS bp FROM part_occurrence
            WHERE plasmid_id = ? AND part_id IN ({placeholders})""",
        (plasmid_id, *part_ids)).fetchone()["bp"] if part_ids else 0
    return ParentMatch(row["id"], row["name"], row["owner"], row["seq_length"],
                       bool(row["circular"]), row["verified"], held, bp)


def analyse_swap(conn: sqlite3.Connection, resolved_modules: list[dict],
                 cfg: Config, require_types: set | None = None,
                 require_family: str = "",
                 prefer_plasmid_id: int | None = None,
                 pin_recipient_id: int | None = None,
                 pin_donor_id: int | None = None) -> SwapAnalysis | None:
    """Plan the shortest edit: one region out of a recipient, one piece in.

    ``resolved_modules`` is the plan's module list, in target order, each
    with ``part_id``, ``part_name`` and ``type``.
    """
    wanted = [m for m in resolved_modules
              if m.get("part_id") and m.get("type") in IDENTITY_TYPES]
    if not wanted:
        return None
    part_ids = [m["part_id"] for m in wanted]

    candidates = rank_parents(conn, part_ids, cfg,
                              require_types=require_types,
                              require_family=require_family)

    # The plasmid being built into is the obvious thing to edit. Preferring
    # it keeps the plan coherent: proposing a swap out of some other plasmid
    # while assembling into this one describes two different experiments.
    if prefer_plasmid_id is not None:
        preferred = _as_parent(conn, prefer_plasmid_id, part_ids)
        # Catalogue identity is too brittle to decide what a plasmid already
        # has. The same protein exists as several parts with slightly
        # different boundaries, so PP072 holding "Voltron-2" does not make it
        # hold the part id resolved for "LR-Voltron2". Ask the sequence.
        if preferred is not None:
            row = conn.execute("SELECT seq, circular FROM plasmid WHERE id = ?",
                               (prefer_plasmid_id,)).fetchone()
            if row is not None:
                host = row["seq"].upper()
                for module in wanted:
                    covered = sequence_coverage(module.get("seq") or "", host)
                    if covered >= VARIANT_COVERAGE:
                        preferred.matched_part_ids.add(module["part_id"])
                        preferred.matched_bp += len(module.get("seq") or "")
        if preferred is not None and (
                not require_types
                or require_types <= {a["feature_type"] for a in conn.execute(
                    "SELECT DISTINCT feature_type FROM annotation "
                    "WHERE plasmid_id = ?", (prefer_plasmid_id,))}):
            candidates = [preferred] + [c for c in candidates
                                        if c.plasmid_id != prefer_plasmid_id]
    if not candidates:
        if require_types or require_family:
            needed = ", ".join(sorted(require_types)) or require_family
            return SwapAnalysis(
                ParentMatch(0, "", "", 0, False, ""), None, SwapRegion(0, 0),
                None, None, None, [], [],
                [f"no plasmid in the inventory has both the parts this design "
                 f"needs and the architecture it needs ({needed}), so there "
                 f"is nothing to edit; this one has to be built from parts"],
                feasible=False)
        return None
    if pin_recipient_id is not None:
        pinned = _as_parent(conn, pin_recipient_id, part_ids)
        if pinned is not None:
            candidates = [pinned] + [c for c in candidates
                                     if c.plasmid_id != pin_recipient_id]
    recipient = candidates[0]
    recipient_options = [c.to_dict() for c in candidates[:6]]

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
    left_occ = (_occurrence(conn, recipient.plasmid_id, before["part_id"],
                            before.get("seq", "")) if before else None)
    right_occ = (_occurrence(conn, recipient.plasmid_id, after["part_id"],
                             after.get("seq", "")) if after else None)

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
    # Donors are searched by sequence, not by catalogue identity. "GFP" and
    # "EGFP" are separate parts differing by three bases, so asking only for
    # the exact part that was resolved would hide every plasmid that holds
    # the same protein under the other name.
    # Both levels. DNA k-mers catch boundary variants of the same gene;
    # protein comparison catches codon-optimised copies that share almost no
    # DNA, and the badly-labelled ones ("ORF frame 2") that are the protein
    # under no useful name at all.
    donor_parts = sorted(set(
        variant_part_ids(conn, new_module.get("seq") or "",
                         new_module["part_id"])
        + same_protein_parts(conn, new_module.get("seq") or "",
                             new_module["part_id"],
                             new_module.get("type", "CDS"))))
    donors = find_donors(conn, donor_parts, require_family,
                         exclude_id=recipient.plasmid_id)
    if pin_donor_id is not None:
        donors = ([d for d in donors if d.plasmid_id == pin_donor_id] or
                  [_as_parent(conn, pin_donor_id, [new_module["part_id"]])]) + \
                 [d for d in donors if d.plasmid_id != pin_donor_id]
        donors = [d for d in donors if d is not None]
    donor = donors[0] if donors else None
    donor_options = [d.to_dict() for d in donors[:6]]
    insert = None
    if donor:
        occ = _occurrence(conn, donor.plasmid_id, new_module["part_id"],
                          new_module.get("seq", ""))
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

    analysis = SwapAnalysis(
        recipient=recipient, donor=donor, remove=remove, insert=insert,
        left_anchor=anchor_left, right_anchor=anchor_right,
        kept=[m["part_name"] for m in kept],
        missing=[m["part_name"] for m in missing],
        notes=notes, feasible=insert is not None)
    # Which donor to use is a lab preference as much as a computation: a
    # dozen plasmids may carry the same protein. The alternatives travel
    # with the plan so the choice stays visible.
    analysis.donor_options = donor_options
    analysis.recipient_options = recipient_options
    return analysis


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
    scheme["donor_options"] = list(getattr(analysis, "donor_options", []))
    scheme["recipient_options"] = list(getattr(analysis, "recipient_options", []))

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


# --------------------------------------------------------------------------
# Protein-level identity
# --------------------------------------------------------------------------

#: Two coding sequences are the same part if their proteins agree this far.
PROTEIN_IDENTITY = 0.95


def protein_of(seq: str) -> str:
    """The protein a coding sequence encodes, read in its best frame.

    Taken as the longest uninterrupted translation across the three forward
    frames, which is what a CDS annotation means even when its boundaries
    are a few bases out.
    """
    from .seqtools import translate

    best = ""
    for frame in range(3):
        for stretch in translate(seq, frame=frame).split("*"):
            if len(stretch) > len(best):
                best = stretch
    return best


def protein_identity(a: str, b: str) -> float:
    """Fraction of the shorter protein shared with the longer one.

    Compared by residue windows rather than aligned properly: the question
    is "are these the same protein under two names", not "what is the exact
    alignment", and a window index answers that without an alignment matrix.
    """
    if not a or not b:
        return 0.0
    short, long_ = (a, b) if len(a) <= len(b) else (b, a)
    k = 8
    if len(short) < k:
        return 1.0 if short in long_ else 0.0
    index = {long_[i:i + k] for i in range(len(long_) - k + 1)}
    total = len(short) - k + 1
    hits = sum(1 for i in range(total) if short[i:i + k] in index)
    return hits / total if total else 0.0


def same_protein_parts(conn: sqlite3.Connection, seq: str,
                       part_id: int | None = None,
                       feature_type: str = "CDS") -> list[int]:
    """Catalogued parts encoding the same protein as `seq`.

    DNA comparison misses codon-optimised versions of a protein entirely:
    two genes can share almost no 40-mers and still make the same thing.
    Asking at the protein level is what groups "GFP", "EGFP" and a
    human-codon-optimised copy of either into one answer.
    """
    target = protein_of(seq or "")
    if len(target) < 30:
        return [part_id] if part_id else []

    ids = [part_id] if part_id else []
    low, high = int(len(seq) * 0.5), int(len(seq) * 2.0)
    rows = conn.execute(
        """SELECT id, seq FROM part
           WHERE feature_type = ? AND seq_length BETWEEN ? AND ?
             AND id != COALESCE(?, -1)""",
        (feature_type, low, high, part_id)).fetchall()
    for row in rows:
        if protein_identity(target, protein_of(row["seq"])) >= PROTEIN_IDENTITY:
            ids.append(row["id"])
    return ids
