"""Stages 1-2: discover files, parse them, normalise plasmids and parts.

The importer is deliberately conservative:

* It reads the share but never writes to it.
* Unannotated sequences are marked ``annotation_incomplete`` rather than
  being silently guessed at.
* Identical sequences are collapsed into one plasmid with multiple alternate
  source files, instead of becoming duplicate templates.
* Sanger reads are stored but kept out of the template pool.
"""

from __future__ import annotations

import hashlib
import re
import sqlite3
import time
from dataclasses import dataclass, field
from pathlib import Path

from .config import Config
from .seqio import (READABLE_EXT, ParseError, Record, classify_role,
                    is_hidden_or_sidecar, read_record, sniff_format)
from .seqtools import revcomp
from .vocab import (classify_label, find_embedded_elements,
                    find_motifs, min_part_length, split_at_elements)


@dataclass
class ImportStats:
    scanned: int = 0
    parsed: int = 0
    failed: int = 0
    plasmids: int = 0
    duplicates: int = 0
    reads: int = 0
    annotations: int = 0
    parts: int = 0
    errors: list = field(default_factory=list)

    def as_dict(self) -> dict:
        d = {k: v for k, v in self.__dict__.items() if k != "errors"}
        d["errors"] = self.errors[:50]
        return d


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def least_rotation(s: str) -> str:
    """Lexicographically smallest rotation of `s` (Booth's algorithm, O(n)).

    Used so that the same circular plasmid saved with different origins
    produces the same checksum. The naive "min of all rotations" is O(n^2)
    and is far too slow on a share full of 10 kb AAV plasmids.
    """
    n = len(s)
    if n < 2:
        return s
    t = s + s
    f = [-1] * len(t)
    k = 0
    for j in range(1, len(t)):
        sj = t[j]
        i = f[j - k - 1]
        while i != -1 and sj != t[k + i + 1]:
            if sj < t[k + i + 1]:
                k = j - i - 1
            i = f[i]
        if sj != t[k + i + 1]:
            if sj < t[k]:
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    return t[k:k + n]


def seq_checksum(seq: str, circular: bool = False) -> str:
    """Checksum that is stable under reverse-complement, and under rotation
    for circular molecules, so the same plasmid saved two ways collapses."""
    seq = seq.upper()
    if not seq:
        return hashlib.sha256(b"").hexdigest()
    candidates = [seq, revcomp(seq)]
    if circular:
        canonical = min(least_rotation(s) for s in candidates)
    else:
        canonical = min(candidates)
    return hashlib.sha256(canonical.encode()).hexdigest()


def infer_owner(rel_path: Path, cfg: Config) -> str:
    """Top-level folder is the person. Normalise 'X Clones'/'X clones'/'X'."""
    parts = rel_path.parts
    if not parts or len(parts) < 2:
        return "(root)"
    top = parts[0]
    if top in cfg.non_person_dirs:
        return top
    # Folder names on the share are inconsistent: "Bill Jia clones",
    # "Alex Cordova", "Veena Venkatachalam Clones (locations given in ...)".
    # Strip any trailing parenthetical note first, then the collection suffix.
    cleaned = re.sub(r"\s*\([^)]*\)\s*$", "", top).strip()
    for suffix in (" Clones", " clones", " CLONES", " Clone", " clone",
                   " plasmids", " Plasmids", " stocks", " Stocks"):
        if cleaned.endswith(suffix):
            cleaned = cleaned[: -len(suffix)]
            break
    return cleaned.strip() or top


def iter_candidate_files(root: Path, cfg: Config):
    """Walk the share, skipping sidecars, hidden dirs and oversized files."""
    root = Path(root)
    stack = [(root, 0)]
    while stack:
        current, depth = stack.pop()
        if depth > cfg.max_depth:
            continue
        try:
            entries = list(current.iterdir())
        except (OSError, PermissionError):
            continue
        for entry in entries:
            name = entry.name
            if is_hidden_or_sidecar(entry):
                continue
            try:
                if entry.is_dir():
                    if any(pat in name for pat in cfg.skip_dir_patterns):
                        continue
                    stack.append((entry, depth + 1))
                elif entry.is_file():
                    if entry.suffix.lower() not in READABLE_EXT:
                        continue
                    yield entry
            except OSError:
                continue


# --------------------------------------------------------------------------
# Annotation + part extraction
# --------------------------------------------------------------------------

def annotations_for(rec: Record, run_motif_scan: bool = True) -> list[dict]:
    """Combine file-borne annotations with de-novo motif hits."""
    out: list[dict] = []
    for f in rec.features:
        ftype = classify_label(f.label, f.type)
        out.append({
            "feature_type": ftype,
            "name": f.label or f.type,
            "original_label": f.label,
            "original_key": f.type,
            "start": f.start, "end": f.end, "strand": f.strand,
            "evidence": "file", "confidence": 1.0,
            "qualifiers": f.qualifiers,
        })

    if run_motif_scan and rec.seq:
        existing = {(a["start"], a["end"], a["strand"]) for a in out}
        for m in find_motifs(rec.seq, circular=rec.circular):
            key = (m["start"], m["end"], m["strand"])
            if key in existing:
                continue
            # Skip motif hits that merely restate a file annotation.
            overlaps = any(
                a["feature_type"] == m["feature_type"]
                and not (a["end"] <= m["start"] or a["start"] >= m["end"])
                for a in out
            )
            if overlaps:
                continue
            out.append({
                "feature_type": m["feature_type"],
                "name": m["name"],
                "original_label": "",
                "original_key": f"motif:{m['level']}",
                "start": m["start"], "end": m["end"], "strand": m["strand"],
                "evidence": "motif",
                "confidence": 0.9 if m["level"] == "protein" else 0.85,
                "qualifiers": {"note": m.get("note", "")},
            })
    return out


#: Feature types worth promoting into the reusable-part catalogue.
PART_TYPES = {"promoter", "CDS", "trafficking_motif", "linker", "p2A", "IRES",
              "lox_site", "AAV_ITR", "polyA", "WPRE", "tag", "enhancer",
              "resistance", "origin", "DIO"}

TRANSLATED_TYPES = {"CDS", "trafficking_motif", "linker", "p2A", "tag"}


def _json_dumps(value):
    import json
    return json.dumps(value)


def extract_parts(conn: sqlite3.Connection, plasmid_id: int, rec: Record,
                  annotations: list[dict]) -> int:
    """Promote annotations into deduplicated reusable parts (Stage 2)."""
    from .seqtools import slice_region

    created = 0
    for a in annotations:
        if a["feature_type"] not in PART_TYPES:
            continue
        length = a["end"] - a["start"]
        if length < min_part_length(a["feature_type"]) or length > 20000:
            # Too short to be a believable instance of this feature type.
            continue
        sub = slice_region(rec.seq, a["start"], a["end"], circular=rec.circular)
        if a["strand"] == -1:
            sub = revcomp(sub)
        if not sub:
            continue
        chk = hashlib.sha256(sub.upper().encode()).hexdigest()
        row = conn.execute("SELECT id FROM part WHERE seq_checksum = ?",
                           (chk,)).fetchone()
        if row:
            part_id = row["id"]
        else:
            embedded = find_embedded_elements(sub, a["feature_type"])
            pieces = split_at_elements(sub, embedded)
            cur = conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                                    seq_length, translated, exact_required,
                                    confidence, contains)
                   VALUES (?,?,?,?,?,?,?,?,?)""",
                (a["feature_type"], a["name"] or a["feature_type"], sub, chk,
                 len(sub), int(a["feature_type"] in TRANSLATED_TYPES),
                 int(a["feature_type"] in ("lox_site", "AAV_ITR", "p2A")),
                 a["confidence"],
                 _json_dumps(embedded) if embedded else None),
            )
            part_id = cur.lastrowid
            created += 1
        conn.execute(
            """INSERT OR IGNORE INTO part_occurrence
               (part_id, plasmid_id, start, end, strand) VALUES (?,?,?,?,?)""",
            (part_id, plasmid_id, a["start"], a["end"], a["strand"]),
        )

        # A fused CDS is catalogued twice: whole, and split into the proteins
        # it actually contains. Both are real answers to "what part is this",
        # so neither replaces the other.
        if row is None and pieces:
            created += _store_segments(conn, plasmid_id, a, sub, pieces,
                                       part_id)
    return created


def _store_segments(conn: sqlite3.Connection, plasmid_id: int, annotation: dict,
                    parent_seq: str, pieces: list[dict],
                    parent_id: int) -> int:
    """Catalogue the pieces a fused CDS splits into.

    A segment whose sequence already exists as a part is linked to that part
    rather than duplicated, which is how the split pieces inherit the names
    the lab already uses for them.
    """
    created = 0
    coding = [p for p in pieces if p["kind"] == "segment"]
    total = len(coding)
    index = 0

    for piece in pieces:
        if piece["kind"] != "segment":
            continue
        index += 1
        seq = piece["seq"]
        if len(seq) < min_part_length("CDS"):
            continue
        chk = hashlib.sha256(seq.upper().encode()).hexdigest()
        existing = conn.execute("SELECT id FROM part WHERE seq_checksum = ?",
                                (chk,)).fetchone()
        if existing:
            part_id = existing["id"]
        else:
            cur = conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length, translated, confidence, derived_from,
                        segment_index, segment_of)
                   VALUES (?,?,?,?,?,1,?,?,?,?)""",
                ("CDS",
                 f"{annotation['name']} [{index} of {total}]",
                 seq, chk, len(seq), 0.8, parent_id, index, total))
            part_id = cur.lastrowid
            created += 1

        # Position of this piece within the plasmid, honouring strand.
        if annotation["strand"] == 1:
            start = annotation["start"] + piece["start"]
            end = annotation["start"] + piece["end"]
        else:
            start = annotation["end"] - piece["end"]
            end = annotation["end"] - piece["start"]
        conn.execute(
            """INSERT OR IGNORE INTO part_occurrence
               (part_id, plasmid_id, start, end, strand) VALUES (?,?,?,?,?)""",
            (part_id, plasmid_id, start, end, annotation["strand"]))
    return created


def infer_topology(rec: Record, annotations: list[dict]) -> tuple[bool, bool]:
    """Guess whether a molecule is really circular.

    ApE in particular writes ``linear`` for plasmids that are plainly
    circular. A record carrying both an origin and a resistance marker is
    almost certainly a plasmid, so restriction planning should treat it as
    circular. Returns ``(circular, was_inferred)``.
    """
    if rec.circular:
        return True, False
    types = {a["feature_type"] for a in annotations}
    if {"origin", "resistance"} <= types and len(rec.seq) >= 1500:
        return True, True
    return False, False


# --------------------------------------------------------------------------
# Main entry point
# --------------------------------------------------------------------------

def import_tree(conn: sqlite3.Connection, cfg: Config, root: str | Path | None = None,
                limit: int | None = None, progress=None,
                motif_scan: bool = True) -> ImportStats:
    root = Path(root or cfg.dna_root)
    stats = ImportStats()
    if not root.exists():
        stats.errors.append(f"DNA root not reachable: {root}")
        return stats

    last_report = time.time()
    for path in iter_candidate_files(root, cfg):
        if limit is not None and stats.scanned >= limit:
            break
        stats.scanned += 1
        if progress and time.time() - last_report > 2.0:
            progress(stats)
            last_report = time.time()

        try:
            rel = path.relative_to(root)
        except ValueError:
            rel = path
        rel_str = str(rel)

        try:
            st = path.stat()
        except OSError as exc:
            stats.failed += 1
            stats.errors.append(f"{rel_str}: {exc}")
            continue

        if st.st_size > cfg.max_file_bytes:
            continue

        existing = conn.execute(
            "SELECT id, file_checksum, mtime FROM source_file WHERE rel_path = ?",
            (rel_str,)).fetchone()
        if existing and existing["mtime"] == st.st_mtime:
            continue  # unchanged since last scan

        owner = infer_owner(rel, cfg)
        fmt = sniff_format(path)

        try:
            rec = read_record(path)
            stats.parsed += 1
        except (ParseError, OSError, ValueError, MemoryError) as exc:
            stats.failed += 1
            stats.errors.append(f"{rel_str}: {type(exc).__name__}: {exc}")
            conn.execute(
                """INSERT INTO source_file(rel_path, owner, filename, ext,
                        format, role, file_checksum, size_bytes, mtime,
                        parse_error)
                   VALUES (?,?,?,?,?,?,?,?,?,?)
                   ON CONFLICT(rel_path) DO UPDATE SET
                        parse_error=excluded.parse_error, mtime=excluded.mtime""",
                (rel_str, owner, path.name, path.suffix.lower(), fmt, "error",
                 None, st.st_size, st.st_mtime, str(exc)[:500]),
            )
            continue

        role, reason = classify_role(path, rec)
        chk_file = sha256_file(path)
        chk_seq = seq_checksum(rec.seq, rec.circular) if rec.seq else ""

        cur = conn.execute(
            """INSERT INTO source_file(rel_path, owner, filename, ext, format,
                    role, role_reason, file_checksum, seq_checksum, size_bytes,
                    mtime, seq_length, has_annotations, parse_warnings)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(rel_path) DO UPDATE SET
                    owner=excluded.owner, format=excluded.format,
                    role=excluded.role, role_reason=excluded.role_reason,
                    file_checksum=excluded.file_checksum,
                    seq_checksum=excluded.seq_checksum,
                    size_bytes=excluded.size_bytes, mtime=excluded.mtime,
                    seq_length=excluded.seq_length,
                    has_annotations=excluded.has_annotations,
                    parse_warnings=excluded.parse_warnings,
                    parse_error=NULL
               RETURNING id""",
            (rel_str, owner, path.name, path.suffix.lower(), fmt, role, reason,
             chk_file, chk_seq, st.st_size, st.st_mtime, len(rec.seq),
             int(bool(rec.features)), "; ".join(rec.warnings)[:1000]),
        )
        source_id = cur.fetchone()[0]

        if role == "read":
            stats.reads += 1
            continue
        if role in ("empty", "fragment") and not rec.features:
            continue

        # --- duplicate detection by sequence --------------------------------
        dup = conn.execute(
            "SELECT id FROM plasmid WHERE seq_checksum = ?", (chk_seq,)).fetchone()
        if dup:
            conn.execute(
                """INSERT OR IGNORE INTO plasmid_alt_source(plasmid_id, source_file_id)
                   VALUES (?,?)""", (dup["id"], source_id))
            stats.duplicates += 1
            continue

        annots = annotations_for(rec, run_motif_scan=motif_scan)
        circular, inferred = infer_topology(rec, annots)
        ann_status = "complete" if rec.features else "incomplete"

        cur = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, topology_inferred, description, annotation_status,
                    source_file_id, available_template)
               VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
            (rec.name or path.stem, owner, rec.seq, len(rec.seq), chk_seq,
             int(circular), int(inferred), rec.description, ann_status,
             source_id, 1),
        )
        plasmid_id = cur.lastrowid
        stats.plasmids += 1

        for a in annots:
            import json as _json
            conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        original_label, original_key, start, end, strand,
                        evidence, confidence, qualifiers)
                   VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
                (plasmid_id, a["feature_type"], a["name"], a["original_label"],
                 a["original_key"], a["start"], a["end"], a["strand"],
                 a["evidence"], a["confidence"], _json.dumps(a["qualifiers"])),
            )
        stats.annotations += len(annots)
        stats.parts += extract_parts(conn, plasmid_id, rec, annots)

        if stats.scanned % 200 == 0:
            conn.commit()

    conn.commit()
    return stats
