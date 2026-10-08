"""SQLite storage.

Design rule from the specification: the backend never modifies the original
DNA files. Each import stores an immutable record of the source (path,
checksum, size, mtime) plus a separately editable normalised record, so a
re-scan can detect drift without ever writing to the share.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

SCHEMA_VERSION = 1

SCHEMA = """
PRAGMA journal_mode=WAL;
PRAGMA foreign_keys=ON;

CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT
);

-- Stage 1: immutable record of a discovered file -------------------------
CREATE TABLE IF NOT EXISTS source_file (
    id            INTEGER PRIMARY KEY,
    rel_path      TEXT UNIQUE NOT NULL,
    owner         TEXT,
    filename      TEXT NOT NULL,
    ext           TEXT,
    format        TEXT,
    role          TEXT,              -- plasmid | read | fragment | empty
    role_reason   TEXT,
    file_checksum TEXT,
    seq_checksum  TEXT,
    size_bytes    INTEGER,
    mtime         REAL,
    imported_at   TEXT DEFAULT CURRENT_TIMESTAMP,
    seq_length    INTEGER,
    has_annotations INTEGER DEFAULT 0,
    parse_warnings TEXT,
    parse_error   TEXT
);
CREATE INDEX IF NOT EXISTS idx_source_owner ON source_file(owner);
CREATE INDEX IF NOT EXISTS idx_source_seqsum ON source_file(seq_checksum);
CREATE INDEX IF NOT EXISTS idx_source_role ON source_file(role);

-- Stage 2: normalised, editable plasmid record ---------------------------
CREATE TABLE IF NOT EXISTS plasmid (
    id            INTEGER PRIMARY KEY,
    name          TEXT NOT NULL,
    owner         TEXT,
    seq           TEXT NOT NULL,
    seq_length    INTEGER NOT NULL,
    seq_checksum  TEXT NOT NULL,
    circular      INTEGER DEFAULT 1,
    topology_inferred INTEGER DEFAULT 0,
    description   TEXT,
    notes         TEXT,
    verified      TEXT DEFAULT 'unverified',  -- unverified|sequence_verified|functional
    last_verified TEXT,
    available_template INTEGER DEFAULT 1,
    annotation_status TEXT DEFAULT 'incomplete',
    source_file_id INTEGER REFERENCES source_file(id) ON DELETE SET NULL,
    created_at    TEXT DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX IF NOT EXISTS idx_plasmid_name ON plasmid(name);
CREATE INDEX IF NOT EXISTS idx_plasmid_checksum ON plasmid(seq_checksum);
CREATE INDEX IF NOT EXISTS idx_plasmid_owner ON plasmid(owner);

-- Alternate source files holding an identical sequence -------------------
CREATE TABLE IF NOT EXISTS plasmid_alt_source (
    plasmid_id     INTEGER REFERENCES plasmid(id) ON DELETE CASCADE,
    source_file_id INTEGER REFERENCES source_file(id) ON DELETE CASCADE,
    PRIMARY KEY (plasmid_id, source_file_id)
);

CREATE TABLE IF NOT EXISTS annotation (
    id           INTEGER PRIMARY KEY,
    plasmid_id   INTEGER NOT NULL REFERENCES plasmid(id) ON DELETE CASCADE,
    feature_type TEXT NOT NULL,
    name         TEXT,
    original_label TEXT,
    original_key   TEXT,
    start        INTEGER NOT NULL,   -- 0-based inclusive
    end          INTEGER NOT NULL,   -- 0-based exclusive
    strand       INTEGER NOT NULL DEFAULT 1,
    evidence     TEXT DEFAULT 'file',   -- file | motif | manual
    confidence   REAL DEFAULT 1.0,
    qualifiers   TEXT
);
CREATE INDEX IF NOT EXISTS idx_annotation_plasmid ON annotation(plasmid_id);
CREATE INDEX IF NOT EXISTS idx_annotation_type ON annotation(feature_type);

-- Reusable biological part, deduplicated by exact sequence ---------------
CREATE TABLE IF NOT EXISTS part (
    id           INTEGER PRIMARY KEY,
    feature_type TEXT NOT NULL,
    name         TEXT NOT NULL,
    seq          TEXT NOT NULL,
    seq_checksum TEXT NOT NULL UNIQUE,
    seq_length   INTEGER NOT NULL,
    translated   INTEGER DEFAULT 0,
    exact_required INTEGER DEFAULT 0,
    variants_ok  INTEGER DEFAULT 1,
    confidence   REAL DEFAULT 1.0,
    notes        TEXT,
    contains     TEXT          -- JSON: elements embedded inside this part
);
CREATE INDEX IF NOT EXISTS idx_part_type ON part(feature_type);
CREATE INDEX IF NOT EXISTS idx_part_name ON part(name);

CREATE TABLE IF NOT EXISTS part_occurrence (
    id         INTEGER PRIMARY KEY,
    part_id    INTEGER NOT NULL REFERENCES part(id) ON DELETE CASCADE,
    plasmid_id INTEGER NOT NULL REFERENCES plasmid(id) ON DELETE CASCADE,
    start      INTEGER NOT NULL,
    end        INTEGER NOT NULL,
    strand     INTEGER NOT NULL DEFAULT 1,
    UNIQUE(part_id, plasmid_id, start, end, strand)
);
CREATE INDEX IF NOT EXISTS idx_occ_part ON part_occurrence(part_id);
CREATE INDEX IF NOT EXISTS idx_occ_plasmid ON part_occurrence(plasmid_id);

-- Stage 4/11: construct specs, designs and outcomes ----------------------
CREATE TABLE IF NOT EXISTS construct (
    id         INTEGER PRIMARY KEY,
    name       TEXT NOT NULL,
    owner      TEXT,
    spec       TEXT NOT NULL,          -- JSON
    created_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS design (
    id           INTEGER PRIMARY KEY,
    construct_id INTEGER REFERENCES construct(id) ON DELETE CASCADE,
    version      INTEGER NOT NULL DEFAULT 1,
    owner        TEXT,
    plan         TEXT NOT NULL,        -- JSON: ranked routes + workflow
    final_seq    TEXT,
    validation   TEXT,                 -- JSON warnings
    created_at   TEXT DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX IF NOT EXISTS idx_design_construct ON design(construct_id);

CREATE TABLE IF NOT EXISTS outcome (
    id          INTEGER PRIMARY KEY,
    design_id   INTEGER REFERENCES design(id) ON DELETE CASCADE,
    route_id    TEXT,
    step_index  INTEGER,
    operation   TEXT,
    result      TEXT,                  -- success | failure | partial
    notes       TEXT,
    recorded_by TEXT,
    recorded_at TEXT DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX IF NOT EXISTS idx_outcome_design ON outcome(design_id);
"""


def connect(db_path: str | Path) -> sqlite3.Connection:
    """Open a connection, creating the schema only the first time.

    Re-running the schema script on every connection would take a write lock,
    which means a reader (the web server opening a connection per thread)
    would block behind a long-running import. Readers therefore only check
    that the schema exists.
    """
    path = Path(db_path)
    if str(path) != ":memory:":
        path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(path), timeout=30.0)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA busy_timeout = 30000")
    conn.execute("PRAGMA foreign_keys = ON")

    already = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='meta'"
    ).fetchone()
    if already:
        _migrate(conn)
        return conn

    conn.executescript(SCHEMA)
    conn.execute(
        "INSERT OR REPLACE INTO meta(key, value) VALUES ('schema_version', ?)",
        (str(SCHEMA_VERSION),),
    )
    conn.commit()
    return conn


#: Columns added after the first release, applied on open so existing
#: databases keep working without a re-import.
_ADDED_COLUMNS = (("part", "contains", "TEXT"),)


def _migrate(conn: sqlite3.Connection) -> None:
    for table, column, decl in _ADDED_COLUMNS:
        cols = {r["name"] for r in conn.execute(f"PRAGMA table_info({table})")}
        if column not in cols:
            conn.execute(f"ALTER TABLE {table} ADD COLUMN {column} {decl}")
            conn.commit()


def set_meta(conn: sqlite3.Connection, key: str, value) -> None:
    conn.execute("INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)",
                 (key, json.dumps(value)))
    conn.commit()


def get_meta(conn: sqlite3.Connection, key: str, default=None):
    row = conn.execute("SELECT value FROM meta WHERE key = ?", (key,)).fetchone()
    if row is None:
        return default
    try:
        return json.loads(row["value"])
    except (json.JSONDecodeError, TypeError):
        return row["value"]
