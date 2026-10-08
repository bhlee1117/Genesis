"""Sequence file readers: GenBank/ApE, FASTA, SnapGene, plain text.

Grounded in what the Cohen Lab share actually contains:

* ``.ape`` and ``.str`` are ApE files, which are GenBank underneath.
* Many old files use CR-only or CRLF line endings, so every reader
  normalises newlines before parsing.
* ``.seq`` files are usually Sanger reads (FASTA with an ``.ab1`` header and
  long N runs), not plasmids, so they are classified as reads.
* ``._*`` AppleDouble sidecars must be ignored.
"""

from __future__ import annotations

import re
import struct
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

from .seqtools import clean_seq

# Extensions we attempt to read, mapped to a parser family.
GENBANK_EXT = {".gb", ".gbk", ".genbank", ".ape", ".str", ".ap"}
FASTA_EXT = {".fa", ".fasta", ".fas", ".fsa", ".seq", ".txt", ".nucl"}
SNAPGENE_EXT = {".dna"}
READABLE_EXT = GENBANK_EXT | FASTA_EXT | SNAPGENE_EXT


class ParseError(Exception):
    pass


@dataclass
class Feature:
    type: str              # original GenBank/SnapGene feature key
    start: int             # 0-based inclusive
    end: int               # 0-based exclusive
    strand: int            # +1 or -1
    label: str = ""        # human readable name
    qualifiers: dict = field(default_factory=dict)
    segments: list = field(default_factory=list)  # for joined locations

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Record:
    name: str
    seq: str
    circular: bool = False
    features: list[Feature] = field(default_factory=list)
    description: str = ""
    source_format: str = ""
    warnings: list[str] = field(default_factory=list)
    notes: str = ""

    @property
    def length(self) -> int:
        return len(self.seq)


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def _normalize_text(raw: bytes) -> str:
    """Decode and normalise CR / CRLF / LF line endings."""
    for enc in ("utf-8", "latin-1"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    else:  # pragma: no cover - latin-1 never fails
        text = raw.decode("utf-8", errors="replace")
    return text.replace("\r\n", "\n").replace("\r", "\n")


def is_hidden_or_sidecar(path: Path) -> bool:
    return path.name.startswith("._") or path.name.startswith(".")


def sniff_format(path: Path, head: bytes | None = None) -> str:
    """Identify the file family from magic bytes first, extension second."""
    if head is None:
        try:
            with open(path, "rb") as fh:
                head = fh.read(512)
        except OSError:
            return "unknown"

    if len(head) > 9 and head[0] == 0x09 and head[5:13] == b"SnapGene":
        return "snapgene"

    text = _normalize_text(head).lstrip()
    if text.startswith("LOCUS"):
        return "genbank"
    if text.startswith(">"):
        return "fasta"

    ext = path.suffix.lower()
    if ext in SNAPGENE_EXT:
        return "snapgene"
    if ext in GENBANK_EXT:
        return "genbank"
    if ext in FASTA_EXT:
        return "fasta"
    return "unknown"


# --------------------------------------------------------------------------
# GenBank / ApE
# --------------------------------------------------------------------------

_LOCUS_RE = re.compile(
    r"^LOCUS\s+(?P<name>\S+)\s+(?P<length>\d+)\s*bp(?P<rest>.*)$", re.IGNORECASE
)
_LOC_NUM_RE = re.compile(r"[<>]?(\d+)")


def parse_location(loc: str) -> tuple[list[tuple[int, int]], int]:
    """Parse a GenBank location string.

    Returns ``(segments, strand)`` with 0-based half-open segments. Handles
    ``complement``, ``join``, ``order`` and the ``<``/``>`` partial markers.
    """
    loc = loc.strip().replace(" ", "")
    strand = 1
    if loc.lower().startswith("complement("):
        strand = -1
        loc = loc[len("complement("):-1]
    if loc.lower().startswith(("join(", "order(")):
        inner = loc[loc.index("(") + 1:-1]
        segments: list[tuple[int, int]] = []
        depth = 0
        buf = ""
        for ch in inner:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            if ch == "," and depth == 0:
                segments.extend(parse_location(buf)[0])
                buf = ""
            else:
                buf += ch
        if buf:
            segments.extend(parse_location(buf)[0])
        return segments, strand

    nums = _LOC_NUM_RE.findall(loc)
    if not nums:
        raise ParseError(f"unparsable location: {loc!r}")
    if len(nums) == 1:
        start = int(nums[0]) - 1
        return [(start, start + 1)], strand
    return [(int(nums[0]) - 1, int(nums[1]))], strand


def parse_genbank(raw: bytes, name_hint: str = "") -> Record:
    text = _normalize_text(raw)
    lines = text.split("\n")
    warnings: list[str] = []

    name = name_hint
    circular = False
    description = ""

    # --- header -----------------------------------------------------------
    for line in lines:
        m = _LOCUS_RE.match(line)
        if m:
            name = m.group("name")
            rest = m.group("rest").lower()
            circular = "circular" in rest
            break
    else:
        warnings.append("no LOCUS line; header fields unavailable")

    for line in lines:
        if line.startswith("DEFINITION"):
            description = line[len("DEFINITION"):].strip().strip(".")
            break

    # --- sequence ---------------------------------------------------------
    seq_parts: list[str] = []
    in_origin = False
    for line in lines:
        upper = line.upper()
        if upper.startswith("ORIGIN"):
            in_origin = True
            continue
        if in_origin:
            if line.startswith("//"):
                break
            seq_parts.append(line)
    seq = clean_seq("".join(seq_parts))

    # --- features ---------------------------------------------------------
    features: list[Feature] = []
    in_features = False
    cur_key = None
    cur_loc = ""
    cur_quals: list[str] = []

    def flush():
        nonlocal cur_key, cur_loc, cur_quals
        if not cur_key:
            return
        try:
            segments, strand = parse_location(cur_loc)
        except ParseError as exc:
            warnings.append(str(exc))
            cur_key, cur_loc, cur_quals = None, "", []
            return
        quals: dict[str, str] = {}
        for q in cur_quals:
            if "=" in q:
                k, v = q.split("=", 1)
                quals[k.strip().lstrip("/")] = v.strip().strip('"')
            elif q.strip():
                quals[q.strip().lstrip("/")] = "true"
        label = (quals.get("label") or quals.get("gene") or quals.get("product")
                 or quals.get("note") or quals.get("standard_name") or "")
        start = min(s for s, _ in segments)
        end = max(e for _, e in segments)
        features.append(Feature(
            type=cur_key, start=start, end=end, strand=strand,
            label=label.strip(), qualifiers=quals,
            segments=segments if len(segments) > 1 else [],
        ))
        cur_key, cur_loc, cur_quals = None, "", []

    for line in lines:
        upper = line.upper()
        if upper.startswith("FEATURES"):
            in_features = True
            continue
        if not in_features:
            continue
        if upper.startswith("ORIGIN") or line.startswith("//"):
            flush()
            break
        if not line.startswith(" "):
            flush()
            in_features = False
            continue
        stripped = line.strip()
        if not stripped:
            continue
        # A feature key starts in column 6 (0-based 5) and is followed by a
        # location; qualifier continuation lines start with '/' or wrap text.
        if len(line) > 5 and line[5] != " " and not stripped.startswith("/"):
            flush()
            parts = stripped.split(None, 1)
            cur_key = parts[0]
            cur_loc = parts[1] if len(parts) > 1 else ""
        elif stripped.startswith("/"):
            cur_quals.append(stripped)
        elif cur_quals:
            # continuation of the previous qualifier value
            cur_quals[-1] = cur_quals[-1].rstrip('"') + stripped.strip('"')
        elif cur_key:
            cur_loc += stripped
    flush()

    # ApE stores some features in a COMMENT block rather than FEATURES.
    if "Features not parsed" in text:
        warnings.append("file contains an ApE 'Features not parsed' block")

    if not seq:
        warnings.append("no sequence found in ORIGIN block")

    return Record(name=name or name_hint, seq=seq, circular=circular,
                  features=features, description=description,
                  source_format="genbank", warnings=warnings)


# --------------------------------------------------------------------------
# FASTA / plain text
# --------------------------------------------------------------------------

def parse_fasta(raw: bytes, name_hint: str = "") -> Record:
    text = _normalize_text(raw)
    lines = [ln for ln in text.split("\n")]
    name = name_hint
    desc = ""
    seq_lines: list[str] = []
    seen_header = False
    for line in lines:
        if line.startswith(">"):
            if seen_header and seq_lines:
                break  # only the first record
            header = line[1:].strip()
            parts = header.split(None, 1)
            name = parts[0] if parts else name_hint
            desc = parts[1] if len(parts) > 1 else ""
            seen_header = True
        else:
            seq_lines.append(line)
    seq = clean_seq("".join(seq_lines))
    warnings = []
    if not seq:
        warnings.append("no sequence characters found")
    non_acgt = sum(1 for c in seq if c not in "ACGT")
    if seq and non_acgt / len(seq) > 0.05:
        warnings.append(f"{non_acgt} ambiguous bases ({non_acgt / len(seq):.0%})")
    return Record(name=name or name_hint, seq=seq, circular=False,
                  description=desc, source_format="fasta", warnings=warnings)


# --------------------------------------------------------------------------
# SnapGene .dna
# --------------------------------------------------------------------------

def _iter_snapgene_blocks(data: bytes):
    i = 0
    while i + 5 <= len(data):
        btype = data[i]
        (blen,) = struct.unpack(">I", data[i + 1:i + 5])
        i += 5
        if i + blen > len(data):
            break
        yield btype, data[i:i + blen]
        i += blen


def parse_snapgene(raw: bytes, name_hint: str = "") -> Record:
    seq = ""
    circular = False
    features: list[Feature] = []
    notes = ""
    description = ""
    warnings: list[str] = []

    for btype, payload in _iter_snapgene_blocks(raw):
        if btype == 0:  # DNA sequence; first byte is a flags byte
            if not payload:
                continue
            circular = bool(payload[0] & 0x01)
            seq = clean_seq(payload[1:].decode("ascii", errors="replace"))
        elif btype == 10:  # features, as XML
            try:
                root = ET.fromstring(payload.decode("utf-8", errors="replace"))
            except ET.ParseError as exc:
                warnings.append(f"feature XML unreadable: {exc}")
                continue
            for fe in root.findall("Feature"):
                ftype = fe.get("type", "misc_feature")
                label = fe.get("name", "")
                directionality = fe.get("directionality", "1")
                strand = -1 if directionality == "2" else 1
                spans = []
                for seg in fe.findall("Segment"):
                    rng = seg.get("range", "")
                    if "-" in rng:
                        a, b = rng.split("-", 1)
                        try:
                            spans.append((int(a) - 1, int(b)))
                        except ValueError:
                            continue
                if not spans:
                    continue
                quals = {"label": label}
                for q in fe.findall("Q"):
                    qn = q.get("name", "")
                    v = q.find("V")
                    if v is not None:
                        quals[qn] = v.get("text") or v.get("predef") or v.get("int") or ""
                features.append(Feature(
                    type=ftype, start=min(s for s, _ in spans),
                    end=max(e for _, e in spans), strand=strand,
                    label=label, qualifiers=quals,
                    segments=spans if len(spans) > 1 else [],
                ))
        elif btype == 6:  # notes XML
            try:
                root = ET.fromstring(payload.decode("utf-8", errors="replace"))
                bits = []
                for tag in ("Description", "Comments", "ReferencesTitle"):
                    el = root.find(tag)
                    if el is not None and el.text:
                        bits.append(el.text.strip())
                notes = "\n".join(bits)
                if bits:
                    description = bits[0][:200]
            except ET.ParseError:
                pass

    if not seq:
        raise ParseError("SnapGene file contained no DNA block")

    return Record(name=name_hint, seq=seq, circular=circular, features=features,
                  description=description, source_format="snapgene",
                  warnings=warnings, notes=notes)


# --------------------------------------------------------------------------
# Dispatch + role classification
# --------------------------------------------------------------------------

#: LOCUS names that editors write when the author never named the record.
#: They carry no information, so the filename is a better label.
PLACEHOLDER_NAMES = {
    "new_dna", "newdna", "untitled", "unnamed", "sequence", "seq",
    "exported", "export", "consensus", ".", "-", "null", "none",
}


def read_record(path: Path) -> Record:
    raw = Path(path).read_bytes()
    fmt = sniff_format(Path(path), raw[:512])
    stem = Path(path).stem
    if fmt == "genbank":
        rec = parse_genbank(raw, stem)
    elif fmt == "fasta":
        rec = parse_fasta(raw, stem)
    elif fmt == "snapgene":
        rec = parse_snapgene(raw, stem)
    else:
        raise ParseError(f"unrecognised format for {path}")

    # ApE and SnapGene default to "New_DNA"/"Untitled" when the author never
    # set a name, and GenBank truncates LOCUS names to 16 characters. The
    # filename is almost always the more informative label in this share.
    normalised = re.sub(r"[^a-z0-9]", "", (rec.name or "").lower())
    if not rec.name or normalised in PLACEHOLDER_NAMES:
        rec.warnings.append(
            f"record had no meaningful name ({rec.name!r}); using the filename")
        rec.name = stem
    elif len(rec.name) >= 16 and stem.lower().startswith(rec.name.lower()[:12]):
        # LOCUS was truncated; the filename holds the full name.
        rec.name = stem
    return rec


def classify_role(path: Path, rec: Record) -> tuple[str, str]:
    """Decide whether a record is a plasmid, a sequencing read, or a fragment.

    The lab share mixes finished plasmid maps with thousands of raw Sanger
    reads. Treating a read as a cloning template would be a real hazard, so
    reads are separated out explicitly.
    """
    seq = rec.seq
    n = len(seq)
    if n == 0:
        return "empty", "no sequence"

    ambiguous = sum(1 for c in seq if c not in "ACGT") / n

    # Sanger reads: ~300-1200 bp, N-rich ends, often named after an .ab1 well.
    looks_like_ab1 = bool(re.search(r"\.ab1\b", rec.description or "", re.I)) or \
        bool(re.search(r"\.ab1\b", rec.name or "", re.I))
    if looks_like_ab1:
        return "read", "header references an .ab1 trace"
    if rec.source_format == "fasta" and n < 1500 and ambiguous > 0.01:
        return "read", f"short unannotated sequence with {ambiguous:.1%} ambiguous bases"
    if path.suffix.lower() == ".seq" and n < 1500 and not rec.features:
        return "read", "short .seq file with no annotations"

    if rec.circular or n >= 1500:
        return "plasmid", ""
    return "fragment", "short linear sequence"
