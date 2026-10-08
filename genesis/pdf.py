"""A small PDF writer, so protocol export needs no third-party package.

Only the base-14 PDF fonts are used (Helvetica, Helvetica-Bold, Courier),
which every reader provides, so nothing has to be embedded. The widths
tables below are the standard AFM metrics; they are what makes line
wrapping and column alignment come out right rather than ragged.

The public surface is a simple flowing-layout API - ``heading``,
``paragraph``, ``bullet``, ``table``, ``mono`` - that paginates itself.
"""

from __future__ import annotations

import zlib
from dataclasses import dataclass, field

A4 = (595.28, 841.89)
LETTER = (612.0, 792.0)

# Standard AFM widths (per 1000 units) for ASCII 32..126.
_HELV = [
    278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278,
    278, 556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584,
    584, 556, 1015, 667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556,
    833, 722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 278,
    278, 278, 469, 556, 333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222,
    500, 222, 833, 556, 556, 556, 556, 333, 500, 278, 556, 500, 722, 500, 500,
    500, 334, 260, 334, 584,
]
_HELV_BOLD = [
    278, 333, 474, 556, 556, 889, 722, 238, 333, 333, 389, 584, 278, 333, 278,
    278, 556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 333, 333, 584, 584,
    584, 611, 975, 722, 722, 722, 722, 667, 611, 778, 722, 278, 556, 722, 611,
    833, 722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 333,
    278, 333, 584, 556, 333, 556, 611, 556, 611, 556, 333, 611, 611, 278, 278,
    556, 278, 889, 611, 611, 611, 611, 389, 556, 333, 611, 556, 778, 556, 556,
    500, 389, 280, 389, 584,
]

FONTS = {
    "F1": ("Helvetica", _HELV),
    "F2": ("Helvetica-Bold", _HELV_BOLD),
    "F3": ("Courier", None),          # monospace: every glyph is 600
    "F4": ("Courier-Bold", None),
}


def text_width(s: str, font: str, size: float) -> float:
    """Width of `s` in points, using the real font metrics."""
    widths = FONTS[font][1]
    if widths is None:
        return len(s) * 600 / 1000.0 * size
    total = 0
    for ch in s:
        code = ord(ch)
        total += widths[code - 32] if 32 <= code <= 126 else 556
    return total / 1000.0 * size


def wrap(s: str, font: str, size: float, max_width: float) -> list[str]:
    """Greedy word wrap that measures with the actual glyph widths."""
    if not s:
        return [""]
    lines: list[str] = []
    for paragraph in s.split("\n"):
        words = paragraph.split(" ")
        current = ""
        for word in words:
            trial = word if not current else current + " " + word
            if text_width(trial, font, size) <= max_width or not current:
                current = trial
                # a single word longer than the line must be hard-split
                while text_width(current, font, size) > max_width and len(current) > 1:
                    cut = len(current)
                    while cut > 1 and text_width(current[:cut], font, size) > max_width:
                        cut -= 1
                    lines.append(current[:cut])
                    current = current[cut:]
            else:
                lines.append(current)
                current = word
        lines.append(current)
    return lines


def _escape(s: str) -> bytes:
    """Encode to WinAnsi and escape the three PDF string metacharacters."""
    out = []
    for ch in s:
        code = ord(ch)
        if ch in "()\\":
            out.append("\\" + ch)
        elif 32 <= code <= 126:
            out.append(ch)
        elif code in (0x2019, 0x2018):
            out.append("'")
        elif code in (0x201C, 0x201D):
            out.append('"')
        elif code in (0x2013, 0x2014):
            out.append("-")
        elif code == 0x00B0:
            out.append("\\260")      # degree sign in WinAnsi
        elif code == 0x00B5:
            out.append("\\265")      # micro sign
        elif code == 0x2192:
            out.append("->")
        elif 160 <= code <= 255:
            out.append(f"\\{code:03o}")
        else:
            out.append("?")
    return "".join(out).encode("latin-1", "replace")


@dataclass
class _Page:
    ops: list = field(default_factory=list)


class PDF:
    """A paginating document builder."""

    def __init__(self, title: str = "", page_size=A4,
                 margin_left: float = 56, margin_right: float = 48,
                 margin_top: float = 56, margin_bottom: float = 56,
                 footer: str = ""):
        self.title = title
        self.width, self.height = page_size
        self.ml, self.mr = margin_left, margin_right
        self.mt, self.mb = margin_top, margin_bottom
        self.footer = footer
        self.pages: list[_Page] = []
        self._new_page()

    # -- geometry ------------------------------------------------------
    @property
    def content_width(self) -> float:
        return self.width - self.ml - self.mr

    def _new_page(self):
        self.pages.append(_Page())
        self.y = self.height - self.mt

    def _ensure(self, needed: float):
        if self.y - needed < self.mb:
            self._new_page()

    # -- primitives ----------------------------------------------------
    def _draw(self, s: str, x: float, y: float, font: str, size: float,
              gray: float | None = None):
        ops = self.pages[-1].ops
        if gray is not None:
            ops.append(f"q {gray:.2f} g")
        ops.append(f"BT /{font} {size:g} Tf 1 0 0 1 {x:.2f} {y:.2f} Tm")
        ops.append(b"(" + _escape(s) + b") Tj")
        ops.append("ET")
        if gray is not None:
            ops.append("Q")

    def rule(self, gray: float = 0.75, pad: float = 6):
        self._ensure(pad * 2)
        self.y -= pad
        self.pages[-1].ops.append(
            f"q {gray:.2f} G 0.6 w {self.ml:.2f} {self.y:.2f} m "
            f"{self.width - self.mr:.2f} {self.y:.2f} l S Q")
        self.y -= pad

    def fill_rect(self, x, y, w, h, gray=0.93, rgb=None):
        colour = (f"{rgb[0]:.3f} {rgb[1]:.3f} {rgb[2]:.3f} rg" if rgb
                  else f"{gray:.2f} g")
        self.pages[-1].ops.append(
            f"q {colour} {x:.2f} {y:.2f} {w:.2f} {h:.2f} re f Q")

    def stroke_rect(self, x, y, w, h, rgb=(0.2, 0.2, 0.2), width=0.6):
        self.pages[-1].ops.append(
            f"q {rgb[0]:.3f} {rgb[1]:.3f} {rgb[2]:.3f} RG {width} w "
            f"{x:.2f} {y:.2f} {w:.2f} {h:.2f} re S Q")

    def line(self, x1, y1, x2, y2, rgb=(0.3, 0.3, 0.3), width=0.6,
             dash: str = ""):
        d = f"[{dash}] 0 d " if dash else ""
        self.pages[-1].ops.append(
            f"q {rgb[0]:.3f} {rgb[1]:.3f} {rgb[2]:.3f} RG {width} w {d}"
            f"{x1:.2f} {y1:.2f} m {x2:.2f} {y2:.2f} l S Q")

    def tri(self, x, y, size=4, rgb=(0.8, 0.1, 0.1), up=True):
        """A small solid triangle, used as a cut or insert marker."""
        dy = size if up else -size
        self.pages[-1].ops.append(
            f"q {rgb[0]:.3f} {rgb[1]:.3f} {rgb[2]:.3f} rg "
            f"{x - size:.2f} {y:.2f} m {x + size:.2f} {y:.2f} l "
            f"{x:.2f} {y + dy:.2f} l f Q")

    def text_at(self, s: str, x: float, y: float, font: str = "F1",
                size: float = 8, rgb=None, centre: bool = False):
        if centre:
            x -= text_width(s, font, size) / 2
        ops = self.pages[-1].ops
        if rgb:
            ops.append(f"q {rgb[0]:.3f} {rgb[1]:.3f} {rgb[2]:.3f} rg")
        ops.append(f"BT /{font} {size:g} Tf 1 0 0 1 {x:.2f} {y:.2f} Tm")
        ops.append(b"(" + _escape(s) + b") Tj")
        ops.append("ET")
        if rgb:
            ops.append("Q")

    def spacer(self, h: float = 8):
        self._ensure(h)
        self.y -= h

    def page_break(self):
        self._new_page()

    # -- flowing content ----------------------------------------------
    def heading(self, text: str, level: int = 1):
        size = {1: 17, 2: 12.5, 3: 10.5}.get(level, 11)
        lead = size * 1.25
        space_before = {1: 0, 2: 14, 3: 10}.get(level, 10)
        self._ensure(lead + space_before + 6)
        self.y -= space_before
        for line in wrap(text, "F2", size, self.content_width):
            self._ensure(lead)
            self.y -= lead
            self._draw(line, self.ml, self.y, "F2", size)
        if level <= 2:
            self.rule(0.8, 4)
        else:
            self.y -= 3

    def paragraph(self, text: str, size: float = 9.5, gray: float | None = None,
                  font: str = "F1", indent: float = 0, lead_factor: float = 1.32):
        lead = size * lead_factor
        for line in wrap(text, font, size, self.content_width - indent):
            self._ensure(lead)
            self.y -= lead
            self._draw(line, self.ml + indent, self.y, font, size, gray)
        self.y -= 2

    def bullet(self, text: str, size: float = 9.5, marker: str = "•",
               indent: float = 12, gray: float | None = None):
        lead = size * 1.32
        body = self.content_width - indent - 10
        lines = wrap(text, "F1", size, body)
        for i, line in enumerate(lines):
            self._ensure(lead)
            self.y -= lead
            if i == 0:
                self._draw("-" if marker == "•" else marker,
                           self.ml + indent, self.y, "F1", size, gray)
            self._draw(line, self.ml + indent + 10, self.y, "F1", size, gray)
        self.y -= 1

    def mono(self, text: str, size: float = 8.2, indent: float = 12,
             gray: float | None = None):
        lead = size * 1.3
        for line in wrap(text, "F3", size, self.content_width - indent):
            self._ensure(lead)
            self.y -= lead
            self._draw(line, self.ml + indent, self.y, "F3", size, gray)

    def step(self, number: str, title: str, body: str = "",
             badge: str = ""):
        """A numbered protocol step with an optional operation-code badge."""
        size = 10
        lead = size * 1.3
        self._ensure(lead * 2 + 8)
        self.y -= 12

        num_w = 22
        x = self.ml + num_w
        badge_w = (text_width(badge, "F4", 7.6) + 8) if badge else 0
        text_x = x + (badge_w + 6 if badge else 0)

        lines = wrap(title, "F2", size, self.width - self.mr - text_x)

        # First row carries the number, the badge and the first title line.
        self.y -= lead
        self._draw(number, self.ml, self.y + 2, "F2", size)
        if badge:
            self.fill_rect(x - 2, self.y - 1, badge_w, lead, 0.90)
            self._draw(badge, x + 2, self.y + 2.5, "F4", 7.6, 0.15)
        self._draw(lines[0], text_x, self.y + 2, "F2", size)

        # Continuation lines align under the title, clear of the number.
        for line in lines[1:]:
            self._ensure(lead)
            self.y -= lead
            self._draw(line, x, self.y + 2, "F2", size)

        if body:
            self.y -= 3
            self.paragraph(body, 9, indent=num_w)
        self.y -= 2

    def table(self, headers: list[str], rows: list[list[str]],
              widths: list[float] | None = None, size: float = 8.4,
              mono_cols: set[int] | None = None, zebra: bool = True):
        """A simple table. `widths` are fractions of the content width."""
        mono_cols = mono_cols or set()
        n = len(headers)
        if widths is None:
            widths = [1.0 / n] * n
        cols = [w * self.content_width for w in widths]
        lead = size * 1.3
        pad = 3.0

        def render_row(cells, font_default, gray=None, bg=None):
            wrapped = []
            for i, cell in enumerate(cells):
                f = "F3" if i in mono_cols else font_default
                wrapped.append(wrap(str(cell), f, size, cols[i] - 2 * pad))
            height = max(len(w) for w in wrapped) * lead + 2
            self._ensure(height)
            if bg is not None:
                self.fill_rect(self.ml, self.y - height, self.content_width,
                               height, bg)
            top = self.y
            for i, lines in enumerate(wrapped):
                x = self.ml + sum(cols[:i]) + pad
                yy = top
                f = "F3" if i in mono_cols else font_default
                for line in lines:
                    yy -= lead
                    self._draw(line, x, yy + 2, f, size, gray)
            self.y -= height

        self._ensure(lead * 3)
        render_row(headers, "F2", None, 0.88)
        for idx, row in enumerate(rows):
            render_row(row, "F1", None, 0.965 if (zebra and idx % 2) else None)
        self.y -= 4

    def key_values(self, pairs: list[tuple[str, str]], size: float = 9.5,
                   label_width: float = 110):
        lead = size * 1.35
        for label, value in pairs:
            lines = wrap(str(value), "F1", size,
                         self.content_width - label_width)
            self._ensure(lead * len(lines))
            first = True
            for line in lines:
                self.y -= lead
                if first:
                    self._draw(label, self.ml, self.y, "F2", size)
                    first = False
                self._draw(line, self.ml + label_width, self.y, "F1", size)
        self.y -= 2

    # -- serialisation --------------------------------------------------
    def _page_content(self, page: _Page, page_no: int, total: int) -> bytes:
        parts: list[bytes] = []
        for op in page.ops:
            parts.append(op if isinstance(op, bytes) else op.encode("latin-1"))
        # footer
        label = self.footer or self.title
        foot = f"{label}".strip()
        if foot:
            parts.append(b"BT /F1 7.5 Tf 0.45 g 1 0 0 1 "
                         + f"{self.ml:.2f} {self.mb - 22:.2f}".encode()
                         + b" Tm (" + _escape(foot) + b") Tj ET")
        pg = f"{page_no} / {total}"
        w = text_width(pg, "F1", 7.5)
        parts.append(b"BT /F1 7.5 Tf 0.45 g 1 0 0 1 "
                     + f"{self.width - self.mr - w:.2f} {self.mb - 22:.2f}".encode()
                     + b" Tm (" + _escape(pg) + b") Tj ET")
        return b"\n".join(parts)

    def to_bytes(self, compress: bool = True) -> bytes:
        objects: list[bytes] = []

        def add(body: bytes) -> int:
            objects.append(body)
            return len(objects)          # 1-based object number

        total = len(self.pages)
        font_ids = {}
        for key, (base, _) in FONTS.items():
            font_ids[key] = add(
                f"<< /Type /Font /Subtype /Type1 /BaseFont /{base} "
                f"/Encoding /WinAnsiEncoding >>".encode())

        resources = ("<< /Font << " +
                     " ".join(f"/{k} {v} 0 R" for k, v in font_ids.items()) +
                     " >> >>")

        pages_obj = add(b"")            # placeholder, patched below
        page_ids: list[int] = []
        for i, page in enumerate(self.pages, start=1):
            raw = self._page_content(page, i, total)
            if compress:
                stream = zlib.compress(raw)
                meta = (f"<< /Length {len(stream)} /Filter /FlateDecode >>")
            else:
                stream = raw
                meta = f"<< /Length {len(stream)} >>"
            content_id = add(meta.encode() + b"\nstream\n" + stream + b"\nendstream")
            page_ids.append(add(
                f"<< /Type /Page /Parent {pages_obj} 0 R "
                f"/MediaBox [0 0 {self.width:.2f} {self.height:.2f}] "
                f"/Resources {resources} /Contents {content_id} 0 R >>".encode()))

        kids = " ".join(f"{pid} 0 R" for pid in page_ids)
        objects[pages_obj - 1] = (
            f"<< /Type /Pages /Kids [{kids}] /Count {len(page_ids)} >>".encode())

        info_id = add(f"<< /Title ({_escape(self.title).decode('latin-1')}) "
                      f"/Producer (Genesis) >>".encode("latin-1"))
        catalog_id = add(f"<< /Type /Catalog /Pages {pages_obj} 0 R >>".encode())

        out = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
        offsets = [0]
        for num, body in enumerate(objects, start=1):
            offsets.append(len(out))
            out += f"{num} 0 obj\n".encode() + body + b"\nendobj\n"

        xref_at = len(out)
        out += f"xref\n0 {len(objects) + 1}\n".encode()
        out += b"0000000000 65535 f \n"
        for off in offsets[1:]:
            out += f"{off:010d} 00000 n \n".encode()
        out += (f"trailer\n<< /Size {len(objects) + 1} /Root {catalog_id} 0 R "
                f"/Info {info_id} 0 R >>\nstartxref\n{xref_at}\n%%EOF\n").encode()
        return bytes(out)

    def save(self, path) -> None:
        from pathlib import Path
        Path(path).write_bytes(self.to_bytes())


# --------------------------------------------------------------------------
# Plasmid scheme drawing
# --------------------------------------------------------------------------

#: Muted fills per feature type, readable in both print and on screen.
FEATURE_COLOURS = {
    "promoter": (0.22, 0.55, 0.28), "enhancer": (0.35, 0.62, 0.40),
    "CDS": (0.20, 0.40, 0.70), "trafficking_motif": (0.48, 0.35, 0.68),
    "p2A": (0.85, 0.55, 0.15), "IRES": (0.80, 0.60, 0.25),
    "linker": (0.55, 0.55, 0.58), "tag": (0.30, 0.60, 0.65),
    "polyA": (0.40, 0.45, 0.70), "WPRE": (0.35, 0.50, 0.55),
    "lox_site": (0.78, 0.22, 0.28), "DIO": (0.70, 0.30, 0.45),
    "AAV_ITR": (0.45, 0.45, 0.45), "origin": (0.60, 0.60, 0.62),
    "resistance": (0.50, 0.52, 0.55),
}
CUT_RED = (0.80, 0.14, 0.14)
TAKE_GREEN = (0.13, 0.55, 0.25)

#: Colours for translation landmarks drawn above a track.
CODON_COLOURS = {
    "start": (0.13, 0.55, 0.25),           # ATG, green
    "stop": (0.25, 0.30, 0.60),            # a proper terminal stop, blue
    "premature_stop": (0.80, 0.14, 0.14),  # in-frame stop inside the fusion
    "missing_start": (0.80, 0.45, 0.05),   # should have been an ATG
    "missing_stop": (0.80, 0.45, 0.05),    # should have been a stop
}


def _feature_rgb(ftype: str):
    return FEATURE_COLOURS.get(ftype, (0.55, 0.57, 0.60))


def plasmid_track(doc: "PDF", title: str, subtitle: str, length: int,
                  features: list, highlight: dict | None = None,
                  highlight_kind: str = "cut", height: float = 15,
                  codon_marks: list | None = None):
    """Draw one parent plasmid as a scaled linear map.

    `highlight` marks either the region to cut out of a recipient (red,
    struck through) or the region to take from a donor (green). Labels are
    only drawn inside blocks wide enough to hold them; everything else is
    called out underneath. Everything that annotates the highlight is placed
    *below* the track so nothing can collide with the title above it.
    """
    codon_marks = codon_marks or []
    # Translation landmarks are drawn above the track, so they never collide
    # with the cut/take annotations below it. Two rows, staggered, because
    # a start and its stop can fall close together on a scaled map.
    mark_rows = 2 if len(codon_marks) > 1 else (1 if codon_marks else 0)
    mark_height = mark_rows * 9 + (4 if mark_rows else 0)

    doc._ensure(height + 74 + mark_height)
    width = doc.content_width
    x0 = doc.ml

    # --- title row -----------------------------------------------------
    doc.y -= 12
    doc.text_at(title, x0, doc.y, "F2", 9.5)
    if subtitle:
        doc.text_at(subtitle, x0 + text_width(title, "F2", 9.5) + 8, doc.y,
                    "F1", 8, rgb=(0.45, 0.45, 0.48))

    # --- translation landmarks, above the track -------------------------
    doc.y -= mark_height
    marks_base = doc.y
    scale0 = (width / length) if length else 0.0
    for idx, mark in enumerate(sorted(codon_marks, key=lambda m: m["start"])):
        mx = x0 + mark["start"] * scale0
        rgb = CODON_COLOURS.get(mark["kind"], (0.4, 0.4, 0.4))
        row = idx % max(mark_rows, 1)
        my = marks_base + row * 9
        label = mark.get("label", "")
        lw = text_width(label, "F2", 6)
        lx = min(max(mx, x0 + lw / 2), x0 + width - lw / 2)
        doc.text_at(label, lx, my + 2, "F2", 6, rgb=rgb, centre=True)
        # a hairline from the label down to the exact base it marks
        doc.line(mx, my, mx, marks_base - 4 - row * 0, rgb=rgb, width=0.5)

    # --- track ---------------------------------------------------------
    doc.y -= height + 10
    track_y = doc.y
    scale = (width / length) if length else 0.0

    doc.fill_rect(x0, track_y + height / 2 - 1.2, width, 2.4,
                  rgb=(0.80, 0.80, 0.82))

    overflow: list[str] = []
    for f in features:
        fx = x0 + f["start"] * scale
        fw = max((f["end"] - f["start"]) * scale, 1.2)
        doc.fill_rect(fx, track_y, fw, height, rgb=_feature_rgb(f["feature_type"]))
        if f.get("strand") == -1:
            # a small notch on the left shows the feature reads right to left
            doc.tri(fx + 3, track_y + height / 2 - 2.5, 2.5, (1, 1, 1), up=True)
        label = f["name"] or f["feature_type"]
        if text_width(label, "F2", 6.6) + 6 <= fw:
            doc.text_at(label, fx + fw / 2, track_y + height / 2 - 2.3, "F2",
                        6.6, rgb=(1, 1, 1), centre=True)
        else:
            overflow.append(label)

    for mark in codon_marks:
        mx = x0 + mark["start"] * scale
        rgb = CODON_COLOURS.get(mark["kind"], (0.4, 0.4, 0.4))
        doc.tri(mx, track_y + height + 1, 3.0, rgb, up=False)

    # --- highlight, annotated entirely below the track ------------------
    doc.y = track_y - 9
    if highlight and highlight.get("length"):
        hx = x0 + highlight["start"] * scale
        hw = max(highlight["length"] * scale, 1.5)
        rgb = CUT_RED if highlight_kind == "cut" else TAKE_GREEN

        doc.stroke_rect(hx - 1, track_y - 3, hw + 2, height + 6, rgb=rgb,
                        width=1.3)
        # markers point up at the two cut boundaries
        doc.tri(hx, track_y - 7.5, 3.4, rgb, up=True)
        doc.tri(hx + hw, track_y - 7.5, 3.4, rgb, up=True)
        if highlight_kind == "cut":
            doc.line(hx, track_y + height / 2, hx + hw, track_y + height / 2,
                     rgb=rgb, width=1.1)

        doc.text_at(f"{highlight['start']:,}", hx, track_y - 16, "F1", 6.4,
                    rgb=rgb, centre=True)
        doc.text_at(f"{highlight['end']:,}", hx + hw, track_y - 16, "F1", 6.4,
                    rgb=rgb, centre=True)
        caption = (f"cut out {highlight['length']:,} bp"
                   if highlight_kind == "cut"
                   else f"take {highlight['length']:,} bp")
        if highlight.get("label"):
            caption += f"  ({highlight['label']})"
        cx = min(max(hx + hw / 2, x0 + text_width(caption, "F2", 7) / 2),
                 x0 + width - text_width(caption, "F2", 7) / 2)
        doc.text_at(caption, cx, track_y - 25, "F2", 7, rgb=rgb, centre=True)
        doc.y = track_y - 32
    else:
        doc.y = track_y - 14

    # --- scale ---------------------------------------------------------
    doc.text_at("0", x0, doc.y, "F1", 6.4, rgb=(0.5, 0.5, 0.5))
    end_label = f"{length:,} bp"
    doc.text_at(end_label, x0 + width - text_width(end_label, "F1", 6.4),
                doc.y, "F1", 6.4, rgb=(0.5, 0.5, 0.5))
    doc.y -= 8
    if overflow:
        doc.paragraph("also: " + ", ".join(overflow[:10]), 7, gray=0.5)
    doc.y -= 4


# --------------------------------------------------------------------------
# Circular plasmid maps
# --------------------------------------------------------------------------

import math as _math


def _polar(cx, cy, radius, fraction):
    """Point on a circle. Zero is at the top and angles run clockwise, which
    is how plasmid maps are always drawn."""
    angle = -_math.pi / 2 + 2 * _math.pi * fraction
    return cx + radius * _math.cos(angle), cy + radius * _math.sin(angle)


def _arc_path(doc, cx, cy, radius, start_f, end_f, width, rgb):
    """Stroke an arc as a run of short segments.

    PDF has Bezier curves but no arc operator, and approximating an arc with
    Beziers for a shape this small buys nothing over flattening it.
    """
    span = (end_f - start_f) % 1.0 or 1.0
    steps = max(int(span * 180), 2)
    pts = [_polar(cx, cy, radius, start_f + span * i / steps)
           for i in range(steps + 1)]
    path = f"{pts[0][0]:.2f} {pts[0][1]:.2f} m " + " ".join(
        f"{x:.2f} {y:.2f} l" for x, y in pts[1:])
    doc.pages[-1].ops.append(
        f"q {rgb[0]:.3f} {rgb[1]:.3f} {rgb[2]:.3f} RG {width} w 1 J "
        f"{path} S Q")


def circular_map(doc: "PDF", title: str, subtitle: str, length: int,
                 features: list, radius: float = 108,
                 codon_marks: list | None = None,
                 highlight: dict | None = None,
                 highlight_kind: str = "cut"):
    """Draw a plasmid the way a plasmid is drawn: as a circle.

    A linear track is fine for comparing two parents side by side, but it
    misrepresents the molecule: it implies the ends are ends, when in a
    plasmid they are the same point. A construct whose cassette wraps the
    origin looks broken on a linear map and obvious on a circular one.
    """
    if not length:
        return
    needed = radius * 2 + 86
    doc._ensure(needed)

    doc.y -= 14
    doc.text_at(title, doc.ml, doc.y, "F2", 9.5)
    if subtitle:
        doc.text_at(subtitle, doc.ml + text_width(title, "F2", 9.5) + 8,
                    doc.y, "F1", 8, rgb=(0.45, 0.45, 0.48))

    cx = doc.ml + doc.content_width / 2
    cy = doc.y - radius - 16
    doc.y = cy - radius - 26

    # backbone circle
    _arc_path(doc, cx, cy, radius, 0.0, 1.0, 1.2, (0.78, 0.78, 0.80))

    # origin tick at twelve o'clock
    ox, oy = _polar(cx, cy, radius + 7, 0.0)
    doc.line(cx, cy + radius - 4, ox, oy, rgb=(0.55, 0.55, 0.58), width=0.7)
    doc.text_at("1", cx, cy + radius + 10, "F1", 6.4, rgb=(0.5, 0.5, 0.5),
                centre=True)

    labelled = []
    for f in features:
        start_f = (f["start"] % length) / length
        end_f = (f["end"] % length) / length or 1.0
        if f["end"] - f["start"] >= length:
            start_f, end_f = 0.0, 1.0
        inner = f.get("strand", 1) == -1
        r = radius - (9 if inner else 0)
        _arc_path(doc, cx, cy, r, start_f, end_f, 7.5,
                  _feature_rgb(f["feature_type"]))
        mid = start_f + (((end_f - start_f) % 1.0) or 1.0) / 2
        labelled.append((mid % 1.0, f, r))

    # labels outside the circle, nudged apart so they stay readable
    labelled.sort(key=lambda item: item[0])
    last_y = {}
    for mid, f, r in labelled:
        if f["end"] - f["start"] < length * 0.012:
            continue
        lx, ly = _polar(cx, cy, radius + 13, mid)
        right = lx >= cx
        name = (f["name"] or f["feature_type"])[:26]
        w = text_width(name, "F1", 6.6)
        side = round(ly / 9)
        while side in last_y:
            side -= 1
        last_y[side] = True
        ly = side * 9
        doc.line(*_polar(cx, cy, r + 4, mid), lx, ly, rgb=(0.6, 0.6, 0.62),
                 width=0.4)
        doc.text_at(name, lx + (3 if right else -3 - w), ly - 2, "F1", 6.6,
                    rgb=_feature_rgb(f["feature_type"]))

    # The region being cut out of a recipient, or taken from a donor, drawn
    # as a ring outside the features with its own boundary markers.
    if highlight and highlight.get("length"):
        rgb = CUT_RED if highlight_kind == "cut" else TAKE_GREEN
        hs = (highlight["start"] % length) / length
        he = (highlight["end"] % length) / length
        if he <= hs:
            he += 1.0
        _arc_path(doc, cx, cy, radius + 7, hs, he, 2.4, rgb)
        for edge in (hs, he % 1.0):
            a = _polar(cx, cy, radius + 2, edge)
            b = _polar(cx, cy, radius + 12, edge)
            doc.line(a[0], a[1], b[0], b[1], rgb=rgb, width=1.4)
        caption = (f"cut out {highlight['length']:,} bp"
                   if highlight_kind == "cut"
                   else f"take {highlight['length']:,} bp")
        if highlight.get("label"):
            caption += f" ({highlight['label'][:28]})"
        doc.text_at(caption, cx, cy - 24, "F2", 7.5, rgb=rgb, centre=True)

    for mark in codon_marks or []:
        mf = (mark["start"] % length) / length
        rgb = CODON_COLOURS.get(mark["kind"], (0.4, 0.4, 0.4))
        x1, y1 = _polar(cx, cy, radius - 14, mf)
        x2, y2 = _polar(cx, cy, radius + 3, mf)
        doc.line(x1, y1, x2, y2, rgb=rgb, width=1.0)

    doc.text_at(f"{length:,} bp", cx, cy + 2, "F2", 9, centre=True)
    doc.text_at("circular", cx, cy - 10, "F1", 7, rgb=(0.5, 0.5, 0.52),
                centre=True)
    doc.y -= 6
