"""A small PDF writer, so protocol export needs no third-party package.

Only the base-14 PDF fonts are used (Helvetica, Courier, Times),
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

_TIMES = [
    250, 333, 408, 500, 500, 833, 778, 180, 333, 333, 500, 564, 250, 333, 250,
    278, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 278, 278, 564, 564,
    564, 444, 921, 722, 667, 667, 722, 611, 556, 722, 722, 333, 389, 722, 611,
    889, 722, 722, 556, 722, 667, 556, 611, 722, 722, 944, 722, 722, 611, 333,
    278, 333, 469, 500, 333, 444, 500, 444, 500, 444, 333, 500, 500, 278, 278,
    500, 278, 778, 500, 500, 500, 500, 333, 389, 278, 500, 500, 722, 500, 500,
    444, 480, 200, 480, 541,
]
_TIMES_BOLD = [
    250, 333, 555, 500, 500, 1000, 833, 278, 333, 333, 500, 570, 250, 333, 250,
    278, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 333, 333, 570, 570,
    570, 500, 930, 722, 667, 722, 722, 667, 611, 778, 778, 389, 500, 778, 667,
    944, 722, 778, 611, 778, 722, 556, 667, 722, 722, 1000, 722, 722, 667, 333,
    278, 333, 581, 500, 333, 500, 556, 444, 556, 444, 333, 500, 556, 278, 333,
    556, 278, 833, 556, 500, 556, 556, 444, 389, 333, 556, 500, 722, 500, 500,
    444, 394, 220, 394, 520,
]
_TIMES_ITALIC = [
    250, 333, 420, 500, 500, 833, 778, 214, 333, 333, 500, 675, 250, 333, 250,
    278, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 333, 333, 675, 675,
    675, 500, 920, 611, 611, 667, 722, 611, 611, 722, 722, 333, 444, 667, 556,
    833, 667, 722, 611, 722, 611, 500, 556, 722, 611, 833, 611, 556, 556, 389,
    278, 389, 422, 500, 333, 500, 500, 444, 500, 444, 278, 500, 500, 278, 278,
    444, 278, 722, 500, 500, 500, 500, 389, 389, 278, 500, 444, 667, 444, 444,
    389, 400, 275, 400, 541,
]

FONTS = {
    "F1": ("Helvetica", _HELV),
    "F2": ("Helvetica-Bold", _HELV_BOLD),
    "F3": ("Courier", None),          # monospace: every glyph is 600
    "F4": ("Courier-Bold", None),
    # A vendor protocol is quoted material in a document otherwise set in
    # a UI face, and setting it in a serif marks it as such - as well as
    # fitting more of a reagent name on a line than Courier can.
    "F5": ("Times-Roman", _TIMES),
    "F6": ("Times-Bold", _TIMES_BOLD),
    "F7": ("Times-Italic", _TIMES_ITALIC),
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


#: Greek letters that turn up in gene and promoter names, spelled out.
#: The base-14 fonts are WinAnsi only, so the alternative is a question
#: mark in the middle of a name.
_GREEK = {
    "\u03b1": "alpha", "\u03b2": "beta", "\u03b3": "gamma",
    "\u03b4": "delta", "\u03b5": "epsilon", "\u03ba": "kappa",
    "\u03bb": "lambda", "\u03bc": "u", "\u03c3": "sigma",
    "\u03c6": "phi", "\u03c8": "psi", "\u03c9": "omega",
    "\u0391": "Alpha", "\u0392": "Beta", "\u0393": "Gamma",
    "\u0394": "Delta", "\u039b": "Lambda", "\u03a9": "Omega",
}


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
        elif ch in _GREEK:
            # "chicken beta-actin" is a name somebody has to read; WinAnsi
            # has no beta, and "chicken ?-actin" is not a name.
            out.append(_GREEK[ch])
        elif 160 <= code <= 255:
            out.append(f"\\{code:03o}")
        else:
            out.append("?")
    return "".join(out).encode("latin-1", "replace")


@dataclass
class _Page:
    ops: list = field(default_factory=list)
    #: (x0, y0, x1, y1, uri) rectangles that open a URL when clicked.
    links: list = field(default_factory=list)


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

    def fill_tri(self, p1, p2, p3, rgb=(0.2, 0.2, 0.2)):
        """A solid triangle at any orientation, from three points."""
        self.pages[-1].ops.append(
            f"q {rgb[0]:.3f} {rgb[1]:.3f} {rgb[2]:.3f} rg "
            f"{p1[0]:.2f} {p1[1]:.2f} m {p2[0]:.2f} {p2[1]:.2f} l "
            f"{p3[0]:.2f} {p3[1]:.2f} l f Q")

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

    def link(self, label: str, uri: str, size: float = 8,
             indent: float = 0, font: str = "F1"):
        """A line of text that opens a URL.

        A protocol that cites a vendor page should be able to go there:
        printing the address and leaving the reader to type it is the
        kind of friction that stops anybody checking the source.
        """
        lead = size * 1.32
        self._ensure(lead)
        self.y -= lead
        x = self.ml + indent
        self._draw(label, x, self.y, font, size, gray=None)
        w = text_width(label, font, size)
        self.pages[-1].ops.append(
            f"q 0.12 0.35 0.65 RG 0.5 w {x:.2f} {self.y - 1.5:.2f} m "
            f"{x + w:.2f} {self.y - 1.5:.2f} l S Q")
        self.pages[-1].links.append(
            (x, self.y - 3, x + w, self.y + size * 0.78, uri))
        self.y -= 2

    def boxed(self, draw, pad: float = 7, rgb=(0.80, 0.82, 0.85),
              fill=None):
        """Run `draw` and put a rule around whatever it printed.

        The border has to go underneath the text, and content operators
        are appended in order, so the rectangle is spliced in at the
        position the block started from rather than drawn after it.
        """
        self.y -= pad
        page = self.pages[-1]
        at = len(page.ops)
        top = self.y
        draw()
        if self.pages[-1] is not page:
            # the block paginated; a box around half of it would be a lie
            self.y -= pad
            return
        self.y -= pad
        h = top - self.y
        x0, x1 = self.ml - pad + 2, self.width - self.mr + pad - 2
        box = []
        if fill:
            box.append(f"q {fill[0]:.3f} {fill[1]:.3f} {fill[2]:.3f} rg "
                       f"{x0:.2f} {self.y:.2f} {x1 - x0:.2f} {h:.2f} re f Q")
        box.append(f"q {rgb[0]:.3f} {rgb[1]:.3f} {rgb[2]:.3f} RG 0.7 w "
                   f"{x0:.2f} {self.y:.2f} {x1 - x0:.2f} {h:.2f} re S Q")
        page.ops[at:at] = box
        self.y -= 4

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
                  font: str = "F1", indent: float = 0,
                  lead_factor: float = 1.32, rgb=None):
        lead = size * lead_factor
        for line in wrap(text, font, size, self.content_width - indent):
            self._ensure(lead)
            self.y -= lead
            self._mark(line, self.ml + indent, self.y, font, size, gray, rgb)
        self.y -= 2

    def bullet(self, text: str, size: float = 9.5, marker: str = "•",
               indent: float = 12, gray: float | None = None,
               font: str = "F1", rgb=None):
        lead = size * 1.32
        body = self.content_width - indent - 10
        lines = wrap(text, font, size, body)
        for i, line in enumerate(lines):
            self._ensure(lead)
            self.y -= lead
            if i == 0:
                self._mark("-" if marker == "•" else marker,
                           self.ml + indent, self.y, font, size, gray, rgb)
            self._mark(line, self.ml + indent + 10, self.y, font, size,
                       gray, rgb)
        self.y -= 1

    def _mark(self, s, x, y, font, size, gray, rgb):
        """Draw in a colour when one is given, else in grey."""
        if rgb is None:
            self._draw(s, x, y, font, size, gray)
        else:
            self.text_at(s, x, y, font, size, rgb=rgb)

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
            annots = ""
            if page.links:
                ids = []
                for x0, y0, x1, y1, uri in page.links:
                    ids.append(add(
                        f"<< /Type /Annot /Subtype /Link "
                        f"/Rect [{x0:.2f} {y0:.2f} {x1:.2f} {y1:.2f}] "
                        f"/Border [0 0 0] /A << /S /URI /URI ("
                        .encode("latin-1") + _escape(uri) + b") >> >>"))
                annots = ("/Annots [" +
                          " ".join(f"{i} 0 R" for i in ids) + "] ")
            page_ids.append(add(
                f"<< /Type /Page /Parent {pages_obj} 0 R "
                f"/MediaBox [0 0 {self.width:.2f} {self.height:.2f}] "
                f"/Resources {resources} {annots}"
                f"/Contents {content_id} 0 R >>".encode()))

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
    "Kozak": (0.20, 0.62, 0.52),
    "polyA": (0.40, 0.45, 0.70), "WPRE": (0.35, 0.50, 0.55),
    "lox_site": (0.78, 0.22, 0.28), "DIO": (0.70, 0.30, 0.45),
    "AAV_ITR": (0.45, 0.45, 0.45), "origin": (0.60, 0.60, 0.62),
    "resistance": (0.50, 0.52, 0.55),
}
CUT_RED = (0.80, 0.14, 0.14)
#: For the one or two lines in a protocol that cost a week when missed.
WARN_RED = (0.72, 0.10, 0.10)
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
        # The ends of the protein, named. On a linear track the left edge
        # is the amino end only when the feature reads forwards; a reverse
        # strand coding sequence starts at its right edge.
        ends_fit = fw >= 16 and f["feature_type"] in _PROTEIN_ENDS
        if ends_fit:
            reverse = f.get("strand", 1) == -1
            n_x = fx + fw - 4 if reverse else fx + 4
            c_x = fx + 4 if reverse else fx + fw - 4
            for text, tx in (("N", n_x), ("C", c_x)):
                doc.text_at(text, tx, track_y + height / 2 - 2.3, "F2", 5.6,
                            rgb=(1, 1, 1), centre=True)

        label = f["name"] or f["feature_type"]
        room = fw - (18 if ends_fit else 0)
        if text_width(label, "F2", 6.6) + 6 <= room:
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


#: Feature types whose two ends are the ends of a protein. A 2A peptide
#: or a trafficking motif is part of somebody else's protein, so giving it
#: its own N and C would say there are more proteins here than there are.
_PROTEIN_ENDS = {"CDS"}

#: Below this an annotation is too narrow to carry two labels without
#: them sitting on top of each other.
_ENDS_MIN_SPAN = 0.035


def _terminus_ends(f: dict) -> tuple:
    """(where the N-terminus is, where the C-terminus is), as fractions.

    Translation starts at the low coordinate of a forward feature and at
    the high coordinate of a reverse one, so which physical end is the
    amino terminus depends on the strand and not on the drawing.
    """
    start_f, end_f = f["_start_f"], f["_end_f"]
    if f.get("strand", 1) == -1:
        return end_f, start_f
    return start_f, end_f


def _arrow_head(doc, cx, cy, radius, fraction, strand, rgb):
    """A chevron on a feature arc, pointing the way it is translated.

    A coding sequence on a plasmid map is a block, and a block says
    nothing about which way the ribosome reads it. On a circular map that
    matters more than on a linear one, because an inverted payload - a
    Cre-lox cassette, say - looks exactly like an upright one.
    """
    step = 0.0085 * (1 if strand == 1 else -1)
    tip = _polar(cx, cy, radius, fraction + step)
    out = _polar(cx, cy, radius + 5.2, fraction)
    inn = _polar(cx, cy, radius - 5.2, fraction)
    doc.fill_tri(tip, out, inn, rgb)


#: What each codon mark is called on the map.
CODON_GLYPH = {
    "start": "ATG",
    "stop": "stop",
    "premature_stop": "early stop",
    "missing_start": "no ATG",
    "missing_stop": "no stop",
}

#: Feature types that are read by a ribosome, and so get a direction.
_TRANSLATED = {"CDS", "p2A", "tag", "trafficking_motif", "linker"}


def circular_map(doc: "PDF", title: str, subtitle: str, length: int,
                 features: list, radius: float = 108,
                 codon_marks: list | None = None,
                 highlight: dict | None = None,
                 highlight_kind: str = "cut",
                 primer_marks: list | None = None):
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
        rgb_f = _feature_rgb(f["feature_type"])
        _arc_path(doc, cx, cy, r, start_f, end_f, 7.5, rgb_f)
        span_f = ((end_f - start_f) % 1.0) or 1.0
        # Which way it is read, drawn on the arc itself.
        if f["feature_type"] in _TRANSLATED and span_f > 0.015:
            strand = f.get("strand", 1)
            at = (end_f if strand == 1 else start_f) % 1.0
            _arrow_head(doc, cx, cy, r, at, strand, rgb_f)
        # Which end of the block is the start of the protein and which is
        # the end of it. On a circle the two look alike, and on a reverse
        # strand feature the answer is the opposite of the way it is drawn.
        if f["feature_type"] in _PROTEIN_ENDS and span_f > _ENDS_MIN_SPAN:
            f["_start_f"], f["_end_f"] = start_f, end_f
            n_at, c_at = _terminus_ends(f)
            for label, where in (("N", n_at), ("C", c_at)):
                lx, ly = _polar(cx, cy, r - 8.5, where % 1.0)
                doc.text_at(label, lx, ly - 1.8, "F2", 5.4, rgb=rgb_f,
                            centre=True)
        mid = start_f + span_f / 2
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

    shown = []
    placed: list = []
    for mark in codon_marks or []:
        mf = (mark["start"] % length) / length
        rgb = CODON_COLOURS.get(mark["kind"], (0.4, 0.4, 0.4))
        x1, y1 = _polar(cx, cy, radius - 14, mf)
        x2, y2 = _polar(cx, cy, radius + 3, mf)
        doc.line(x1, y1, x2, y2, rgb=rgb, width=1.0)
        # Name it. An unlabelled tick on a ring is not information.
        glyph = CODON_GLYPH.get(mark["kind"], "")
        if glyph:
            # Two codon marks a few bases apart land on the same spot and
            # overprint each other, which is how "no stop" and "no ATG"
            # came out as one unreadable word. Step the label inwards
            # until it has room.
            depth = 20
            while any(abs(mf - prev) * length < 260 and abs(depth - d) < 9
                      for prev, d in placed):
                depth += 9
            placed.append((mf, depth))
            gx, gy = _polar(cx, cy, radius - depth, mf)
            doc.text_at(glyph, gx, gy - 2, "F2", 5.6, rgb=rgb, centre=True)
        shown.append(mark["kind"])

    # Where the primers sit, drawn as arrows reading the way they extend.
    for pm in primer_marks or []:
        pf = (pm["start"] % length) / length
        ef = (pm["end"] % length) / length
        rgb = (0.55, 0.20, 0.60)
        _arc_path(doc, cx, cy, radius + 13, min(pf, ef), max(pf, ef), 1.6, rgb)
        _arrow_head(doc, cx, cy, radius + 13, ef,
                    pm.get("strand", 1), rgb)
        tx, ty = _polar(cx, cy, radius + 19, (pf + ef) / 2)
        doc.text_at(str(pm.get("name", ""))[:18], tx, ty - 2, "F1", 5.6,
                    rgb=rgb, centre=True)

    # Size only. A plasmid name runs off both sides of a circle this size
    # and overprints the ring; the heading above already carries it.
    doc.text_at(f"{length:,} bp", cx, cy + 2, "F2", 9, centre=True)
    doc.text_at("circular", cx, cy - 10, "F1", 7, rgb=(0.5, 0.5, 0.52),
                centre=True)

    key = []
    if any(f["feature_type"] in _TRANSLATED for f in features):
        key.append("arrowheads point the way each coding sequence is read")
    if any(f["feature_type"] in _PROTEIN_ENDS for f in features):
        key.append("N and C mark the amino and carboxy ends of each protein")
    for kind in ("start", "stop", "premature_stop",
                 "missing_start", "missing_stop"):
        if kind in shown:
            key.append(CODON_GLYPH[kind] + {
                "start": " = start codon", "stop": " = stop codon",
                "premature_stop": " = in-frame stop inside the fusion",
                "missing_start": " = no start codon where one belongs",
                "missing_stop": " = no stop codon where one belongs",
            }[kind])
    if primer_marks:
        key.append("purple arcs are primers, pointing 5' to 3'")
    if key:
        doc.y -= 2
        doc.text_at("; ".join(key), doc.ml, doc.y, "F1", 6.2,
                    rgb=(0.45, 0.45, 0.48))
        doc.y -= 8
    doc.y -= 6


# --------------------------------------------------------------------------
# The cloning figure
# --------------------------------------------------------------------------
# Until now the PDF said where each fragment was cut from by drawing a
# full plasmid track per source and then tabulating where each landed.
# That is the data, not the picture: it never showed the reaction, so
# "two PCRs and a HiFi" and "a digest, then a ligation" looked identical
# on the page. This is the figure a cloning paper draws - one lane per
# piece, one column per bench stage - and it is the same drawing the web
# app shows, from the same lane model in genesis.figure.

#: Geometry, in points. The whole figure has to fit the A4 content width
#: (491 pt) at its widest, which is two stage columns.
_FIG = {
    "r": 30,            # parent circle
    "pr": 34,           # product circle
    "row": 110,
    "parent": 100,      # the parent column, its labels included
    "stage": 54,
    "frag": 86,
    "bracket": 10,
    "label": 50,        # how wide a name beside a circle may be
}

_FIG_GREY = (0.52, 0.54, 0.57)
_FIG_INK = (0.16, 0.18, 0.21)
_FIG_PRIMER = (0.13, 0.55, 0.25)


def _fig_label(doc, placed, cx, cy, r, fraction, text, rgb, size=5.6):
    """A name beside a circle, nudged clear of whatever is already there.

    Two cut sites a few hundred bases apart sit at nearly the same angle
    and their names overprint into a smudge, and a label at three
    o'clock lands on the stage arrow's caption. Clashing labels are
    pushed away from the circle - never across it - and a nudged label
    keeps a leader line back to the arc it names, so moving it does not
    make it ambiguous.
    """
    if not text:
        return
    at = fraction % 1.0
    anchor = _polar(cx, cy, r + 2, at)
    # At twelve and six o'clock a label set beside the point lies along
    # the tangent and so runs straight across the ring. There it is
    # centred on the pole and pushed clear instead.
    pole = min(at, 1 - at) < 0.06 or abs(at - 0.5) < 0.06
    x, y = _polar(cx, cy, r + (13 if pole else 6), at)
    right = at < 0.5
    while len(text) > 3 and text_width(text, "F1", size) > _FIG["label"]:
        text = text[:-1]
    w = text_width(text, "F1", size)
    # PDF y grows upwards, so "away from the circle" in the lower half
    # means a smaller y. The JavaScript port has the opposite sign.
    step = -7 if 0.25 <= at <= 0.75 else 7
    y0 = y
    for _ in range(9):
        if not any((p[1] == right or pole) and abs(p[0] - y) < 6.5
                   for p in placed):
            break
        y += step
    placed.append((y, right))
    if pole:
        placed.append((y, not right))
    tx = (x - w / 2) if pole else x + (2 if right else -2 - w)
    if abs(y - y0) > 3:
        doc.line(anchor[0], anchor[1], tx + (0 if right else w), y - 2,
                 rgb=(0.72, 0.73, 0.75), width=0.3)
    doc.text_at(text, tx, y - 2, "F1", size, rgb=rgb)


def _fig_primer(doc, cx, cy, r, fraction, forward, placed):
    """A primer, drawn the way a cloning figure draws one: a short arc on
    the backbone with a head on the end it extends towards."""
    w = 0.03 * (1 if forward else -1)
    _arc_path(doc, cx, cy, r + 4.5, min(fraction, fraction + w) % 1.0,
              max(fraction, fraction + w) % 1.0, 1.3, _FIG_PRIMER)
    tip = _polar(cx, cy, r + 4.5, fraction + w)
    doc.fill_tri(tip, _polar(cx, cy, r + 2, fraction + w * 0.55),
                 _polar(cx, cy, r + 7, fraction + w * 0.55), _FIG_PRIMER)
    _fig_label(doc, placed, cx, cy, r + 3, fraction + w * 0.5,
               "F" if forward else "R", _FIG_PRIMER, 5.4)


def _fig_cut(doc, cx, cy, r, fraction, name, placed):
    """An enzyme cut: a tick through the backbone, with the enzyme named."""
    a = _polar(cx, cy, r - 5, fraction)
    b = _polar(cx, cy, r + 5, fraction)
    doc.line(a[0], a[1], b[0], b[1], rgb=CUT_RED, width=1.3)
    _fig_label(doc, placed, cx, cy, r + 2, fraction, name, CUT_RED)


def _fig_caption(doc, cx, cy, r, name, size_text):
    """The plasmid's name under its circle, its size inside it.

    Inventory names run to thirty characters; inside a 60 pt circle they
    either overflow the ring or have to be cut so short they stop
    identifying the plasmid. Under it there is a whole column of room.
    """
    doc.text_at(size_text, cx, cy - 2, "F1", 6.2, rgb=_FIG_GREY, centre=True)
    text = str(name or "template")
    # Shrink to fit the column before breaking the name in half; splitting
    # "paav-ef1a-dio-asap7y-wpre_" after 24 characters left "e_" on a line
    # of its own.
    size = 6.6
    while size > 5.2 and text_width(text, "F2", size) > _FIG["parent"] - 4:
        size -= 0.3
    while len(text) > 6 and text_width(text, "F2", size) > _FIG["parent"] - 4:
        text = (text[:-4] if text.endswith("...") else text[:-1]) + "..."
    doc.text_at(text, cx, cy - r - 16, "F2", size, rgb=_FIG_INK, centre=True)


def _within(s, e, lo, hi) -> bool:
    """Is the arc s..e inside the arc lo..hi, on a circle?"""
    def rel(x):
        return (x - lo) % 1.0
    span = (hi - lo) % 1.0 or 1.0
    return rel(s) <= span and rel(e) <= span and rel(s) <= rel(e)


def _fig_parent(doc, lane, cx, cy):
    """One parent plasmid: what it is, what comes out, and how."""
    src = lane.get("src") or {}
    span = src.get("template_length") or 0
    r = _FIG["r"]
    _arc_path(doc, cx, cy, r, 0.0, 1.0, 1.0, (0.76, 0.77, 0.79))

    f0, f1 = src.get("from_start"), src.get("from_end")
    if not span or not isinstance(f0, int) or not isinstance(f1, int):
        # Nobody simulated this route, so there is no geometry to draw.
        # An empty ring with the right name beats an invented arc.
        _fig_caption(doc, cx, cy, r, src.get("template"),
                     f"{span:,} bp" if span else "size unknown")
        return

    from .figure import figure_features
    #: The stage arrow's verb and size sit level with the circle on the
    #: right. Claiming that band up front moves a label out of it
    #: instead of printing "EF-1-alpha pro" across "1,305 bp".
    placed: list = [(cy + 6, True), (cy, True), (cy - 6, True),
                    (cy - 12, True),
                    # and the plasmid's own name, under the circle
                    (cy - r - 16, True), (cy - r - 16, False),
                    (cy - r - 10, True), (cy - r - 10, False)]
    take_s, take_e = (f0 % span) / span, (f1 % span) / span
    for a in figure_features(src):
        s, e = (a["start"] % span) / span, (a["end"] % span) / span
        _arc_path(doc, cx, cy, r, s, e, 3.6, _feature_rgb(a["feature_type"]))
        # Only the promoter that ends up in the product earns a name. A
        # donor's own promoter is not part of what is being taken out of
        # it, and labelling it crowded the figure with a fact about the
        # parent that the construct never inherits.
        if a["feature_type"] == "promoter" and _within(s, e, take_s, take_e):
            mid = s + (((e - s) % 1.0) or 1.0) / 2
            _fig_label(doc, placed, cx, cy, r, mid,
                       (a.get("name") or "promoter")[:14],
                       _feature_rgb("promoter"))

    s, e = take_s, take_e
    _arc_path(doc, cx, cy, r, s, e, 4.6, lane["colour"])

    cuts = lane.get("cuts_parent") or []
    if cuts:
        second = cuts[1] if len(cuts) > 1 else cuts[0]
        # Two sites a few bases apart are one tick at this scale, and two
        # names stacked on it read as a smudge rather than as two
        # enzymes. Say so once.
        if min((s - e) % 1.0, (e - s) % 1.0) < 0.02:
            _fig_cut(doc, cx, cy, r, s, f"{cuts[0]}/{second}", placed)
        else:
            _fig_cut(doc, cx, cy, r, s, cuts[0], placed)
            _fig_cut(doc, cx, cy, r, e, second, placed)
    else:
        _fig_primer(doc, cx, cy, r, s, True, placed)
        _fig_primer(doc, cx, cy, r, e, False, placed)

    _fig_caption(doc, cx, cy, r, src.get("template"), f"{span:,} bp")


def _fig_order(doc, lane, cx, cy):
    """A piece nobody has a plasmid for: it gets ordered."""
    w, hgt = 74, 34
    x0, y0, x1, y1 = cx - w / 2, cy - hgt / 2, cx + w / 2, cy + hgt / 2
    for a, b in (((x0, y0), (x1, y0)), ((x1, y0), (x1, y1)),
                 ((x1, y1), (x0, y1)), ((x0, y1), (x0, y0))):
        doc.line(a[0], a[1], b[0], b[1], rgb=_FIG_GREY, width=0.8,
                 dash="3 2")
    doc.text_at("synthesise", cx, cy + 3, "F2", 6.6, rgb=_FIG_INK, centre=True)
    bp = (lane.get("src") or {}).get("bp") or 0
    doc.text_at(f"{bp:,} bp gBlock", cx, cy - 7, "F1", 6.2, rgb=_FIG_GREY,
                centre=True)


def _fig_arrow(doc, x, y, w, verb, sub):
    """A stage: the operation, an arrow, and what comes out of it."""
    x1, x2 = x + 7, x + w - 7
    doc.line(x1, y, x2 - 4, y, rgb=_FIG_GREY, width=0.8)
    doc.fill_tri((x2, y), (x2 - 5, y + 2.6), (x2 - 5, y - 2.6), _FIG_GREY)
    doc.text_at(verb, (x1 + x2) / 2, y + 6, "F2", 7, rgb=_FIG_INK, centre=True)
    if sub and text_width(sub, "F1", 6.2) <= w - 14:
        doc.text_at(sub, (x1 + x2) / 2, y - 12, "F1", 6.2, rgb=_FIG_GREY,
                    centre=True)


def _fig_fragment(doc, lane, index, joins, x, y, w):
    """The linear piece in the tube, with whatever its ends are.

    Both halves of a join wear the same colour and the same letter, so
    which end goes into which is something the reader sees rather than
    works out; two grey caps both labelled "overlap" said only that an
    overlap existed.
    """
    doc.fill_rect(x, y - 3, w, 6, rgb=lane["colour"])
    sticky = lane.get("end_kind") == "sticky"
    for i, side in enumerate(("left", "right")):
        pair = joins.get((index, side))
        fallback = (lane.get("ends") or [None, None])[i] \
            if len(lane.get("ends") or []) > i else None
        if not pair and not fallback:
            continue
        label = pair[0] if pair else fallback
        colour = pair[1] if pair else (CUT_RED if sticky
                                       else (0.62, 0.64, 0.67))
        ex = x if i == 0 else x + w
        doc.fill_rect(ex - 1.8, y - 5.5, 3.6, 11, rgb=colour)
        doc.text_at(label, ex, y - 15, "F1", 5.8, rgb=colour, centre=True)
    doc.text_at(str(lane.get("piece") or "")[:26], x + w / 2, y + 7, "F2", 7,
                rgb=_FIG_INK, centre=True)
    bp = (lane.get("src") or {}).get("bp") or 0
    doc.text_at(f"{bp:,} bp", x + w / 2, y - 15, "F1", 6.2, rgb=_FIG_GREY,
                centre=True)


def _fig_product(doc, lanes, construct, product_bp, cx, cy):
    """The new plasmid, each lane's contribution in that lane's colour."""
    span = product_bp or (construct or {}).get("length") or 1
    r = _FIG["pr"]
    _arc_path(doc, cx, cy, r, 0.0, 1.0, 1.0, (0.76, 0.77, 0.79))

    for lane in lanes:
        src = lane.get("src") or {}
        if not isinstance(src.get("at_start"), int):
            continue
        _arc_path(doc, cx, cy, r, src["at_start"] / span,
                  src["at_end"] / span, 4.2, lane["colour"])

    from .figure import product_features
    named = product_features((construct or {}).get("layout") or [], span)
    placed: list = [(cy - r - 16, True), (cy - r - 16, False),
                    (cy - r - 10, True), (cy - r - 10, False)]
    for a in named[:3]:
        s, e = (a["start"] % span) / span, (a["end"] % span) / span
        _arc_path(doc, cx, cy, r - 7, s, e, 3.6, _feature_rgb(a["type"]))
        mid = s + (((e - s) % 1.0) or 1.0) / 2
        _fig_label(doc, placed, cx, cy, r, mid,
                   (a.get("name") or a["type"])[:14], _feature_rgb(a["type"]))

    _fig_caption(doc, cx, cy, r, (construct or {}).get("name")
                 or "new plasmid", f"{span:,} bp")


def figure_height(lanes: list) -> float:
    """How much vertical room the figure wants, so a caller can keep the
    heading with it. ``_ensure`` inside the figure only protects the
    drawing; the heading above it had already gone on the previous page.
    """
    if not lanes:
        return 0.0
    return max(len(lanes) * _FIG["row"], 2 * _FIG["pr"] + 40) + 30


def cloning_figure(doc: "PDF", lanes: list, join: str,
                   construct: dict | None = None, product_bp: int = 0):
    """Draw the whole route: parents, stages, pieces, and the join."""
    if not lanes:
        return
    max_stages = max(len(l.get("stages") or []) for l in lanes) or 1
    x_parent = doc.ml
    x_stage = x_parent + _FIG["parent"]
    x_frag = x_stage + max_stages * _FIG["stage"]
    x_bracket = x_frag + _FIG["frag"] + 6
    x_join = x_bracket + _FIG["bracket"]
    x_prod = x_join + _FIG["stage"] + _FIG["pr"]

    rows = len(lanes)
    height = max(rows * _FIG["row"], 2 * _FIG["pr"] + 40)
    doc._ensure(height + 20)
    doc.y -= 10
    top_y = doc.y

    from .figure import junction_colours
    joins = junction_colours(lanes)
    centres = []
    for i, lane in enumerate(lanes):
        cy = top_y - _FIG["row"] / 2 - i * _FIG["row"]
        centres.append(cy)
        cx = x_parent + _FIG["parent"] / 2
        if lane.get("origin") == "order":
            _fig_order(doc, lane, cx, cy)
        else:
            _fig_parent(doc, lane, cx, cy)
        for k, stage in enumerate(lane.get("stages") or []):
            _fig_arrow(doc, x_stage + k * _FIG["stage"], cy, _FIG["stage"],
                       stage.get("verb", ""), stage.get("sub", ""))
        _fig_fragment(doc, lane, i, joins, x_frag, cy, _FIG["frag"])

    # One bracket gathering the pieces into the single reaction.
    if rows > 1:
        top, bot = centres[0], centres[-1]
        doc.line(x_bracket, top, x_bracket + 6, top, rgb=_FIG_GREY, width=0.8)
        doc.line(x_bracket + 6, top, x_bracket + 6, bot, rgb=_FIG_GREY,
                 width=0.8)
        doc.line(x_bracket + 6, bot, x_bracket, bot, rgb=_FIG_GREY, width=0.8)
    mid_y = (centres[0] + centres[-1]) / 2
    _fig_arrow(doc, x_join, mid_y, _FIG["stage"], join,
               f"{rows} piece" + ("" if rows == 1 else "s"))
    _fig_product(doc, lanes, construct, product_bp, x_prod, mid_y)

    doc.y = top_y - height - 6
    key = []
    if any(not (l.get("cuts_parent") or l.get("origin") == "order")
           for l in lanes):
        key.append("green arrows on a circle are primers")
    if any(l.get("cuts_parent") or l.get("cuts_fragment") for l in lanes):
        key.append("red ticks are enzyme cuts")
    if any(l.get("origin") == "order" for l in lanes):
        key.append("a dashed box is DNA nobody has, so it is ordered")
    key.append("each colour follows one piece into the product")
    if joins:
        key.append("matching end colours are compatible overhangs"
                   if lanes[0].get("end_kind") == "sticky"
                   else "the two ends sharing a colour and letter are one "
                        "overlap")
    doc.text_at("; ".join(key), doc.ml, doc.y, "F1", 6.2, rgb=(0.45, 0.45, 0.48))
    doc.y -= 10
