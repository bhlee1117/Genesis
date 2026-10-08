/* Genesis web client.
 *
 * Vanilla ES2020, no build step. Everything the server returns is lab-file
 * derived text, so it is only ever put into the DOM as a text node (see h()).
 */
'use strict';

/* ------------------------------------------------------------------ */
/* Controlled vocabulary                                               */
/* ------------------------------------------------------------------ */

/* The server owns the vocabulary and it has changed under us before (DIO
 * became Cre-lox), so the list below is only the offline fallback: init()
 * replaces it with GET /api/vocab before the first row is drawn. */
var FEATURE_TYPES = [
  'backbone', 'promoter', 'enhancer', 'CDS', 'trafficking_motif', 'linker',
  'p2A', 'IRES', 'Cre-lox', 'lox_site', 'AAV_ITR', 'resistance', 'origin',
  'tag', 'polyA', 'WPRE', 'restriction_site', 'primer_binding_site',
  'terminator', 'intron', 'unknown'
];

/* Mid-saturation hues, legible on both light and dark backgrounds. */
/* The lab builds AAVs by default; Auto remains one click away. */
var DEFAULT_VECTOR = 'AAV';

var TYPE_COLORS = {
  backbone: '#8a8f98', promoter: '#2f7f4f', enhancer: '#4f9e6a',
  CDS: '#2f6fb0', trafficking_motif: '#7a5fb0', linker: '#9aa0a8',
  p2A: '#b07a2f', IRES: '#b0652f', 'Cre-lox': '#a0398c', lox_site: '#c03a5b',
  AAV_ITR: '#b02f2f', resistance: '#6b8f2f', origin: '#2f8f8f',
  tag: '#8f6f2f', polyA: '#4a5fa8', WPRE: '#5f7f9e', restriction_site: '#8f4f6f',
  primer_binding_site: '#7f8f4f', terminator: '#6f5f8f', intron: '#5f8f7f',
  unknown: '#8a8f98',
  /* Parts indexed before the rename still carry the old type name. */
  DIO: '#a0398c'
};

function colorFor(type) {
  return TYPE_COLORS[type] || TYPE_COLORS.unknown;
}

/* ---- Cre-lox ------------------------------------------------------- */

/* The arrangements genesis/construct.py RECOMBINASE_MODES offers. `short` is
 * what the group header says at a glance; `long` is the tooltip. */
var CRE_LOX_TYPE = 'Cre-lox';

var CRE_LOX_MODES = [
  { key: 'DIO', title: 'double-floxed inverted orientation',
    short: 'payload inverted, Cre switches ON',
    long: 'DIO: the payload is stored inverted between lox2272 and loxP ' +
          'pairs, so nothing is expressed until Cre flips it — Cre switches ' +
          'expression ON.' },
  { key: 'FLEX', title: 'flip-excision switch',
    short: 'payload inverted, Cre switches ON',
    long: 'FLEX: the flip-excision arrangement of the same lox2272 / loxP ' +
          'pairs. The payload is stored inverted, so Cre switches ' +
          'expression ON.' },
  { key: 'FAS', title: 'loxP / loxFAS switch',
    short: 'payload inverted, Cre switches ON (loxFAS)',
    long: 'FAS: as DIO, but the second site pair is loxFAS instead of ' +
          'lox2272. The payload is stored inverted, so Cre switches ' +
          'expression ON.' },
  { key: 'DO', title: 'double-floxed, Cre-OFF',
    short: 'payload forward, Cre switches OFF',
    long: 'DO: the payload is stored forward and is expressed from the ' +
          'start; Cre inverts it away — Cre switches expression OFF.' }
];

var DEFAULT_CRE_LOX_MODE = 'DIO';

/* Module types a Cre-lox cassette encloses; it ends at the first module that
 * is not one of these (a polyA or a WPRE stays outside).
 * Mirrors CASSETTE_PAYLOAD_TYPES in genesis/construct.py. */
var CASSETTE_PAYLOAD_TYPES = ['CDS', 'trafficking_motif', 'linker', 'p2A',
                              'tag', 'IRES'];

function isCreLox(type) { return type === CRE_LOX_TYPE || type === 'DIO'; }

function creLoxMode(key) {
  var want = String(key || DEFAULT_CRE_LOX_MODE).toUpperCase();
  for (var i = 0; i < CRE_LOX_MODES.length; i++) {
    if (CRE_LOX_MODES[i].key === want) return CRE_LOX_MODES[i];
  }
  return CRE_LOX_MODES[0];
}

/** "DIO — payload inverted, Cre switches ON" */
function creLoxLabel(key) {
  var m = creLoxMode(key);
  return m.key + ' — ' + m.short;
}

/* ------------------------------------------------------------------ */
/* DOM helpers                                                         */
/* ------------------------------------------------------------------ */

/** h('div', {class:'x', onclick:fn}, 'text', node, [nodes]) */
function h(tag, attrs, ...kids) {
  var node = document.createElement(tag);
  applyAttrs(node, attrs);
  appendKids(node, kids);
  return node;
}

var SVG_NS = 'http://www.w3.org/2000/svg';

function svg(tag, attrs, ...kids) {
  var node = document.createElementNS(SVG_NS, tag);
  if (attrs) {
    Object.keys(attrs).forEach(function (k) {
      if (attrs[k] === null || attrs[k] === undefined) return;
      node.setAttribute(k, String(attrs[k]));
    });
  }
  appendKids(node, kids);
  return node;
}

function applyAttrs(node, attrs) {
  if (!attrs) return;
  Object.keys(attrs).forEach(function (k) {
    var v = attrs[k];
    if (v === null || v === undefined || v === false) return;
    if (k.indexOf('on') === 0 && typeof v === 'function') {
      node.addEventListener(k.slice(2), v);
    } else if (k === 'class') {
      node.className = v;
    } else if (k === 'dataset') {
      Object.keys(v).forEach(function (d) { node.dataset[d] = v[d]; });
    } else if (k === 'value' || k === 'checked' || k === 'disabled' || k === 'hidden') {
      node[k] = v;
    } else {
      node.setAttribute(k, String(v));
    }
  });
}

function appendKids(node, kids) {
  kids.forEach(function (kid) {
    if (kid === null || kid === undefined || kid === false) return;
    if (Array.isArray(kid)) { appendKids(node, kid); return; }
    // Strings become text nodes: server data is never parsed as markup.
    node.appendChild(kid instanceof Node ? kid : document.createTextNode(String(kid)));
  });
}

function clear(node) { while (node.firstChild) node.removeChild(node.firstChild); }

function $(id) { return document.getElementById(id); }

/** Build a <table> from rows + column definitions. */
function table(cols, rows, opts) {
  opts = opts || {};
  var head = h('tr', null, cols.map(function (c) {
    return h('th', { class: c.numeric ? 'num' : null }, c.label);
  }));
  var body = rows.map(function (row, i) {
    var tr = h('tr', {
      class: opts.rowClass ? opts.rowClass(row) : null,
      onclick: opts.onRow ? function () { opts.onRow(row, i); } : null
    }, cols.map(function (c) {
      var v = c.cell ? c.cell(row) : row[c.key];
      return h('td', { class: c.numeric ? 'num' : null },
        v === null || v === undefined ? '' : v);
    }));
    if (opts.onRow) tr.classList.add('clickable');
    return tr;
  });
  return h('table', { class: 'grid' },
    h('thead', null, head),
    h('tbody', null, body.length ? body : h('tr', null,
      h('td', { colspan: cols.length, class: 'muted' }, 'No rows.'))));
}

/* ------------------------------------------------------------------ */
/* Fetch helpers                                                       */
/* ------------------------------------------------------------------ */

/** JSON fetch. Throws Error carrying the server's `error` field. */
async function api(path, options) {
  var res;
  try {
    res = await fetch(path, options);
  } catch (netErr) {
    throw new Error('cannot reach the Genesis server (' + netErr.message + ')');
  }
  var body = null;
  var text = await res.text();
  if (text) {
    try { body = JSON.parse(text); } catch (e) { body = null; }
  }
  if (!res.ok) {
    var msg = (body && body.error) || text || ('HTTP ' + res.status);
    throw new Error(msg);
  }
  if (body === null) throw new Error('server returned a non-JSON response');
  return body;
}

function postJSON(path, payload) {
  return api(path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(payload)
  });
}

function qs(params) {
  var parts = [];
  Object.keys(params).forEach(function (k) {
    if (params[k] !== '' && params[k] !== null && params[k] !== undefined) {
      parts.push(encodeURIComponent(k) + '=' + encodeURIComponent(params[k]));
    }
  });
  return parts.length ? '?' + parts.join('&') : '';
}

function setStatus(node, text, kind) {
  clear(node);
  node.className = 'status' + (kind ? ' status-' + kind : '');
  if (text) node.appendChild(document.createTextNode(text));
}

function errorBox(message) {
  return h('div', { class: 'errorbox' },
    h('strong', null, 'Request failed'), h('span', null, message));
}

/* ------------------------------------------------------------------ */
/* Modal dialogs                                                       */
/* ------------------------------------------------------------------ */

/* One dialog at a time, built into #modal-host. Used for the two places
 * where a stray click must not be enough on its own: naming a Save as…, and
 * confirming a delete. Both are deliberately *not* window.confirm/prompt,
 * because the text has to spell out exactly what is about to be lost and a
 * failed save has to be able to re-ask without losing what was typed. */

var openModal = null;   // { close } for the dialog currently on screen

/**
 * opts: { title, body (Node|[Node]), confirmLabel, danger, onConfirm(close),
 *         focus (element to focus), confirmDisabled }
 * onConfirm receives a `close` callback and may keep the dialog open (for
 * instance to show a server error and let the user try again).
 */
function modal(opts) {
  if (openModal) openModal.close();

  var host = $('modal-host');
  var errorLine = h('p', { class: 'modal-error', hidden: true });
  var busy = false;

  var confirmBtn = h('button', {
    type: 'button',
    class: opts.danger ? 'danger-solid' : 'primary',
    onclick: function () { fire(); }
  }, opts.confirmLabel || 'Confirm');

  var cancelBtn = h('button', {
    type: 'button', class: 'secondary', onclick: function () { close(); }
  }, opts.cancelLabel || 'Cancel');

  var card = h('div', {
    class: 'modal-card', role: 'dialog', 'aria-modal': 'true',
    'aria-label': opts.title || 'Confirm',
    // A click inside must not reach the backdrop handler below.
    onclick: function (e) { e.stopPropagation(); }
  },
    h('h3', { class: 'modal-title' }, opts.title || ''),
    h('div', { class: 'modal-body' }, opts.body || null),
    errorLine,
    h('div', { class: 'modal-foot' }, cancelBtn, confirmBtn));

  var backdrop = h('div', { class: 'modal-backdrop', onclick: function () { close(); } }, card);

  var previous = document.activeElement;

  function onKey(e) {
    if (e.key === 'Escape') { e.preventDefault(); close(); return; }
    if (e.key !== 'Tab') return;
    // Keep focus inside the dialog.
    var focusable = card.querySelectorAll(
      'button:not(:disabled), input:not(:disabled), select:not(:disabled), textarea:not(:disabled)');
    if (!focusable.length) return;
    var first = focusable[0];
    var last = focusable[focusable.length - 1];
    if (e.shiftKey && document.activeElement === first) {
      e.preventDefault(); last.focus();
    } else if (!e.shiftKey && document.activeElement === last) {
      e.preventDefault(); first.focus();
    }
  }

  function close() {
    if (openModal !== api_) return;
    document.removeEventListener('keydown', onKey, true);
    clear(host);
    host.hidden = true;
    openModal = null;
    if (previous && previous.focus) previous.focus();
  }

  function setError(message) {
    clear(errorLine);
    errorLine.hidden = !message;
    if (message) errorLine.appendChild(document.createTextNode(message));
  }

  function setBusy(on) {
    busy = on;
    confirmBtn.disabled = on;
    cancelBtn.disabled = on;
  }

  function fire() {
    if (busy || confirmBtn.disabled) return;
    var result = opts.onConfirm ? opts.onConfirm(api_) : null;
    if (result && typeof result.then === 'function') {
      setBusy(true);
      result.then(function () { setBusy(false); },
                  function () { setBusy(false); });
    }
  }

  var api_ = {
    close: close, setError: setError, setBusy: setBusy,
    confirmButton: confirmBtn, card: card
  };

  clear(host);
  host.hidden = false;
  host.appendChild(backdrop);
  document.addEventListener('keydown', onKey, true);
  openModal = api_;
  (opts.focus || confirmBtn).focus();
  if (opts.focus && opts.focus.select) opts.focus.select();
  return api_;
}

/** POST a plan and save the returned attachment to disk. */
async function downloadExport(kind, payload, fallbackName) {
  var res = await fetch('/api/export/' + kind, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(payload)
  });
  if (!res.ok) {
    var msg = 'HTTP ' + res.status;
    try {
      var j = JSON.parse(await res.text());
      if (j && j.error) msg = j.error;
    } catch (e) { /* keep the status line */ }
    throw new Error(msg);
  }
  var disp = res.headers.get('Content-Disposition') || '';
  var match = /filename\*?=(?:UTF-8'')?"?([^";]+)"?/i.exec(disp);
  var name = match ? decodeURIComponent(match[1]) : fallbackName;
  var blob = await res.blob();
  var url = URL.createObjectURL(blob);
  var a = h('a', { href: url, download: name });
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(function () { URL.revokeObjectURL(url); }, 2000);
}

/* ------------------------------------------------------------------ */
/* Formatting                                                          */
/* ------------------------------------------------------------------ */

function bp(n) {
  if (n === null || n === undefined || n === '') return '';
  return Number(n).toLocaleString() + ' bp';
}

function hours(minutes) {
  if (!minutes && minutes !== 0) return '';
  return (Number(minutes) / 60).toFixed(1) + ' h';
}

function pct(fraction) {
  if (fraction === null || fraction === undefined) return '';
  return (Number(fraction) * 100).toFixed(1) + '%';
}

function num(value, digits) {
  if (value === null || value === undefined || value === '') return '';
  return Number(value).toFixed(digits);
}

function yesNo(v) { return v ? 'yes' : 'no'; }

function plural(n, word) { return n + ' ' + word + (n === 1 ? '' : 's'); }

/* Cell renderers shared by several tables. */

function typeTag(type) {
  return h('span', { class: 'typetag' },
    h('span', { class: 'swatch', style: 'background:' + colorFor(type) }), type);
}

function verifiedPill(value) {
  return h('span', { class: 'pill pill-' + (value === 'unverified' ? 'warn' : 'ok') },
    value || 'unverified');
}

function strandSign(strand) { return Number(strand) === -1 ? '−' : '+'; }

function warnCell(list) {
  return (list && list.length)
    ? h('span', { class: 'inline-warn' }, list.join('; ')) : '';
}

/* ------------------------------------------------------------------ */
/* Shared: linear feature map (inline SVG)                             */
/* ------------------------------------------------------------------ */

var MAP_AXIS = 76;       // y of the axis line when no label band is needed
var MAP_HEIGHT = 132;
var MAP_PAD = 12;
var FEAT_H = 18;         // height of a feature block
var FEAT_RISE = 44;      // axis -> top of the forward block
var FEAT_DROP = 26;      // axis -> top of the reverse block

/* Codon markers: one staggered label row is MARK_ROW_H tall, and a marker
 * that cannot claim one of MARK_ROWS rows keeps its tick but drops its label. */
var MARK_ROW_H = 13;
var MARK_ROWS = 2;
var MARK_FONT = 10;      // px; kept in sync with .mark-label in style.css

/**
 * How each `codon_markers[].kind` is drawn and described.
 *
 * `cls` picks the colour token in style.css. The severities mirror
 * genesis/validate.py: a premature stop is an error there, a missing start or
 * stop is a warning, so the map and the validation panel agree.
 */
var CODON_KINDS = {
  start:          { cls: 'mark-ok',   legend: 'start codon',    defect: false },
  stop:           { cls: 'mark-stop', legend: 'stop codon',     defect: false },
  premature_stop: { cls: 'mark-err',  legend: 'premature stop', defect: true },
  missing_start:  { cls: 'mark-warn', legend: 'missing start',  defect: true },
  missing_stop:   { cls: 'mark-warn', legend: 'missing stop',   defect: true }
};

/** Unknown kinds from a newer server still draw, in the neutral colour. */
function codonKind(kind) {
  return CODON_KINDS[kind] ||
    { cls: 'mark-stop', legend: String(kind || 'codon'), defect: false };
}

function typeOf(f) { return f.feature_type || f.type || 'unknown'; }

/** Rough advance width of a label; deliberately generous so labels never touch. */
function textWidth(s) { return String(s).length * (MARK_FONT * 0.6) + 6; }

/**
 * Place markers for one strand into staggered label rows.
 *
 * Defects claim rows first, so when a clean marker and a flagged one compete
 * for the same spot it is the clean one that loses its label. A marker that
 * fits nowhere gets `row = -1`: the tick and tooltip stay, the text goes.
 * Returns the number of rows actually used.
 */
function layoutMarkerRows(marks, width) {
  var rows = [];
  for (var r = 0; r < MARK_ROWS; r++) rows.push([]);

  var order = marks.map(function (m, i) { return i; });
  order.sort(function (a, b) {
    var da = marks[a].defect ? 0 : 1;
    var db = marks[b].defect ? 0 : 1;
    if (da !== db) return da - db;
    return marks[a].px - marks[b].px;
  });

  var used = 0;
  order.forEach(function (i) {
    var m = marks[i];
    var half = m.halfWidth;
    // Keep the text inside the drawing area; the leader slants back to the base.
    var lo = MAP_PAD + half;
    var hi = width - MAP_PAD - half;
    m.textX = hi < lo ? width / 2 : Math.max(lo, Math.min(m.px, hi));

    var a = m.textX - half;
    var b = m.textX + half;
    for (var r = 0; r < MARK_ROWS; r++) {
      var clash = rows[r].some(function (iv) { return a < iv[1] && b > iv[0]; });
      if (clash) continue;
      rows[r].push([a, b]);
      m.row = r;
      if (r + 1 > used) used = r + 1;
      return;
    }
    m.row = -1;
  });
  return used;
}

/** Tooltip text for one marker: label, module, position, strand. */
function markerTitle(m) {
  var pos = (m.end > m.start) ? (m.start + '–' + m.end) : String(m.start);
  return m.label + '\n' +
    'in ' + (m.module || '(unknown module)') + '\n' +
    'bp ' + pos + ' · strand ' + strandSign(m.strand) +
    (m.run ? '\nreading frame: ' + m.run : '');
}

/** Draw one marker: hit area, tick across the feature block, leader, label. */
function drawMarker(root, m, geom) {
  var g = svg('g', { class: 'mark ' + m.cls });
  g.appendChild(svg('title', null, markerTitle(m)));

  var tickFrom = geom.tickFrom;   // y at the edge of the block facing the labels
  var tickTo = geom.tickTo;       // y at the far edge of the block
  var baseline = null;
  if (m.row >= 0) {
    baseline = geom.up
      ? tickFrom - 6 - m.row * MARK_ROW_H
      : tickFrom + 13 + m.row * MARK_ROW_H;
  }

  // Invisible, comfortably wide hit area so the 2 px tick is easy to hover.
  var hy = [tickFrom, tickTo];
  if (baseline !== null) hy.push(geom.up ? baseline - MARK_FONT : baseline + 4);
  var top = Math.min.apply(null, hy);
  var bot = Math.max.apply(null, hy);
  g.appendChild(svg('rect', {
    x: Math.min(m.px, m.textX) - 6, y: top,
    width: Math.abs(m.textX - m.px) + 12, height: Math.max(bot - top, 4),
    class: 'mark-hit'
  }));

  g.appendChild(svg('line', {
    x1: m.px, y1: tickFrom, x2: m.px, y2: tickTo, class: 'mark-tick'
  }));
  g.appendChild(svg('circle', {
    cx: m.px, cy: tickFrom, r: m.defect ? 2.8 : 2, class: 'mark-head'
  }));

  if (baseline !== null) {
    g.appendChild(svg('line', {
      x1: m.textX, y1: baseline + (geom.up ? 3 : -8),
      x2: m.px, y2: tickFrom, class: 'mark-leader'
    }));
    g.appendChild(svg('text', {
      x: m.textX, y: baseline, class: 'mark-label', 'text-anchor': 'middle'
    }, m.label));
  }
  root.appendChild(g);
}

/** Draw the map at an explicit pixel width, so nothing is scaled or squashed. */
function buildMapSvg(width, span, list, opts) {
  var inner = Math.max(width - MAP_PAD * 2, 10);
  var x = function (pos) {
    return MAP_PAD + (Math.max(0, Math.min(pos, span)) / span) * inner;
  };

  /* Codon markers follow the feature convention: a forward marker sits above
   * the forward track, a reverse one below the reverse track, so an inverted
   * DIO payload and its ATG/stop always read as belonging together. */
  var marks = (opts.markers || []).map(function (m) {
    var meta = codonKind(m.kind);
    return {
      label: String(m.label || m.codon || m.kind || ''),
      module: m.module, run: m.run, strand: m.strand,
      start: m.start, end: (typeof m.end === 'number' ? m.end : m.start),
      cls: meta.cls, defect: meta.defect,
      reverse: Number(m.strand) === -1,
      px: x(m.start), halfWidth: 0, textX: 0, row: -1
    };
  });
  marks.forEach(function (m) { m.halfWidth = textWidth(m.label) / 2 + 3; });

  var fwd = marks.filter(function (m) { return !m.reverse; });
  var rev = marks.filter(function (m) { return m.reverse; });
  var fwdRows = layoutMarkerRows(fwd, width);
  var revRows = layoutMarkerRows(rev, width);
  if (opts.info) {
    opts.info.hidden = marks.filter(function (m) { return m.row < 0; }).length;
  }

  // With no markers these collapse to the original 76 / 132 geometry exactly.
  var axis = MAP_AXIS;
  if (fwdRows) axis = Math.max(axis, 66 + (fwdRows - 1) * MARK_ROW_H);
  var revBottom = axis + FEAT_DROP + FEAT_H;
  var height = revBottom +
    (revRows ? 17 + (revRows - 1) * MARK_ROW_H : MAP_HEIGHT - MAP_AXIS - FEAT_DROP - FEAT_H);

  var root = svg('svg', {
    width: width, height: height,
    viewBox: '0 0 ' + width + ' ' + height,
    class: 'featuremap', role: 'img'
  });
  root.appendChild(svg('title', null,
    (opts.label || 'Feature map') + ' — ' + span + ' bp'));

  // Axis plus evenly spaced ticks; fewer ticks when the panel is narrow.
  root.appendChild(svg('line', {
    x1: MAP_PAD, y1: axis, x2: width - MAP_PAD, y2: axis, class: 'axis-line'
  }));
  var ticks = width < 380 ? 4 : (width < 700 ? 6 : 10);
  for (var t = 0; t <= ticks; t++) {
    var pos = Math.round(span * t / ticks);
    var tx = x(pos);
    root.appendChild(svg('line', {
      x1: tx, y1: axis, x2: tx, y2: axis + 5, class: 'axis-tick'
    }));
    root.appendChild(svg('text', {
      x: tx, y: axis + 18, class: 'axis-label',
      'text-anchor': t === 0 ? 'start' : (t === ticks ? 'end' : 'middle')
    }, String(pos)));
  }

  list.forEach(function (f) {
    var type = typeOf(f);
    var reverse = Number(f.strand) === -1;
    var start = Math.min(f.start, f.end);
    var end = Math.max(f.start, f.end);
    var bx = x(start);
    var bw = Math.max(x(end) - bx, 1.5);
    var by = reverse ? axis + FEAT_DROP : axis - FEAT_RISE;

    var rect = svg('rect', {
      x: bx, y: by, width: bw, height: FEAT_H, rx: 2,
      fill: colorFor(type), class: 'feat' + (reverse ? ' feat-rev' : '')
    });
    rect.appendChild(svg('title', null,
      (f.name || '(unnamed)') + '\n' + type + '\n' +
      start + '–' + end + ' (' + (end - start) + ' bp), strand ' +
      (reverse ? '−' : '+') +
      (f.evidence ? '\nevidence: ' + f.evidence : '')));
    root.appendChild(rect);

    // Label only blocks wide enough to hold readable text (~6 px per char).
    if (bw > 44 && f.name) {
      root.appendChild(svg('text', {
        x: bx + bw / 2, y: by + 13, class: 'feat-label', 'text-anchor': 'middle'
      }, String(f.name).slice(0, Math.floor(bw / 6))));
    }
  });

  // Markers last, so their ticks read on top of the feature blocks.
  var fwdTop = axis - FEAT_RISE;
  fwd.forEach(function (m) {
    drawMarker(root, m, { up: true, tickFrom: fwdTop, tickTo: fwdTop + FEAT_H });
  });
  rev.forEach(function (m) {
    drawMarker(root, m, { up: false, tickFrom: revBottom, tickTo: revBottom - FEAT_H });
  });
  return root;
}

/**
 * items:        [{name, feature_type|type, start, end, strand}]
 * opts.markers: [{kind, start, end, strand, codon, module, run, label}]
 *
 * Forward features sit above the axis, reverse features below it, and codon
 * markers sit just outside their own strand's track. `markers` is optional:
 * designs saved before the server returned `codon_markers` draw exactly as
 * they always did, with no extra row and no empty legend.
 */
/* --------------------------------------------------------------------------
 * Circular plasmid map
 *
 * A plasmid is a circle, and a linear track implies the ends are ends. The
 * difference is not cosmetic: a cassette that wraps the coordinate origin -
 * which real pAAV backbones do constantly - looks broken on a linear map and
 * obvious on a circular one.
 *
 * Forward features sit on the outer ring, reverse features on an inner ring,
 * and anything translated carries an arrowhead at its 3' end so the direction
 * of translation can be read off the picture.
 * ------------------------------------------------------------------------ */

/* Types that are translated, and so have a direction worth drawing. */
var CODING_TYPES = {
  CDS: 1, p2A: 1, trafficking_motif: 1, tag: 1, linker: 1, IRES: 1
};

function polar(cx, cy, r, f) {
  var a = -Math.PI / 2 + 2 * Math.PI * f;           /* 0 at twelve o'clock */
  return [cx + r * Math.cos(a), cy + r * Math.sin(a)];
}

function arcPath(cx, cy, r, from, to) {
  var span = to - from;
  if (span <= 0) span += 1;
  if (span >= 0.999) span = 0.999;                  /* keep it an arc */
  var a = polar(cx, cy, r, from), b = polar(cx, cy, r, from + span);
  var large = span > 0.5 ? 1 : 0;
  return 'M ' + a[0].toFixed(2) + ' ' + a[1].toFixed(2) +
         ' A ' + r + ' ' + r + ' 0 ' + large + ' 1 ' +
         b[0].toFixed(2) + ' ' + b[1].toFixed(2);
}

/* A triangle at `f`, pointing the way the feature is read. */
function arrowHead(cx, cy, r, f, forward, fill) {
  var tipF = f + (forward ? 0.004 : -0.004);
  var tip = polar(cx, cy, r, tipF);
  var a = polar(cx, cy, r + 4.5, f);
  var b = polar(cx, cy, r - 4.5, f);
  return svg('polygon', {
    points: [tip, a, b].map(function (p) {
      return p[0].toFixed(2) + ',' + p[1].toFixed(2);
    }).join(' '),
    fill: fill
  });
}

function circularMap(length, items, opts) {
  opts = opts || {};
  var span = Math.max(Number(length) || 0, 1);
  var list = (items || []).filter(function (f) {
    return f && typeof f.start === 'number' && typeof f.end === 'number' &&
      f.type !== 'primer_binding_site' && f.type !== 'unknown';
  });

  var host = h('div', { class: 'map-circular' });
  var captionHost = h('div');
  var caption = null;
  var lastWidth = 0;

  function draw() {
    var w = Math.round(host.clientWidth) || 560;
    if (w === lastWidth) return;
    lastWidth = w;
    caption = null;
    clear(host);

    var size = Math.max(Math.min(w, 560), 320);
    var cx = size / 2, cy = size / 2;
    var radius = size / 2 - 86;                     /* room for labels */
    var node = svg('svg', {
      width: '100%', viewBox: '0 0 ' + size + ' ' + size,
      class: 'circ-svg', role: 'img',
      'aria-label': (opts.label || 'construct') + ', ' + span + ' base pairs, circular'
    });

    node.appendChild(svg('circle', {
      cx: cx, cy: cy, r: radius, fill: 'none',
      class: 'circ-backbone', 'stroke-width': 2
    }));

    /* position 1, at the top */
    node.appendChild(svg('line', {
      x1: cx, y1: cy - radius - 9, x2: cx, y2: cy - radius + 5,
      class: 'circ-origin'
    }));
    node.appendChild(svg('text', {
      x: cx, y: cy - radius - 13, 'text-anchor': 'middle', class: 'circ-tick'
    }, '1'));

    var placed = [];
    list.forEach(function (f) {
      var from = (f.start % span) / span;
      var to = (f.end % span) / span;
      if (f.end - f.start >= span) { from = 0; to = 0.999; }
      if (to <= from) to += 1;
      var reverse = f.strand === -1;
      var r = radius + (reverse ? -11 : 0);
      var colour = colorFor(f.type);

      var g = svg('g', { class: 'circ-feature' });
      g.appendChild(svg('title', null,
        f.name + ' · ' + f.type + ' · ' + bp(f.end - f.start) +
        ' · ' + (reverse ? 'reverse' : 'forward') + ' strand · ' +
        (f.start + 1) + '–' + f.end));
      g.appendChild(svg('path', {
        d: arcPath(cx, cy, r, from, to), fill: 'none',
        stroke: colour, 'stroke-width': 8, 'stroke-linecap': 'butt'
      }));
      /* direction of translation, at the 3' end of the coding sequence */
      if (CODING_TYPES[f.type]) {
        g.appendChild(arrowHead(cx, cy, r, reverse ? from : to % 1,
                                !reverse, colour));
      }
      node.appendChild(g);
      placed.push({ mid: ((from + to) / 2) % 1, f: f, r: r, colour: colour });
    });

    /* labels outside the circle, nudged apart so none overprint */
    placed.sort(function (a, b) { return a.mid - b.mid; });
    var usedLeft = [], usedRight = [];
    placed.forEach(function (item) {
      if (item.f.end - item.f.start < span * 0.01) return;
      var anchor = polar(cx, cy, item.r + 5, item.mid);
      var out = polar(cx, cy, radius + 16, item.mid);
      var right = out[0] >= cx;
      var used = right ? usedRight : usedLeft;
      var y = out[1];
      while (used.some(function (v) { return Math.abs(v - y) < 11; })) {
        y += right ? 11 : -11;
      }
      used.push(y);
      var x = right ? cx + radius + 18 : cx - radius - 18;
      node.appendChild(svg('polyline', {
        points: anchor[0].toFixed(1) + ',' + anchor[1].toFixed(1) + ' ' +
                out[0].toFixed(1) + ',' + out[1].toFixed(1) + ' ' +
                x.toFixed(1) + ',' + y.toFixed(1),
        class: 'circ-leader', fill: 'none'
      }));
      node.appendChild(svg('text', {
        x: right ? x + 4 : x - 4, y: y + 3,
        'text-anchor': right ? 'start' : 'end',
        class: 'circ-label', fill: item.colour
      }, item.f.name || item.f.type));
    });

    /* The region being cut from a recipient, or taken from a donor. */
    var hl = opts.highlight;
    if (hl && hl.length) {
      var kind = opts.highlightKind === 'take' ? 'take' : 'cut';
      var hs = (hl.start % span) / span;
      var he = (hl.end % span) / span;
      if (he <= hs) he += 1;
      node.appendChild(svg('path', {
        d: arcPath(cx, cy, radius + 9, hs, he), fill: 'none',
        class: 'circ-highlight hl-' + kind, 'stroke-width': 3
      }));
      [hs, he % 1].forEach(function (edge) {
        var a = polar(cx, cy, radius + 3, edge), b = polar(cx, cy, radius + 15, edge);
        node.appendChild(svg('line', {
          x1: a[0].toFixed(1), y1: a[1].toFixed(1),
          x2: b[0].toFixed(1), y2: b[1].toFixed(1),
          class: 'circ-highlight hl-' + kind
        }));
      });
      /* The caption goes under the map, not inside it: element names make
         it wider than the circle and it would cross the ring. */
      caption = h('p', { class: 'circ-caption hl-' + kind },
        (kind === 'cut' ? 'cut out ' : 'take ') + bp(hl.length) +
        (hl.label ? ' \u2014 ' + hl.label : ''));
    }

    (opts.markers || []).forEach(function (m) {
      if (typeof m.start !== 'number') return;
      var f = (m.start % span) / span;
      var a = polar(cx, cy, radius - 16, f), b = polar(cx, cy, radius + 5, f);
      var g = svg('g', { class: 'circ-codon ' + ('mark-' + m.kind) });
      g.appendChild(svg('title', null,
        m.label + ' · ' + m.module + ' · bp ' + (m.start + 1)));
      g.appendChild(svg('line', {
        x1: a[0].toFixed(1), y1: a[1].toFixed(1),
        x2: b[0].toFixed(1), y2: b[1].toFixed(1)
      }));
      node.appendChild(g);
    });

    /* Only the size goes in the middle. A plasmid name is long enough to
       run off both sides of the circle and overprint the ring, and the card
       header already carries it. */
    node.appendChild(svg('text', {
      x: cx, y: cy + 2, 'text-anchor': 'middle', class: 'circ-size'
    }, bp(span)));
    node.appendChild(svg('text', {
      x: cx, y: cy + 16, 'text-anchor': 'middle', class: 'circ-topo'
    }, 'circular'));

    host.appendChild(node);
    clear(captionHost);
    if (caption) captionHost.appendChild(caption);
  }

  // Draw once up front: the element is not in the document yet, so a
  // ResizeObserver has nothing to report and would leave the map blank
  // until something happened to resize it.
  draw();
  if (typeof ResizeObserver === 'function') new ResizeObserver(draw).observe(host);
  else requestAnimationFrame(draw);

  var note = h('p', { class: 'map-note muted small' },
    'Outer ring reads forward, inner ring reverse. Arrowheads point the way ' +
    'a coding sequence is translated. Position 1 is at the top.');
  return h('div', null, host, captionHost, note);
}

function featureMap(length, items, opts) {
  opts = opts || {};
  var span = Math.max(Number(length) || 0, 1);
  var list = (items || []).filter(function (f) {
    return f && typeof f.start === 'number' && typeof f.end === 'number';
  });
  var marks = (opts.markers || []).filter(function (m) {
    return m && typeof m.start === 'number';
  });

  // Row assignment depends on pixel width, so how many labels had to be
  // dropped is only known after a draw; the footnote is refreshed from here.
  var info = { hidden: 0 };
  var mapOpts = {
    label: opts.label, circular: opts.circular, markers: marks, info: info
  };

  var host = h('div', { class: 'map-svg' });
  var hiddenNote = h('span');
  var lastWidth = 0;
  function draw() {
    var w = Math.round(host.clientWidth) || 600;
    if (w === lastWidth) return;
    lastWidth = w;
    clear(host);
    host.appendChild(buildMapSvg(w, span, list, mapOpts));
    clear(hiddenNote);
    if (info.hidden) {
      hiddenNote.appendChild(document.createTextNode(
        ' · ' + info.hidden + ' label' + (info.hidden === 1 ? '' : 's') +
        ' hidden (hover the tick)'));
    }
  }
  if (typeof ResizeObserver === 'function') new ResizeObserver(draw).observe(host);
  else requestAnimationFrame(draw);

  var legendTypes = [];
  list.forEach(function (f) {
    if (legendTypes.indexOf(typeOf(f)) < 0) legendTypes.push(typeOf(f));
  });
  var legendKinds = [];
  marks.forEach(function (m) {
    if (legendKinds.indexOf(m.kind) < 0) legendKinds.push(m.kind);
  });

  var defects = marks.filter(function (m) { return codonKind(m.kind).defect; });

  return h('div', { class: 'map-wrap' },
    host,
    h('div', { class: 'map-foot muted' },
      bp(span) + (opts.circular === undefined ? '' :
        (opts.circular ? ' · circular' : ' · linear')) +
      ' · ' + list.length + ' feature' + (list.length === 1 ? '' : 's') +
      ' · forward above, reverse below' +
      (marks.length
        ? ' · ' + marks.length + ' codon marker' + (marks.length === 1 ? '' : 's') +
          (defects.length ? ', ' + defects.length + ' flagged' : '')
        : ''),
      hiddenNote),
    (legendTypes.length || legendKinds.length)
      ? h('div', { class: 'legend' },
          legendTypes.map(function (t) {
            return h('span', { class: 'legend-item' },
              h('span', { class: 'swatch', style: 'background:' + colorFor(t) }), t);
          }),
          legendKinds.map(function (k) {
            var meta = codonKind(k);
            return h('span', {
              class: 'legend-item legend-codon ' + meta.cls +
                (meta.defect ? ' legend-defect' : '')
            }, h('span', { class: 'legend-tick' }), meta.legend);
          }))
      : null);
}

/* ------------------------------------------------------------------ */
/* Application state                                                   */
/* ------------------------------------------------------------------ */

var state = {
  tab: 'design',
  modules: [newModule('promoter')],
  plan: null,
  spec: null,
  routeId: null,
  appliedSuggestion: null,  // suggested_name the user (or we) already accepted
  /* The saved design the plan on screen is tied to, if any: set when a design
   * is loaded from the Designs tab and when a save returns a row. Overwrite
   * is only possible while this is known, because the server replaces a
   * version by id and there is nothing to replace without one. */
  saved: null               // { design_id, construct_id, version, name }
};

/** A builder row. `part_id` pins an exact inventory part when one was picked. */
function newModule(type) {
  return {
    type: type, name: '', part_id: null, form: 'all',
    // Only meaningful on a Cre-lox row; kept on the row either way so that
    // switching the type away and back does not forget the arrangement.
    mode: isCreLox(type) ? DEFAULT_CRE_LOX_MODE : ''
  };
}

/* ------------------------------------------------------------------ */
/* Controlled vocabulary, from the server                              */
/* ------------------------------------------------------------------ */

/**
 * Replace the built-in feature-type list with the server's own.
 *
 * The vocabulary is the server's to define — it renamed DIO to Cre-lox — so
 * the type selects are filled from it rather than from a list in here that
 * can silently drift. A failed request leaves the fallback in place: the
 * builder still works, it just offers the types this file knows about.
 */
async function loadVocab() {
  try {
    var v = await api('/api/vocab');
    if (Array.isArray(v.feature_types) && v.feature_types.length) {
      FEATURE_TYPES = v.feature_types.slice();
    }
  } catch (err) {
    /* Offline or an older server: keep the fallback list. */
  }
}

/* ------------------------------------------------------------------ */
/* Header stats                                                        */
/* ------------------------------------------------------------------ */

async function loadStats() {
  var node = $('stats');
  clear(node);
  node.appendChild(h('span', { class: 'stat-loading' }, 'loading stats…'));
  try {
    var s = await api('/api/stats');
    clear(node);
    [['plasmid', s.plasmids], ['part', s.parts], ['owner', s.owners],
     ['design', s.designs]].forEach(function (pair) {
      var n = pair[1] || 0;
      node.appendChild(h('span', { class: 'stat' },
        h('b', null, String(n)), ' ', pair[0] + (n === 1 ? '' : 's')));
    });
    populatePartTypes(s.by_type || []);
  } catch (err) {
    clear(node);
    node.appendChild(h('span', { class: 'stat-error' }, err.message));
  }
}

/* ------------------------------------------------------------------ */
/* Tab 1: Inventory                                                    */
/* ------------------------------------------------------------------ */

async function loadOwners() {
  var sel = $('inv-owner');
  try {
    var owners = await api('/api/owners');
    owners.forEach(function (o) {
      if (!o.owner) return;
      sel.appendChild(h('option', { value: o.owner }, o.owner + ' (' + o.n + ')'));
    });
  } catch (err) {
    sel.appendChild(h('option', { value: '', disabled: true }, 'owners unavailable'));
  }
}

async function renderInventory() {
  var out = $('inv-table');
  var status = $('inv-status');
  setStatus(status, 'searching…', 'busy');
  clear(out);
  out.appendChild(h('p', { class: 'muted' }, 'Loading plasmids…'));
  try {
    var rows = await api('/api/plasmids' + qs({
      q: $('inv-q').value.trim(), owner: $('inv-owner').value, limit: 300
    }));
    clear(out);
    out.appendChild(table([
      { label: 'Name', key: 'name' },
      { label: 'Owner', key: 'owner' },
      { label: 'Length', numeric: true, cell: function (r) { return bp(r.seq_length); } },
      { label: 'Topology', cell: function (r) { return r.circular ? 'circular' : 'linear'; } },
      { label: 'Features', numeric: true, key: 'n_features' },
      { label: 'Verified', cell: function (r) { return verifiedPill(r.verified); } },
      { label: 'Annotation', cell: function (r) {
          return h('span', { class: 'pill' }, r.annotation_status || ''); } }
    ], rows, { onRow: function (r) { showPlasmid(r.id); } }));
    setStatus(status, plural(rows.length, 'plasmid'));
  } catch (err) {
    clear(out);
    out.appendChild(errorBox(err.message));
    setStatus(status, 'failed', 'error');
  }
}

async function showPlasmid(id) {
  var side = $('inv-detail');
  clear(side);
  side.appendChild(h('p', { class: 'muted' }, 'Loading plasmid ' + id + '…'));
  try {
    var p = await api('/api/plasmid/' + id);
    clear(side);
    var meta = [
      ['Owner', p.owner], ['Length', bp(p.seq_length)],
      ['Topology', p.circular ? 'circular' : 'linear'],
      ['Verified', p.verified], ['Annotation', p.annotation_status],
      ['Template available', yesNo(p.available_template)],
      ['Last verified', p.last_verified], ['Added', p.created_at],
      ['Source file', p.source ? p.source.rel_path : '(none)'],
      ['Source format', p.source ? p.source.format : ''],
      ['Description', p.description], ['Notes', p.notes]
    ];
    side.appendChild(h('h3', null, p.name || '(unnamed)'));
    side.appendChild(h('dl', { class: 'meta' }, meta.map(function (m) {
      if (!m[1] && m[1] !== 0) return null;
      return [h('dt', null, m[0]), h('dd', null, String(m[1]))];
    })));

    if (p.alt_sources && p.alt_sources.length) {
      side.appendChild(h('p', { class: 'muted small' },
        'Also found in: ' + p.alt_sources.map(function (a) { return a.rel_path; }).join(', ')));
    }

    side.appendChild(h('h4', null, 'Feature map'));
    side.appendChild(featureMap(p.seq_length, p.annotations || [],
      { label: p.name, circular: !!p.circular }));

    side.appendChild(h('h4', null, 'Annotations (' + (p.annotations || []).length + ')'));
    side.appendChild(h('div', { class: 'scroll-y mid' }, table([
      { label: 'Type', cell: function (a) { return typeTag(a.feature_type); } },
      { label: 'Name', key: 'name' },
      { label: 'Start', numeric: true, key: 'start' },
      { label: 'End', numeric: true, key: 'end' },
      { label: 'Strand', numeric: true, cell: function (a) { return strandSign(a.strand); } },
      { label: 'Evidence', key: 'evidence' }
    ], p.annotations || [])));
  } catch (err) {
    clear(side);
    side.appendChild(errorBox(err.message));
  }
}

/* ------------------------------------------------------------------ */
/* Tab 2: Parts                                                        */
/* ------------------------------------------------------------------ */

function populatePartTypes(byType) {
  var sel = $('part-type');
  if (sel.dataset.filled) return;
  sel.dataset.filled = '1';
  var counts = {};
  byType.forEach(function (r) { counts[r.feature_type] = r.n; });
  // Parts indexed under a type the vocabulary has since renamed (DIO before
  // it became Cre-lox) are still in the inventory, so list those too rather
  // than leaving them unreachable from this filter.
  var types = FEATURE_TYPES.concat(Object.keys(counts).filter(function (t) {
    return t && FEATURE_TYPES.indexOf(t) < 0;
  }));
  types.forEach(function (t) {
    var n = counts[t];
    sel.appendChild(h('option', { value: t },
      t + (n ? ' (' + n + ')' : '')));
  });
}

async function renderParts() {
  var out = $('part-table');
  var status = $('part-status');
  setStatus(status, 'searching…', 'busy');
  clear(out);
  out.appendChild(h('p', { class: 'muted' }, 'Loading parts…'));
  try {
    var rows = await api('/api/parts' + qs({
      type: $('part-type').value, q: $('part-q').value.trim()
    }));
    clear(out);
    out.appendChild(table([
      { label: 'Name', key: 'name' },
      { label: 'Type', cell: function (r) { return typeTag(r.feature_type); } },
      { label: 'Length', numeric: true, cell: function (r) { return bp(r.seq_length); } },
      { label: 'Sources', numeric: true, key: 'n_sources' }
    ], rows, { onRow: function (r) { showPartSources(r); } }));
    setStatus(status, plural(rows.length, 'part'));
  } catch (err) {
    clear(out);
    out.appendChild(errorBox(err.message));
    setStatus(status, 'failed', 'error');
  }
}

async function showPartSources(part) {
  var side = $('part-detail');
  clear(side);
  side.appendChild(h('h3', null, part.name || '(unnamed)'));
  side.appendChild(h('p', { class: 'muted small' },
    part.feature_type + ' · ' + bp(part.seq_length)));
  var holder = h('div', null, h('p', { class: 'muted' }, 'Loading sources…'));
  side.appendChild(holder);
  try {
    var rows = await api('/api/part/' + part.id + '/sources');
    clear(holder);
    holder.appendChild(h('h4', null, 'Source plasmids (' + rows.length + ')'));
    holder.appendChild(table([
      { label: 'Plasmid', key: 'name' },
      { label: 'Owner', key: 'owner' },
      { label: 'Verified', cell: function (r) { return verifiedPill(r.verified); } },
      { label: 'Length', numeric: true, cell: function (r) { return bp(r.seq_length); } },
      { label: 'Start', numeric: true, key: 'start' },
      { label: 'End', numeric: true, key: 'end' },
      { label: 'Strand', numeric: true, cell: function (r) { return strandSign(r.strand); } }
    ], rows, { onRow: function (r) { selectTab('inventory'); showPlasmid(r.id); } }));
  } catch (err) {
    clear(holder);
    holder.appendChild(errorBox(err.message));
  }
}

/* ------------------------------------------------------------------ */
/* Tab 3: Design - module builder                                      */
/* ------------------------------------------------------------------ */

/* The inventory does not change while the page is open, so repeated lookups
   for the same (type, query) are answered from memory. */
var PART_CACHE = {};
var COMBO_SEQ = 0;
var COMBO_LIMIT = 50;      // rows shown; the endpoint itself caps at 300

function fetchParts(type, query, form) {
  var key = type + '\u0000' + query.toLowerCase() + '\u0000' + form;
  if (!PART_CACHE[key]) {
    PART_CACHE[key] = api('/api/parts' + qs({ type: type, q: query, form: form }))
      .catch(function (err) { delete PART_CACHE[key]; throw err; });
  }
  return PART_CACHE[key];
}

var SEGMENT_CACHE = {};

function fetchSegments(partId) {
  if (!SEGMENT_CACHE[partId]) {
    SEGMENT_CACHE[partId] = api('/api/part/' + partId + '/segments')
      .catch(function (err) { delete SEGMENT_CACHE[partId]; throw err; });
  }
  return SEGMENT_CACHE[partId];
}

/**
 * Split a generated segment name into its stem and its "[1 of 2]" suffix.
 *
 * The suffix is the whole point of the row — it says which protein of the
 * fusion this is — so it is rendered as its own element that never wraps or
 * truncates; the long stem is what gives way when the list is narrow.
 */
function splitSegmentName(name) {
  var m = /^(.*?)\s*(\[\s*\d+\s+of\s+\d+\s*\])\s*$/.exec(String(name || ''));
  return m ? { stem: m[1], suffix: m[2] } : { stem: String(name || ''), suffix: '' };
}

/** How a /api/parts row relates to a fusion: 'fusion', 'piece' or ''. */
function partForm(row) {
  if (Number(row.n_segments) > 0) return 'fusion';
  if (row.derived_from) return 'piece';
  return '';
}

/**
 * Typeahead over /api/parts for one module row.
 *
 * `mod` is the row's model object, so the picker keeps working after rows are
 * reordered or removed (no index lookups). Choosing a row writes both the
 * name and `part_id`; typing by hand drops `part_id` again, because the text
 * no longer necessarily names that part.
 *
 * Many lab CDS annotations are fusions, and the inventory catalogues both the
 * fusion and the proteins it splits into. Both are legitimate things to pin,
 * so the list marks which is which, a fusion row can be opened to show its
 * own proteins, and a `form` filter narrows the list to one or the other.
 */
function partPicker(mod) {
  var listId = 'combo-list-' + (++COMBO_SEQ);

  var input = h('input', {
    type: 'text', value: mod.name, autocomplete: 'off', spellcheck: 'false',
    placeholder: 'part name — type to search the inventory (optional)',
    role: 'combobox', 'aria-autocomplete': 'list', 'aria-expanded': 'false',
    'aria-controls': listId
  });
  var list = h('div', { class: 'combo-list', id: listId, role: 'listbox', hidden: true });
  var node = h('div', { class: 'combo' }, input, list);
  var chip = h('span', { class: 'pin-chip', hidden: true });

  var formSel = h('select', {
    class: 'form-filter',
    title: 'Show whole fusions, the individual proteins they split into, or both',
    'aria-label': 'Part form filter',
    onchange: function (e) {
      mod.form = e.target.value;
      refresh();
      // Only re-open the list if the user is still in the text box; changing
      // the filter must not pop a dropdown over whatever else is on screen.
      if (document.activeElement === input && input.value.trim()) search();
    }
  },
    h('option', { value: 'all' }, 'all forms'),
    h('option', { value: 'whole' }, 'whole only'),
    h('option', { value: 'split' }, 'pieces only'));
  formSel.value = mod.form || 'all';

  var items = [];      // rows currently shown; see rowFor() for the shape
  var active = -1;     // highlighted row, -1 for none
  var timer = null;
  var reqId = 0;
  var view = 'parts';  // 'parts' | 'segments'
  var segmentsOf = null;

  function renderChip() {
    clear(chip);
    if (!mod.part_id) { chip.hidden = true; return; }
    chip.hidden = false;
    var extra = mod.pin_note ? ' ' + mod.pin_note : '';
    chip.appendChild(h('span', {
      class: 'pin-label',
      title: 'This module is pinned to inventory part ' + mod.part_id +
             (mod.pin_note ? ' (' + mod.pin_note + ' of its fusion)' : '') +
             '. The server will use exactly this part instead of guessing.'
    }, 'pinned #' + mod.part_id + extra));
    chip.appendChild(h('button', {
      type: 'button', class: 'pin-clear', title: 'Unpin — match this name loosely again',
      onclick: function () {
        mod.part_id = null; mod.pin_note = ''; renderChip(); input.focus();
      }
    }, '×'));
  }

  function close() {
    list.hidden = true;
    input.setAttribute('aria-expanded', 'false');
    input.removeAttribute('aria-activedescendant');
    active = -1;
    view = 'parts';
    segmentsOf = null;
  }

  function paint() {
    var kids = list.querySelectorAll('.combo-item');
    Array.prototype.forEach.call(kids, function (el, i) {
      el.classList.toggle('active', i === active);
      el.setAttribute('aria-selected', i === active ? 'true' : 'false');
    });
    var el = active >= 0 ? kids[active] : null;
    if (!el) { input.removeAttribute('aria-activedescendant'); return; }
    input.setAttribute('aria-activedescendant', el.id);
    if (el.offsetTop < list.scrollTop) {
      list.scrollTop = el.offsetTop;
    } else if (el.offsetTop + el.offsetHeight > list.scrollTop + list.clientHeight) {
      list.scrollTop = el.offsetTop + el.offsetHeight - list.clientHeight;
    }
  }

  /** Pin one row — a whole part, a fusion, or one segment of a fusion. */
  function choose(i) {
    var row = items[i];
    if (!row) return;
    if (row.kind === 'back') { showParts(); return; }
    mod.name = row.name || '';
    mod.part_id = row.id;
    mod.pin_note = row.note || '';
    input.value = mod.name;
    renderChip();
    close();
  }

  function show(node_) {
    clear(list);
    list.appendChild(node_);
    list.hidden = false;
    input.setAttribute('aria-expanded', 'true');
  }

  /** Render `items` into the dropdown. Each item supplies its own cells. */
  function paintList() {
    clear(list);
    items.forEach(function (r, i) {
      var el = h('div', {
        class: 'combo-item' + (r.kind === 'back' ? ' combo-back' : '') +
               (r.kind === 'segment' ? ' combo-segment' : ''),
        role: 'option', id: listId + '-' + i, 'aria-selected': 'false',
        // mousedown fires before blur, and preventDefault keeps focus in the
        // input, so the click always lands on a list that is still open.
        onmousedown: function (e) { e.preventDefault(); choose(i); },
        onmousemove: function () { if (active !== i) { active = i; paint(); } }
      }, r.cells());
      list.appendChild(el);
    });
    list.scrollTop = 0;
    list.hidden = false;
    input.setAttribute('aria-expanded', 'true');
  }

  function segBadge(kind, text, title) {
    return h('span', { class: 'seg-badge seg-badge-' + kind, title: title }, text);
  }

  function nameCell(name) {
    var parts = splitSegmentName(name);
    return h('span', { class: 'combo-name' },
      h('span', { class: 'combo-stem' }, parts.stem || '(unnamed)'),
      parts.suffix
        ? h('span', {
            class: 'seg-suffix',
            title: 'this is one protein of a fusion CDS'
          }, parts.suffix)
        : null);
  }

  function sourcesCell(n) {
    n = Number(n || 0);
    return h('span', {
      class: 'combo-sources' + (n ? '' : ' zero'),
      title: n + ' plasmid' + (n === 1 ? '' : 's') + ' in the inventory carry this part'
    }, n + ' src');
  }

  /** One /api/parts row. */
  function rowFor(r) {
    var form = partForm(r);
    return {
      kind: 'part', id: r.id, name: r.name || '', source: r,
      note: form === 'piece'
        ? '[' + r.segment_index + ' of ' + r.segment_of + ']' : '',
      cells: function () {
        return [
          nameCell(r.name),
          form === 'fusion'
            ? segBadge('fusion', 'fusion · ' + r.n_segments,
                'a fusion CDS catalogued as ' + r.n_segments +
                ' separate proteins — open it to pin one of them')
            : null,
          form === 'piece'
            ? segBadge('piece',
                'piece ' + r.segment_index + ' of ' + r.segment_of,
                'protein ' + r.segment_index + ' of ' + r.segment_of +
                ' of a fusion CDS')
            : null,
          h('span', { class: 'combo-len mono' }, bp(r.seq_length)),
          sourcesCell(r.n_sources),
          form === 'fusion'
            ? h('button', {
                type: 'button', class: 'seg-open',
                title: 'Show the ' + r.n_segments + ' proteins this fusion splits into',
                tabindex: '-1',
                onmousedown: function (e) {
                  e.preventDefault(); e.stopPropagation(); showSegments(r);
                }
              }, 'split ▸')
            : null
        ];
      }
    };
  }

  function showParts() {
    view = 'parts';
    segmentsOf = null;
    search();
  }

  /** Swap the dropdown for the proteins `fusion` splits into. */
  function showSegments(fusion) {
    var mine = ++reqId;
    view = 'segments';
    segmentsOf = fusion;
    show(h('div', { class: 'combo-empty muted' }, 'Loading the parts of this fusion…'));
    fetchSegments(fusion.id).then(function (data) {
      if (mine !== reqId) return;
      var whole = data.whole || {};
      var segs = data.segments || [];
      items = [{
        kind: 'back', cells: function () {
          return [h('span', { class: 'combo-name' }, '◂ back to the search results')];
        }
      }, {
        kind: 'part', id: whole.id, name: whole.name || '', note: '',
        cells: function () {
          return [nameCell(whole.name),
            segBadge('fusion', 'whole fusion', 'pin the fusion as one part'),
            h('span', { class: 'combo-len mono' }, bp(whole.seq_length))];
        }
      }].concat(segs.map(function (s) {
        return {
          kind: 'segment', id: s.id, name: s.name || '',
          note: '[' + s.segment_index + ' of ' + s.segment_of + ']',
          cells: function () {
            return [nameCell(s.name),
              segBadge('piece', 'piece ' + s.segment_index + ' of ' + s.segment_of,
                'pin only this protein'),
              h('span', { class: 'combo-len mono' }, bp(s.seq_length))];
          }
        };
      }));
      active = segs.length ? 2 : 1;
      paintList();
      if ((data.separators || []).length) {
        list.appendChild(h('div', { class: 'combo-foot muted' },
          'split at ' + data.separators.map(function (s) {
            return (s.name || s.feature_type || '?');
          }).join(', ')));
      }
      paint();
    }).catch(function (err) {
      if (mine !== reqId) return;
      items = [];
      show(h('div', { class: 'combo-empty inline-warn' },
        'could not load the parts of this fusion: ' + err.message));
    });
  }

  function draw(rows, query) {
    items = rows.map(rowFor);
    active = -1;
    if (!rows.length) {
      var filter = (mod.form === 'whole') ? ' whole' :
        (mod.form === 'split' ? ' fusion-piece' : '');
      show(h('div', { class: 'combo-empty muted' },
        'No' + filter + ' ' + mod.type + ' part matching “' + query + '”'));
      return;
    }
    paintList();
  }

  function search() {
    var query = input.value.trim();
    if (query.length < 1) { items = []; close(); return; }
    var mine = ++reqId;
    var type = mod.type;
    var form = mod.form || 'all';
    view = 'parts';
    segmentsOf = null;
    fetchParts(type, query, form).then(function (rows) {
      // Ignore a slow response that a later keystroke or type change replaced.
      if (mine !== reqId || type !== mod.type) return;
      draw(rows.slice(0, COMBO_LIMIT), query);
    }).catch(function (err) {
      if (mine !== reqId) return;
      items = [];
      show(h('div', { class: 'combo-empty inline-warn' }, 'search failed: ' + err.message));
    });
  }

  function schedule() {
    if (timer) clearTimeout(timer);
    timer = setTimeout(search, 180);
  }

  function refresh() {
    if (timer) clearTimeout(timer);
    reqId++;
    items = [];
    close();
  }

  input.addEventListener('input', function () {
    mod.name = input.value;
    if (mod.part_id) { mod.part_id = null; mod.pin_note = ''; renderChip(); }
    schedule();
  });
  input.addEventListener('focus', function () {
    if (input.value.trim()) search();
  });
  input.addEventListener('blur', close);
  input.addEventListener('keydown', function (e) {
    if (e.key === 'ArrowDown') {
      e.preventDefault();
      if (list.hidden) { search(); return; }
      if (items.length) { active = (active + 1) % items.length; paint(); }
    } else if (e.key === 'ArrowUp') {
      e.preventDefault();
      if (list.hidden || !items.length) return;
      active = (active <= 0 ? items.length : active) - 1;
      paint();
    } else if (e.key === 'ArrowRight') {
      // Only hijack → when the caret is already at the end of the text, so
      // it still moves the cursor while the user is editing.
      if (list.hidden || view !== 'parts' || active < 0) return;
      if (input.selectionStart !== input.value.length) return;
      var row = items[active];
      if (!row || row.kind !== 'part') return;
      var source = row.source;
      if (!source || !(Number(source.n_segments) > 0)) return;
      e.preventDefault();
      showSegments(source);
    } else if (e.key === 'ArrowLeft') {
      if (!list.hidden && view === 'segments') { e.preventDefault(); showParts(); }
    } else if (e.key === 'Enter') {
      if (!list.hidden && active >= 0) { e.preventDefault(); choose(active); }
    } else if (e.key === 'Escape') {
      if (!list.hidden && view === 'segments') {
        e.preventDefault(); e.stopPropagation(); showParts(); return;
      }
      if (!list.hidden) { e.preventDefault(); e.stopPropagation(); close(); }
    }
  });

  renderChip();
  return {
    node: node,
    chip: chip,
    formFilter: formSel,
    // The feature_type select owns the `type` filter: when it changes the old
    // suggestions are for the wrong table slice and must be thrown away.
    refresh: function () {
      refresh();
      if (document.activeElement === input && input.value.trim()) search();
    }
  };
}

/**
 * Which rows each Cre-lox cassette encloses.
 *
 * A Cre-lox module is an architecture, not a part: the server wraps the
 * coding modules that follow it and stops at the first module that is not
 * one of them. This recomputes that span from the flat list, so the bracket
 * drawn in the builder is derived from the same rule and never has to be
 * kept in step by hand.
 *
 * Returns [{ head, members: [index…] }] in row order.
 */
function cassetteGroups(mods) {
  var groups = [];
  for (var i = 0; i < mods.length; i++) {
    if (!isCreLox(mods[i].type)) continue;
    var members = [];
    var j = i + 1;
    while (j < mods.length &&
           CASSETTE_PAYLOAD_TYPES.indexOf(mods[j].type) >= 0) {
      members.push(j);
      j++;
    }
    groups.push({ head: i, members: members });
    // A cassette cannot open inside another one's payload: carry on from
    // the first row this cassette did not take.
    i = j - 1;
  }
  return groups;
}

/** The mode selector shown on a Cre-lox row. */
function modeSelect(mod, i, onPicked) {
  var sel = h('select', {
    class: 'mode-select',
    'aria-label': 'Module ' + (i + 1) + ' Cre-lox arrangement',
    onchange: function (e) {
      mod.mode = e.target.value;
      sel.title = creLoxMode(mod.mode).long;
      onPicked();
    }
  }, CRE_LOX_MODES.map(function (m) {
    return h('option', { value: m.key, title: m.long },
      m.key + ' — ' + m.title);
  }));
  sel.value = creLoxMode(mod.mode).key;
  sel.title = creLoxMode(mod.mode).long;
  return sel;
}

/**
 * One builder row, without any of the cassette decoration around it.
 *
 * `focus` is 'type' or 'mode' when this row is being redrawn because that
 * control was just used: changing either one re-groups the whole builder, so
 * the control has to be given its focus back on the new row.
 */
function moduleRow(mod, i, focus) {
  var picker = partPicker(mod);
  var typeSel = h('select', {
    'aria-label': 'Module ' + (i + 1) + ' feature type',
    onchange: function (e) {
      mod.type = e.target.value;
      if (isCreLox(mod.type) && !mod.mode) mod.mode = DEFAULT_CRE_LOX_MODE;
      // The type decides whether this row opens a cassette and whether the
      // rows under it are enclosed, so the whole builder is redrawn; focus
      // is put back on the select that was just used.
      renderModules(i, 'type');
    }
  }, FEATURE_TYPES.map(function (t) {
    return h('option', { value: t, selected: t === mod.type }, t);
  }));
  // A spec loaded with a type this vocabulary no longer lists would
  // otherwise silently become the first option.
  if (FEATURE_TYPES.indexOf(mod.type) < 0 && mod.type) {
    typeSel.appendChild(h('option', { value: mod.type }, mod.type + ' (legacy)'));
  }
  typeSel.value = mod.type;
  if (focus === 'type') requestAnimationFrame(function () { typeSel.focus(); });

  var modeSel = isCreLox(mod.type)
    ? modeSelect(mod, i, function () { renderModules(i, 'mode'); })
    : null;
  if (modeSel && focus === 'mode') {
    requestAnimationFrame(function () { modeSel.focus(); });
  }

  return h('div', { class: 'module-row' + (isCreLox(mod.type) ? ' is-crelox' : '') },
    h('span', { class: 'module-index' }, String(i + 1)),
    typeSel,
    modeSel,
    picker.node,
    picker.formFilter,
    picker.chip,
    h('div', { class: 'row-buttons' },
      h('button', {
        type: 'button', title: 'Move up', disabled: i === 0,
        onclick: function () { swapModules(i, i - 1); }
      }, '↑'),
      h('button', {
        type: 'button', title: 'Move down', disabled: i === state.modules.length - 1,
        onclick: function () { swapModules(i, i + 1); }
      }, '↓'),
      h('button', {
        type: 'button', class: 'danger', title: 'Remove',
        onclick: function () { state.modules.splice(i, 1); renderModules(); }
      }, '×')));
}

/**
 * Draw the builder, bracketing each Cre-lox cassette around what it encloses.
 *
 * The bracket is presentation only: `state.modules` stays the flat, ordered
 * list the server is sent, and the row numbers keep counting straight
 * through the group.
 */
function renderModules(focusIndex, focusWhat) {
  var host = $('modules');
  clear(host);
  if (!state.modules.length) {
    host.appendChild(h('p', { class: 'muted' }, 'No modules yet. Add one below.'));
    return;
  }

  var groups = cassetteGroups(state.modules);
  var headOf = {};
  groups.forEach(function (g) { headOf[g.head] = g; });

  var i = 0;
  while (i < state.modules.length) {
    var group = headOf[i];
    if (!group) {
      host.appendChild(moduleRow(state.modules[i], i,
        i === focusIndex ? (focusWhat || 'type') : null));
      i++;
      continue;
    }

    var mod = state.modules[i];
    var meta = creLoxMode(mod.mode);
    var members = group.members;

    var box = h('div', {
      class: 'cassette' + (members.length ? '' : ' cassette-empty'),
      role: 'group',
      'aria-label': 'Cre-lox cassette, ' + meta.key + ' arrangement, enclosing ' +
        plural(members.length, 'module')
    }, moduleRow(mod, i, i === focusIndex ? (focusWhat || 'type') : null));

    box.appendChild(h('div', { class: 'cassette-label', title: meta.long },
      h('span', { class: 'cassette-mode' }, meta.key),
      h('span', { class: 'cassette-meaning' }, '— ' + meta.short),
      h('span', { class: 'cassette-count muted small' },
        members.length
          ? 'encloses ' + plural(members.length, 'module')
          : 'encloses nothing yet')));

    if (members.length) {
      box.appendChild(h('div', { class: 'cassette-body' },
        members.map(function (idx) {
          return moduleRow(state.modules[idx], idx,
            idx === focusIndex ? (focusWhat || 'type') : null);
        })));
      box.appendChild(h('div', { class: 'cassette-foot muted small' },
        'end of cassette · ' +
        (state.modules[members[members.length - 1] + 1]
          ? 'the rows below stay outside it'
          : 'nothing follows it')));
      i = members[members.length - 1] + 1;
    } else {
      box.appendChild(h('p', { class: 'cassette-warn inline-warn small' },
        'A Cre-lox module has to enclose something. Put a CDS (or a tag, ' +
        'linker, p2A, trafficking motif or IRES) directly after it, or the ' +
        'server drops the cassette.'));
      i++;
    }

    host.appendChild(box);
  }
}

function swapModules(a, b) {
  if (b < 0 || b >= state.modules.length) return;
  var tmp = state.modules[a];
  state.modules[a] = state.modules[b];
  state.modules[b] = tmp;
  renderModules();
}

function usingYaml() {
  var checked = document.querySelector('input[name="spec-mode"]:checked');
  return checked && checked.value === 'yaml';
}

/** Build the `spec` payload: a YAML string, or the object form. */
/* Parents the user chose explicitly. Null means "let Genesis recommend",
   which is the default and must stay reachable once something is pinned. */
var pinnedParents = { recipient: null, donor: null };

/* A pinned parent has to survive whichever way the spec is being written.
   In YAML mode collectSpec hands back the user's own text, so the pins are
   merged into it rather than silently dropped. */
function withPins(spec) {
  if (typeof spec !== 'string') return spec;
  var lines = spec.split('\n').filter(function (line) {
    return !/^\s*(recipient|donor)_plasmid_id\s*:/.test(line);
  });
  var extra = [];
  if (pinnedParents.recipient) {
    extra.push('recipient_plasmid_id: ' + pinnedParents.recipient);
  }
  if (pinnedParents.donor) {
    extra.push('donor_plasmid_id: ' + pinnedParents.donor);
  }
  if (!extra.length) return lines.join('\n');
  // Top-level keys must come before the modules list in this YAML dialect.
  var at = lines.findIndex(function (line) { return /^modules\s*:/.test(line); });
  if (at < 0) at = lines.length;
  return lines.slice(0, at).concat(extra, lines.slice(at)).join('\n');
}

function collectSpec() {
  if (usingYaml()) return withPins($('spec-yaml').value);
  var backbone = $('spec-backbone').value.trim();
  var spec = {
    name: $('spec-name').value.trim() || 'untitled',
    owner: $('spec-owner').value.trim(),
    circular: $('spec-circular').checked,
    // "" means Auto: the server infers the family from the backbone.
    vector: $('spec-vector').value,
    modules: state.modules.map(function (m) {
      var out = { type: m.type, name: (m.name || '').trim() };
      // Only send a pin when the row actually carries one; an absent part_id
      // leaves the server free to pick a candidate as it always did.
      if (m.part_id) out.part_id = m.part_id;
      // `mode` means nothing on any other type, and the server falls back to
      // reading the module's *name* as a mode when it is absent, so a
      // Cre-lox row always states it.
      if (isCreLox(m.type)) out.mode = creLoxMode(m.mode).key;
      return out;
    })
  };
  // A bare integer is treated as a plasmid id, anything else as a name.
  if (/^\d+$/.test(backbone)) spec.backbone_plasmid_id = parseInt(backbone, 10);
  else spec.backbone = backbone;
  if (pinnedParents.recipient) spec.recipient_plasmid_id = pinnedParents.recipient;
  if (pinnedParents.donor) spec.donor_plasmid_id = pinnedParents.donor;
  return spec;
}

/** Fill the builder from a stored spec, so a loaded design can be edited. */
function applySpecToBuilder(spec) {
  if (!spec || typeof spec !== 'object') return false;
  $('spec-name').value = spec.name || '';
  $('spec-backbone').value = spec.backbone_plasmid_id
    ? String(spec.backbone_plasmid_id) : (spec.backbone || '');
  $('spec-owner').value = spec.owner || '';
  $('spec-circular').checked = spec.circular !== false;
  // Designs saved before the vector field existed have no `vector` key; those
  // fall back to Auto, which is what they were planned with.
  var vec = $('spec-vector');
  var wanted = spec.vector || '';
  vec.value = Array.prototype.some.call(vec.options, function (o) {
    return o.value === wanted;
  }) ? wanted : '';
  state.modules = (spec.modules || []).map(function (m) {
    var type = m.type || 'unknown';
    return {
      type: type,
      name: m.name || '',
      part_id: m.part_id || null,
      pin_note: '',
      form: 'all',
      // A Cre-lox row saved without a mode is the DIO arrangement, which is
      // what the server assumes for it too.
      mode: isCreLox(type) ? creLoxMode(m.mode).key : (m.mode || '')
    };
  });
  if (!state.modules.length) state.modules = [newModule('promoter')];
  // The builder is the only view that can show a parsed spec.
  var builderRadio = document.querySelector('input[name="spec-mode"][value="builder"]');
  if (builderRadio && !builderRadio.checked) {
    builderRadio.checked = true;
    $('builder-pane').hidden = false;
    $('yaml-pane').hidden = true;
  }
  renderModules();
  return true;
}

/**
 * Enable/disable the three save buttons and describe what they would do.
 *
 * Overwrite needs a design_id — the server replaces a stored version by id —
 * so it stays disabled, with a tooltip saying why, until the plan on screen
 * is tied to a saved design (loaded from the Designs tab, or just saved).
 */
function refreshSaveButtons() {
  var version = $('save-version');
  var overwrite = $('save-overwrite');
  var saveAs = $('save-as');
  var link = $('save-link');

  var hasPlan = !!(state.plan && state.spec);
  version.disabled = !hasPlan;
  saveAs.disabled = !hasPlan;

  if (!hasPlan) {
    var why = state.plan
      ? 'This plan was loaded without its specification, so there is nothing to re-save. Generate a plan first.'
      : 'Generate a plan first.';
    version.title = why;
    saveAs.title = why;
  } else {
    version.title = 'Save a new version under this construct name';
    saveAs.title = 'Store this plan under a different construct name, leaving the original alone';
  }

  var saved = state.saved;
  overwrite.disabled = !(hasPlan && saved && saved.design_id);
  overwrite.title = overwrite.disabled
    ? (hasPlan
        ? 'Overwrite replaces a stored version in place, so it needs a saved ' +
          'design to replace. This plan has not been saved yet and was not ' +
          'loaded from the Designs tab — use Save, which adds version 1.'
        : 'Generate a plan first.')
    : 'Replace version ' + saved.version + ' of “' + saved.name +
      '” (design #' + saved.design_id + ') in place, without adding a version';

  clear(link);
  if (hasPlan && saved && saved.design_id) {
    link.appendChild(document.createTextNode(
      'Linked to design #' + saved.design_id + ' — version ' + saved.version +
      ' of “' + saved.name + '”. Overwrite replaces that version; Save adds a new one.'));
  }
}

/* ------------------------------------------------------------------ */
/* Tab 3: Design - clearing the specification                          */
/* ------------------------------------------------------------------ */

/** The construct fields, as they are on a freshly loaded page. */
function fieldsArePristine() {
  var name = $('spec-name').value.trim().toLowerCase();
  return (name === '' || name === 'untitled') &&
         !$('spec-backbone').value.trim() &&
         !$('spec-owner').value.trim() &&
         !$('spec-vector').value &&
         $('spec-circular').checked;
}

function modulesArePristine() {
  if (state.modules.length !== 1) return false;
  var m = state.modules[0];
  return m.type === 'promoter' && !(m.name || '').trim() && !m.part_id;
}

/**
 * What a Clear would throw away, in words, or [] when there is nothing.
 *
 * Clearing is one click and cannot be undone, so the confirmation names the
 * work by what it is — a plan, pinned parts, typed YAML — rather than asking
 * "are you sure?" about an unspecified something.
 */
function unsavedWork() {
  var bits = [];

  if (state.plan) {
    bits.push(state.saved && state.saved.design_id
      ? 'the plan on screen (saved as version ' + state.saved.version +
        ' of “' + state.saved.name + '”, so it can be re-opened from the ' +
        'Designs tab)'
      : 'the generated plan on screen, which has not been saved');
  }

  if (!modulesArePristine()) {
    var pinned = state.modules.filter(function (m) { return !!m.part_id; }).length;
    bits.push(plural(state.modules.length, 'module row') +
      (pinned ? ' (' + pinned + ' pinned to an exact inventory part)' : ''));
  }

  if (!fieldsArePristine()) bits.push('the construct details that were typed in');
  if ($('spec-yaml').value.trim()) bits.push('the YAML specification');

  return bits;
}

/** Reset the whole Design tab to the state the page opens in. */
function clearSpec() {
  $('spec-name').value = 'untitled';
  $('spec-backbone').value = '';
  $('spec-owner').value = '';
  $('spec-vector').value = DEFAULT_VECTOR;
  $('spec-circular').checked = true;
  $('spec-yaml').value = '';

  // A fresh row: no name, and no `part_id`, so nothing stays pinned.
  state.modules = [newModule('promoter')];
  // Nor do the chosen parents survive: a cleared specification has no
  // recipient or donor to be pinned to.
  pinnedParents.recipient = null;
  pinnedParents.donor = null;
  state.plan = null;
  state.spec = null;
  state.routeId = null;
  state.saved = null;
  state.appliedSuggestion = null;

  clear($('plan-out'));
  renderModules();
  refreshSaveButtons();   // also clears the saved-design link
  setStatus($('design-status'), 'specification cleared');
}

function clearSpecWithConfirm() {
  var losing = unsavedWork();
  if (!losing.length) { clearSpec(); return; }

  modal({
    title: 'Clear the construct specification?',
    confirmLabel: 'Clear',
    danger: true,
    body: [
      h('p', null, 'This resets the Design tab to an empty specification. ' +
        'It discards:'),
      h('ul', { class: 'modal-list' }, losing.map(function (b) {
        return h('li', null, b);
      })),
      h('p', { class: 'muted small' },
        'Saved designs are not touched — only what is on screen.')
    ],
    onConfirm: function (api_) { api_.close(); clearSpec(); }
  });
}

async function generatePlan() {
  var status = $('design-status');
  var out = $('plan-out');
  setStatus(status, 'running stages 4–10…', 'busy');
  clear(out);
  out.appendChild(h('div', { class: 'card' },
    h('p', { class: 'muted' }, 'Resolving modules, validating and ranking routes…')));
  try {
    var spec = collectSpec();
    var plan = await postJSON('/api/design', { spec: spec });
    state.spec = spec;
    state.plan = plan;
    state.appliedSuggestion = null;   // a fresh plan may suggest a new name
    state.routeId = (plan.routes && plan.routes.length) ? plan.routes[0].id : null;
    // `state.saved` is deliberately kept: regenerating after loading design #5
    // is exactly the case Overwrite exists for.
    refreshSaveButtons();
    renderPlan();
    setStatus(status, 'plan ready', 'ok');
  } catch (err) {
    clear(out);
    out.appendChild(errorBox(err.message));
    setStatus(status, 'failed', 'error');
  }
}

/** Name of the construct this plan would be saved under. */
function specName(spec) {
  if (typeof spec === 'string') {
    var m = /^\s*name:\s*(.+)$/m.exec(spec);
    return m ? m[1].trim().replace(/^["']|["']$/g, '') : 'untitled';
  }
  return (spec && spec.name) || 'untitled';
}

/** Return a copy of the spec with a different construct name. */
function specWithName(spec, name) {
  if (typeof spec === 'string') {
    return /^\s*name:/m.test(spec)
      ? spec.replace(/^\s*name:.*$/m, 'name: ' + name)
      : 'name: ' + name + '\n' + spec;
  }
  var copy = {};
  Object.keys(spec).forEach(function (k) { copy[k] = spec[k]; });
  copy.name = name;
  return copy;
}

/**
 * POST /api/save. `mode` is version | overwrite | save_as.
 *
 * Throws on failure so callers that need to react (Save as… re-asking after a
 * 409) can, rather than only seeing a status line.
 */
async function postSave(spec, mode) {
  var payload = { spec: spec, plan: state.plan, mode: mode };
  if (mode === 'overwrite') payload.design_id = state.saved.design_id;
  var res = await postJSON('/api/save', payload);
  state.spec = spec;
  state.saved = {
    design_id: res.design_id,
    construct_id: res.construct_id,
    version: res.version,
    name: specName(spec)
  };
  refreshSaveButtons();
  loadStats();   // the design count in the header is now stale
  if (LOADED.designs) renderDesigns();
  return res;
}

function savedMessage(res, name) {
  if (res.mode === 'overwrite') {
    return 'overwrote version ' + res.version + ' of “' + name +
           '” (design #' + res.design_id + ')';
  }
  return 'saved version ' + res.version + ' of “' + name +
         '” (design #' + res.design_id + ')';
}

async function saveDesign(mode) {
  var status = $('design-status');
  if (!state.plan) { setStatus(status, 'generate a plan first', 'error'); return; }
  if (!state.spec) {
    setStatus(status,
      'this plan was loaded without its specification; regenerate it to save',
      'error');
    return;
  }
  if (mode === 'overwrite' && !(state.saved && state.saved.design_id)) {
    setStatus(status, 'there is no saved design to overwrite', 'error');
    return;
  }
  setStatus(status, mode === 'overwrite' ? 'overwriting…' : 'saving…', 'busy');
  try {
    var res = await postSave(state.spec, mode);
    setStatus(status, savedMessage(res, specName(state.spec)), 'ok');
  } catch (err) {
    setStatus(status, 'save failed: ' + err.message, 'error');
  }
}

/**
 * Save as…: ask for a name, and keep asking if the server rejects it.
 *
 * The server answers 409 when the name is taken, which is the whole reason
 * this dialog stays open on failure: the message names the clash and the
 * typed name is still there to edit.
 */
function saveAsDesign() {
  var status = $('design-status');
  if (!state.plan || !state.spec) {
    setStatus(status, 'generate a plan first', 'error');
    return;
  }
  var current = specName(state.spec);
  var field = h('input', {
    type: 'text', class: 'modal-input', value: current + ' copy',
    spellcheck: 'false', autocomplete: 'off', 'aria-label': 'New construct name'
  });

  var dlg = modal({
    title: 'Save as a new construct',
    confirmLabel: 'Save as',
    focus: field,
    body: [
      h('p', null,
        'Store this plan under a new construct name as version 1. ' +
        '“' + current + '” is left untouched.'),
      field,
      h('p', { class: 'muted small' },
        'The name has to be free — the server refuses to merge into an ' +
        'existing construct’s history.')
    ],
    onConfirm: function (api_) {
      var name = field.value.trim();
      if (!name) { api_.setError('Give the new construct a name.'); field.focus(); return; }
      api_.setError('');
      setStatus(status, 'saving as “' + name + '”…', 'busy');
      return postSave(specWithName(state.spec, name), 'save_as').then(function (res) {
        // Reflect the new name in the builder, since that is what is now on
        // screen and what a following Save/Overwrite would act on.
        if (!usingYaml()) $('spec-name').value = name;
        api_.close();
        setStatus(status, savedMessage(res, name), 'ok');
      }, function (err) {
        // 409 and anything else: stay open, say why, let them pick again.
        api_.setError(err.message);
        setStatus(status, 'save as failed', 'error');
        field.focus();
        field.select();
      });
    }
  });

  field.addEventListener('keydown', function (e) {
    if (e.key === 'Enter') { e.preventDefault(); dlg.confirmButton.click(); }
  });
}

/* ------------------------------------------------------------------ */
/* Tab 3: Design - plan rendering                                      */
/* ------------------------------------------------------------------ */

/**
 * Offer the name the planner derived from the resolved architecture.
 *
 * The user's own name is never overwritten silently: the only exception is a
 * field that is empty or still carries the placeholder, which nobody chose.
 */
function suggestedName(construct) {
  var suggestion = String(construct.suggested_name || '').trim();
  if (!suggestion) return null;

  var field = $('spec-name');
  var current = field.value.trim();

  function applied() {
    return h('p', { class: 'suggest-name' },
      h('span', { class: 'muted small' }, 'Construct name set from the architecture:'),
      h('span', { class: 'mono small suggest-value' }, suggestion));
  }

  if (current === suggestion) {
    // renderPlan also runs when a route is clicked, so the confirmation is
    // kept for as long as this plan is on screen rather than flickering away.
    return state.appliedSuggestion === suggestion ? applied() : null;
  }

  if (!current || current.toLowerCase() === 'untitled') {
    field.value = suggestion;
    state.appliedSuggestion = suggestion;
    return applied();
  }

  return h('p', { class: 'suggest-name' },
    h('span', { class: 'muted small' }, 'Suggested:'),
    h('span', { class: 'mono small suggest-value' }, suggestion),
    h('button', {
      type: 'button', class: 'linkish',
      title: 'Put this name in the Construct name field',
      onclick: function () {
        field.value = suggestion;
        state.appliedSuggestion = suggestion;
        renderPlan();
      }
    }, 'use'));
}

/* --------------------------------------------------------------------------
 * Cloning scheme: the two parent plasmids, drawn as plasmids.
 *
 * A swap is an edit to a molecule that already exists, so the thing worth
 * looking at is the recipient with the region to cut marked on it, and the
 * donor with the region to take. Both are circles for the same reason the
 * product is.
 * ------------------------------------------------------------------------ */

/* A list of candidate parents, with the reason each was offered, and a way
 * to choose one. Which plasmid to start from is a lab judgement as much as
 * a computation, so the recommendation is a default rather than a verdict
 * and "let Genesis recommend" stays reachable once something is pinned. */
function parentPicker(role, options, chosenId) {
  options = options || [];
  if (options.length < 2 && !pinnedParents[role]) return null;

  var box = h('div', { class: 'parent-picker' });
  box.appendChild(h('div', { class: 'picker-head muted small' },
    options.length + ' plasmid' + (options.length === 1 ? '' : 's') +
    ' could serve as the ' + role + ' \u2014 ranked, best first'));

  var list = h('div', { class: 'picker-list', role: 'radiogroup',
                        'aria-label': role + ' options' });

  function row(option) {
    var id = option ? option.id : null;
    var isPinned = pinnedParents[role] === id;
    var isActive = id === null
      ? !pinnedParents[role]
      : (isPinned || (!pinnedParents[role] && id === chosenId));

    var item = h('button', {
      type: 'button', role: 'radio',
      'aria-checked': isActive ? 'true' : 'false',
      class: 'picker-row' + (isActive ? ' active' : '') +
             (isPinned ? ' pinned' : ''),
      onclick: function () {
        pinnedParents[role] = id;
        generatePlan();
      }
    });

    if (id === null) {
      item.appendChild(h('span', { class: 'picker-name' },
        'Let Genesis recommend'));
      item.appendChild(h('span', { class: 'picker-why muted small' },
        'clears the pin and re-ranks'));
      return item;
    }

    item.appendChild(h('span', { class: 'picker-name' }, option.name));
    if (isPinned) item.appendChild(h('span', { class: 'pill pill-accent' }, 'pinned'));
    else if (isActive) item.appendChild(h('span', { class: 'pill' }, 'recommended'));
    item.appendChild(h('span', { class: 'picker-why muted small' },
      (option.reasons || []).join(' \u00b7 ')));
    if (option.owner) {
      item.appendChild(h('span', { class: 'picker-owner muted small' },
        option.owner));
    }
    return item;
  }

  options.forEach(function (option) { list.appendChild(row(option)); });
  if (pinnedParents[role]) list.appendChild(row(null));
  box.appendChild(list);
  return box;
}

function parentCard(role, parent, kind, options) {
  if (!parent) {
    if (!(options || []).length) return null;
    var empty = h('div', { class: 'parent-card' },
      h('div', { class: 'parent-head' },
        h('span', { class: 'parent-role ' + role }, role.toUpperCase()),
        h('span', { class: 'muted' }, 'none chosen')));
    var pickEmpty = parentPicker(role, options, null);
    if (pickEmpty) empty.appendChild(pickEmpty);
    return empty;
  }
  var head = h('div', { class: 'parent-head' },
    h('span', { class: 'parent-role ' + role }, role.toUpperCase()),
    h('strong', null, parent.name),
    h('span', { class: 'muted small' },
      bp(parent.length) + (parent.owner ? ' · ' + parent.owner : '')));
  var card = h('div', { class: 'parent-card' }, head);
  card.appendChild(circularMap(parent.length, (parent.features || []).map(
    function (f) {
      return { type: f.feature_type, name: f.name, start: f.start,
               end: f.end, strand: f.strand };
    }), {
      label: parent.name, highlight: parent.highlight, highlightKind: kind
    }));
  var picker = parentPicker(role, options, parent.id);
  if (picker) card.appendChild(picker);
  return card;
}

function schemeCard(plan) {
  var scheme = plan.scheme || {};
  if (!scheme.recipient) return null;

  var card = h('div', { class: 'card scheme-card' },
    h('div', { class: 'card-head' },
      h('h2', null, 'Cloning scheme'),
      h('span', { class: 'muted small' },
        'two parents: cut a region out of the recipient, put a region from ' +
        'the donor in its place')));

  var grid = h('div', { class: 'parent-grid' });
  var rec = parentCard('recipient', scheme.recipient, 'cut',
                       scheme.recipient_options);
  var don = parentCard('donor', scheme.donor, 'take', scheme.donor_options);
  if (rec) grid.appendChild(rec);
  if (don) grid.appendChild(don);
  card.appendChild(grid);

  var anchors = [];
  if (scheme.left_anchor) {
    anchors.push(scheme.left_anchor.name + ' (ends ' +
                 bp(scheme.left_anchor.end) + ')');
  }
  if (scheme.right_anchor) {
    anchors.push(scheme.right_anchor.name + ' (starts ' +
                 bp(scheme.right_anchor.start) + ')');
  }
  if (anchors.length) {
    card.appendChild(h('p', { class: 'muted small' },
      'The cut is bounded by ' + anchors.join(' and ') +
      (anchors.length === 1 ? ', which is kept' : ', which are kept') +
      ' and is what the primers anneal to or the enzymes cut between.'));
  }

  return card;
}

function renderPlan() {
  var out = $('plan-out');
  clear(out);
  var plan = state.plan;
  if (!plan) return;

  var c = plan.construct || {};

  /* --- construct summary + map ------------------------------------ */
  var summary = h('div', { class: 'card' },
    h('div', { class: 'card-head' },
      h('h2', null, c.name || '(unnamed construct)'),
      h('span', { class: 'muted small' },
        bp(c.length) + ' · ' + (c.circular ? 'circular' : 'linear') +
        (plan.backbone ? ' · backbone: ' + plan.backbone.name +
          ' (' + bp(plan.backbone.length) + ')' : ' · no backbone'))),
    suggestedName(c),
    creLoxSummary(plan));

  if (plan.modules && plan.modules.length) {
    summary.appendChild(h('h4', null, 'Modules'));
    summary.appendChild(table([
      { label: 'Type', cell: function (m) { return typeTag(m.type); } },
      { label: 'Requested', key: 'requested' },
      { label: 'Resolved part', key: 'part_name' },
      { label: 'Length', numeric: true, cell: function (m) { return bp(m.length); } },
      { label: 'Sources', numeric: true, cell: function (m) { return (m.sources || []).length; } },
      { label: 'Notes', cell: function (m) { return warnCell(m.warnings); } }
    ], plan.modules));
  }

  if (plan.fragments && plan.fragments.length) {
    summary.appendChild(h('h4', null, 'Planned fragments'));
    summary.appendChild(table([
      { label: 'Fragment', key: 'name' },
      { label: 'Length', numeric: true, cell: function (f) { return bp(f.length); } },
      { label: 'Source', key: 'source' },
      { label: 'Synthesis', cell: function (f) {
          return f.needs_synthesis
            ? h('span', { class: 'pill pill-warn' }, 'required') : 'no'; } }
    ], plan.fragments));
  }
  out.appendChild(summary);

  var scheme = schemeCard(plan);
  if (scheme) out.appendChild(scheme);

  /* --- the designed molecule itself -------------------------------- */
  out.appendChild(renderConstruct(c));

  /* --- read-back verification -------------------------------------- */
  if (plan.verification) out.appendChild(renderVerification(plan.verification));

  /* --- vector family ----------------------------------------------- */
  if (plan.vector) out.appendChild(renderVector(plan.vector, plan));

  /* --- validation -------------------------------------------------- */
  out.appendChild(renderValidation(plan.validation || {}));

  /* --- routes ------------------------------------------------------ */
  out.appendChild(renderRoutes(plan.routes || []));

  /* --- selected route detail --------------------------------------- */
  var route = selectedRoute();
  if (route) out.appendChild(renderRouteDetail(route));

  /* --- export ------------------------------------------------------ */
  out.appendChild(renderExports());
}

/* ------------------------------------------------------------------ */
/* Plan: the Cre-lox arrangement that was asked for                    */
/* ------------------------------------------------------------------ */

/**
 * One line naming each cassette in the specification and what it does.
 *
 * `spec_order` is the order as typed, before the cassette was expanded, so
 * it is the only place in the plan that still carries the mode.
 */
function creLoxSummary(plan) {
  var order = (plan.construct && plan.construct.spec_order) || [];
  var cassettes = order.filter(function (m) { return isCreLox(m.type); });
  if (!cassettes.length) return null;

  return h('p', { class: 'crelox-summary' }, cassettes.map(function (m) {
    var meta = creLoxMode(m.mode);
    return h('span', { class: 'crelox-chip', title: meta.long },
      h('span', { class: 'crelox-key' }, 'Cre-lox ' + meta.key),
      h('span', { class: 'muted small' }, meta.short));
  }));
}

/* ------------------------------------------------------------------ */
/* Plan: the full designed construct                                   */
/* ------------------------------------------------------------------ */

var SEQ_LINE = 60;        // bases per rendered line
var SEQ_GROUP = 10;       // bases per visual group within a line

/** Fraction of G + C, counted without allocating a copy of the sequence. */
function gcFraction(seq) {
  var gc = 0;
  for (var i = 0; i < seq.length; i++) {
    var ch = seq.charCodeAt(i);
    // G g C c
    if (ch === 71 || ch === 103 || ch === 67 || ch === 99) gc++;
  }
  return seq.length ? gc / seq.length : 0;
}

/**
 * The numbered sequence listing.
 *
 * Built as one string and handed to a single text node: a 15 kb construct is
 * 250 lines, and one text node costs a fraction of what 250 elements and
 * their layout would. Called lazily, the first time the block is opened.
 */
function sequenceText(seq) {
  var width = String(seq.length).length;
  var lines = [];
  for (var i = 0; i < seq.length; i += SEQ_LINE) {
    var chunk = seq.slice(i, i + SEQ_LINE);
    var groups = [];
    for (var j = 0; j < chunk.length; j += SEQ_GROUP) {
      groups.push(chunk.slice(j, j + SEQ_GROUP));
    }
    lines.push(String(i + 1).padStart(width, ' ') + '  ' + groups.join(' '));
  }
  return lines.join('\n');
}

/** Put `text` on the clipboard, with a fallback for browsers without the API. */
async function copyToClipboard(text) {
  if (navigator.clipboard && navigator.clipboard.writeText) {
    try {
      await navigator.clipboard.writeText(text);
      return;
    } catch (err) { /* fall through to the textarea below */ }
  }
  var ta = h('textarea', {
    style: 'position:fixed;top:-1000px;left:-1000px;opacity:0', 'aria-hidden': 'true'
  });
  ta.value = text;
  document.body.appendChild(ta);
  ta.select();
  var ok = false;
  try { ok = document.execCommand('copy'); } catch (err) { ok = false; }
  ta.remove();
  if (!ok) throw new Error('this browser would not let the page copy');
}

/** A button that copies whatever `get()` returns, and says that it did. */
function copyButton(label, get, status, what) {
  return h('button', {
    type: 'button', class: 'secondary',
    onclick: async function (e) {
      var btn = e.currentTarget;
      btn.disabled = true;
      try {
        var text = get();
        await copyToClipboard(text);
        setStatus(status, what + ' copied (' +
          Number(text.length).toLocaleString() + ' characters)', 'ok');
      } catch (err) {
        setStatus(status, 'copy failed: ' + err.message, 'error');
      } finally {
        btn.disabled = false;
      }
    }
  }, label);
}

/** FASTA, so the copied sequence can be pasted straight into another tool. */
function fastaText(name, seq) {
  return '>' + (String(name || 'construct').replace(/\s+/g, '_')) + ' ' +
    seq.length + ' bp\n' + sequenceWrapped(seq) + '\n';
}

function sequenceWrapped(seq) {
  var lines = [];
  for (var i = 0; i < seq.length; i += SEQ_LINE) lines.push(seq.slice(i, i + SEQ_LINE));
  return lines.join('\n');
}

/**
 * Everything that was designed: the feature map, every element in it, and the
 * sequence itself.
 *
 * The element table and the sequence listing are both collapsed, and the
 * sequence is only built when it is first opened, so a 15 kb construct costs
 * nothing until somebody asks to read it.
 */
function renderConstruct(c) {
  var seq = typeof c.sequence === 'string' ? c.sequence : '';
  var layout = (c.layout || []).filter(function (f) {
    return f && typeof f.start === 'number';
  });
  var status = h('span', { class: 'status' });

  var card = h('div', { class: 'card construct-full' },
    h('div', { class: 'card-head' },
      h('h2', null, 'Full construct'),
      h('span', { class: 'muted small' },
        seq ? 'every base and every element of the designed molecule'
            : 'this plan carries no sequence')));

  // A circular plasmid is shown as one by default; the linear track is
  // still there for reading coordinates off, behind a toggle.
  var mapHost = h('div', { class: 'map-host' });
  function drawMap(mode) {
    clear(mapHost);
    mapHost.appendChild(mode === 'circular'
      ? circularMap(c.length, layout,
          { label: c.name, markers: c.codon_markers })
      : featureMap(c.length, layout,
          { label: c.name, circular: !!c.circular,
            markers: c.codon_markers }));
  }
  if (c.circular) {
    var mode = 'circular';
    var toggle = h('div', { class: 'map-toggle' });
    ['circular', 'linear'].forEach(function (option) {
      var b = h('button', {
        type: 'button',
        class: option === mode ? 'active' : '',
        onclick: function () {
          mode = option;
          [].forEach.call(toggle.children, function (child) {
            child.className = child.textContent === mode ? 'active' : '';
          });
          drawMap(mode);
        }
      }, option);
      toggle.appendChild(b);
    });
    card.appendChild(toggle);
  }
  card.appendChild(mapHost);
  drawMap(c.circular ? 'circular' : 'linear');

  if (!seq) {
    card.appendChild(h('p', { class: 'muted' },
      'The sequence is not part of this plan, so only the map can be shown.'));
    return card;
  }

  card.appendChild(h('dl', { class: 'vector-facts construct-facts' },
    statPair('length', bp(seq.length)),
    statPair('topology', c.circular ? 'circular' : 'linear'),
    statPair('GC', pct(gcFraction(seq))),
    statPair('elements', String(layout.length))));

  card.appendChild(h('div', { class: 'button-row' },
    copyButton('Copy sequence', function () { return seq; }, status, 'sequence'),
    copyButton('Copy FASTA', function () { return fastaText(c.name, seq); },
      status, 'FASTA'),
    status));

  /* --- every element ------------------------------------------------ */
  card.appendChild(h('details', { class: 'construct-elements' },
    h('summary', null, 'Every element (' + layout.length + ')'),
    table([
      { label: 'Element', key: 'name' },
      { label: 'Type', cell: function (f) { return typeTag(typeOf(f)); } },
      { label: 'Start', numeric: true, cell: function (f) { return f.start + 1; } },
      { label: 'End', numeric: true, key: 'end' },
      { label: 'Length', numeric: true, cell: function (f) {
          return bp(Number(f.end) - Number(f.start)); } },
      { label: 'Strand', numeric: true, cell: function (f) {
          return strandSign(f.strand); } }
    ], layout)));

  /* --- the sequence itself ------------------------------------------ */
  var pre = h('pre', { class: 'seq-block mono' });
  var built = false;
  var details = h('details', {
    class: 'construct-sequence',
    ontoggle: function () {
      if (built || !details.open) return;
      built = true;
      clear(pre);
      pre.appendChild(document.createTextNode(sequenceText(seq)));
    }
  },
    h('summary', null, 'Sequence — ' + bp(seq.length) + ', ' + SEQ_LINE +
      ' bases per line'),
    h('p', { class: 'muted small' },
      'Positions are 1-based and count from the start of the molecule' +
      (c.circular ? ', which for a circular plasmid is where the backbone ' +
        'file began' : '') + '.'),
    pre);
  card.appendChild(details);

  return card;
}

/* ------------------------------------------------------------------ */
/* Plan: read-back verification                                        */
/* ------------------------------------------------------------------ */

/**
 * The server's independent read-back of the finished sequence.
 *
 * `passed` is the headline and is shown either way, because a design that
 * quietly failed its own check is the one thing nobody should have to scroll
 * for. Note that `problems` can be non-empty on a pass: a part that occurs
 * more than once is reported but is not by itself a failure.
 */
function renderVerification(v) {
  var passed = !!v.passed;
  var problems = v.problems || [];
  var checks = v.checks || [];
  var observed = v.observed || [];
  var rec = v.recombined;

  var card = h('div', {
    class: 'card verification ' + (passed ? 'verify-pass' : 'verify-fail')
  },
    h('div', { class: 'card-head' },
      h('h2', null, 'Verification'),
      h('span', { class: 'muted small' },
        'the assembled sequence read back on its own, without the layout ' +
        'that built it')));

  card.appendChild(h('div', { class: 'verdict' },
    h('span', { class: 'verdict-pill ' + (passed ? 'pill-ok' : 'pill-bad') },
      passed ? 'PASSED' : 'FAILED'),
    h('span', { class: 'verdict-text' },
      v.summary || (passed ? 'the construct reads back as specified'
                           : 'the construct does not read back as specified'))));

  card.appendChild(h('dl', { class: 'vector-facts' },
    statPair('parts found',
      checks.filter(function (c) { return c.found; }).length + ' of ' + checks.length),
    statPair('order', v.order_ok
      ? 'as specified'
      : 'NOT as specified — the forward-strand parts occur in a different order'),
    statPair('re-annotated', plural(observed.length, 'feature') +
      ' found by scanning the finished molecule')));

  if (problems.length) {
    card.appendChild(h('h4', { class: 'lvl lvl-' + (passed ? 'warning' : 'error') },
      plural(problems.length, 'problem') +
      (passed ? ' reported (none of them fails the check)' : '')));
    card.appendChild(h('ul', {
      class: 'issues issues-' + (passed ? 'warning' : 'error')
    }, problems.map(function (p) {
      return h('li', null, h('span', { class: 'msg' }, p));
    })));
  }

  /* --- after recombination ------------------------------------------ */
  if (rec) {
    var payload = rec.payload_modules || [];
    var forward = rec.now_forward || [];
    var switched = !!rec.all_switched;
    var still = payload.filter(function (p) { return forward.indexOf(p) < 0; });

    // The server only counts a module as payload when it was stored
    // *inverted*, so a Cre-OFF (DO) cassette legitimately has none: there is
    // nothing to turn on, and calling that a failure would be wrong.
    var kind = payload.length ? (switched ? 'is-ok' : 'is-bad') : 'is-neutral';
    var verdict = payload.length
      ? (switched ? 'payload switched' : 'payload did not switch')
      : 'nothing stored inverted';
    var sentence;
    if (!payload.length) {
      sentence = 'Simulating Cre gives a ' + bp(rec.length) + ' molecule. ' +
        'No module was stored inverted, so there is nothing for Cre to turn ' +
        'on — which is what a Cre-OFF (DO) arrangement looks like.';
    } else if (switched) {
      sentence = 'Simulating Cre gives a ' + bp(rec.length) + ' molecule in ' +
        'which every inverted payload module now reads on the forward ' +
        'strand, so the cassette switches as designed.';
    } else {
      sentence = 'Simulating Cre gives a ' + bp(rec.length) + ' molecule in ' +
        'which ' + still.join(', ') + ' still cannot be read on the forward ' +
        'strand.';
    }

    card.appendChild(h('div', { class: 'recombined ' + kind },
      h('h4', null, 'After Cre recombination'),
      h('p', { class: 'recombined-verdict' },
        h('span', {
          class: payload.length
            ? 'pill pill-' + (switched ? 'ok' : 'bad') : 'pill'
        }, verdict),
        sentence),
      h('dl', { class: 'vector-facts' },
        statPair('stored inverted', payload.join(' + ') || 'nothing'),
        statPair('forward after Cre', forward.join(' + ') || 'none'),
        statPair('recombined length', bp(rec.length)))));
  }

  /* --- per-part read-back ------------------------------------------- */
  card.appendChild(h('h4', null, 'Parts read back (' + checks.length + ')'));
  card.appendChild(table([
    { label: 'Part', key: 'name' },
    { label: 'Type', cell: function (r) { return typeTag(r.feature_type); } },
    { label: 'Expected', numeric: true, cell: function (r) { return bp(r.expected_bp); } },
    { label: 'Found', cell: function (r) {
        return h('span', { class: 'pill pill-' + (r.found ? 'ok' : 'bad') },
          r.found ? 'found' : 'missing'); } },
    // Two lox sites of each type, or two ITRs, are what makes the construct
    // work: the server only calls a repeat a problem when it is one, and
    // says so in `detail`, so only those are flagged here.
    { label: 'Copies', numeric: true, cell: function (r) {
        return h('span', {
          class: (r.occurrences > 1 && r.detail) ? 'inline-warn' : null,
          title: r.occurrences > 1 && !r.detail
            ? 'more than one copy is expected for this element' : null
        }, String(r.occurrences)); } },
    { label: 'Strand', numeric: true, cell: function (r) {
        return r.found ? strandSign(r.strand) : ''; } },
    { label: 'Position', numeric: true, cell: function (r) {
        return (r.start === null || r.start === undefined)
          ? '' : Number(r.start + 1).toLocaleString(); } },
    { label: 'Protein', cell: function (r) {
        if (r.protein_ok === null || r.protein_ok === undefined) return '';
        return h('span', { class: 'pill pill-' + (r.protein_ok ? 'ok' : 'bad') },
          r.protein_ok ? 'no internal stop' : 'internal stop'); } },
    { label: 'Detail', cell: function (r) {
        return r.detail ? h('span', { class: 'inline-warn' }, r.detail) : ''; } }
  ], checks, {
    rowClass: function (r) { return r.found ? null : 'row-bad'; }
  }));

  /* --- independent re-annotation ------------------------------------ */
  card.appendChild(h('details', { class: 'observed' },
    h('summary', null, 'What a fresh scan of the finished molecule finds (' +
      observed.length + ')'),
    h('p', { class: 'muted small' },
      'The same motif scanner that annotates imported files, run over the ' +
      'assembled sequence with no knowledge of how it was built.'),
    table([
      { label: 'Feature', key: 'name' },
      { label: 'Type', cell: function (r) { return typeTag(r.feature_type); } },
      { label: 'Start', numeric: true, cell: function (r) { return r.start + 1; } },
      { label: 'End', numeric: true, key: 'end' },
      { label: 'Strand', numeric: true, cell: function (r) { return strandSign(r.strand); } },
      { label: 'Level', key: 'level' }
    ], observed)));

  return card;
}

/* ------------------------------------------------------------------ */
/* Plan: vector family                                                 */
/* ------------------------------------------------------------------ */

var VECTOR_LABELS = {
  AAV: 'AAV', lentiviral: 'lentiviral', plasmid: 'plain plasmid'
};

function vectorLabel(key) { return VECTOR_LABELS[key] || key || 'unknown'; }

/**
 * One line of detection evidence.
 *
 * The server reports a name match as the regular expression it matched
 * (`the name matches '\bpAAV'`). The pattern is still the honest answer, so
 * it is shown rather than reworded — but set in monospace, so a lab reader
 * sees at a glance that the backslashes belong to a pattern and not to the
 * plasmid's name.
 */
function evidenceText(line) {
  var m = /^(.*?)('[^']*')(.*)$/.exec(String(line));
  if (!m) return String(line);
  return h('span', null, m[1], h('span', { class: 'mono' }, m[2]), m[3]);
}

/** A row of element names, e.g. the required features that are absent. */
function elementChips(names, kind) {
  return h('span', { class: 'elements' }, (names || []).map(function (n) {
    return h('span', { class: 'element element-' + kind }, n);
  }));
}

/**
 * What family this construct is, why, and whether it carries the full set.
 *
 * The validation panel already spells out each missing element as a
 * `packaging` issue, so this is deliberately a summary: the family, the
 * evidence behind it, and the element names grouped by severity. A conflict
 * between AAV and lentiviral evidence is shown first and loudest, because it
 * means the family underneath everything else may be wrong.
 */
function renderVector(vector, plan) {
  var done = vector.completeness || {};
  var spec = vector.spec || {};
  var complete = !!done.complete;
  var missingRequired = (done.missing_required || []).concat(done.missing_named || []);
  var missingExpected = done.missing_expected || [];

  var card = h('div', {
    class: 'card vector' + (vector.warning ? ' has-conflict' :
      (complete ? ' is-complete' : (missingRequired.length ? ' has-gaps' : '')))
  },
    h('div', { class: 'card-head' },
      h('h2', null, 'Vector family'),
      h('span', { class: 'muted small' },
        vector.requested
          ? 'requested: ' + vectorLabel(vector.requested)
          : 'family inferred from the backbone')));

  if (vector.warning) {
    card.appendChild(h('div', { class: 'conflict' },
      h('strong', null, 'Conflicting vector evidence'),
      h('span', null, vector.warning)));
  }

  card.appendChild(h('dl', { class: 'vector-facts' },
    statPair('checked as', spec.label || vectorLabel(vector.family)),
    statPair('detected', vectorLabel(vector.detected)),
    vector.requested && vector.requested !== vector.detected
      ? statPair('mismatch',
          'you asked for ' + vectorLabel(vector.requested) +
          ', the molecule looks like ' + vectorLabel(vector.detected))
      : null));

  var evidence = vector.detected_evidence || [];
  card.appendChild(h('p', { class: 'vector-evidence' + (evidence.length ? '' : ' muted') },
    h('span', { class: 'muted small' }, 'evidence: '),
    evidence.length
      ? evidence.map(function (e, i) {
          return [i ? ' · ' : null, evidenceText(e)];
        })
      : 'none — nothing in the backbone name or the annotations names a viral family'));

  // A requested family filters the backbone search, so asking for AAV and
  // naming a lentiviral backbone leaves the design with no backbone at all.
  if (vector.requested && !plan.backbone) {
    card.appendChild(h('p', { class: 'vector-note inline-warn' },
      'No backbone was loaded: when a vector family is requested, only a ' +
      'backbone of that family is accepted. Give a plasmid id to use one ' +
      'anyway, or set Vector to Auto.'));
  }

  if (complete) {
    card.appendChild(h('p', { class: 'vector-verdict ok' },
      h('span', { class: 'pill pill-ok' }, 'complete'),
      'This construct carries the full ' +
      (spec.label || vectorLabel(vector.family)) + ' element set.'));
  } else {
    card.appendChild(h('p', { class: 'vector-verdict' },
      h('span', { class: 'pill pill-' + (missingRequired.length ? 'bad' : 'warn') },
        'incomplete'),
      'Some of the ' + (spec.label || vectorLabel(vector.family)) +
      ' element set is absent; each one is listed under ' +
      'packaging in Validation below.'));
    if (missingRequired.length) {
      card.appendChild(h('div', { class: 'vector-missing' },
        h('span', { class: 'lvl lvl-error' },
          plural(missingRequired.length, 'required element') + ' missing'),
        elementChips(missingRequired, 'required')));
    }
    if (missingExpected.length) {
      card.appendChild(h('div', { class: 'vector-missing' },
        h('span', { class: 'lvl lvl-warning' },
          plural(missingExpected.length, 'expected element') + ' missing'),
        elementChips(missingExpected, 'expected')));
    }
  }

  card.appendChild(h('details', { class: 'vector-set' },
    h('summary', null, 'The full ' + (spec.label || vectorLabel(vector.family)) +
      ' element set'),
    h('dl', { class: 'vector-facts' },
      statPair('required', (spec.required || []).join(', ') || '—'),
      statPair('expected', (spec.expected || []).join(', ') || '—')),
    spec.note ? h('p', { class: 'muted small' }, spec.note) : null));

  return card;
}

function selectedRoute() {
  if (!state.plan || !state.plan.routes) return null;
  var found = state.plan.routes.filter(function (r) { return r.id === state.routeId; });
  return found.length ? found[0] : (state.plan.routes[0] || null);
}

function renderValidation(validation) {
  var issues = validation.issues || [];
  var groups = { error: [], warning: [], info: [] };
  issues.forEach(function (i) {
    (groups[i.level] || groups.info).push(i);
  });

  var card = h('div', { class: 'card validation' + (groups.error.length ? ' has-errors' : '') },
    h('div', { class: 'card-head' },
      h('h2', null, 'Validation'),
      h('div', { class: 'counts' },
        h('span', { class: 'count count-error' + (groups.error.length ? '' : ' zero') },
          plural(validation.errors || 0, 'error')),
        h('span', { class: 'count count-warning' + (groups.warning.length ? '' : ' zero') },
          plural(validation.warnings || 0, 'warning')),
        h('span', { class: 'count count-info' },
          String(validation.info || 0) + ' info'))));

  if (!issues.length) {
    card.appendChild(h('p', { class: 'muted' }, 'No issues reported.'));
    return card;
  }

  [['error', 'Errors'], ['warning', 'Warnings'], ['info', 'Information']].forEach(function (g) {
    var level = g[0];
    var list = groups[level];
    if (!list.length) return;
    card.appendChild(h('h4', { class: 'lvl lvl-' + level }, g[1] + ' (' + list.length + ')'));
    card.appendChild(h('ul', { class: 'issues issues-' + level }, list.map(function (i) {
      return h('li', null,
        h('span', { class: 'cat' }, i.category || '?'),
        h('span', { class: 'msg' }, i.message || ''),
        i.location ? h('span', { class: 'loc mono' },
          Array.isArray(i.location) ? i.location.join('–') : String(i.location)) : null);
    })));
  });
  return card;
}

function renderRoutes(routes) {
  var card = h('div', { class: 'card' },
    h('div', { class: 'card-head' },
      h('h2', null, 'Routes'),
      h('span', { class: 'muted small' }, 'ranked best first')));

  if (!routes.length) {
    card.appendChild(h('p', { class: 'muted' }, 'No routes were generated.'));
    return card;
  }

  card.appendChild(h('div', { class: 'routes' }, routes.map(function (r, i) {
    var selected = r.id === state.routeId;
    var cls = 'route' + (selected ? ' selected' : '') + (r.feasible ? '' : ' infeasible');
    return h('div', {
      class: cls, tabindex: '0',
      onclick: function () { state.routeId = r.id; renderPlan(); },
      onkeydown: function (e) {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); state.routeId = r.id; renderPlan(); }
      }
    },
      h('div', { class: 'route-head' },
        h('span', { class: 'rank' }, '#' + (i + 1)),
        h('h3', null, r.title || r.id),
        r.feasible ? null : h('span', { class: 'pill pill-bad' }, 'not feasible')),
      h('dl', { class: 'route-stats' },
        statPair('cost', num(r.cost, 2)),
        statPair('steps', r.n_steps),
        statPair('PCRs', r.n_pcr),
        statPair('assemblies', r.n_assembly),
        statPair('hands-on', hours(r.hands_on_minutes)),
        statPair('primers', (r.primers || []).length)),
      r.feasible
        ? ((r.warnings || []).length
            ? h('p', { class: 'route-note inline-warn' }, plural(r.warnings.length, 'warning'))
            : h('p', { class: 'route-note muted' }, 'no warnings'))
        : h('p', { class: 'route-note reason' },
            (r.warnings && r.warnings.length) ? r.warnings[0] : 'no reason given'));
  })));
  return card;
}

function statPair(label, value) {
  return [h('dt', null, label),
          h('dd', null, value === null || value === undefined ? '' : String(value))];
}

function renderRouteDetail(route) {
  var card = h('div', { class: 'card' },
    h('div', { class: 'card-head' },
      h('h2', null, route.title || route.id),
      h('span', { class: 'muted small' },
        'strategy: ' + (route.strategy || '?') + ' · cost ' + num(route.cost, 2) +
        ' · ' + hours(route.hands_on_minutes) + ' hands-on')));

  if (!route.feasible) {
    card.appendChild(h('div', { class: 'errorbox' },
      h('strong', null, 'This route is not feasible'),
      h('span', null, (route.warnings && route.warnings[0]) || 'no reason given')));
  }

  if (route.warnings && route.warnings.length) {
    card.appendChild(h('details', { class: 'warn-details' },
      h('summary', null, plural(route.warnings.length, 'route warning')),
      h('ul', { class: 'issues issues-warning' }, route.warnings.map(function (w) {
        return h('li', null, h('span', { class: 'msg' }, w));
      }))));
  }

  /* --- steps ------------------------------------------------------- */
  card.appendChild(h('h4', null, 'Steps (' + (route.steps || []).length + ')'));
  if (!(route.steps || []).length) {
    card.appendChild(h('p', { class: 'muted' }, 'No steps.'));
  } else {
    card.appendChild(h('ol', { class: 'steps' }, route.steps.map(function (s) {
      return h('li', null,
        h('div', { class: 'step-head' },
          h('span', { class: 'opcode mono' }, s.operation || ''),
          h('span', { class: 'step-desc' }, s.description || ''),
          h('span', { class: 'step-min muted' }, hours(s.minutes))),
        s.output ? h('div', { class: 'step-io muted small' }, '→ ' + s.output) : null,
        s.notes ? h('div', { class: 'step-note' }, s.notes) : null);
    })));
  }

  /* --- primers ----------------------------------------------------- */
  var primers = route.primers || [];
  card.appendChild(h('h4', null, 'Primers (' + primers.length + ')'));
  card.appendChild(table([
    { label: 'Name', key: 'name' },
    { label: 'Dir', cell: function (p) { return p.direction; } },
    { label: 'Length', numeric: true, key: 'length' },
    { label: 'Anneal Tm', numeric: true, cell: function (p) { return num(p.anneal_tm, 1); } },
    { label: 'GC', numeric: true, cell: function (p) { return pct(p.gc); } },
    { label: 'Sequence (tail / anneal)', cell: function (p) {
        return h('span', { class: 'mono seq' },
          p.tail ? h('span', { class: 'tail' }, String(p.tail).toLowerCase()) : null,
          h('span', { class: 'anneal' }, String(p.anneal || '').toUpperCase())); } },
    { label: 'Warnings', cell: function (p) { return warnCell(p.warnings); } }
  ], primers));

  /* --- junctions --------------------------------------------------- */
  var junctions = route.junctions || [];
  var side = function (name, type) {
    return h('span', null, name,
      h('span', { class: 'muted small' }, ' (' + (type || '') + ')'));
  };
  card.appendChild(h('h4', null, 'Junctions (' + junctions.length + ')'));
  card.appendChild(table([
    { label: '#', numeric: true, key: 'index' },
    { label: 'Left', cell: function (j) { return side(j.left_name, j.left_type); } },
    { label: 'Right', cell: function (j) { return side(j.right_name, j.right_type); } },
    { label: 'Method', key: 'method' },
    { label: 'Overlap', numeric: true, key: 'overlap_length' },
    { label: 'Frame shift', numeric: true, key: 'frame_shift' },
    { label: 'Added residues', cell: function (j) {
        return h('span', { class: 'mono' }, j.added_residues || ''); } },
    { label: 'Sequence', cell: function (j) {
        return h('span', { class: 'mono seq' }, String(j.sequence || '').toUpperCase()); } },
    { label: 'Warnings', cell: function (j) { return warnCell(j.warnings); } }
  ], junctions));

  return card;
}

function renderExports() {
  var status = h('span', { class: 'status' });
  var plan = state.plan;
  var route = selectedRoute();
  var base = ((plan && plan.construct && plan.construct.name) || 'construct')
    .replace(/\s+/g, '_');

  function exportButton(label, kind, fallback) {
    return h('button', {
      type: 'button', class: 'secondary',
      onclick: async function (e) {
        var btn = e.currentTarget;
        btn.disabled = true;
        setStatus(status, 'exporting ' + kind + '…', 'busy');
        try {
          await downloadExport(kind,
            { plan: plan, route_id: route ? route.id : null }, fallback);
          setStatus(status, kind + ' downloaded', 'ok');
        } catch (err) {
          setStatus(status, 'export failed: ' + err.message, 'error');
        } finally {
          btn.disabled = false;
        }
      }
    }, label);
  }

  return h('div', { class: 'card exports' },
    h('div', { class: 'card-head' },
      h('h2', null, 'Export'),
      h('span', { class: 'muted small' },
        route ? 'route: ' + (route.title || route.id) : 'no route selected')),
    h('div', { class: 'button-row' },
      exportButton('Protocol (.md)', 'protocol', base + '_protocol.md'),
      exportButton('PDF protocol (.pdf)', 'pdf', base + '_protocol.pdf'),
      exportButton('Primers (.csv)', 'primers', base + '_primers.csv'),
      exportButton('GenBank (.gb)', 'genbank', base + '.gb'),
      status));
}

/* ------------------------------------------------------------------ */
/* Tab 4: Designs                                                      */
/* ------------------------------------------------------------------ */

/* The listing is capped at 100 rows server-side, so a version count taken
 * from it is a count of what is listed. At the cap the confirmation says
 * "at least", because an irreversible action must not overstate what it
 * knows. */
var DESIGNS_LIMIT = 100;

/** Group the listing by construct: { construct_id: [rows, newest first] }. */
function groupByConstruct(rows) {
  var by = {};
  rows.forEach(function (r) {
    var key = String(r.construct_id);
    (by[key] = by[key] || []).push(r);
  });
  return by;
}

async function renderDesigns() {
  var out = $('designs-table');
  var status = $('designs-status');
  setStatus(status, 'loading…', 'busy');
  clear(out);
  out.appendChild(h('p', { class: 'muted' }, 'Loading saved designs…'));
  try {
    var rows = await api('/api/designs');
    var by = groupByConstruct(rows);
    var capped = rows.length >= DESIGNS_LIMIT;
    clear(out);
    out.appendChild(table([
      { label: 'Design', numeric: true, key: 'id' },
      { label: 'Construct', key: 'name' },
      { label: 'Version', numeric: true, key: 'version' },
      { label: 'Owner', key: 'owner' },
      { label: 'Created', key: 'created_at' },
      {
        label: 'Delete',
        cell: function (r) {
          var siblings = by[String(r.construct_id)] || [r];
          return h('span', {
            class: 'row-delete',
            // The row itself opens the design; the controls must not.
            onclick: function (e) { e.stopPropagation(); }
          },
            h('button', {
              type: 'button', class: 'danger',
              title: 'Delete version ' + r.version + ' of “' + r.name + '”',
              onclick: function () { confirmDeleteDesign(r, siblings, capped); }
            }, 'version'),
            h('button', {
              type: 'button', class: 'danger',
              title: 'Delete “' + r.name + '” and every version of it',
              onclick: function () { confirmDeleteConstruct(r, siblings, capped); }
            }, 'construct' + (siblings.length > 1
              ? ' (' + siblings.length + ')' : '')));
        }
      }
    ], rows, { onRow: function (r) { openDesign(r.id); } }));
    setStatus(status, plural(rows.length, 'design'));
  } catch (err) {
    clear(out);
    out.appendChild(errorBox(err.message));
    setStatus(status, 'failed', 'error');
  }
}

/** "versions 1, 2 and 3" — the exact list, so nothing is deleted vaguely. */
function versionList(rows) {
  var vs = rows.map(function (r) { return r.version; })
    .sort(function (a, b) { return a - b; });
  if (vs.length === 1) return 'version ' + vs[0];
  return 'versions ' + vs.slice(0, -1).join(', ') + ' and ' + vs[vs.length - 1];
}

/** Run a delete and refresh everything the deleted rows could be showing. */
async function runDelete(payload) {
  var result = await postJSON('/api/delete-design', payload);
  // A plan on screen that pointed at a deleted design is no longer saved.
  if (state.saved) {
    var gone = (payload.design_id && state.saved.design_id === payload.design_id) ||
      (payload.construct_id && state.saved.construct_id === payload.construct_id);
    if (gone) { state.saved = null; refreshSaveButtons(); }
  }
  loadStats();
  // Awaited, so the caller's "deleted …" line is written after the reload's
  // own "N designs" status rather than being overwritten by it.
  await renderDesigns();
  return result;
}

function confirmDeleteDesign(row, siblings, capped) {
  var status = $('designs-status');
  var last = siblings.length <= 1;
  var others = siblings.length - 1;

  modal({
    title: 'Delete one version',
    danger: true,
    confirmLabel: 'Delete version ' + row.version,
    body: [
      h('p', null,
        'This permanently deletes ',
        h('b', null, 'version ' + row.version + ' of “' + row.name + '”'),
        ' (design #' + row.id + ') and any experimental outcomes recorded ' +
        'against it.'),
      h('p', { class: last ? 'modal-stress' : '' }, last
        ? 'It is the only version, so the construct “' + row.name +
          '” is removed along with it.'
        : (capped
            ? 'At least ' + plural(others, 'other version') +
              ' of “' + row.name + '” remain' + (others === 1 ? 's' : '') + '.'
            : plural(others, 'other version') + ' of “' + row.name +
              '” remain' + (others === 1 ? 's' : '') + ': ' +
              versionList(siblings.filter(function (s) { return s.id !== row.id; })) + '.')),
      h('p', { class: 'modal-stress' }, 'This cannot be undone.')
    ],
    onConfirm: function (api_) {
      setStatus(status, 'deleting design ' + row.id + '…', 'busy');
      return runDelete({ design_id: row.id }).then(function (res) {
        api_.close();
        setStatus(status,
          'deleted version ' + res.version + ' of “' + res.name + '”' +
          (res.versions_left
            ? ' (' + plural(res.versions_left, 'version') + ' left)'
            : ' (the construct is gone too)'), 'ok');
      }, function (err) {
        api_.setError(err.message);
        setStatus(status, 'delete failed', 'error');
      });
    }
  });
}

function confirmDeleteConstruct(row, siblings, capped) {
  var status = $('designs-status');
  var n = siblings.length;

  modal({
    title: 'Delete a whole construct',
    danger: true,
    confirmLabel: 'Delete “' + row.name + '” and ' +
      (capped ? 'every version' : plural(n, 'version')),
    body: [
      h('p', null,
        'This permanently deletes the construct ',
        h('b', null, '“' + row.name + '”'),
        ' and ',
        h('b', null, capped
          ? 'every version of it (at least ' + n + ' listed here)'
          : 'all ' + plural(n, 'version') + ' of it'),
        ', together with any experimental outcomes recorded against them.'),
      capped ? null : h('p', null,
        'Going: ' + versionList(siblings) + ' — design' + (n === 1 ? ' ' : 's ') +
        siblings.map(function (s) { return '#' + s.id; }).join(', ') + '.'),
      h('p', { class: 'modal-stress' }, 'This cannot be undone.')
    ],
    onConfirm: function (api_) {
      setStatus(status, 'deleting construct ' + row.construct_id + '…', 'busy');
      return runDelete({ construct_id: row.construct_id }).then(function (res) {
        api_.close();
        setStatus(status,
          'deleted “' + res.name + '” and ' +
          plural(res.versions_deleted, 'version'), 'ok');
      }, function (err) {
        api_.setError(err.message);
        setStatus(status, 'delete failed', 'error');
      });
    }
  });
}

async function openDesign(id) {
  var status = $('designs-status');
  setStatus(status, 'loading design ' + id + '…', 'busy');
  try {
    var d = await api('/api/design/' + id);
    state.plan = d.plan;
    state.appliedSuggestion = null;
    // The stored row carries the spec as well as the plan, so the design can
    // be edited, re-saved as a new version, or overwritten in place.
    state.spec = d.spec || null;
    state.saved = {
      design_id: d.id,
      construct_id: d.construct_id,
      version: d.version,
      name: d.construct_name || specName(d.spec || {})
    };
    state.routeId = (d.plan && d.plan.routes && d.plan.routes.length)
      ? d.plan.routes[0].id : null;
    var filled = applySpecToBuilder(d.spec);
    refreshSaveButtons();
    selectTab('design');
    renderPlan();
    setStatus($('design-status'),
      'loaded design ' + d.id + ' (version ' + d.version + ' of “' +
      state.saved.name + '”)' +
      (filled ? '' : ' — no stored specification, so it cannot be re-saved'),
      filled ? 'ok' : 'error');
    setStatus(status, '');
  } catch (err) {
    setStatus(status, 'failed: ' + err.message, 'error');
  }
}

/* ------------------------------------------------------------------ */
/* Tabs and wiring                                                     */
/* ------------------------------------------------------------------ */

var LOADED = {};

function selectTab(name) {
  state.tab = name;
  var buttons = document.querySelectorAll('#tabs .tab');
  Array.prototype.forEach.call(buttons, function (b) {
    b.classList.toggle('active', b.dataset.tab === name);
  });
  ['inventory', 'parts', 'design', 'designs'].forEach(function (t) {
    $('panel-' + t).hidden = (t !== name);
  });

  // Each list tab loads once on first view; the buttons reload on demand.
  if (name === 'inventory' && !LOADED.inventory) { LOADED.inventory = true; renderInventory(); }
  if (name === 'parts' && !LOADED.parts) { LOADED.parts = true; renderParts(); }
  if (name === 'designs') { LOADED.designs = true; renderDesigns(); }
}

async function init() {
  document.querySelectorAll('#tabs .tab').forEach(function (btn) {
    btn.addEventListener('click', function () { selectTab(btn.dataset.tab); });
  });

  $('inv-go').addEventListener('click', renderInventory);
  $('inv-q').addEventListener('keydown', function (e) {
    if (e.key === 'Enter') renderInventory();
  });
  $('inv-owner').addEventListener('change', renderInventory);

  $('part-go').addEventListener('click', renderParts);
  $('part-q').addEventListener('keydown', function (e) {
    if (e.key === 'Enter') renderParts();
  });
  $('part-type').addEventListener('change', renderParts);

  $('add-module').addEventListener('click', function () {
    state.modules.push(newModule('CDS'));
    renderModules();
  });
  $('generate').addEventListener('click', generatePlan);
  $('clear-spec').addEventListener('click', clearSpecWithConfirm);
  $('save-version').addEventListener('click', function () { saveDesign('version'); });
  $('save-overwrite').addEventListener('click', function () { saveDesign('overwrite'); });
  $('save-as').addEventListener('click', saveAsDesign);
  $('designs-go').addEventListener('click', renderDesigns);

  document.querySelectorAll('input[name="spec-mode"]').forEach(function (radio) {
    radio.addEventListener('change', function () {
      var yaml = usingYaml();
      $('builder-pane').hidden = yaml;
      $('yaml-pane').hidden = !yaml;
      if (yaml && !$('spec-yaml').value.trim()) $('spec-yaml').value = specToYaml();
    });
  });

  // The type selects and the cassette rule are both driven by the server's
  // vocabulary, so it is fetched before the first row is drawn. Everything
  // after this point is independent and loads in parallel.
  await loadVocab();

  renderModules();
  refreshSaveButtons();
  loadStats();
  loadOwners();
  selectTab('design');
}

/** Seed the YAML textarea from whatever the builder currently holds. */
function specToYaml() {
  var backbone = $('spec-backbone').value.trim();
  var lines = ['name: ' + ($('spec-name').value.trim() || 'untitled')];
  if (backbone) {
    lines.push((/^\d+$/.test(backbone) ? 'backbone_plasmid_id: ' : 'backbone: ') + backbone);
  }
  var owner = $('spec-owner').value.trim();
  if (owner) lines.push('owner: ' + owner);
  var vector = $('spec-vector').value;
  if (vector) lines.push('vector: ' + vector);
  lines.push('circular: ' + ($('spec-circular').checked ? 'true' : 'false'));
  lines.push('modules:');
  state.modules.forEach(function (m) {
    lines.push('  - type: ' + m.type);
    if (isCreLox(m.type)) lines.push('    mode: ' + creLoxMode(m.mode).key);
    if ((m.name || '').trim()) lines.push('    name: ' + m.name.trim());
    if (m.part_id) lines.push('    part_id: ' + m.part_id);
  });
  if (pinnedParents.recipient) {
    lines.splice(1, 0, 'recipient_plasmid_id: ' + pinnedParents.recipient);
  }
  if (pinnedParents.donor) {
    lines.splice(1, 0, 'donor_plasmid_id: ' + pinnedParents.donor);
  }
  return lines.join('\n') + '\n';
}

document.addEventListener('DOMContentLoaded', init);
