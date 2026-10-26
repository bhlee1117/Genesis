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
      /* An item may carry its own colour: the simulation colours a stretch
         by which reaction produced it, not by what kind of feature it is. */
      fill: f.colour || colorFor(type),
      class: 'feat' + (reverse ? ' feat-rev' : '')
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

/* Only these have two ends of their own. A 2A peptide or a trafficking
   motif belongs to somebody else's protein, so labelling it N and C
   would claim there are more proteins on the map than there are. */
var PROTEIN_END_TYPES = { CDS: 1 };

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
      typeOf(f) !== 'primer_binding_site' && typeOf(f) !== 'unknown';
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
      var colour = f.colour || colorFor(typeOf(f));

      var g = svg('g', { class: 'circ-feature' });
      g.appendChild(svg('title', null,
        f.name + ' · ' + typeOf(f) + ' · ' + bp(f.end - f.start) +
        ' · ' + (reverse ? 'reverse' : 'forward') + ' strand · ' +
        (f.start + 1) + '–' + f.end));
      g.appendChild(svg('path', {
        d: arcPath(cx, cy, r, from, to), fill: 'none',
        stroke: colour, 'stroke-width': 8, 'stroke-linecap': 'butt'
      }));
      /* direction of translation, at the 3' end of the coding sequence */
      if (CODING_TYPES[typeOf(f)]) {
        g.appendChild(arrowHead(cx, cy, r, reverse ? from : to % 1,
                                !reverse, colour));
      }
      /* Which end is the start of the protein and which is the end of it.
         A reverse-strand coding sequence begins at the end that is drawn
         second, so the labels follow the strand and not the drawing. */
      if (PROTEIN_END_TYPES[typeOf(f)] && (to - from) > 0.035) {
        var nAt = (reverse ? to : from) % 1;
        var cAt = (reverse ? from : to) % 1;
        [['N', nAt], ['C', cAt]].forEach(function (pair) {
          var at = polar(cx, cy, r - 9, pair[1]);
          var t = svg('text', {
            x: at[0].toFixed(1), y: (at[1] + 2).toFixed(1),
            'text-anchor': 'middle', class: 'circ-terminus', fill: colour
          });
          t.textContent = pair[0];
          g.appendChild(t);
        });
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
      }, item.f.name || typeOf(item.f)));
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
    'Outer ring forward, inner ring reverse; arrowheads show which way a '
    + 'CDS is read. Position 1 at the top.');
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

    side.appendChild(h('h4', null, 'Map'));
    side.appendChild(mapWithToggle(p.seq_length, p.annotations || [],
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

/* The standard genetic code, as the one table it is. */
var CODONS = (function () {
  var aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  var b = 'TCAG', table = {}, i = 0;
  for (var x = 0; x < 4; x++)
    for (var y = 0; y < 4; y++)
      for (var z = 0; z < 4; z++)
        table[b[x] + b[y] + b[z]] = aas[i++];
  return table;
})();

function translateDna(seq) {
  var out = '';
  seq = String(seq || '').toUpperCase();
  for (var i = 0; i + 2 < seq.length; i += 3) {
    out += CODONS[seq.substr(i, 3)] || 'X';
  }
  return out;
}

/** A sequence in numbered rows, the way a sequence is read.

   One unbroken string is technically the same information and is
   unusable: nobody can tell base 400 from base 430 in it. Sixty to a
   line, grouped in tens, with the coordinate of the first base. */
function seqBlock(seq, opts) {
  opts = opts || {};
  var per = opts.per || 60, group = opts.group || 10;
  var pre = h('pre', { class: 'seq-block' });
  var lines = [];
  for (var i = 0; i < seq.length; i += per) {
    var chunk = seq.substr(i, per), parts = [];
    for (var j = 0; j < chunk.length; j += group) {
      parts.push(chunk.substr(j, group));
    }
    lines.push(String(i + 1).padStart(6, ' ') + '  ' + parts.join(' '));
  }
  pre.textContent = lines.join('\n');
  return pre;
}

function seqSection(title, seq, note, extra) {
  var open = false;
  var body = h('div', { hidden: 'hidden' });
  var head = h('button', {
    type: 'button', class: 'link-button',
    onclick: function () {
      open = !open;
      body.hidden = !open;
      head.textContent = (open ? 'hide ' : 'show ') + title;
    }
  }, 'show ' + title);
  body.appendChild(h('p', { class: 'muted small' }, note));
  body.appendChild(seqBlock(seq, extra));
  return h('div', { class: 'seq-wrap' }, head, body);
}

async function showPartSources(part) {
  var side = $('part-detail');
  clear(side);
  side.appendChild(h('h3', null, part.name || '(unnamed)'));
  side.appendChild(h('p', { class: 'muted small' },
    part.feature_type + ' · ' + bp(part.seq_length)));
  var seqHolder = h('div', null);
  side.appendChild(seqHolder);
  var holder = h('div', null, h('p', { class: 'muted' }, 'Loading sources…'));
  side.appendChild(holder);

  /* The sequence itself - which is, after all, what a part is. */
  api('/api/part/' + part.id).then(function (full) {
    var dna = String(full.seq || '').toUpperCase();
    if (!dna) return;
    clear(seqHolder);
    seqHolder.appendChild(seqSection(
      'DNA (' + bp(dna.length) + ")", dna,
      "5' to 3', numbered from the first base of the part."));
    if (part.feature_type === 'CDS' || full.feature_type === 'CDS') {
      var aa = translateDna(dna);
      var trailing = aa.length && aa[aa.length - 1] === '*';
      seqHolder.appendChild(seqSection(
        'protein (' + (aa.length - (trailing ? 1 : 0)) + ' aa)', aa,
        'Translated in frame 1, N to C. '
        + (trailing ? 'The terminal * is the stop codon.'
                    : 'No terminal stop codon in this part.'),
        { per: 60, group: 10 }));
    }
  }).catch(function () { /* the sources below are the point */ });
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

var FAMILY_CACHE = {};

/* One entry per protein rather than per annotation. A name search for
   "CheRiff" returns the same molecule seventeen times and hides the
   soma-targeted one behind an identical label. */
function fetchFamilies(type, query) {
  var key = type + '\u0000' + query.toLowerCase();
  if (!FAMILY_CACHE[key]) {
    FAMILY_CACHE[key] = api('/api/parts' + qs({
      type: type, q: query, group: 'protein', limit: 14
    })).catch(function (err) { delete FAMILY_CACHE[key]; throw err; });
  }
  return FAMILY_CACHE[key];
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
    h('option', { value: 'protein' }, 'by protein'),
    h('option', { value: 'all' }, 'all forms'),
    h('option', { value: 'whole' }, 'whole only'),
    h('option', { value: 'split' }, 'pieces only'));
  /* Grouping is the useful default for a coding sequence, where the same
     protein is catalogued many times over; for everything else the plain
     list is what people expect. */
  formSel.value = mod.form || (mod.type === 'CDS' ? 'protein' : 'all');
  mod.form = formSel.value;

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
        mod.part_id = null; mod.pin_note = ''; mod.part_bp = null;
        renderChip(); renderSpecPreview(); input.focus();
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
  /* Move the highlight, skipping family headings: they name a group, they
     are not something to pin, so landing on one would be a dead stop. */
  function step(from, delta) {
    var n = items.length;
    if (!n) return -1;
    var at = from;
    for (var k = 0; k < n; k++) {
      at = ((at + delta) % n + n) % n;
      if (items[at].kind !== 'head') return at;
    }
    return from;
  }

  function choose(i) {
    var row = items[i];
    if (!row) return;
    if (row.kind === 'head') return;   // a family heading names, it does not pin
    if (row.kind === 'back') { showParts(); return; }
    mod.name = row.name || '';
    mod.part_id = row.id;
    mod.pin_note = row.note || '';
    /* The sketch was drawn from typical lengths for each kind of element.
       Once a real part is chosen its real length is known, so the drawing
       stops being a guess for that module. */
    var src = row.source || {};
    mod.part_bp = Number(src.bp || src.seq_length) || null;
    input.value = mod.name;
    renderChip();
    close();
    renderSpecPreview();
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
               (r.kind === 'segment' ? ' combo-segment' : '') +
               (r.kind === 'head' ? ' combo-head' : '') +
               (r.kind === 'part' && r.source && r.source.decorations
                  ? ' combo-variant' : ''),
        role: r.kind === 'head' ? 'presentation' : 'option',
        id: listId + '-' + i, 'aria-selected': 'false',
        title: r.title || null,
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
      title: n + ' plasmid' + (n === 1 ? '' : 's') +
             ' in the inventory carry this part — it is how many tubes you '
             + 'could amplify it from, not how good it is'
    }, n + (n === 1 ? ' plasmid' : ' plasmids'));
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

  /** A heading row naming the protein a group of variants share. */
  function familyHead(fam) {
    return {
      kind: 'head', id: null, name: fam.core.name,
      cells: function () {
        var hosts = (fam.hosts || []).slice(0, 2).join(', ');
        return [
          h('span', { class: 'fam-core' }, fam.core.name),
          h('span', { class: 'fam-aa mono' }, fam.core.aa + ' aa'),
          /* Several families can carry one name - QuasAr2 is four distinct
             proteins here - so say how this one differs and where it is. */
          fam.differs
            ? h('span', { class: 'fam-differs', title:
                'a different protein from the most-used one of this name' },
                fam.differs)
            : null,
          hosts
            ? h('span', { class: 'fam-hosts', title:
                'plasmids carrying this exact protein: ' +
                (fam.hosts || []).join(', ') }, 'in ' + hosts)
            : null,
          h('span', {
            class: 'fam-count',
            title: fam.n_rows + ' catalogued row(s) encode this protein, '
              + 'across ' + fam.sources + ' plasmid occurrence(s)'
          }, fam.n_variants + (fam.n_variants === 1 ? ' form' : ' forms'))
        ];
      }
    };
  }

  /** One member of a family, named for what is fused to the core. */
  function variantRow(fam, v) {
    var alsoTitle = v.also && v.also.length
      ? '\nalso catalogued as: ' + v.also.join(', ') : '';
    /* Which tubes this exact protein is in. "20 plasmids" answers how
       many; the only question after that is which, and the answer was
       nowhere on the page. */
    var hosts = v.hosts || [];
    var more = Math.max(0, (v.n_hosts || hosts.length) - hosts.length);
    var where = hosts.length
      ? 'In ' + (v.n_hosts || hosts.length) + ' plasmid'
        + ((v.n_hosts || hosts.length) === 1 ? '' : 's') + ':\n  '
        + hosts.join('\n  ') + (more ? '\n  …and ' + more + ' more' : '')
      : 'No plasmid in the inventory carries this.';
    return {
      kind: 'part', id: v.part_id, name: fam.core.name, source: v,
      note: v.label === 'plain' ? '' : v.label,
      title: where + alsoTitle,
      cells: function () {
        return [
          h('span', {
            class: 'var-label' + (v.is_core ? ' plain' : ''),
            title: where
          }, v.label === 'plain' ? fam.core.name : v.label),
          v.n_rows > 1
            ? h('span', {
                class: 'var-rows',
                title: v.n_rows + ' catalogued rows are this same protein, '
                  + 'differing only in where the annotation starts or ends'
                  + alsoTitle
              }, '×' + v.n_rows)
            : null,
          h('span', { class: 'combo-len mono' }, bp(v.bp)),
          h('span', {
            class: 'combo-sources' + (v.n_sources ? '' : ' zero'),
            title: where
          }, v.n_sources + (v.n_sources === 1 ? ' plasmid' : ' plasmids')),
          /* A visible hint that there is something to hover. */
          hosts.length
            ? h('span', { class: 'var-where', title: where },
                hosts[0].length > 22 ? hosts[0].slice(0, 21) + '…' : hosts[0])
            : null
        ];
      }
    };
  }

  function drawFamilies(groups, query) {
    items = [];
    groups.forEach(function (fam) {
      items.push(familyHead(fam));
      fam.variants.forEach(function (v) { items.push(variantRow(fam, v)); });
    });
    active = -1;
    if (!items.length) {
      show(h('div', { class: 'combo-empty muted' },
        'No protein matching “' + query + '”'));
      return;
    }
    paintList();
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
    if (form === 'protein' && type === 'CDS') {
      fetchFamilies(type, query).then(function (data) {
        if (mine !== reqId || type !== mod.type) return;
        drawFamilies((data && data.groups) || [], query);
      }).catch(function (err) {
        if (mine !== reqId) return;
        items = [];
        show(h('div', { class: 'combo-empty inline-warn' },
          'search failed: ' + err.message));
      });
      return;
    }
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
    if (mod.part_id) {
      mod.part_id = null; mod.pin_note = ''; mod.part_bp = null; renderChip();
    }
    /* Redraw as it is typed, not only when the row list is rebuilt. */
    renderSpecPreview();
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
      if (items.length) { active = step(active, 1); paint(); }
    } else if (e.key === 'ArrowUp') {
      e.preventDefault();
      if (list.hidden || !items.length) return;
      active = step(active, -1);
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

  var custom = customSeqEditor(mod, i);

  return h('div', { class: 'module-wrap' },
    h('div', { class: 'module-row' + (isCreLox(mod.type) ? ' is-crelox' : '')
                      + (mod.seq ? ' is-literal' : '') },
    h('span', { class: 'module-index' }, String(i + 1)),
    typeSel,
    modeSel,
    mod.seq ? h('span', { class: 'literal-note' },
                'own sequence, ' + bp(cleanDna(mod.seq).length)) : picker.node,
    mod.seq ? null : picker.formFilter,
    mod.seq ? null : picker.chip,
    h('div', { class: 'row-buttons' },
      h('button', {
        type: 'button', class: mod.seq ? 'is-on' : '',
        title: mod.seq ? 'Edit the sequence' : 'Paste your own sequence',
        onclick: function () { custom.toggle(); }
      }, 'DNA'),
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
      }, '×'))),
    custom.node);
}

/**
 * Draw the builder, bracketing each Cre-lox cassette around what it encloses.
 *
 * The bracket is presentation only: `state.modules` stays the flat, ordered
 * list the server is sent, and the row numbers keep counting straight
 * through the group.
 */

/* ------------------------------------------------------------------ */
/* Specification: a live sketch of what is being described             */
/* ------------------------------------------------------------------ */
/* The module list says what goes in and in what order, but not what the
   plasmid looks like, and you only found that out by designing. This
   draws it while you type. It is a sketch: nothing has been resolved, so
   the lengths are typical ones for each kind of element, and the vector
   is whatever the backbone will bring. */

var NOMINAL_BP = {
  promoter: 800, enhancer: 400, CDS: 1000, p2A: 66, linker: 45,
  tag: 30, trafficking_motif: 60, IRES: 580, polyA: 480, WPRE: 590,
  lox_site: 34, AAV_ITR: 141, origin: 590, resistance: 860,
  Kozak: 6, intron: 1000, terminator: 220
};

/* What a Cre-lox cassette puts around its payload. */
var LOX_OF = { DIO: ['lox2272', 'loxP'], FLEX: ['lox2272', 'loxP'],
               FAS: ['loxFAS', 'loxP'], DO: ['lox2272', 'loxP'] };

function nominalBp(mod) {
  // A sequence that was typed in is not nominal; it is the length.
  if (mod.seq) return cleanDna(mod.seq).length || 1;
  if (mod.part_bp) return Number(mod.part_bp);
  return NOMINAL_BP[mod.type] || 300;
}

/** Only the bases, out of whatever was pasted.

   People paste FASTA, they paste a GenBank ORIGIN block with its line
   numbers, they paste something with a space every ten. It is all the
   same sequence. */
function cleanDna(text) {
  return String(text || '')
    .replace(/^>.*$/gm, '')
    .replace(/[^A-Za-z]/g, '')
    .toUpperCase();
}

/** A box for a sequence the inventory does not have.

   Some parts are not in any plasmid here yet - a tag out of a paper, a
   linker somebody designed this morning. Genesis has always taken a
   literal sequence through the API and had nowhere to type one. */
function customSeqEditor(mod, i) {
  var box = h('div', { class: 'literal-box', hidden: 'hidden' });
  var status = h('div', { class: 'literal-status muted small' }, '');
  var area = h('textarea', {
    rows: '4', spellcheck: 'false',
    placeholder: "Paste DNA 5' to 3'. FASTA headers, numbers and spaces "
                 + "are ignored.",
    oninput: function (e) { review(e.target.value); }
  });
  area.value = mod.seq || '';

  function review(raw) {
    var dna = cleanDna(raw);
    var bad = dna.replace(/[ACGTUNRYKMSWBDHV]/g, '');
    var notes = [], wrong = false;
    if (!dna) {
      notes.push('empty - the row goes back to searching the inventory');
    } else {
      notes.push(bp(dna.length));
      if (bad) {
        wrong = true;
        notes.push('not DNA: ' + Array.from(new Set(bad.split(''))).join(' '));
      }
      if (mod.type === 'CDS') {
        if (dna.length % 3) {
          wrong = true;
          notes.push('NOT a multiple of 3 - this shifts the frame');
        } else {
          notes.push('a whole number of codons');
          var prot = translateDna(dna);
          var inner = prot.slice(0, -1).indexOf('*');
          if (inner >= 0) {
            wrong = true;
            notes.push('stop codon at residue ' + (inner + 1)
                       + ', inside the sequence');
          }
        }
      }
    }
    clear(status);
    status.appendChild(h('span', { class: wrong ? 'bad' : '' },
                         notes.join(' \u00b7 ')));
    return dna;
  }

  box.appendChild(area);
  box.appendChild(status);
  box.appendChild(h('div', { class: 'button-row' },
    h('button', {
      type: 'button', class: 'primary',
      onclick: function () {
        mod.seq = cleanDna(area.value);
        if (mod.seq) mod.part_id = null;   // a literal is not a pinned part
        renderModules();
        renderSpecPreview();
      }
    }, 'Use this sequence'),
    h('button', {
      type: 'button',
      onclick: function () {
        mod.seq = '';
        renderModules();
        renderSpecPreview();
      }
    }, 'Clear')));
  review(area.value);

  return {
    node: box,
    toggle: function () {
      box.hidden = !box.hidden;
      if (!box.hidden) area.focus();
    }
  };
}

/** Build the sketch's feature list from the modules as typed. */
function specPreviewItems() {
  var items = [];
  var at = 0;
  function put(type, name, bp_, strand) {
    items.push({ feature_type: type, name: name, start: at,
                 end: at + bp_, strand: strand || 1 });
    at += bp_;
  }

  /* The vector the backbone brings, drawn muted so the payload reads as
     the thing being designed. */
  put('AAV_ITR', "5' ITR", NOMINAL_BP.AAV_ITR);

  var groups = cassetteGroups(state.modules);
  var inCassette = {};
  var lastOf = {};
  groups.forEach(function (g) {
    (g.members || []).forEach(function (k) { inCassette[k] = g; });
    if ((g.members || []).length) {
      lastOf[g.members[g.members.length - 1]] = g;
    }
  });

  state.modules.forEach(function (mod, i) {
    if (!isCreLox(mod.type) && !(mod.name || '').trim() && !mod.part_id) return;
    if (isCreLox(mod.type)) {
      var pair = LOX_OF[(mod.mode || 'DIO').toUpperCase()] || LOX_OF.DIO;
      put('lox_site', pair[0], NOMINAL_BP.lox_site);
      put('lox_site', pair[1], NOMINAL_BP.lox_site);
      return;
    }
    var g = inCassette[i];
    put(mod.type, mod.name || mod.type, nominalBp(mod), g ? -1 : 1);
    /* Close the cassette after its last enclosed module. */
    if (lastOf[i]) {
      g = lastOf[i];
      var pair2 = LOX_OF[((state.modules[g.head] || {}).mode || 'DIO')
        .toUpperCase()] || LOX_OF.DIO;
      put('lox_site', pair2[0], NOMINAL_BP.lox_site, -1);
      put('lox_site', pair2[1], NOMINAL_BP.lox_site, -1);
    }
  });

  put('WPRE', 'WPRE', NOMINAL_BP.WPRE);
  put('polyA', 'polyA', NOMINAL_BP.polyA);
  put('AAV_ITR', "3' ITR", NOMINAL_BP.AAV_ITR);
  put('origin', 'ori', NOMINAL_BP.origin);
  put('resistance', 'AmpR', NOMINAL_BP.resistance);
  return { items: items, length: at };
}

function renderSpecPreview() {
  var host = $('spec-preview');
  if (!host) return;
  clear(host);
  /* An empty row is not a module yet: drawing a 4.6 kb plasmid from one
     blank promoter box claimed more than the user had said. */
  var real = state.modules.filter(function (m) {
    return m.type && !isCreLox(m.type) && ((m.name || '').trim() || m.part_id);
  }).length;
  if (!real) {
    host.appendChild(h('p', { class: 'muted small' },
      'Add a module and the plasmid it describes will be sketched here.'));
    return;
  }
  var built = specPreviewItems();
  host.appendChild(h('div', { class: 'card-head' },
    h('h3', null, 'Sketch'),
    h('span', { class: 'muted small' },
      '~' + bp(built.length) + ' · nothing resolved yet')));
  host.appendChild(circularMap(built.length, built.items, {
    label: ($('spec-name') && $('spec-name').value) || 'construct'
  }));
  host.appendChild(h('p', { class: 'muted small' },
    'Typical lengths, not real parts \u2014 Scan source replaces them.'));
}

function renderModules(focusIndex, focusWhat) {
  var host = $('modules');
  clear(host);
  if (!state.modules.length) {
    host.appendChild(h('p', { class: 'muted' }, 'No modules yet. Add one below.'));
    renderSpecPreview();
    return;
  }
  renderSpecPreview();

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

/* Which plasmid to take each part from, as {part_id: plasmid_id}. Chosen in
   the Candidate parents card; the planner otherwise ranks the sources. */
var pinnedSources = {};

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
      // A sequence typed in by hand is the part: the server takes it as
      // given and orders it rather than hunting the inventory for it.
      if (m.seq) out.seq = m.seq;
      return out;
    })
  };
  // A bare integer is treated as a plasmid id, anything else as a name.
  if (/^\d+$/.test(backbone)) spec.backbone_plasmid_id = parseInt(backbone, 10);
  else spec.backbone = backbone;
  if (pinnedParents.recipient) spec.recipient_plasmid_id = pinnedParents.recipient;
  if (pinnedParents.donor) spec.donor_plasmid_id = pinnedParents.donor;
  if (Object.keys(pinnedSources).length) spec.source_plasmid_ids = pinnedSources;
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
  pinnedSources = {};
  state.plan = null;
  state.scan = null;
  state.spec = null;
  state.routeId = null;
  state.saved = null;
  state.appliedSuggestion = null;

  resetSteps();
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

/**
 * Resolve the specification and show what it could be built from.
 *
 * Designing primers is the expensive half and every one of them depends on
 * which parents are used, so the parents come first: this stops at the
 * candidates, and the protocol is written once they are settled.
 */
async function scanSources() {
  var status = $('design-status');
  var out = $('plan-out');
  setStatus(status, 'resolving modules...', 'busy');
  clear(out);
  out.appendChild(h('div', { class: 'card' },
    h('p', { class: 'muted' }, 'Looking for what these parts are in...')));
  try {
    var spec = collectSpec();
    var scan = await postJSON('/api/scan', { spec: spec });
    state.spec = spec;
    state.scan = scan;
    state.plan = null;              /* the old protocol is not this spec's */
    state.routeId = null;
    refreshSaveButtons();
    renderScan();
    setStatus(status,
      'scanned - choose the parents, then generate the protocol', 'ok');
  } catch (err) {
    clear(out);
    out.appendChild(errorBox(err.message));
    setStatus(status, 'scan failed', 'bad');
  }
}


/* ------------------------------------------------------------------ */
/* Design steps                                                        */
/* ------------------------------------------------------------------ */
/* The design is a sequence - say what you want, see it, look for parts,
   choose parents, pick a route, read the scheme - and it was laid out as
   one column several screens tall. It is now one pane per step, so the
   thing you are working on is the thing on screen. */

var STEPS = [
  { id: 'spec',    n: 1, label: 'Specification',
    hint: 'what to build' },
  { id: 'target',  n: 2, label: 'Target construct',
    hint: 'the molecule asked for' },
  /* Scanning the source and reading the parents it turned up are one
     act, not two: the candidates ARE the result of the scan, and
     splitting them meant clicking across to see what the click did. */
  { id: 'scan',    n: 3, label: 'Scan source',
    hint: 'what the inventory has' },
  /* "Select" was a step of its own, which put a click between choosing
     a parent and seeing what you had chosen - on a page where the
     choosing happens. It is a summary at the top of Scan source now. */
  { id: 'routes',  n: 4, label: 'Routes',
    hint: 'how to build it' },
  { id: 'scheme',  n: 5, label: 'Protocol',
    hint: 'what to do at the bench' }
];

var stepState = { current: 'spec', ready: { spec: true } };

function stepPane(id) { return $('step-' + id); }

/** Mark a step as having something to show. */
function stepReady(id, yes) {
  stepState.ready[id] = yes !== false;
  paintSteps();
}

function showStep(id) {
  if (!stepState.ready[id]) return;
  stepState.current = id;
  STEPS.forEach(function (st) {
    var pane = stepPane(st.id);
    if (pane) pane.hidden = st.id !== id;
  });
  paintSteps();
  var panel = $('panel-design');
  if (panel && panel.scrollIntoView) {
    panel.scrollIntoView({ block: 'start', behavior: 'auto' });
  }
}

function paintSteps() {
  var nav = $('steps');
  if (!nav) return;
  clear(nav);
  /* A step with nothing to show is not a step. Greying them out left the
     rail cluttered with places you cannot go; the sequence is clearer
     when it only contains what exists. */
  var live = STEPS.filter(function (st) { return !!stepState.ready[st.id]; });
  /* Numbered by where they are in the rail, not by their place in the
     full sequence: hiding the empty ones left it reading 1, 3, 4, 5. */
  live.forEach(function (st, i) {
    var active = stepState.current === st.id;
    var btn = h('button', {
      type: 'button',
      class: 'step' + (active ? ' active' : ''),
      'aria-current': active ? 'step' : null,
      title: st.hint,
      onclick: function () { showStep(st.id); }
    },
      h('span', { class: 'step-n' }, String(i + 1)),
      h('span', { class: 'step-label' }, st.label),
      h('span', { class: 'step-hint' }, st.hint));
    nav.appendChild(btn);
    if (i < live.length - 1) {
      nav.appendChild(h('span', { class: 'step-sep', 'aria-hidden': 'true' }, ''));
    }
  });
}

/** Put a built section into a step, and light the step up. */
function fillStep(id, nodes) {
  var pane = stepPane(id);
  if (!pane) return false;
  clear(pane);
  var list = [].concat(nodes).filter(Boolean);
  if (!list.length) { stepReady(id, false); return false; }
  list.forEach(function (n) { pane.appendChild(n); });
  stepReady(id, true);
  return true;
}

/** Clear every step but the specification — a new run starts over. */
function resetSteps() {
  STEPS.forEach(function (st) {
    if (st.id === 'spec') return;
    var pane = stepPane(st.id);
    if (pane) clear(pane);
    stepState.ready[st.id] = false;
  });
  showStep('spec');
}

/** A fingerprint of the specification, to tell a stale plan from a live one. */
function specSignature() {
  return JSON.stringify((state.modules || []).map(function (m) {
    return [m.type, m.name, m.part_id, m.mode];
  })) + '|' + (($('spec-backbone') || {}).value || '');
}

function renderScan() {
  var scan = state.scan;
  if (!scan) return;
  /* Scanning used to clear the rail, so going back for another look at
     the parents threw away the routes - six steps became three and the
     way forward vanished. The plan only goes stale when the thing it was
     built from changes. */
  if (state.planSignature && state.planSignature !== specSignature()) {
    resetSteps();
  }

  var card = h('div', { class: 'card' },
    h('div', { class: 'card-head' },
      h('h2', null, scan.construct.name || '(unnamed construct)'),
      h('span', { class: 'muted small' },
        scan.backbone
          ? 'would start from ' + scan.backbone.name +
            ' (' + bp(scan.backbone.length) + ')'
          : 'no backbone chosen yet')));
  card.appendChild(h('p', { class: 'muted small' },
    'Nothing has been designed yet. Settle the parents below, then ' +
    'generate the protocol.'));
  if ((scan.modules || []).length) {
    card.appendChild(table([
      { label: 'Type', cell: function (m) { return typeTag(m.type); } },
      { label: 'Requested', key: 'requested' },
      { label: 'Resolved part', key: 'part_name' },
      { label: 'Length', numeric: true,
        cell: function (m) { return bp(m.length); } },
      { label: 'Found in', numeric: true,
        cell: function (m) { return (m.sources || []).length; } },
      { label: 'Notes', cell: function (m) { return warnCell(m.warnings); } }
    ], scan.modules));
  }
  fillStep('scan', [selectionSummary(scan), card, renderCandidates(scan),
                    nextFromScan()]);
  showStep('scan');
}

/** The move out of Routes: the protocol for the route now chosen. */
function nextFromRoutes() {
  return h('div', { class: 'card step-next' },
    h('button', {
      type: 'button', class: 'primary',
      onclick: function () { showStep('scheme'); }
    }, 'Move on to Protocol \u2192'),
    h('span', { class: 'muted small' },
      'The numbered protocol, the primers to order and the checks.'));
}

/** The move out of Scan source: straight to the routes. */
function nextFromScan() {
  var haveRoutes = !!(state.plan && (state.plan.routes || []).length);
  return h('div', { class: 'card step-next' },
    haveRoutes
      ? h('button', {
          type: 'button', class: 'primary',
          onclick: function () { showStep('routes'); }
        }, 'Move on to Routes \u2192')
      : h('button', {
          type: 'button', class: 'primary',
          title: 'Design with the parents chosen here and rank the ways '
               + 'to build it',
          onclick: function () { generatePlan(); }
        }, 'Move on to Routes \u2192'),
    h('span', { class: 'muted small' },
      haveRoutes
        ? 'Already designed \u2014 this goes back to the ranked routes.'
        : 'Runs the planner with whatever is chosen above.'));
}

/* ------------------------------------------------------------------ */
/* Step 5: the choices currently in force                              */
/* ------------------------------------------------------------------ */
/* Candidate parents is a long list to read; this is the short answer to
   "so what is it going to use?", which is the thing worth checking before
   spending a protocol on it. */

function selectionSummary(data) {
  if (!data) return null;
  var rows = [];
  var scheme = data.scheme || {};
  var rec = (scheme.recipient || {}).name || (data.backbone || {}).name;
  if (rec) rows.push(['Backbone / recipient', rec, pinnedParents.recipient]);
  var don = (scheme.donor || {}).name;
  if (don) rows.push(['Donor', don, pinnedParents.donor]);

  var card = h('div', { class: 'card' },
    h('div', { class: 'card-head' }, h('h2', null, 'Selected parents')),
    h('p', { class: 'muted small' },
      'What the plan will use. A value marked ' +
      '“pinned” was chosen by hand on the previous step; the rest were ' +
      'picked by the ranking and can be overridden there.'));

  if (!rows.length && !Object.keys(pinnedSources).length) {
    card.appendChild(h('p', { class: 'muted' },
      'Nothing chosen yet — the planner will pick for you.'));
    return card;
  }
  rows.forEach(function (r) {
    card.appendChild(h('div', { class: 'sel-row' },
      h('span', { class: 'sel-what' }, r[0]),
      h('span', { class: 'sel-name' }, r[1]),
      r[2] ? h('span', { class: 'pill pill-ok' }, 'pinned') : null));
  });

  var pins = Object.keys(pinnedSources || {});
  if (pins.length) {
    card.appendChild(h('h4', null, 'Pinned part sources'));
    pins.forEach(function (k) {
      card.appendChild(h('div', { class: 'sel-row' },
        h('span', { class: 'sel-what' }, 'part #' + k),
        h('span', { class: 'sel-name' }, 'plasmid #' + pinnedSources[k]),
        h('span', { class: 'pill pill-ok' }, 'pinned')));
    });
  }

  return card;
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

/** Say it before the protocol: the lab may already have this.

   Two weeks of cloning to arrive at a plasmid that is in the freezer
   is the most expensive thing this tool can fail to mention, and the
   inventory is the only place that knows. */
function alreadyMade(plan) {
  var hits = plan.already_made || [];
  if (!hits.length) return null;
  var made = hits.filter(function (h) {
    return h.exact || h.coverage >= 0.98;
  });
  var close = hits.filter(function (h) { return made.indexOf(h) < 0; });

  var card = h('div', { class: 'card' + (made.length ? ' card-alarm' : '') },
    h('div', { class: 'card-head' },
      h('h2', null, made.length ? 'The lab already has this construct'
                                : 'Close relatives already in the inventory'),
      h('span', { class: 'muted small' },
        made.length ? 'check one before cloning it again'
                    : 'different constructs, but worth a look')));

  function row(x, strong) {
    return h('div', { class: 'made-row' },
      h('span', { class: strong ? 'made-name' : 'made-name muted' }, x.name),
      h('span', { class: 'muted small' },
        bp(x.length) + ' · ' + (x.owner || 'unknown owner')
        + (x.verified ? ' · ' + String(x.verified).replace(/_/g, ' ') : '')),
      h('span', { class: 'muted small' },
        x.exact ? 'the same molecule, base for base'
                : Math.round(x.coverage * 100) + '% shared'));
  }
  made.forEach(function (x) { card.appendChild(row(x, true)); });
  if (made.length) {
    card.appendChild(h('p', { class: 'muted small' },
      'Sequence one of these before building it again. The plan below '
      + 'is only worth running if it turns out not to be what the name '
      + 'says.'));
  }
  if (close.length) {
    if (made.length) {
      card.appendChild(h('p', { class: 'muted small' },
        'Close relatives, for comparison:'));
    }
    close.forEach(function (x) { card.appendChild(row(x, false)); });
  }
  return card;
}

function renderPlan() {
  resetSteps();
  var plan = state.plan;
  if (!plan) return;
  var scanCards = [alreadyMade(plan), selectionSummary(plan), scanRecap(plan),
                   renderCandidates(plan), nextFromScan()].filter(Boolean);

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
  /* --- 2. the molecule that was asked for -------------------------- */
  fillStep('target', [summary, renderConstruct(c),
                      plan.vector ? renderVector(plan.vector, plan) : null]);

  /* --- 3. where the parts come from, and what to start from -------- */
  fillStep('scan', scanCards);

  /* --- 6. how to build it ------------------------------------------ */
  var route = selectedRoute();
  fillStep('routes', [renderRoutes(plan.routes || []),
                      renderCloningFlow(route, plan),
                      nextFromRoutes()]);

  /* --- 7. the plan at the bench, and whether it holds up ------------ */
  /* schemeCard drew the two parents and the product - which is what
     "How this route builds it" now does on Routes, for every route
     rather than only a swap. Two drawings of one thing, on two steps,
     is one too many. */
  fillStep('scheme', [route ? renderRouteDetail(route) : null,
                      renderSimulation(route, plan),
                      plan.verification
                        ? renderVerification(plan.verification) : null,
                      renderValidation(plan.validation || {}),
                      renderExports()]);

  state.planSignature = specSignature();

  /* Land on the routes: with a plan in hand that is the decision to
     make, and every earlier step stays one click away. */
  showStep(stepState.ready.routes ? 'routes' : 'target');
}

/** Step 3 once a plan exists: where each module resolved to. */
function scanRecap(plan) {
  if (!(plan.modules || []).length) return null;
  var card = h('div', { class: 'card' },
    h('div', { class: 'card-head' }, h('h2', null, 'Sources found'),
      h('span', { class: 'muted small' },
        plan.modules.length + ' module(s) resolved')));
  card.appendChild(table([
    { label: 'Type', cell: function (m) { return typeTag(m.type); } },
    { label: 'Requested', key: 'requested' },
    { label: 'Resolved part', key: 'part_name' },
    { label: 'Length', numeric: true, cell: function (m) { return bp(m.length); } },
    { label: 'Found in', numeric: true,
      cell: function (m) { return (m.sources || []).length; } },
    { label: 'Notes', cell: function (m) { return warnCell(m.warnings); } }
  ], plan.modules));
  return card;
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
/**
 * A map with a circular/linear toggle, defaulting to the molecule's own
 * topology.
 *
 * A plasmid is a circle and a linear track misrepresents it - a cassette
 * that wraps the origin looks broken on one and obvious on the other -
 * but the track is still the easier thing to read coordinates off. So
 * both, with the honest one first.
 */
function mapWithToggle(length, items, opts) {
  opts = opts || {};
  var host = h('div');
  var mapHost = h('div', { class: 'map-host' });
  var mode = opts.circular ? 'circular' : 'linear';

  function draw() {
    clear(mapHost);
    mapHost.appendChild(mode === 'circular'
      ? circularMap(length, items, {
          label: opts.label, markers: opts.markers,
          highlight: opts.highlight, highlightKind: opts.highlightKind })
      : featureMap(length, items, {
          label: opts.label, circular: !!opts.circular,
          markers: opts.markers }));
  }

  if (opts.circular) {
    var toggle = h('div', { class: 'map-toggle' });
    ['circular', 'linear'].forEach(function (option) {
      toggle.appendChild(h('button', {
        type: 'button',
        class: option === mode ? 'active' : '',
        onclick: function () {
          mode = option;
          [].forEach.call(toggle.children, function (child) {
            child.className = child.textContent === mode ? 'active' : '';
          });
          draw();
        }
      }, option));
    });
    host.appendChild(toggle);
  }
  host.appendChild(mapHost);
  draw();
  return host;
}

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

  card.appendChild(mapWithToggle(c.length, layout, {
    label: c.name, markers: c.codon_markers, circular: !!c.circular }));

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
/* ------------------------------------------------------------------ */
/* Plan: choosing the parents before anything is designed              */
/* ------------------------------------------------------------------ */

/**
 * The request read as an architecture, with the parents on offer.
 *
 * The planner picks the plasmid that already carries the longest leading
 * run, and takes whatever is left from as few donors as it can. Both are
 * only ever the best guess from what is catalogued: somebody may know that
 * one of these is on their own bench, or sequenced, or simply the one the
 * lab uses. Choosing here re-plans from that choice.
 */

/* ------------------------------------------------------------------ */
/* Scan source: the plasmid you are pointing at                        */
/* ------------------------------------------------------------------ */
/* Choosing a parent off a list of names is choosing blind. This draws
   whichever one is selected, with the parts the request asked for picked
   out of it, so "already has CAG" is something you can see rather than
   something you take on trust. */

var PLASMID_CACHE = {};

function fetchPlasmid(id) {
  if (!PLASMID_CACHE[id]) {
    PLASMID_CACHE[id] = api('/api/plasmid/' + id)
      .catch(function (e) { delete PLASMID_CACHE[id]; throw e; });
  }
  return PLASMID_CACHE[id];
}

/** Does this annotation answer one of the names the request asked for? */
function wantedBy(ann, wants) {
  var nm = String(ann.name || '').toLowerCase();
  if (!nm) return false;
  return wants.some(function (w) {
    w = String(w || '').toLowerCase();
    return w && (nm === w || nm.indexOf(w) >= 0 || w.indexOf(nm) >= 0);
  });
}

function showCandidateMap(host, plasmidId, wants, role) {
  clear(host);
  if (!plasmidId) {
    host.appendChild(h('p', { class: 'muted small' },
      'Pick a plasmid to see what is in it.'));
    return;
  }
  host.appendChild(h('p', { class: 'muted small' }, 'Loading\u2026'));
  fetchPlasmid(plasmidId).then(function (p) {
    clear(host);
    var anns = (p.annotations || []).filter(function (a) {
      return a && typeof a.start === 'number' && a.end > a.start;
    });
    var hits = anns.filter(function (a) { return wantedBy(a, wants); });
    host.appendChild(h('div', { class: 'card-head' },
      h('h4', null, p.name || ('plasmid #' + plasmidId)),
      h('span', { class: 'muted small' },
        bp(p.seq_length) + ' \u00b7 ' + (role || 'parent'))));
    host.appendChild(mapWithToggle(p.seq_length, anns.map(function (a) {
      return {
        feature_type: a.feature_type, name: a.name,
        start: a.start, end: a.end, strand: a.strand,
        /* The thing being looked for, lit up against the rest. */
        colour: wantedBy(a, wants) ? '#f0a32a' : null
      };
    }), { label: p.name, circular: !!p.circular }));
    if (hits.length) {
      host.appendChild(h('p', { class: 'muted small' },
        'In amber: ' + hits.map(function (a) {
          return a.name + ' (' + (a.start + 1) + '\u2013' + a.end + ')';
        }).join(', ') + '.'));
    } else {
      host.appendChild(h('p', { class: 'muted small' },
        'None of the requested parts is annotated in this plasmid by name '
        + '\u2014 it may still carry them under another label.'));
    }
  }).catch(function (e) {
    clear(host);
    host.appendChild(h('p', { class: 'inline-warn' },
      'could not load that plasmid: ' + e.message));
  });
}

function renderCandidates(plan) {
  var cand = plan.candidates || {};
  var recipients = cand.recipients || [];
  var donors = cand.donors || [];
  if (!recipients.length && !donors.length) return null;

  var chosenBackbone = (plan.backbone && plan.backbone.id) || null;

  var card = h('div', { class: 'card cand-card' },
    h('div', { class: 'card-head' },
      h('h2', null, 'Candidate parents'),
      h('span', { class: 'muted small' },
        'start from the one that already has the most of this')));
  card.appendChild(h('p', { class: 'muted small' },
    'Read as an architecture: ' + (cand.modules || []).join(' – ') + '.'));

  var wants = (cand.modules || []).slice();
  /* One panel per parent rather than one shared one. A two-donor plan
     has three plasmids in it, and a single panel could only ever show
     the last thing clicked - which is the one you are least likely to
     be comparing against. */
  var maps = h('aside', { class: 'cand-maps' });
  var panels = {};
  function panelFor(key, role) {
    if (!panels[key]) {
      panels[key] = h('div', { class: 'cand-map card' });
      maps.appendChild(panels[key]);
    }
    return panels[key];
  }
  var lists = h('div', { class: 'cand-lists' });
  card.appendChild(h('div', { class: 'cand-split' }, lists, maps));

  /* Draw what is already chosen, so nothing is blank before a click. */
  showCandidateMap(panelFor('recipient'), chosenBackbone, wants, 'recipient');
  (cand.donors || []).forEach(function (d, di) {
    var ids = partIdsFor(plan, d.covers || []);
    var picked = (d.options || []).filter(function (o) {
      return ids.length && ids.every(function (id) {
        return pinnedSources[id] === o.plasmid_id;
      });
    })[0] || (d.options || [])[0];
    if (picked) {
      showCandidateMap(panelFor('donor' + di), picked.plasmid_id,
                       d.covers || [],
                       'donor ' + (di + 1) + ' \u2014 ' +
                       (d.covers || []).join(' + '));
    }
  });

  if (recipients.length) {
    lists.appendChild(h('h4', null, 'Start from'));
    /* Four lines of prose explaining one column, above the table it
       explains. The score itself already carries the explanation as a
       tooltip on every cell; this is for whoever wants it spelled out. */
    lists.appendChild(disclosure('how the score is read',
      'The score counts how many of these a plasmid really carries, in '
      + 'order from the start: 2/2 has both and only the vector has to be '
      + 'opened, 1/2 has the first and the rest must be built. A plasmid '
      + 'marked variant holds something close rather than the part itself, '
      + 'which is a reason to start from it and not a reason to take the '
      + 'part out of it.'));
    /* When every candidate carries the same thing the column separates
       nothing, and six identical cells read as data. */
    var varies = recipients.some(function (r) {
      return (r.carries || []).join(',')
           !== (recipients[0].carries || []).join(',');
    });
    /* Six rows of radio buttons to pick one thing is a list pretending
       to be a control. The chosen one shows as a row - with its score,
       its owner and its map - and the rest live in the select until
       somebody wants to change it. */
    lists.appendChild(candidatePicker(
      'Start from', recipients, chosenBackbone,
      function (r) {
        return r.exact_run + '/' + (cand.modules || []).length + '  '
             + r.name + '  ' + bp(r.length)
             + (r.owner ? ' - ' + r.owner : '')
             + (r.exact ? '' : '  (variant)');
      },
      function (r) {
        $('spec-backbone').value = String(r.plasmid_id);
        showCandidateMap(panelFor('recipient'), r.plasmid_id, wants,
                         'recipient');
        replan();
      }));
    var list = h('div', { class: 'cand-list cand-rest'
                          + (varies ? '' : ' cand-nohas') });
    list.appendChild(h('div', { class: 'cand-row cand-head' },
      h('span', null, ''),
      h('span', { class: 'cand-score' }, 'has'),
      h('span', null, 'plasmid'),
      h('span', null, 'size · owner'),
      varies ? h('span', null, 'already has')
             : h('span', { class: 'muted small cand-head-note',
                           title: 'every candidate here already has '
                                + ((recipients[0].carries || []).join(', ')
                                   || 'none of the parts') },
                 'all carry '
                 + ((recipients[0].carries || []).join(', ') || 'nothing'))));
    recipients.forEach(function (r) {
      var on = r.plasmid_id === chosenBackbone;
      var row = h('label', { class: 'cand-row' + (on ? ' cand-on' : '') },
        h('input', {
          type: 'radio', name: 'cand-recipient',
          checked: on ? 'checked' : null,
          onchange: function () {
            $('spec-backbone').value = String(r.plasmid_id);
            showCandidateMap(panelFor('recipient'), r.plasmid_id, wants,
                             'recipient');
            replan();
          }
        }),
        h('span', {
          class: 'cand-score' + (r.exact ? '' : ' cand-score-soft'),
          title: 'matches the first ' + r.run + ' of ' +
                 (cand.modules || []).length + ' modules in order; ' +
                 r.exact_run + ' of those ' +
                 (r.exact_run === 1 ? 'is the part itself'
                                    : 'are the parts themselves') +
                 ', the rest are relatives that would be replaced'
        }, r.exact_run + '/' + (cand.modules || []).length),
        h('span', { class: 'cand-name' }, r.name,
          r.exact ? null : h('span', { class: 'pill pill-soft',
                                       title: 'a relative, not the part ' +
                                              'itself - this is what gets ' +
                                              'swapped out' }, 'variant')),
        h('span', { class: 'muted small' },
          bp(r.length) + (r.owner ? ' · ' + r.owner : '') +
          (r.has_cassette ? ' · has lox cassette' : '')),
        varies ? h('span', { class: 'muted small cand-has' },
                     'already has ' + (r.carries || []).join(', '))
               : h('span', { class: 'cand-has' }, ''));
      list.appendChild(row);
    });
    lists.appendChild(foldAway(
      'compare all ' + recipients.length + ' side by side', list));
  }

  if (donors.length) {
    var n = cand.n_donors || 0;
    lists.appendChild(h('h4', null,
      'Then take the rest from ' + n + ' donor' + (n === 1 ? '' : 's')));
    donors.forEach(function (d, di) {
      var box = h('div', { class: 'cand-donor' },
        h('div', null, h('strong', null, (d.covers || []).join(' – '))));
      if (d.note) {
        box.appendChild(h('p', { class: 'issue issue-warning' }, d.note));
      }
      var opts = d.options || [];
      var ids0 = partIdsFor(plan, d.covers || []);
      var picked = opts.filter(function (o) {
        return ids0.length && ids0.every(function (id) {
          return pinnedSources[id] === o.plasmid_id;
        });
      })[0];
      box.appendChild(candidatePicker(
        'use', opts, picked ? picked.plasmid_id : null,
        function (o) {
          return o.name + '  ' + bp(o.length)
               + (o.owner ? ' - ' + o.owner : '')
               + (o.exact ? '  (exact)' : '  (variant)');
        },
        function (o) {
          ids0.forEach(function (id) { pinnedSources[id] = o.plasmid_id; });
          showCandidateMap(panelFor('donor' + di), o.plasmid_id,
                           d.covers || [],
                           'donor ' + (di + 1) + ' \u2014 '
                           + (d.covers || []).join(' + '));
          replan();
        }));
      var rest = h('div', { class: 'cand-list cand-rest' });
      opts.forEach(function (o) {
        var ids = partIdsFor(plan, d.covers || []);
        var on = ids.length && ids.every(function (id) {
          return pinnedSources[id] === o.plasmid_id;
        });
        rest.appendChild(h('label', { class: 'cand-row' + (on ? ' cand-on' : '') },
          h('input', {
            type: 'radio', name: 'cand-donor-' + (d.covers || []).join('-'),
            checked: on ? 'checked' : null,
            onchange: function () {
              ids.forEach(function (id) { pinnedSources[id] = o.plasmid_id; });
              /* Show the donor being pointed at, with the part it is
                 wanted for picked out of it. */
              showCandidateMap(panelFor('donor' + di), o.plasmid_id,
                               d.covers || [],
                               'donor ' + (di + 1) + ' \u2014 ' +
                               (d.covers || []).join(' + '));
              replan();
            }
          }),
          /* keeps a donor row on the same columns as a recipient row */
          h('span', { class: 'cand-score' }, ''),
          h('span', { class: 'cand-name' }, o.name,
            h('span', {
              class: 'pill ' + (o.exact ? 'pill-ok' : 'pill-soft'),
              title: o.exact
                ? 'carries the part that was asked for'
                : 'carries a relative; taking it would clone that instead'
            }, o.exact ? 'exact' : 'variant')),
          h('span', { class: 'muted small' },
            bp(o.length) + (o.owner ? ' · ' + o.owner : '')),
          h('span', { class: 'muted small cand-has' },
            (o.carries || []).length ? 'has ' + o.carries.join(', ') : '')));
      });
      if (opts.length > 1) {
        box.appendChild(foldAway(
          'compare all ' + opts.length + ' side by side', rest));
      }
      lists.appendChild(box);
    });
  }

  if (Object.keys(pinnedSources).length) {
    lists.appendChild(h('p', { class: 'muted small' },
      h('button', {
        type: 'button', class: 'secondary',
        onclick: function () { pinnedSources = {}; replan(); }
      }, 'Clear donor choices')));
  }
  return card;
}

/**
 * Redo whichever half the user is looking at.
 *
 * Changing a parent while only the scan is on screen should not spend a
 * minute designing primers nobody asked for yet; changing one after the
 * protocol exists should rebuild it.
 */
function replan() {
  if (state.plan) generatePlan(); else scanSources();
}

/** The part ids behind a donor segment, so pinning one pins those parts. */
function partIdsFor(plan, names) {
  var wanted = {};
  names.forEach(function (n) { wanted[n] = true; });
  return (plan.modules || []).filter(function (m) {
    return wanted[m.part_name] && m.part_id;
  }).map(function (m) { return m.part_id; });
}

/* ------------------------------------------------------------------ */
/* Plan: the protocol run on paper                                     */
/* ------------------------------------------------------------------ */

/** A colour per fragment, so the same piece is the same colour everywhere. */
function fragmentHue(i) {
  return 'hsl(' + ((i * 67) % 360) + ' 62% 46%)';
}

/**
 * The simulation, drawn rather than only tabulated.
 *
 * Each source plasmid is shown as the ring it is, with the stretch that
 * gets amplified marked on it, and the construct is shown with the same
 * colours, so "this piece comes from there and lands here" is one glance
 * rather than three cross-references.
 */
function renderSimulation(route, plan) {
  var sim = route && route.simulation;
  if (!sim || !sim.simulated) return null;

  var ok = sim.passed;
  var card = h('div', { class: 'card sim-card' },
    h('div', { class: 'card-head' },
      h('h2', null, 'Cloning simulation'),
      h('span', { class: 'pill ' + (ok ? 'pill-ok' : 'pill-bad') },
        ok ? (sim.identical ? 'builds the designed sequence'
                            : 'all parts present') : 'does not build it')));

  card.appendChild(h('p', { class: 'muted small' },
    'The primers below were amplified off the named templates and the ' +
    'products joined on their shared ends. What follows is the molecule ' +
    'those steps make, not the one the design describes.'));
  card.appendChild(h('p', null, h('strong', null, sim.summary || '')));

  (sim.problems || []).forEach(function (p) {
    card.appendChild(h('p', { class: 'issue issue-error' }, p));
  });
  (sim.warnings || []).forEach(function (w) {
    card.appendChild(h('p', { class: 'issue issue-warning' }, w));
  });

  if ((sim.modules || []).length) {
    card.appendChild(h('h4', null, 'Requested parts, looked for in the product'));
    card.appendChild(table([
      { label: 'Part', key: 'name' },
      { label: 'Length', numeric: true,
        cell: function (m) { return bp(m.bp); } },
      { label: 'In the product', cell: function (m) {
          return h('span', { class: 'pill ' + (m.present ? 'pill-ok' : 'pill-bad') },
                   m.present ? 'yes' : 'no'); } },
      { label: 'Covered', numeric: true, cell: function (m) {
          return Math.round((m.coverage || 0) * 100) + '%'; } }
    ], sim.modules));
  }

  var sources = sim.sources || [];
  if (!sources.length) return card;

  /* --- the schematic ------------------------------------------------ */
  card.appendChild(h('h4', null, 'Where each fragment is cut from'));
  card.appendChild(h('p', { class: 'muted small' },
    'The marked arc is what the primers amplify. Its colour is the same ' +
    'on the assembled construct below.'));

  var grid = h('div', { class: 'sim-grid' });
  sources.forEach(function (src, i) {
    var hue = fragmentHue(i);
    var box = h('div', { class: 'sim-source' });
    box.appendChild(h('div', { class: 'sim-source-head' },
      h('span', { class: 'sim-swatch', style: 'background:' + hue }),
      h('strong', null, src.name),
      h('span', { class: 'muted small' }, bp(src.bp))));
    box.appendChild(h('div', { class: 'muted small' },
      'from ' + src.template +
      (src.template_length ? ' (' + bp(src.template_length) + ')' : '')));
    var hi = null;
    if (typeof src.from_start === 'number' && typeof src.from_end === 'number') {
      hi = { start: src.from_start, end: src.from_end,
             label: src.name, length: src.bp };
    }
    box.appendChild(circularMap(src.template_length,
      (src.features || []).map(function (f) {
        return { type: f.feature_type, name: f.name, start: f.start,
                 end: f.end, strand: f.strand };
      }), { label: src.template, highlight: hi, highlightKind: 'take',
            circular: src.template_circular !== false }));
    grid.appendChild(box);
  });
  card.appendChild(grid);

  /* --- and where they land ------------------------------------------ */
  var placed = sources.filter(function (s) {
    return typeof s.at_start === 'number';
  });
  if (placed.length && sim.product_bp) {
    card.appendChild(h('h4', null, 'How they go together'));
    card.appendChild(h('p', { class: 'muted small' },
      'The simulated product, coloured by which reaction each stretch ' +
      'came out of. An arc that wraps past the top crosses position 1.'));
    card.appendChild(featureMap(sim.product_bp, placed.map(function (s, i) {
      return { type: 'fragment', name: s.name, start: s.at_start,
               end: s.at_end < s.at_start ? sim.product_bp : s.at_end,
               strand: s.flipped ? -1 : 1,
               colour: fragmentHue(sources.indexOf(s)) };
    }), { label: 'assembled product', circular: true }));

    card.appendChild(table([
      { label: 'Fragment', cell: function (s) {
          return h('span', null,
            h('span', { class: 'sim-swatch',
                        style: 'background:' + fragmentHue(sources.indexOf(s)) }),
            ' ' + s.name); } },
      { label: 'Cut from', cell: function (s) {
          return s.template + '  ' + bp(s.from_start) + '-' + bp(s.from_end); } },
      { label: 'Length', numeric: true, cell: function (s) { return bp(s.bp); } },
      { label: 'Lands at', cell: function (s) {
          return bp(s.at_start) + '-' + bp(s.at_end) +
                 (s.flipped ? ' (reversed)' : ''); } }
    ], placed));
  }
  return card;
}

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


/* ------------------------------------------------------------------ */
/* Routes: the cloning drawn the way a paper would draw it             */
/* ------------------------------------------------------------------ */
/* The first version put a full plasmid map next to every parent, which
   draws every annotation the file happens to carry. The question a
   cloning figure answers is narrower: which plasmid, which piece out of
   it, by what reaction, and what joins to what. One lane per piece, one
   column per bench stage, and nothing else on the page. */

var FIG = {
  r: 56,               /* parent circle */
  pr: 68,              /* product circle */
  rowH: 196,
  colParent: 250,      /* width of the parent column, labels included */
  colStage: 96,
  colFrag: 168,
  bracket: 18,
  pad: 14
};

/* One colour per lane, carried from the parent through to the arc it
   becomes in the product, so "where did this part come from" is answered
   by looking rather than by reading. */
var LANE_COLOURS = ['#2f7bd4', '#d4572f', '#2f9e5e', '#8e4fc9', '#b58900'];

/** Drop anything sitting entirely inside something else kept.

   Real files carry nested annotations - the chicken beta-actin promoter
   inside CAG - and auto-annotation leaves junk inside them too: YD030
   has a 609 bp reverse-strand "Voltron2" buried in its CAG promoter
   that has nothing to do with the real Voltron-2 further round the
   plasmid. At figure scale a feature inside another feature is a detail
   of that feature, not a part of the construct. The full map on the
   Protocol step still shows everything. */
function outermost(items, startKey, endKey) {
  var sk = startKey || 'start', ek = endKey || 'end';
  return items.filter(function (a) {
    var s = a[sk], e = a[ek];
    if (typeof s !== 'number' || typeof e !== 'number' || e <= s) return true;
    return !items.some(function (b) {
      return b !== a && typeof b[sk] === 'number' && typeof b[ek] === 'number'
        && b[ek] > b[sk] && b[sk] <= s && e <= b[ek]
        && (b[ek] - b[sk]) > (e - s);
    });
  });
}

/* One colour per junction, so the two halves of a join match by eye.
   Deliberately away from LANE_COLOURS: a lane colour says "this piece",
   a junction colour says "these two ends meet". */
var JUNCTION_COLOURS = ['#d9a60d', '#bd339a', '#0d9ea8', '#73798a'];
var JUNCTION_LETTERS = 'ABCDEF';

/** {laneIndex|side: [label, colour]} for every end that has a partner.

   Sticky ends pair by enzyme - a HindIII end only goes into a HindIII
   end - so the enzyme name is the junction's identity. Homology arms
   have no name of their own, so they are lettered by the order the
   pieces sit in the product: piece k's right end meets piece k+1's
   left end, all the way round the circle. */
function junctionColours(lanes) {
  var ids = [];
  var sticky = lanes.filter(function (l) { return l.endKind === 'sticky'; });
  if (sticky.length && sticky.length === lanes.length) {
    lanes.forEach(function (lane, i) {
      ['left', 'right'].forEach(function (side, k) {
        var name = (lane.ends || [])[k];
        if (name) ids.push([i, side, name]);
      });
    });
  } else {
    var placed = lanes.map(function (l, i) { return [i, l]; })
      .filter(function (t) {
        return typeof ((t[1].src || {}).at_start) === 'number';
      });
    if (placed.length !== lanes.length || lanes.length < 2) return {};
    placed.sort(function (a, b) {
      return a[1].src.at_start - b[1].src.at_start;
    });
    for (var k = 0; k < placed.length; k++) {
      var label = JUNCTION_LETTERS[k % JUNCTION_LETTERS.length];
      ids.push([placed[k][0], 'right', label]);
      ids.push([placed[(k + 1) % placed.length][0], 'left', label]);
    }
  }
  var order = [];
  ids.forEach(function (t) {
    if (order.indexOf(t[2]) < 0) order.push(t[2]);
  });
  var out = {};
  ids.forEach(function (t) {
    out[t[0] + '|' + t[1]] = [t[2],
      JUNCTION_COLOURS[order.indexOf(t[2]) % JUNCTION_COLOURS.length]];
  });
  return out;
}

/** The features worth naming on the assembled plasmid. A construct has
   one promoter and one payload worth labelling at this size; the full
   map on the Protocol step carries the rest. */
function productFeatures(layout, span, limit) {
  span = span || 1;
  var keep = (layout || []).filter(function (a) {
    return a && typeof a.start === 'number'
      && ['promoter', 'CDS'].indexOf(a.type) >= 0
      && (a.end - a.start) > span * 0.03;
  });
  keep = outermost(keep);
  keep.sort(function (a, b) { return (b.end - b.start) - (a.end - a.start); });
  return keep.slice(0, limit || 4);
}

/** The handful of features worth naming on a parent. */
function figFeatures(src) {
  var span = src.template_length || 1;
  var keep = (src.features || []).filter(function (a) {
    return a && typeof a.start === 'number' && a.end > a.start
      && (a.end - a.start) > span * 0.05
      && ['CDS', 'promoter', 'resistance', 'origin'].indexOf(a.feature_type) >= 0;
  });
  keep = outermost(keep);
  keep.sort(function (a, b) { return (b.end - b.start) - (a.end - a.start); });
  return keep.slice(0, 3);
}

/** Text just outside a circle, on whichever side of it the angle is.
 *
 *  Two cut sites a few hundred bases apart sit at nearly the same angle,
 *  and their names landed on top of each other - "HindIII" and "EcoRV"
 *  overprinted into an unreadable smudge. Each circle keeps a list of
 *  where it has already written, and a clashing label is nudged clear.
 */
function ringLabel(g, cx, cy, r, f, text, cls) {
  if (!text) return;
  var at = ((f % 1) + 1) % 1;
  var p = polar(cx, cy, r + 9, at);
  var right = at < 0.5;
  var x = p[0] + (right ? 3 : -3), y = p[1] + 3;
  /* Push a clashing label away from the circle, not across it: nudging
     everything downwards wrote the second label over the ring whenever
     both sat near twelve o'clock. */
  var step = (at < 0.25 || at > 0.75) ? -11 : 11;
  var used = g.__labels || (g.__labels = []);
  for (var guard = 0; guard < 6; guard++) {
    var clash = used.some(function (u) {
      return u.right === right && Math.abs(u.y - y) < 11
          && Math.abs(u.x - x) < 70;
    });
    if (!clash) break;
    y += step;
  }
  used.push({ x: x, y: y, right: right });
  var t = svg('text', {
    x: x, y: y, 'text-anchor': right ? 'start' : 'end',
    class: cls || 'fig-tag'
  });
  t.textContent = text;
  g.appendChild(t);
}

/** A primer, drawn as the little arrow a cloning figure uses. */
function figPrimer(g, cx, cy, r, f, forward, name) {
  var w = 0.022 * (forward ? 1 : -1);
  var a = polar(cx, cy, r + 7, f), b = polar(cx, cy, r + 7, f + w * 0.6);
  var tip = polar(cx, cy, r + 7, f + w);
  var s1 = polar(cx, cy, r + 3.5, f + w * 0.6);
  var s2 = polar(cx, cy, r + 10.5, f + w * 0.6);
  g.appendChild(svg('path', {
    d: 'M ' + a[0] + ' ' + a[1] + ' L ' + b[0] + ' ' + b[1],
    class: 'fig-primer', 'stroke-width': 2.4
  }));
  var head = svg('path', {
    d: 'M ' + s1[0] + ' ' + s1[1] + ' L ' + tip[0] + ' ' + tip[1]
       + ' L ' + s2[0] + ' ' + s2[1] + ' Z', class: 'fig-primer-head'
  });
  if (name) {
    var ti = svg('title', null);
    ti.textContent = name;
    head.appendChild(ti);
  }
  g.appendChild(head);
  /* F and R, not the whole primer name: the name belongs on the protocol
     and what the figure has to say is which end and which way. */
  ringLabel(g, cx, cy, r + 4, f + w * 0.5, forward ? 'F' : 'R',
            'fig-prim-tag');
}

/** An enzyme cut: a tick through the backbone with the enzyme named. */
function figCut(g, cx, cy, r, f, name) {
  var a = polar(cx, cy, r - 8, f), b = polar(cx, cy, r + 8, f);
  g.appendChild(svg('line', { x1: a[0], y1: a[1], x2: b[0], y2: b[1],
    class: 'fig-cut', 'stroke-width': 2.4 }));
  if (name) ringLabel(g, cx, cy, r + 5, f, name, 'fig-enz');
}

/** One parent plasmid: what it is, what comes out, and how. */
function figParent(lane, cx, cy) {
  var src = lane.src || {};
  var span = src.template_length || 1;
  var r = FIG.r;
  var g = svg('g', null);
  g.appendChild(svg('circle', { cx: cx, cy: cy, r: r, fill: 'none',
    class: 'fig-ring', 'stroke-width': 2.5 }));

  figFeatures(src).forEach(function (a) {
    var from = (a.start % span) / span, to = (a.end % span) / span;
    g.appendChild(svg('path', {
      d: arcPath(cx, cy, r, from, to), fill: 'none', 'stroke-width': 6,
      stroke: colorFor(a.feature_type), opacity: 0.4
    }));
    /* Only the promoter earns a label; the rest are context. */
    if (a.feature_type === 'promoter') {
      var mid = from + (((to - from) + 1) % 1) / 2;
      ringLabel(g, cx, cy, r, mid, shortName(a.name, 15));
    }
  });

  if (typeof src.from_start !== 'number'
      || typeof src.from_end !== 'number' || !src.template_length) {
    figCaption(g, cx, cy, src.template,
               src.template_length ? bp(src.template_length)
                                   : 'size unknown');
    return g;
  }
  var f0 = ((src.from_start || 0) % span) / span;
  var f1 = ((src.from_end || 0) % span) / span;
  g.appendChild(svg('path', {
    d: arcPath(cx, cy, r, f0, f1), fill: 'none',
    'stroke-width': 8, stroke: lane.colour
  }));
  var cuts = lane.cutsParent || [];
  if (cuts.length) {
    var second = cuts[1] || cuts[0];
    /* Two sites a few bases apart are one tick at this scale, and two
       names stacked on it read as a smudge rather than as two enzymes. */
    if (Math.min((f0 - f1 + 1) % 1, (f1 - f0 + 1) % 1) < 0.02) {
      figCut(g, cx, cy, r, f0, cuts[0] + '/' + second);
    } else {
      figCut(g, cx, cy, r, f0, cuts[0]);
      figCut(g, cx, cy, r, f1, second);
    }
  } else {
    var pn = lane.primers || [];
    figPrimer(g, cx, cy, r, f0, true, pn[0]);
    figPrimer(g, cx, cy, r, f1, false, pn[1]);
  }

  figCaption(g, cx, cy, src.template, bp(src.template_length));
  return g;
}

/** A piece nobody has a plasmid for: it gets ordered. */
function figOrder(lane, cx, cy) {
  var g = svg('g', null);
  var w = 96, hgt = 52;
  g.appendChild(svg('rect', { x: cx - w / 2, y: cy - hgt / 2, width: w,
    height: hgt, rx: 6, class: 'fig-order' }));
  var t = svg('text', { x: cx, y: cy - 2, 'text-anchor': 'middle',
    class: 'fig-name' });
  t.textContent = 'synthesise';
  g.appendChild(t);
  var sub = svg('text', { x: cx, y: cy + 13, 'text-anchor': 'middle',
    class: 'fig-sub' });
  sub.textContent = bp(lane.src && lane.src.bp) + ' gBlock';
  g.appendChild(sub);
  return g;
}

function figStageArrow(x, y, w, verb, sub, colour) {
  var g = svg('g', null);
  var x1 = x + 10, x2 = x + w - 10;
  g.appendChild(svg('line', { x1: x1, y1: y, x2: x2, y2: y,
    class: 'fig-arrow', 'stroke-width': 1.8,
    'marker-end': 'url(#figHead)' }));
  var t = svg('text', { x: (x1 + x2) / 2, y: y - 8, 'text-anchor': 'middle',
    class: 'fig-verb' });
  t.textContent = verb;
  g.appendChild(t);
  if (sub) {
    var s2 = svg('text', { x: (x1 + x2) / 2, y: y + 15,
      'text-anchor': 'middle', class: 'fig-sub' });
    s2.textContent = sub;
    g.appendChild(s2);
  }
  return g;
}

/** The linear piece in the tube, with whatever its ends are. */
function figFragment(lane, index, joins, x, y, w) {
  var g = svg('g', null);
  g.appendChild(svg('rect', { x: x, y: y - 6, width: w, height: 12, rx: 2,
    fill: lane.colour }));
  /* Both halves of a join wear the same colour and the same letter, so
     "which end goes into which" is a thing you see rather than work
     out. Two grey caps both labelled "overlap" said only that an
     overlap existed. */
  ['left', 'right'].forEach(function (side, i) {
    var pair = joins[index + '|' + side];
    var fallback = (lane.ends || [])[i];
    if (!pair && !fallback) return;
    var label = pair ? pair[0] : fallback;
    var colour = pair ? pair[1] : null;
    var ex = i === 0 ? x : x + w;
    g.appendChild(svg('rect', { x: ex - 3.5, y: y - 10, width: 7, height: 20,
      rx: 1.5, fill: colour || null,
      class: colour ? null : (lane.endKind === 'sticky' ? 'fig-end-cut'
                                                        : 'fig-end-hom') }));
    var t = svg('text', { x: ex, y: y + 22, 'text-anchor': 'middle',
      class: 'fig-end-label', fill: colour || null });
    t.textContent = label;
    g.appendChild(t);
  });
  /* The piece is named here rather than beside its parent: around a
     circle the name of a near-whole-plasmid arc has no sensible angle
     to sit at, and it collided with the plasmid's own name. */
  var nm = svg('text', { x: x + w / 2, y: y - 14, 'text-anchor': 'middle',
    class: 'fig-piece' });
  nm.textContent = shortName(lane.piece, 24);
  g.appendChild(nm);
  var lab = svg('text', { x: x + w / 2, y: y + 22, 'text-anchor': 'middle',
    class: 'fig-sub' });
  lab.textContent = bp(lane.src && lane.src.bp);
  /* The ends are named at the ends; the size sits between them. */
  g.appendChild(lab);
  return g;
}

/** The new plasmid: which piece is which, and where they join.

   It used to carry two rings - the lane contributions outside and the
   construct's own promoter and CDS inside - and the labels pointed at
   the inner one. So the arc labelled ASAP7y was not the arc the
   fragment on the left was drawn in, and the two colours contradicted
   each other. One ring now, the pieces, labelled as the pieces they
   are, with the junctions marked where they actually fall. */
function figProduct(lanes, plan, sim, cx, cy, joins) {
  var c = (plan && plan.construct) || {};
  var span = sim.product_bp || c.length || 1;
  var r = FIG.pr;
  var g = svg('g', null);
  g.appendChild(svg('circle', { cx: cx, cy: cy, r: r, fill: 'none',
    class: 'fig-ring', 'stroke-width': 3 }));

  var placed = lanes.filter(function (l) {
    return typeof ((l.src || {}).at_start) === 'number';
  });
  placed.forEach(function (lane, i) {
    var src = lane.src;
    var from = src.at_start / span, to = src.at_end / span;
    g.appendChild(svg('path', {
      d: arcPath(cx, cy, r, from, to), fill: 'none',
      'stroke-width': 9, stroke: lane.colour
    }));
    var mid = from + ((((to - from) % 1) + 1) % 1) / 2;
    ringLabel(g, cx, cy, r, mid, shortName(lane.piece, 16),
              'fig-tag fig-tag-strong');
  });

  /* The overhangs, on the product as well as on the pieces: the join is
     a place on the finished molecule, not only a property of a tube. */
  var order = lanes.map(function (l, i) { return i; }).filter(function (i) {
    return typeof ((lanes[i].src || {}).at_start) === 'number';
  });
  order.forEach(function (i) {
    var pair = joins && joins[i + '|right'];
    if (!pair) return;
    var at = (lanes[i].src.at_end % span) / span;
    var a = polar(cx, cy, r - 7, at), b = polar(cx, cy, r + 7, at);
    g.appendChild(svg('line', { x1: a[0], y1: a[1], x2: b[0], y2: b[1],
      stroke: pair[1], 'stroke-width': 2.6, 'stroke-linecap': 'round' }));
    var t = svg('text', { x: 0, y: 0, class: 'fig-end-label', fill: pair[1] });
    var o = polar(cx, cy, r + 15, at);
    t.setAttribute('x', o[0]);
    t.setAttribute('y', o[1] + 3);
    t.setAttribute('text-anchor', 'middle');
    t.textContent = pair[0];
    g.appendChild(t);
  });

  figCaption(g, cx, cy, c.name || 'new plasmid', bp(span), r);
  return g;
}

/** The plasmid's name under its circle, its size inside it. */
function figCaption(g, cx, cy, name, size, r) {
  var mid = svg('text', { x: cx, y: cy + 4, 'text-anchor': 'middle',
    class: 'fig-sub' });
  mid.textContent = size;
  g.appendChild(mid);
  var t = svg('text', { x: cx, y: cy + (r || FIG.r) + 26,
    'text-anchor': 'middle', class: 'fig-name' });
  t.textContent = shortName(String(name || 'template'), 28);
  var ti = svg('title', null);
  ti.textContent = String(name || '');
  t.appendChild(ti);
  g.appendChild(t);
}

function shortName(name, n) {
  name = String(name || '');
  return name.length > n ? name.slice(0, n - 1) + '…' : name;
}

/* ------------------------------------------------------------------ */

/** Every enzyme this route mentions anywhere. */
function routeEnzymes(route) {
  var all = [];
  (route.steps || []).forEach(function (st) {
    stepEnzymes(st).forEach(function (e) {
      if (all.indexOf(e) < 0) all.push(e);
    });
  });
  return all;
}

/** The enzyme list one RE_DIGEST step names, wherever it put it. */
function stepEnzymes(st) {
  var ez = ((st.details || {}).enzymes || []).slice();
  if (ez.length) return ez;
  /* Some steps carry them in the inputs instead, as "HindIII, EcoRV". */
  (st.inputs || []).forEach(function (i) {
    if (/^[A-Za-z]+[IVX0-9]+(\s*,\s*[A-Za-z]+[IVX0-9]+)*$/.test(i)) {
      ez = String(i).split(/\s*,\s*/);
    }
  });
  return ez;
}

/** Split "Kozak+ASAP7y-F+Kozak+ASAP7y-R" back into its two primers. */
function stepPrimers(st, route) {
  var names = (route.primers || []).map(function (p) { return p.name; });
  var pair = (st.inputs || []).filter(function (i) {
    return i.indexOf('+') >= 0;
  })[0] || '';
  for (var i = 0; i < names.length; i++) {
    for (var j = 0; j < names.length; j++) {
      if (i !== j && names[i] + '+' + names[j] === pair) {
        return [names[i], names[j]];
      }
    }
  }
  var m = /^(.*-F)\+(.*-R)$/.exec(pair);
  return m ? [m[1], m[2]] : [];
}

/** Where a primer sits on its template, when the plan printed it. */
function primerSpan(route, names) {
  var by = {};
  (route.primers || []).forEach(function (p) { by[p.name] = p; });
  var f = by[names[0]], r = by[names[1]];
  if (!f || typeof f.target_start !== 'number') return null;
  return { start: f.target_start,
           end: typeof f.target_end === 'number' ? f.target_end : r
                && r.target_end };
}

/** Walk the route's steps into one lane per piece that enters the tube. */
function routeLanes(route, sim, plan) {
  var pool = (sim.sources || []).filter(Boolean).slice();
  var atlas = (plan && plan.templates) || {};
  function takeSource(template) {
    for (var i = 0; i < pool.length; i++) {
      var s = pool[i];
      var ok = template === null ? s.synthesised : s.template === template;
      if (ok) { pool.splice(i, 1); return s; }
    }
    return null;
  }
  /* An unsimulated route still knows its plasmids and its primers, which
     is everything the drawing needs. Refusing to draw it was how the
     restriction route came up blank. */
  function stand(template, span) {
    var t = atlas[template] || {};
    return { name: '', template: t.name || template, template_id: t.id,
             template_length: t.bp, template_circular: t.circular !== false,
             features: t.features || [],
             from_start: span && span.start, from_end: span && span.end,
             bp: span && span.end > span.start ? span.end - span.start : 0,
             estimated: true };
  }

  var lanes = [], last = null;
  (route.steps || []).forEach(function (st) {
    var op = st.operation;
    if (op === 'PCR') {
      var tmpl = (st.details || {}).template || (st.inputs || [])[0] || '';
      var names = stepPrimers(st, route);
      var src = takeSource(tmpl) || stand(tmpl, primerSpan(route, names));
      if (!src.bp) src.bp = (st.details || {}).length || src.bp;
      last = { origin: 'plasmid', template: tmpl, src: src,
               piece: pieceLabel(st.output), primers: names, colour: '',
               stages: [{ verb: 'PCR', sub: bp(src.bp) }] };
      lanes.push(last);
    } else if (op === 'RE_DIGEST') {
      var ez = stepEnzymes(st);
      /* An input that is not an enzyme list is the plasmid being cut. */
      var host = (st.inputs || []).filter(function (i) {
        return ez.indexOf(i) < 0 && i !== ez.join(', ');
      })[0];
      if (!ez.length) ez = routeEnzymes(route);
      if (host) {
        var hsrc = takeSource(host) || stand(host, null);
        last = { origin: 'plasmid', template: host, src: hsrc,
                 piece: pieceLabel(st.output), colour: '', cutsParent: ez,
                 stages: [{ verb: 'Digest', sub: ez.join(' + ') }] };
        lanes.push(last);
      } else if (last) {
        /* A bare digest cuts whatever came off the bench last - almost
           always the amplicon, which is how a hybrid route gets sticky
           ends onto a PCR product. That is a second stage, not a
           footnote, so it gets its own column. */
        last.stages.push({ verb: 'Digest', sub: ez.join(' + ') });
        last.cutsFragment = ez;
      }
    } else if (op === 'SYNTHESIS') {
      var osrc = takeSource(null)
              || { bp: (st.details || {}).length || 0, synthesised: true };
      last = { origin: 'order', template: null, src: osrc,
               piece: pieceLabel(st.output), colour: '',
               stages: [{ verb: 'Order', sub: bp(osrc.bp) }] };
      lanes.push(last);
    }
  });

  /* What the ends of each piece are, which is what the join is made of. */
  var ops = (route.steps || []).map(function (s) { return s.operation; });
  var homology = ops.indexOf('HI_FI_ASSEMBLY') >= 0
               || ops.indexOf('GIBSON_ASSEMBLY') >= 0;
  lanes.forEach(function (lane, i) {
    lane.colour = LANE_COLOURS[i % LANE_COLOURS.length];
    var ez = lane.cutsFragment || lane.cutsParent;
    if (ez && ez.length) {
      lane.ends = [ez[0], ez[1] || ez[0]];
      lane.endKind = 'sticky';
    } else if (homology) {
      lane.ends = ['overlap', 'overlap'];
      lane.endKind = 'homology';
    } else {
      lane.ends = [];
    }
    if (!lane.src) lane.src = { bp: 0, template: lane.template };
  });
  return lanes;
}

/** "EGFP amplicon with sites" -> "EGFP"; "linearised vector" stays. */
function pieceLabel(out) {
  return String(out || 'fragment')
    .replace(/\s+(amplicon|fragment)(\s+with sites)?$/i, '')
    .replace(/^(purified|clean|digested)\s+/i, '');
}

/** The reaction that puts the pieces together. */
function joinVerb(route) {
  var ops = (route.steps || []).map(function (s) { return s.operation; });
  if (ops.indexOf('HI_FI_ASSEMBLY') >= 0) return 'NEBuilder HiFi';
  if (ops.indexOf('GIBSON_ASSEMBLY') >= 0) return 'Gibson assembly';
  if (ops.indexOf('LIGATION') >= 0) return 'T4 ligase';
  return 'Assemble';
}

function cloningFigure(route, plan) {
  var sim = route.simulation || {};
  var lanes = routeLanes(route, sim, plan);
  if (!lanes.length) return null;

  var maxStages = lanes.reduce(function (n, l) {
    return Math.max(n, l.stages.length);
  }, 1);
  var xParent = FIG.pad;
  var xStage = xParent + FIG.colParent;
  var xFrag = xStage + maxStages * FIG.colStage;
  var xBracket = xFrag + FIG.colFrag + 10;
  var xJoin = xBracket + FIG.bracket;
  var xProd = xJoin + FIG.colStage + FIG.pr;
  var W = xProd + FIG.pr + FIG.pad + 40;
  var rows = lanes.length;
  var H = Math.max(rows * FIG.rowH, 2 * FIG.pr + 56) + 24;

  var root = svg('svg', {
    viewBox: '0 0 ' + Math.round(W) + ' ' + Math.round(H),
    class: 'fig-svg', role: 'img',
    'aria-label': 'how ' + (route.title || 'this route')
                  + ' builds the construct'
  });
  var defs = svg('defs', null);
  var marker = svg('marker', { id: 'figHead', viewBox: '0 0 10 10',
    refX: '9', refY: '5', markerWidth: '6', markerHeight: '6',
    orient: 'auto-start-reverse' });
  marker.appendChild(svg('path', { d: 'M 0 0 L 10 5 L 0 10 z',
    class: 'fig-arrow-head' }));
  defs.appendChild(marker);
  root.appendChild(defs);

  var joins = junctionColours(lanes);
  lanes.forEach(function (lane, i) {
    var cy = i * FIG.rowH + FIG.rowH / 2;
    var cx = xParent + FIG.colParent / 2;
    root.appendChild(lane.origin === 'order' ? figOrder(lane, cx, cy)
                                             : figParent(lane, cx, cy));
    lane.stages.forEach(function (stage, k) {
      root.appendChild(figStageArrow(xStage + k * FIG.colStage, cy,
                                     FIG.colStage, stage.verb, stage.sub));
    });
    root.appendChild(figFragment(lane, i, joins, xFrag, cy, FIG.colFrag));
  });

  /* One bracket gathering the pieces into the single reaction. */
  var top = FIG.rowH / 2, bot = (rows - 1) * FIG.rowH + FIG.rowH / 2;
  if (rows > 1) {
    root.appendChild(svg('path', {
      d: 'M ' + xBracket + ' ' + top + ' H ' + (xBracket + 10)
         + ' V ' + bot + ' H ' + xBracket,
      fill: 'none', class: 'fig-arrow', 'stroke-width': 1.8
    }));
  }
  var mid = (top + bot) / 2;
  root.appendChild(figStageArrow(xJoin, mid, FIG.colStage, joinVerb(route),
                                 rows + ' pieces'));
  root.appendChild(figProduct(lanes, plan, sim, xProd, mid, joins));
  return root;
}

function figLegend(route, lanes) {
  var ops = (route.steps || []).map(function (s) { return s.operation; });
  var bits = [];
  if (ops.indexOf('PCR') >= 0) bits.push('arrows on a circle are primers');
  if (ops.indexOf('RE_DIGEST') >= 0) bits.push('ticks are enzyme cuts');
  if (lanes.some(function (l) { return l.origin === 'order'; })) {
    bits.push('a dashed box is DNA nobody has, so it is ordered');
  }
  bits.push('each colour follows one piece into the product');
  if (Object.keys(junctionColours(lanes)).length) {
    bits.push(lanes[0] && lanes[0].endKind === 'sticky'
      ? 'matching end colours are compatible overhangs'
      : 'the two ends sharing a colour and letter are one overlap');
  }
  return bits.join(' · ');
}

function renderCloningFlow(route, plan) {
  if (!route) return null;
  var sim = route.simulation || {};
  var card = h('div', { class: 'card' },
    h('div', { class: 'card-head' },
      h('h2', null, 'How this route builds it'),
      h('span', { class: 'muted small' }, route.title || '')));

  var lanes = routeLanes(route, sim, plan);
  if (!lanes.length) {
    card.appendChild(h('p', { class: 'muted' },
      'This route has no amplification, digest or synthesis step to draw. '
      + 'The numbered protocol is on the Protocol step.'));
    return card;
  }

  var figure = cloningFigure(route, plan);
  if (figure) card.appendChild(h('div', { class: 'fig-wrap' }, figure));
  card.appendChild(h('p', { class: 'fig-legend muted small' },
    figLegend(route, lanes)));
  if (!sim.simulated) {
    card.appendChild(h('p', { class: 'muted small' },
      'Drawn from the planned steps; this route was not simulated, so the '
      + 'product shown is the design rather than a checked molecule.'));
  }
  return card;
}


/* ------------------------------------------------------------------ */
/* Routes: filter by technique                                         */
/* ------------------------------------------------------------------ */
/* Which way a construct gets built is a lab preference before it is an
   optimisation - some benches have no HiFi master mix, some would rather
   not design primers - so the ranked list is filterable by what the
   route actually does. */

var METHOD_FILTERS = [
  { key: 'pcr',      label: 'PCR',        hint: 'amplifies at least one piece' },
  { key: 'enzyme',   label: 'Enzyme',     hint: 'cuts with restriction enzymes' },
  { key: 'assembly', label: 'Homology',   hint: 'HiFi or Gibson assembly' },
  { key: 'order',    label: 'Synthesis',  hint: 'orders a fragment' }
];

var routeFilter = {};   // key -> true when that technique is required

/** The set of techniques a route uses. */
function routeMethods(route) {
  var ops = (route.steps || []).map(function (s) { return s.operation; });
  var m = {};
  if (ops.indexOf('PCR') >= 0) m.pcr = true;
  if (ops.indexOf('RE_DIGEST') >= 0) m.enzyme = true;
  if (ops.indexOf('HI_FI_ASSEMBLY') >= 0 ||
      ops.indexOf('GIBSON_ASSEMBLY') >= 0) m.assembly = true;
  if (ops.indexOf('SYNTHESIS') >= 0) m.order = true;
  return m;
}

/** Ticked boxes are the techniques on offer: a route need only use one.

   Requiring all of them answered "which route does everything I ticked",
   which is not a question anybody has. The question is "I am willing to
   do PCR or enzymes - show me those".
 */
function routePassesFilter(route) {
  var want = Object.keys(routeFilter).filter(function (k) {
    return routeFilter[k];
  });
  if (!want.length) return true;
  var has = routeMethods(route);
  return want.some(function (k) { return !!has[k]; });
}

function methodChips(route) {
  var has = routeMethods(route);
  var on = METHOD_FILTERS.filter(function (f) { return has[f.key]; });
  if (!on.length) return null;
  return h('div', { class: 'route-methods' }, on.map(function (f) {
    return h('span', { class: 'method-chip method-' + f.key, title: f.hint },
      f.label);
  }));
}

function methodFilterBar(routes) {
  var bar = h('div', { class: 'method-filter' },
    h('span', { class: 'muted small' }, 'show routes using:'));
  METHOD_FILTERS.forEach(function (f) {
    var n = routes.filter(function (r) { return routeMethods(r)[f.key]; }).length;
    var box = h('input', {
      type: 'checkbox', checked: routeFilter[f.key] ? 'checked' : null,
      disabled: n ? null : 'disabled',
      onchange: function (e) {
        routeFilter[f.key] = e.target.checked;
        renderPlan();
      }
    });
    bar.appendChild(h('label', {
      class: 'method-box' + (n ? '' : ' muted'),
      title: f.hint + ' \u2014 ' + n + ' of ' + routes.length + ' route(s)'
    }, box, h('span', null, f.label + ' (' + n + ')')));
  });
  var any = Object.keys(routeFilter).some(function (k) { return routeFilter[k]; });
  if (any) {
    bar.appendChild(h('button', {
      type: 'button', class: 'link-button',
      onclick: function () { routeFilter = {}; renderPlan(); }
    }, 'clear'));
  }
  return bar;
}

function renderRoutes(routes) {
  var card = h('div', { class: 'card' },
    h('div', { class: 'card-head' },
      h('h2', null, 'Routes'),
      h('span', { class: 'muted small' },
        'ranked best first \u2014 click one to choose it')));

  if (!routes.length) {
    card.appendChild(h('p', { class: 'muted' }, 'No routes were generated.'));
    return card;
  }

  card.appendChild(methodFilterBar(routes));
  var shown = routes.filter(routePassesFilter);
  if (!shown.length) {
    card.appendChild(h('p', { class: 'muted' },
      'No route uses any of the techniques ticked above.'));
    return card;
  }
  if (shown.length < routes.length) {
    card.appendChild(h('p', { class: 'muted small' },
      'Showing ' + shown.length + ' of ' + routes.length
      + '; the rest use none of the techniques ticked.'));
  }

  card.appendChild(h('div', { class: 'routes' }, shown.map(function (r, i) {
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
        r.feasible ? null : h('span', { class: 'pill pill-bad' }, 'not feasible'),
        methodChips(r)),
      /* What actually decides between two routes at the bench: how many
         reactions, how many oligos to order, how long it takes. "cost"
         is a dimensionless ranking score that reads like dollars, and
         an infeasible route carries a +1000 sort sentinel, so it is
         named for what it is and the sentinel is not shown at all.
         "steps" and "assemblies" were the same on every route. */
      h('dl', { class: 'route-stats' },
        statPair('PCRs', r.n_pcr),
        statPair('primers', (r.primers || []).length),
        statPair('hands-on', hours(r.hands_on_minutes)),
        r.feasible ? statPair('effort', Math.round(r.cost)) : null),
      r.feasible
        ? ((r.warnings || []).length
            ? h('p', { class: 'route-note inline-warn',
                       title: r.warnings.join('\n') },
                firstLine(r.warnings[0])
                + (r.warnings.length > 1
                    ? '  (+' + (r.warnings.length - 1) + ' more)' : ''))
            : h('p', { class: 'route-note muted' }, 'no warnings'))
        : h('p', { class: 'route-note reason',
                    title: (r.warnings || []).join('\n') },
            (r.warnings && r.warnings.length) ? r.warnings[0]
                                              : 'no reason given'));
  })));
  return card;
}

/** A warning's first clause, for the one line a card has room for. */
function firstLine(text) {
  var t = String(text || '').replace(/^simulation:\s*/, '');
  if (t.length > 72) t = t.slice(0, 71).replace(/[\s,;:]+\S*$/, '') + '\u2026';
  return t;
}

var STAT_HINTS = {
  effort: 'relative planning score \u2014 reactions, oligos and ordered '
        + 'DNA weighted together; lower is cheaper',
  'hands-on': 'active bench time, not counting overnight growth, '
            + 'screening or sequencing'
};

/** A one-line link that opens a paragraph. Used wherever a long
 *  explanation was being read out before anybody asked for it. */
function disclosure(label, text) {
  var body = h('p', { class: 'muted small', hidden: 'hidden' }, text);
  var link = h('button', {
    type: 'button', class: 'link-button',
    onclick: function () {
      body.hidden = !body.hidden;
      link.textContent = body.hidden ? label : 'hide';
    }
  }, label);
  return h('div', { class: 'disclosure' }, link, body);
}

/** The look-before-you-go-on block under a QC step.

   The step people skip is the one that would have told them, so it
   does not read like the rest of the list: what to run, what a good
   result is, and a row per way it can go wrong with the thing to try
   next. Every number in it came from this plan. */
function checkpointBlock(step) {
  var c = (step.details || {}).checkpoint;
  if (!c) return null;
  var box = h('div', { class: 'check-box' });
  (c.how || []).forEach(function (line) {
    box.appendChild(h('p', { class: 'check-how' }, line));
  });
  box.appendChild(h('div', { class: 'check-head' }, 'Go on if you see'));
  var ul = h('ul', { class: 'check-list' });
  (c.expect || []).forEach(function (line) {
    ul.appendChild(h('li', null, line));
  });
  box.appendChild(ul);
  if ((c.troubles || []).length) {
    box.appendChild(h('div', { class: 'check-head check-bad' }, 'If not'));
    box.appendChild(table([
      { label: 'What you see', cell: function (t) { return t[0]; } },
      { label: 'What it means', cell: function (t) { return t[1]; } },
      { label: 'What to do', cell: function (t) { return t[2]; } }
    ], c.troubles));
  }
  return box;
}

/** A link that reveals a node, for a table that is worth having and
   not worth showing until it is asked for. */
function foldAway(label, node) {
  node.hidden = true;
  var link = h('button', {
    type: 'button', class: 'link-button',
    onclick: function () {
      node.hidden = !node.hidden;
      link.textContent = node.hidden ? label : 'hide';
    }
  }, label);
  return h('div', { class: 'fold-away' }, link, node);
}

/** One chosen parent, with the rest folded into a select.

   A radio list is the right control for three things and the wrong
   one for thirty: the picker scrolls, every option is on screen
   whether or not it is a contender, and the one that matters - the
   chosen one - looks exactly like the others. */
function candidatePicker(label, rows, chosenId, describe, onPick) {
  if (!rows || rows.length < 2) return h('div', null);
  var chosen = rows.filter(function (r) {
    return r.plasmid_id === chosenId;
  })[0] || rows[0];

  var select = h('select', {
    onchange: function (e) {
      var pick = rows[Number(e.target.value)];
      if (pick) onPick(pick);
    }
  }, rows.map(function (r, i) {
    return h('option', {
      value: String(i),
      selected: r === chosen ? 'selected' : null
    }, describe(r));
  }));

  return h('div', { class: 'cand-pick' },
    h('label', { class: 'cand-pick-label' }, label),
    select,
    h('span', { class: 'muted small' },
      rows.length + ' candidate' + (rows.length === 1 ? '' : 's')
      + ' \u00b7 ranked best first'));
}

function statPair(label, value) {
  if (label === null) return [];
  return [h('dt', { title: STAT_HINTS[label] || null }, label),
          h('dd', null, value === null || value === undefined ? '' : String(value))];
}

function renderRouteDetail(route) {
  var card = h('div', { class: 'card' },
    h('div', { class: 'card-head' },
      h('h2', null, route.title || route.id),
      h('span', { class: 'muted small' },
        (route.n_pcr || 0) + ' PCR · '
        + plural((route.primers || []).length, 'primer') + ' · '
        + hours(route.hands_on_minutes) + ' hands-on')));

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
        s.notes ? h('div', { class: 'step-note' }, s.notes) : null,
        /* What this plan found, as opposed to what the operation always
           does: how many sites DpnI has to cut here, for instance. The
           PDF has carried it since it was added; the screen did not. */
        (s.details && s.details.note)
          ? h('div', { class: 'step-note' }, s.details.note) : null,
        checkpointBlock(s));
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
      /* The same primers in the four columns IDT's bulk order page
         takes, so ordering is a paste rather than a transcription. */
      exportButton('IDT order (.csv)', 'order', base + '_IDT_order.csv'),
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
  paintSteps();
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
  var more = $('builder-more'), help = $('builder-help');
  if (more && help) {
    more.addEventListener('click', function () {
      help.hidden = !help.hidden;
      more.textContent = help.hidden ? 'what the markers mean' : 'hide';
    });
  }
  $('scan').addEventListener('click', scanSources);
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
