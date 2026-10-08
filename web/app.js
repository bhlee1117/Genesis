/* Genesis web client.
 *
 * Vanilla ES2020, no build step. Everything the server returns is lab-file
 * derived text, so it is only ever put into the DOM as a text node (see h()).
 */
'use strict';

/* ------------------------------------------------------------------ */
/* Controlled vocabulary (mirrors genesis/vocab.py FEATURE_TYPES)       */
/* ------------------------------------------------------------------ */

var FEATURE_TYPES = [
  'backbone', 'promoter', 'enhancer', 'CDS', 'trafficking_motif', 'linker',
  'p2A', 'IRES', 'DIO', 'lox_site', 'AAV_ITR', 'resistance', 'origin',
  'tag', 'polyA', 'WPRE', 'restriction_site', 'primer_binding_site',
  'terminator', 'intron', 'unknown'
];

/* Mid-saturation hues, legible on both light and dark backgrounds. */
var TYPE_COLORS = {
  backbone: '#8a8f98', promoter: '#2f7f4f', enhancer: '#4f9e6a',
  CDS: '#2f6fb0', trafficking_motif: '#7a5fb0', linker: '#9aa0a8',
  p2A: '#b07a2f', IRES: '#b0652f', DIO: '#a0398c', lox_site: '#c03a5b',
  AAV_ITR: '#b02f2f', resistance: '#6b8f2f', origin: '#2f8f8f',
  tag: '#8f6f2f', polyA: '#4a5fa8', WPRE: '#5f7f9e', restriction_site: '#8f4f6f',
  primer_binding_site: '#7f8f4f', terminator: '#6f5f8f', intron: '#5f8f7f',
  unknown: '#8a8f98'
};

function colorFor(type) {
  return TYPE_COLORS[type] || TYPE_COLORS.unknown;
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
  appliedSuggestion: null   // suggested_name the user (or we) already accepted
};

/** A builder row. `part_id` pins an exact inventory part when one was picked. */
function newModule(type) {
  return { type: type, name: '', part_id: null };
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
  FEATURE_TYPES.forEach(function (t) {
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

function fetchParts(type, query) {
  var key = type + '\u0000' + query.toLowerCase();
  if (!PART_CACHE[key]) {
    PART_CACHE[key] = api('/api/parts' + qs({ type: type, q: query }))
      .catch(function (err) { delete PART_CACHE[key]; throw err; });
  }
  return PART_CACHE[key];
}

/**
 * Typeahead over /api/parts for one module row.
 *
 * `mod` is the row's model object, so the picker keeps working after rows are
 * reordered or removed (no index lookups). Choosing a row writes both the
 * name and `part_id`; typing by hand drops `part_id` again, because the text
 * no longer necessarily names that part.
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

  var items = [];      // rows currently shown
  var active = -1;     // highlighted row, -1 for none
  var timer = null;
  var reqId = 0;

  function renderChip() {
    clear(chip);
    if (!mod.part_id) { chip.hidden = true; return; }
    chip.hidden = false;
    chip.appendChild(h('span', {
      class: 'pin-label',
      title: 'This module is pinned to inventory part ' + mod.part_id +
             '. The server will use exactly this part instead of guessing.'
    }, 'pinned #' + mod.part_id));
    chip.appendChild(h('button', {
      type: 'button', class: 'pin-clear', title: 'Unpin — match this name loosely again',
      onclick: function () { mod.part_id = null; renderChip(); input.focus(); }
    }, '×'));
  }

  function close() {
    list.hidden = true;
    input.setAttribute('aria-expanded', 'false');
    input.removeAttribute('aria-activedescendant');
    active = -1;
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

  function choose(i) {
    var row = items[i];
    if (!row) return;
    mod.name = row.name || '';
    mod.part_id = row.id;
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

  function draw(rows, query) {
    items = rows;
    active = -1;
    if (!rows.length) {
      show(h('div', { class: 'combo-empty muted' },
        'No ' + mod.type + ' part matching “' + query + '”'));
      return;
    }
    clear(list);
    rows.forEach(function (r, i) {
      var n = Number(r.n_sources || 0);
      list.appendChild(h('div', {
        class: 'combo-item', role: 'option', id: listId + '-' + i,
        'aria-selected': 'false',
        // mousedown fires before blur, and preventDefault keeps focus in the
        // input, so the click always lands on a list that is still open.
        onmousedown: function (e) { e.preventDefault(); choose(i); },
        onmousemove: function () { if (active !== i) { active = i; paint(); } }
      },
        h('span', { class: 'combo-name' }, r.name || '(unnamed)'),
        h('span', { class: 'combo-len mono' }, bp(r.seq_length)),
        h('span', {
          class: 'combo-sources' + (n ? '' : ' zero'),
          title: n + ' plasmid' + (n === 1 ? '' : 's') + ' in the inventory carry this part'
        }, n + ' src')));
    });
    list.scrollTop = 0;
    list.hidden = false;
    input.setAttribute('aria-expanded', 'true');
  }

  function search() {
    var query = input.value.trim();
    if (query.length < 1) { items = []; close(); return; }
    var mine = ++reqId;
    var type = mod.type;
    fetchParts(type, query).then(function (rows) {
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

  input.addEventListener('input', function () {
    mod.name = input.value;
    if (mod.part_id) { mod.part_id = null; renderChip(); }
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
    } else if (e.key === 'Enter') {
      if (!list.hidden && active >= 0) { e.preventDefault(); choose(active); }
    } else if (e.key === 'Escape') {
      if (!list.hidden) { e.preventDefault(); e.stopPropagation(); close(); }
    }
  });

  renderChip();
  return {
    node: node,
    chip: chip,
    // The feature_type select owns the `type` filter: when it changes the old
    // suggestions are for the wrong table slice and must be thrown away.
    refresh: function () {
      if (timer) clearTimeout(timer);
      reqId++;
      items = [];
      close();
      if (document.activeElement === input && input.value.trim()) search();
    }
  };
}

function renderModules() {
  var host = $('modules');
  clear(host);
  if (!state.modules.length) {
    host.appendChild(h('p', { class: 'muted' }, 'No modules yet. Add one below.'));
  }
  state.modules.forEach(function (mod, i) {
    var picker = partPicker(mod);
    var typeSel = h('select', {
      'aria-label': 'Module ' + (i + 1) + ' feature type',
      onchange: function (e) { mod.type = e.target.value; picker.refresh(); }
    }, FEATURE_TYPES.map(function (t) {
      return h('option', { value: t, selected: t === mod.type }, t);
    }));
    typeSel.value = mod.type;

    host.appendChild(h('div', { class: 'module-row' },
      h('span', { class: 'module-index' }, String(i + 1)),
      typeSel,
      picker.node,
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
        }, '×'))));
  });
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
function collectSpec() {
  if (usingYaml()) return $('spec-yaml').value;
  var backbone = $('spec-backbone').value.trim();
  var spec = {
    name: $('spec-name').value.trim() || 'untitled',
    owner: $('spec-owner').value.trim(),
    circular: $('spec-circular').checked,
    modules: state.modules.map(function (m) {
      var out = { type: m.type, name: (m.name || '').trim() };
      // Only send a pin when the row actually carries one; an absent part_id
      // leaves the server free to pick a candidate as it always did.
      if (m.part_id) out.part_id = m.part_id;
      return out;
    })
  };
  // A bare integer is treated as a plasmid id, anything else as a name.
  if (/^\d+$/.test(backbone)) spec.backbone_plasmid_id = parseInt(backbone, 10);
  else spec.backbone = backbone;
  return spec;
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
    $('save-design').disabled = false;
    renderPlan();
    setStatus(status, 'plan ready', 'ok');
  } catch (err) {
    clear(out);
    out.appendChild(errorBox(err.message));
    setStatus(status, 'failed', 'error');
  }
}

async function saveDesign() {
  var status = $('design-status');
  if (!state.plan) { setStatus(status, 'generate a plan first', 'error'); return; }
  if (!state.spec) {
    setStatus(status, 'this plan was loaded from storage; regenerate to save a new version', 'error');
    return;
  }
  setStatus(status, 'saving…', 'busy');
  try {
    var res = await postJSON('/api/save', { spec: state.spec, plan: state.plan });
    setStatus(status, 'saved as design ' + res.design_id, 'ok');
    loadStats();   // the design count in the header is now stale
  } catch (err) {
    setStatus(status, 'save failed: ' + err.message, 'error');
  }
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
    featureMap(c.length, c.layout || [], {
      label: c.name, circular: !!c.circular, markers: c.codon_markers
    }));

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

async function renderDesigns() {
  var out = $('designs-table');
  var status = $('designs-status');
  setStatus(status, 'loading…', 'busy');
  clear(out);
  out.appendChild(h('p', { class: 'muted' }, 'Loading saved designs…'));
  try {
    var rows = await api('/api/designs');
    clear(out);
    out.appendChild(table([
      { label: 'Design', numeric: true, key: 'id' },
      { label: 'Construct', key: 'name' },
      { label: 'Version', numeric: true, key: 'version' },
      { label: 'Owner', key: 'owner' },
      { label: 'Created', key: 'created_at' }
    ], rows, { onRow: function (r) { openDesign(r.id); } }));
    setStatus(status, plural(rows.length, 'design'));
  } catch (err) {
    clear(out);
    out.appendChild(errorBox(err.message));
    setStatus(status, 'failed', 'error');
  }
}

async function openDesign(id) {
  var status = $('designs-status');
  setStatus(status, 'loading design ' + id + '…', 'busy');
  try {
    var d = await api('/api/design/' + id);
    state.plan = d.plan;
    state.appliedSuggestion = null;
    // The stored row carries the plan, not the spec, so saving a new version
    // from here would be ambiguous; disable it until a plan is regenerated.
    state.spec = null;
    state.routeId = (d.plan && d.plan.routes && d.plan.routes.length)
      ? d.plan.routes[0].id : null;
    $('save-design').disabled = true;
    selectTab('design');
    renderPlan();
    setStatus($('design-status'),
      'loaded design ' + d.id + ' (version ' + d.version + ')', 'ok');
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
  if (name === 'designs') renderDesigns();
}

function init() {
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
  $('save-design').addEventListener('click', saveDesign);
  $('designs-go').addEventListener('click', renderDesigns);

  document.querySelectorAll('input[name="spec-mode"]').forEach(function (radio) {
    radio.addEventListener('change', function () {
      var yaml = usingYaml();
      $('builder-pane').hidden = yaml;
      $('yaml-pane').hidden = !yaml;
      if (yaml && !$('spec-yaml').value.trim()) $('spec-yaml').value = specToYaml();
    });
  });

  renderModules();
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
  lines.push('circular: ' + ($('spec-circular').checked ? 'true' : 'false'));
  lines.push('modules:');
  state.modules.forEach(function (m) {
    lines.push('  - type: ' + m.type);
    if ((m.name || '').trim()) lines.push('    name: ' + m.name.trim());
    if (m.part_id) lines.push('    part_id: ' + m.part_id);
  });
  return lines.join('\n') + '\n';
}

document.addEventListener('DOMContentLoaded', init);
