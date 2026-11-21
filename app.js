/*
 * MeshLab - blueprint workbench UI
 * Reads the UMD solver from window.MeshLab, renders an SVG result canvas and
 * wires the controls, tables, modal animation and file exports.
 * No external runtime dependencies.
 */
(function () {
  'use strict';

  var M = (typeof window !== 'undefined' && window.MeshLab) ||
    (typeof globalThis !== 'undefined' ? globalThis.MeshLab : null);
  if (!M) return;

  var NS = 'http://www.w3.org/2000/svg';
  var $ = function (id) { return document.getElementById(id); };
  var SVG_W = 1000;

  /* Palette roles are kept in step with the validated --vm-1 .. --vm-5 steps
     in style.css. Dark steps are selected independently, not flipped. */
  var RAMP = {
    dark: ['#184f95', '#256abf', '#5598e7', '#9ec5f4', '#cde2fb'],
    light: ['#5598e7', '#256abf', '#184f95', '#0d366b', '#082443']
  };
  var COLOR = {
    dark: { tension: '#3987e5', compression: '#e66767', accent: '#3987e5',
      load: '#199e70', reaction: '#d95926', neutral: '#8898ac' },
    light: { tension: '#2a78d6', compression: '#e34948', accent: '#2a78d6',
      load: '#1baf7a', reaction: '#eb6834', neutral: '#7b8a9c' }
  };
  var STATUS = { good: 'var(--good)', warning: 'var(--warning)', critical: 'var(--critical)' };

  var state = {
    key: 'pratt', kind: 'truss', model: null, solution: null,
    safety: null, modes: null, plate: null, field: null, load: 1,
    view: 'force', scale: 100, ghost: true, reactions: true, grid: true,
    theme: 'dark', modeIndex: 0, phase: 0, playing: false,
    filter: 'all', sortKey: 'id', sortDir: 1, hovered: -1, pinned: -1,
    resolution: 1
  };
  var patchCheckCache = {};
  var inspectionTargets = [], animationFrame = null, previousTime = null;
  var reducedMotion = false;
  if (typeof window.matchMedia === 'function') {
    var motionQuery = window.matchMedia('(prefers-reduced-motion: reduce)');
    reducedMotion = motionQuery.matches;
    if (typeof motionQuery.addEventListener === 'function') {
      motionQuery.addEventListener('change', function (event) {
        reducedMotion = event.matches;
        if (reducedMotion) stopAnimation();
      });
    } else if (typeof motionQuery.addListener === 'function') {
      motionQuery.addListener(function (mq) { reducedMotion = mq.matches; });
    }
  }

  /* ---------- small DOM helpers, all label text goes through textContent ---------- */

  function svgEl(tag, attrs) {
    var e = document.createElementNS(NS, tag);
    if (attrs) Object.keys(attrs).forEach(function (k) { e.setAttribute(k, String(attrs[k])); });
    return e;
  }
  function hEl(tag, cls, text) {
    var e = document.createElement(tag);
    if (cls) e.setAttribute('class', cls);
    if (text !== undefined && text !== null) e.textContent = String(text);
    return e;
  }
  function clear(node) {
    if (!node) return;
    while (node.firstChild) node.removeChild(node.firstChild);
  }
  function appendRow(table, label, value, cls) {
    var tr = hEl('tr');
    var a = hEl('td', null, label);
    var b = hEl('td', cls || null, value);
    tr.appendChild(a); tr.appendChild(b);
    table.appendChild(tr);
  }
  function clamp(v, lo, hi) { return Math.min(hi, Math.max(lo, v)); }
  function fixed(v, d) {
    if (typeof v !== 'number' || !isFinite(v)) return 'n/a';
    if (v !== 0 && (Math.abs(v) >= 1e7 || Math.abs(v) < 1e-4)) return v.toExponential(2);
    return v.toFixed(d === undefined ? 2 : d);
  }
  function signed(v, d) { return (v > 0 ? '+' : '') + fixed(v, d); }
  function niceStep(x) {
    if (!(x > 0)) return 1;
    var e = Math.pow(10, Math.floor(Math.log(x) / Math.LN10));
    var f = x / e;
    return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * e;
  }
  function pal() { return COLOR[state.theme]; }
  function ramp() { return RAMP[state.theme]; }
  function rampColor(t, steps) {
    var s = steps || ramp();
    return s[clamp(Math.floor(t * s.length), 0, s.length - 1)];
  }

  /* ---------- analysis ---------- */

  function compute(key) {
    state.key = key;
    var preset = M.getPresets().filter(function (p) { return p.key === key; })[0];
    state.kind = preset ? preset.kind : 'truss';
    var model = M.getPreset(key);
    if (model.nodes) model.nodes.forEach(function (n) { if (n.load) n.load = n.load.map(function (v) { return v * state.load; }); });
    if (model.plate) Object.keys(model.plate.loads).forEach(function (k) {
      if (typeof model.plate.loads[k] === 'number') model.plate.loads[k] *= state.load;
    });
    state.model = model;
    state.plate = null; state.solution = null; state.field = null;
    state.safety = null; state.modes = null;
    if (state.kind === 'plate') {
      state.plate = M.solvePlateCST(model.plate.mesh, model.plate.material,
        model.plate.bc, model.plate.loads, { method: 'pcg', tolerance: 1e-12 });
      if (state.view === 'modes' || state.view === 'force') { state.view = 'stress'; syncViewSelect(); }
      clear($('mode-select'));
    } else if (state.kind === 'field') {
      var fc = model.field;
      if (fc.spec && state.resolution !== 1) {
        fc.mesh = M.generateDomainMesh(Object.assign({}, fc.spec, {
          nx: Math.min(60, Math.max(2, Math.round(fc.spec.nx * state.resolution))),
          ny: Math.min(60, Math.max(2, Math.round(fc.spec.ny * state.resolution))) }));
      }
      if (fc.heater !== undefined) fc.heater *= state.load;
      else if (fc.hotEdge !== undefined) {
        fc.hotEdge = (fc.coldEdge || 0) + (fc.hotEdge - (fc.coldEdge || 0)) * state.load;
      }
      state.field = solveField(fc.mesh, fc);
      state.view = 'stress'; syncViewSelect();
    } else {
      var sol = M.solveTruss(model);
      state.solution = sol;
      state.safety = M.analyzeSafety(sol.model, sol);
      state.modes = M.solveModes(model, 4);
    }
    state.modeIndex = 0;
    state.phase = 0;
    state.hovered = -1;
    state.pinned = -1;
    state.sortKey = 'id';
    state.sortDir = 1;
  }

  function solveField(mesh, config) {
    var hotSet = new Set((mesh.bcSets && mesh.bcSets.hot) || []);
    var coldSet = new Set((mesh.bcSets && mesh.bcSets.cold) || []);
    return M.solvePoisson(mesh, { conductivity: config.conductivity,
      source: config.heater === undefined ? undefined : function (x, y) {
        return config.heater * Math.exp(-Math.pow((x - mesh.width * config.heaterX) / config.heaterRadius, 2) -
          Math.pow((y - mesh.height * config.heaterY) / config.heaterRadius, 2));
      },
      boundary: function (node, i) {
        if (hotSet.has(i)) return config.hotEdge;
        if (coldSet.has(i)) return config.coldEdge;
        if (node.x === 0) return config.hotEdge;
        if (Math.abs(node.x - mesh.width) < 1e-12) return config.coldEdge;
      }
    });
  }

  function updateAnalysis(key) {
    var snapshot = Object.assign({}, state);
    snapshot.load = state.solvedLoad === undefined ? 1 : state.solvedLoad; snapshot.playing = false;
    stopAnimation(); hideTooltip();
    try {
      compute(key); renderAll();
      state.solvedLoad = state.load;
      $('solve-status').textContent = 'SOLVED';
      $('solve-status').setAttribute('data-state', 'good');
    } catch (error) {
      Object.assign(state, snapshot);
      $('preset-select').value = state.key;
      $('load-range').value = String(state.load); $('load-output').textContent = state.load.toFixed(1) + 'x';
      $('solve-status').textContent = 'CHECK INPUT';
      $('solve-status').setAttribute('data-state', 'critical');
      $('status-message').textContent = error.message;
    }
  }

  function maxAbsMemberForce() {
    if (!state.solution) return 1;
    var s = state.solution.members;
    var mx = 0;
    for (var i = 0; i < s.length; i++) mx = Math.max(mx, Math.abs(s[i].axialForce));
    return mx || 1;
  }
  function maxAbsStress() {
    if (!state.solution) return 1;
    var s = state.solution.members, mx = 0;
    for (var i = 0; i < s.length; i++) mx = Math.max(mx, Math.abs(s[i].stress));
    return mx || 1;
  }

  /* ---------- projection ---------- */

  function rawPoint(node, dx, dy, dz) {
    var x = node.x + (dx || 0);
    var y = node.y + (dy || 0);
    var z = (node.z || 0) + (dz || 0);
    if (state.model.dim === 3) return [x + 0.62 * z, -(y + 0.30 * z)];
    return [x, -y];
  }

  function buildProjection() {
    var pts = [], i;
    var s = state.view === 'deformed' ? state.scale : 0;
    var d = state.model.dim;
    var displ = null;
    if (state.kind === 'plate') {
      var pn = state.plate.nodes;
      for (i = 0; i < pn.length; i++) {
        pts.push(rawPoint(pn[i], 0, 0, 0));
        pts.push(rawPoint(pn[i], pn[i].ux * s, pn[i].uy * s, 0));
      }
    } else {
      displ = state.solution.displacements;
      var ns = state.solution.nodes;
      for (i = 0; i < ns.length; i++) {
        pts.push(rawPoint(ns[i], 0, 0, 0));
        pts.push(rawPoint(ns[i], displ[d * i] * s, displ[d * i + 1] * s,
          d === 3 ? displ[d * i + 2] * s : 0));
        if (state.view === 'modes') {
          var amp = modalAmplitude();
          var savedPhase = state.phase;
          state.phase = 0;
          pts.push(rawPoint(ns[i], modeDisplacement(amp, d * i), modeDisplacement(amp, d * i + 1), d === 3 ? modeDisplacement(amp, d * i + 2) : 0));
          state.phase = 180;
          pts.push(rawPoint(ns[i], modeDisplacement(amp, d * i), modeDisplacement(amp, d * i + 1), d === 3 ? modeDisplacement(amp, d * i + 2) : 0));
          state.phase = savedPhase;
        }
      }
    }
    var minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    pts.forEach(function (p) {
      minX = Math.min(minX, p[0]); maxX = Math.max(maxX, p[0]);
      minY = Math.min(minY, p[1]); maxY = Math.max(maxY, p[1]);
    });
    var w = Math.max(maxX - minX, 1e-6), h = Math.max(maxY - minY, 1e-6);
    var H = clamp(Math.round(SVG_W * h / w), 340, 640);
    var margin = 56;
    var scale = Math.min((SVG_W - 2 * margin) / w, (H - 2 * margin) / h);
    var ox = (SVG_W - w * scale) / 2 - minX * scale;
    var oy = (H - h * scale) / 2 - minY * scale;
    return {
      W: SVG_W, H: H, scale: scale, ox: ox, oy: oy,
      minX: minX, maxX: maxX, minY: minY, maxY: maxY,
      px: function (p) { return [ox + p[0] * scale, oy + p[1] * scale]; }
    };
  }

  function modeDisplacement(i, dof) {
    var md = state.modes;
    if (!md || !md.modes.length) return 0;
    var m = md.modes[state.modeIndex];
    var raw = m.modeShape[dof] || 0;
    var peak = 0;
    for (var j = 0; j < m.modeShape.length; j++) peak = Math.max(peak, Math.abs(m.modeShape[j]));
    if (peak <= 0) return 0;
    return (raw / peak) * Math.cos(state.phase * Math.PI / 180) * i;
  }
  function modalAmplitude() {
    var nodes = state.solution.nodes, xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    nodes.forEach(function (n) { xmin = Math.min(xmin, n.x); xmax = Math.max(xmax, n.x); ymin = Math.min(ymin, n.y); ymax = Math.max(ymax, n.y); });
    return Math.max(xmax - xmin, ymax - ymin, 1) * 0.055 * state.scale / 100;
  }

  function memberEndOffsets(memberIndex) {
    var sol = state.solution;
    var d = state.model.dim;
    var s = state.scale;
    var m = sol.members[memberIndex];
    if (state.view === 'deformed') {
      var u = sol.displacements;
      return [
        [u[d * m.a] * s, u[d * m.a + 1] * s, d === 3 ? u[d * m.a + 2] * s : 0],
        [u[d * m.b] * s, u[d * m.b + 1] * s, d === 3 ? u[d * m.b + 2] * s : 0]
      ];
    }
    if (state.view === 'modes') {
      var amp = modalAmplitude();
      return [
        [modeDisplacement(amp, d * m.a), modeDisplacement(amp, d * m.a + 1),
          d === 3 ? modeDisplacement(amp, d * m.a + 2) : 0],
        [modeDisplacement(amp, d * m.b), modeDisplacement(amp, d * m.b + 1),
          d === 3 ? modeDisplacement(amp, d * m.b + 2) : 0]
      ];
    }
    return [[0, 0, 0], [0, 0, 0]];
  }

  function memberColor(m) {
    if (state.view === 'force') {
      if (m.axialForce > 1e-9) return pal().tension;
      if (m.axialForce < -1e-9) return pal().compression;
      return pal().neutral;
    }
    if (state.view === 'stress') {
      return rampColor(Math.abs(m.stress) / maxAbsStress());
    }
    return pal().accent;
  }

  /* ---------- tooltip ---------- */

  var tt = null;
  function tooltipNode() { if (!tt) tt = $('tooltip'); return tt; }
  function tooltipPos(evt) {
    var node = tooltipNode();
    if (!node) return;
    var wrap = $('canvas-wrap');
    var rect = wrap && wrap.getBoundingClientRect ? wrap.getBoundingClientRect() : { left: 0, top: 0 };
    var x = (evt && typeof evt.clientX === 'number' ? evt.clientX - rect.left : 24) + 14;
    var y = (evt && typeof evt.clientY === 'number' ? evt.clientY - rect.top : 24) + 14;
    x = clamp(x, 8, Math.max(8, (rect.width || 800) - (node.offsetWidth || 240) - 8));
    y = clamp(y, 8, Math.max(8, (rect.height || 500) - (node.offsetHeight || 180) - 8));
    node.setAttribute('style', 'left:' + x + 'px;top:' + y + 'px;');
  }
  function showTooltip(title, rows, evt) {
    var node = tooltipNode();
    if (!node) return;
    clear(node);
    node.appendChild(hEl('div', 'tt-title', title));
    rows.forEach(function (r) {
      var row = hEl('div', 'tt-row');
      var key = hEl('span', 'tt-key');
      key.setAttribute('style', 'background:' + (r.color || pal().accent) + ';');
      row.appendChild(key);
      row.appendChild(hEl('span', 'tt-value', r.value));
      row.appendChild(hEl('span', 'tt-name', r.name));
      node.appendChild(row);
    });
    node.removeAttribute('hidden');
    tooltipPos(evt);
  }
  function hideTooltip() {
    var node = tooltipNode();
    if (node) node.setAttribute('hidden', '');
  }
  function attachTip(el, title, rowsFn) {
    var index = inspectionTargets.length;
    inspectionTargets.push({ element: el, title: title, rows: rowsFn });
    el.setAttribute('tabindex', index === 0 ? '0' : '-1');
    el.setAttribute('role', 'button');
    el.setAttribute('data-inspect', String(index));
    el.setAttribute('aria-pressed', state.pinned === index ? 'true' : 'false');
    function open(evt) {
      showTooltip(title, rowsFn(), evt);
      if (state.pinned < 0) renderInspection(title, rowsFn());
    }
    ['pointermove', 'pointerenter', 'focus'].forEach(function (ev) {
      el.addEventListener(ev, open);
    });
    ['pointerleave', 'blur'].forEach(function (ev) {
      el.addEventListener(ev, hideTooltip);
    });
    function pin() {
      state.pinned = state.pinned === index ? -1 : index;
      inspectionTargets.forEach(function (target, i) { target.element.setAttribute('aria-pressed', i === state.pinned ? 'true' : 'false'); });
      renderInspection(title, rowsFn()); hideTooltip();
    }
    el.addEventListener('click', pin);
    el.addEventListener('keydown', function (event) {
      if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); pin(); }
      else if (/^Arrow(Left|Right|Up|Down)$/.test(event.key) || event.key === 'Home' || event.key === 'End') {
        event.preventDefault();
        var delta = event.key === 'ArrowLeft' || event.key === 'ArrowUp' ? -1 : 1;
        var next = event.key === 'Home' ? 0 : event.key === 'End' ? inspectionTargets.length - 1 :
          (index + delta + inspectionTargets.length) % inspectionTargets.length;
        inspectionTargets.forEach(function (target, i) { target.element.setAttribute('tabindex', i === next ? '0' : '-1'); });
        inspectionTargets[next].element.focus();
      } else if (event.key === 'Escape') {
        state.pinned = -1; hideTooltip();
        inspectionTargets.forEach(function (target) { target.element.setAttribute('aria-pressed', 'false'); });
        resetInspection();
      }
    });
  }
  function renderInspection(title, rows) {
    var node = $('inspector'); clear(node);
    node.appendChild(hEl('span', 'eyebrow', state.pinned >= 0 ? 'PINNED RESULT / ' + title : 'LIVE RESULT / ' + title));
    var values = hEl('div', 'inspection-values');
    rows.forEach(function (row) { var cell = hEl('div'); cell.appendChild(hEl('span', 'inspection-name', row.name)); cell.appendChild(hEl('strong', null, row.value)); values.appendChild(cell); });
    node.appendChild(values);
  }
  function resetInspection() {
    clear($('inspector')); $('inspector').appendChild(hEl('span', 'eyebrow', 'LIVE INSPECTOR'));
    $('inspector').appendChild(hEl('p', null, 'Hover a result or select it to pin its measurements.'));
  }

  /* ---------- SVG primitives ---------- */

  function drawGrid(g, pr) {
    if (!state.grid) return;
    var step = niceStep((pr.maxX - pr.minX) / 8);
    var x0 = Math.ceil(pr.minX / step) * step;
    var y0 = Math.ceil(pr.minY / step) * step;
    var grid = svgEl('g', { 'class': 'svg-grid' });
    var p, t;
    for (var x = x0; x <= pr.maxX + 1e-9; x += step) {
      p = pr.px([x, pr.minY]); var q = pr.px([x, pr.maxY]);
      grid.appendChild(svgEl('line', { x1: fixed(p[0], 2), y1: fixed(p[1], 2),
        x2: fixed(q[0], 2), y2: fixed(q[1], 2) }));
      t = svgEl('text', { 'class': 'svg-tick-text', x: fixed(p[0] + 4, 2),
        y: fixed(pr.H - 12, 2), 'text-anchor': 'start' });
      t.textContent = fixed(x, step < 1 ? 2 : 0) + ' m';
      grid.appendChild(t);
    }
    for (var y = y0; y <= pr.maxY + 1e-9; y += step) {
      p = pr.px([pr.minX, y]); var q2 = pr.px([pr.maxX, y]);
      grid.appendChild(svgEl('line', { x1: fixed(p[0], 2), y1: fixed(p[1], 2),
        x2: fixed(q2[0], 2), y2: fixed(q2[1], 2) }));
      t = svgEl('text', { 'class': 'svg-tick-text', x: 10,
        y: fixed(p[1] - 4, 2) });
      t.textContent = fixed(-y, step < 1 ? 2 : 0) + ' m';
      grid.appendChild(t);
    }
    g.appendChild(grid);
  }

  function drawArrow(g, at, vec, color, label, into, textClass, lenPx) {
    var mag = Math.sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
    if (!(mag > 1e-12)) return;
    var len = lenPx === undefined ? 34 : lenPx;
    var ux = vec[0] / mag, uy = -vec[1] / mag;
    var tail = into ? [at[0] - ux * len, at[1] - uy * len] : [at[0], at[1]];
    var tip = into ? [at[0], at[1]] : [at[0] + ux * len, at[1] + uy * len];
    var line = svgEl('line', { x1: fixed(tail[0], 2), y1: fixed(tail[1], 2),
      x2: fixed(tip[0], 2), y2: fixed(tip[1], 2), stroke: color,
      'stroke-width': 2, 'stroke-linecap': 'round' });
    line.setAttribute('class', textClass === 'svg-field-flux' ? 'svg-field-flux' : textClass === 'svg-reaction-text' ? 'svg-reaction' : 'svg-load');
    g.appendChild(line);
    var bx = tip[0] - ux * 9, by = tip[1] - uy * 9;
    var px = -uy * 4.5, py = ux * 4.5;
    var pts = fixed(tip[0], 2) + ',' + fixed(tip[1], 2) + ' ' +
      fixed(bx + px, 2) + ',' + fixed(by + py, 2) + ' ' +
      fixed(bx - px, 2) + ',' + fixed(by - py, 2);
    var head = svgEl('polygon', { points: pts, fill: color });
    head.setAttribute('class', textClass === 'svg-field-flux' ? 'svg-field-flux-head' : textClass === 'svg-reaction-text' ? 'svg-reaction-head' : 'svg-load-head');
    g.appendChild(head);
    if (label) {
      var tx = svgEl('text', { 'class': textClass || 'svg-label',
        x: fixed(tail[0] + 4, 2), y: fixed(tail[1] - 5, 2) });
      tx.textContent = label;
      g.appendChild(tx);
    }
  }

  function drawSupport(g, p, support) {
    if (!support || support === 'free') return;
    var s = String(support).toLowerCase();
    var x = p[0], y = p[1];
    if (s === 'fixed') {
      g.appendChild(svgEl('rect', { x: fixed(x - 7, 2), y: fixed(y - 7, 2),
        width: 14, height: 14, 'class': 'svg-node-support' }));
      return;
    }
    var tri = fixed(x, 2) + ',' + fixed(y, 2) + ' ' +
      fixed(x - 8, 2) + ',' + fixed(y + 12, 2) + ' ' +
      fixed(x + 8, 2) + ',' + fixed(y + 12, 2);
    g.appendChild(svgEl('polygon', { points: tri, fill: 'none',
      'class': 'svg-support-line' }));
    g.appendChild(svgEl('line', { x1: fixed(x - 11, 2), y1: fixed(y + 13, 2),
      x2: fixed(x + 11, 2), y2: fixed(y + 13, 2), 'class': 'svg-support-line' }));
    if (s.indexOf('roller') === 0) {
      g.appendChild(svgEl('circle', { cx: fixed(x - 5, 2), cy: fixed(y + 16, 2),
        r: 2, 'class': 'svg-node-support' }));
      g.appendChild(svgEl('circle', { cx: fixed(x + 5, 2), cy: fixed(y + 16, 2),
        r: 2, 'class': 'svg-node-support' }));
    }
  }

  /* ---------- truss canvas ---------- */

  function renderTrussCanvas(svg, pr) {
    var sol = state.solution, model = state.model;
    var d = model.dim;
    var ghostG = svgEl('g');
    var memberG = svgEl('g');
    var overlayG = svgEl('g');
    var hitG = svgEl('g');

    var maxF = maxAbsMemberForce();
    var safetyById = {};
    state.safety.members.forEach(function (r) { safetyById[r.id] = r; });

    sol.members.forEach(function (m, i) {
      var na = model.nodes[m.a], nb = model.nodes[m.b];
      var off = memberEndOffsets(i);
      var p0 = pr.px(rawPoint(na, 0, 0, 0));
      var p1 = pr.px(rawPoint(nb, 0, 0, 0));
      var q0 = pr.px(rawPoint(na, off[0][0], off[0][1], off[0][2]));
      var q1 = pr.px(rawPoint(nb, off[1][0], off[1][1], off[1][2]));
      if (state.ghost) {
        ghostG.appendChild(svgEl('line', { x1: fixed(p0[0], 2), y1: fixed(p0[1], 2),
          x2: fixed(p1[0], 2), y2: fixed(p1[1], 2), 'class': 'svg-ghost' }));
      }
      var width = state.view === 'force'
        ? 2.2 + 3.4 * Math.min(1, Math.abs(m.axialForce) / maxF)
        : 3.4;
      var line = svgEl('line', { x1: fixed(q0[0], 2), y1: fixed(q0[1], 2),
        x2: fixed(q1[0], 2), y2: fixed(q1[1], 2),
        stroke: memberColor(m), 'stroke-width': fixed(width, 2),
        'stroke-linecap': 'round' });
      line.setAttribute('class', 'svg-member' + (state.hovered === i ? ' hovered' : '') + (safetyById[m.id] && safetyById[m.id].critical ? ' is-critical' : ''));
      memberG.appendChild(line);

      var hit = svgEl('line', { x1: fixed(q0[0], 2), y1: fixed(q0[1], 2),
        x2: fixed(q1[0], 2), y2: fixed(q1[1], 2),
        'class': 'svg-member-hit', tabindex: '0', role: 'button',
        'aria-label': m.label + ' axial force ' + fixed(m.axialForce / 1000, 2) + ' kN' });
      (function (idx, mm, sr) {
        attachTip(hit, mm.label, function () {
          return [
            { value: fixed(mm.axialForce / 1000, 2) + ' kN', name: 'Axial N',
              color: mm.axialForce >= 0 ? pal().tension : pal().compression },
            { value: fixed(mm.stress / 1e6, 2) + ' MPa', name: 'Stress' },
            { value: fixed(mm.length, 3) + ' m', name: 'Length' },
            { value: sr && isFinite(sr.fosGoverning) ? fixed(sr.fosGoverning, 2) : 'n/a',
              name: 'FoS governing' },
            { value: Math.abs(mm.axialForce) <= 1e-9 ? 'ZERO FORCE' : mm.axialForce > 0 ? 'TENSION' : 'COMPRESSION', name: 'State' },
            { value: sr ? sr.governingMode.toUpperCase() : 'n/a', name: 'Governing check' }
          ];
        });
        hit.addEventListener('focus', function () { state.hovered = idx; });
      })(i, m, safetyById[m.id]);
      hitG.appendChild(hit);
    });

    var nodeG = svgEl('g');
    var supportG = svgEl('g');
    var arrowG = svgEl('g');
    var maxR = 1, maxL = 1;
    sol.nodes.forEach(function (n) {
      maxR = Math.max(maxR, Math.abs(n.reaction[0]), Math.abs(n.reaction[1]));
      maxL = Math.max(maxL, Math.abs(n.load[0]), Math.abs(n.load[1]));
    });
    sol.nodes.forEach(function (n, i) {
      var off = [0, 0, 0];
      if (state.view === 'deformed') {
        off = [sol.displacements[d * i] * state.scale,
          sol.displacements[d * i + 1] * state.scale,
          d === 3 ? sol.displacements[d * i + 2] * state.scale : 0];
      } else if (state.view === 'modes') {
        var amp = modalAmplitude();
        off = [modeDisplacement(amp, d * i), modeDisplacement(amp, d * i + 1),
          d === 3 ? modeDisplacement(amp, d * i + 2) : 0];
      }
      var p = pr.px(rawPoint(n, off[0], off[1], off[2]));
      var base = pr.px(rawPoint(n, 0, 0, 0));
      if (n.support && n.support !== 'free') {
        drawSupport(supportG, base, n.support);
        nodeG.appendChild(svgEl('circle', { cx: fixed(base[0], 2), cy: fixed(base[1], 2),
          r: 4, 'class': 'svg-node-support' }));
      } else {
        nodeG.appendChild(svgEl('circle', { cx: fixed(p[0], 2), cy: fixed(p[1], 2),
          r: 3.2, 'class': 'svg-node' }));
      }
      var lx = n.load[0] + (d === 3 ? 0.62 * n.load[2] : 0), ly = n.load[1] + (d === 3 ? 0.30 * n.load[2] : 0);
      if (Math.abs(lx) + Math.abs(ly) > 1e-9) {
        drawArrow(arrowG, base, [lx / maxL, ly / maxL], pal().load,
          fixed(Math.sqrt(lx * lx + ly * ly) / 1000, 1) + ' kN', true, 'svg-label');
      }
      var rx = n.reaction[0] + (d === 3 ? 0.62 * n.reaction[2] : 0), ry = n.reaction[1] + (d === 3 ? 0.30 * n.reaction[2] : 0);
      if (state.reactions && (Math.abs(rx) + Math.abs(ry)) > 1e-6) {
        drawArrow(arrowG, base, [rx / maxR, ry / maxR],
          pal().reaction, null, false, 'svg-reaction-text');
      }
    });

    svg.appendChild(ghostG);
    svg.appendChild(memberG);
    svg.appendChild(overlayG);
    svg.appendChild(supportG);
    svg.appendChild(nodeG);
    svg.appendChild(arrowG);
    svg.appendChild(hitG);
  }

  /* ---------- plate canvas ---------- */

  function plateRowsData() {
    return state.plate.elements;
  }

  function renderPlateCanvas(svg, pr) {
    var p = state.plate;
    var mesh = p.mesh;
    var lo = p.minVonMises, hi = p.maxVonMises;
    var denom = Math.max(hi - lo, 1e-9);
    var triG = svgEl('g');
    var edgeG = svgEl('g');
    var hitG = svgEl('g');
    var s = state.scale;
    var deformed = state.view === 'deformed';
    var boundaryEdges = new Map();
    p.elements.forEach(function (e) {
      var ns = e.nodes.map(function (k) { return p.nodes[k]; });
      var pts = ns.map(function (n) {
        var q = deformed ? rawPoint(n, n.ux * s, n.uy * s, 0) : rawPoint(n, 0, 0, 0);
        return pr.px(q);
      });
      var pointStr = pts.map(function (q) { return fixed(q[0], 2) + ',' + fixed(q[1], 2); }).join(' ');
      var fill = rampColor((e.vonMises - lo) / denom);
      var poly = svgEl('polygon', { points: pointStr, fill: fill,
        'class': 'svg-plate-edge' });
      triG.appendChild(poly);
      if (deformed) {
        var outline = svgEl('polygon', { points: pointStr, fill: 'none',
          stroke: pal().accent, 'stroke-width': 0.8, opacity: 0.75 });
        edgeG.appendChild(outline);
      }
      var hit = svgEl('polygon', { points: pointStr,
        tabindex: '0', role: 'button',
        'aria-label': 'Element ' + (e.id + 1) + ' von Mises ' + fixed(e.vonMises / 1e6, 2) + ' MPa' });
      hit.setAttribute('class', 'svg-plate-hit');
      (function (ee) {
        attachTip(hit, 'CST element ' + (ee.id + 1), function () {
          return [
            { value: fixed(ee.vonMises / 1e6, 2) + ' MPa', name: 'von Mises',
              color: rampColor((ee.vonMises - lo) / denom) },
            { value: fixed(ee.stress.sigmaX / 1e6, 2) + ' MPa', name: 'sigma x' },
            { value: fixed(ee.stress.tauXY / 1e6, 2) + ' MPa', name: 'tau xy' },
            { value: fixed(ee.sigma1 / 1e6, 2) + ' MPa', name: 'sigma 1' },
            { value: fixed(ee.area, 5) + ' m2', name: 'Area' }
          ];
        });
      })(e);
      hitG.appendChild(hit);
      if (state.ghost && deformed) for (var edge = 0; edge < 3; edge++) {
        var a = e.nodes[edge], b = e.nodes[(edge + 1) % 3], key = Math.min(a, b) + ':' + Math.max(a, b);
        if (boundaryEdges.has(key)) boundaryEdges.delete(key); else boundaryEdges.set(key, [a, b]);
      }
    });
    svg.appendChild(triG);
    svg.appendChild(edgeG);
    boundaryEdges.forEach(function (edge) {
      var a = pr.px(rawPoint(p.nodes[edge[0]], 0, 0, 0)), b = pr.px(rawPoint(p.nodes[edge[1]], 0, 0, 0));
      svg.appendChild(svgEl('line', { x1: fixed(a[0], 2), y1: fixed(a[1], 2), x2: fixed(b[0], 2), y2: fixed(b[1], 2), 'class': 'svg-ghost' }));
    });
    if (mesh.holeRadius > 0) {
      var holeCenter = pr.px(rawPoint({ x: mesh.width / 2, y: mesh.height / 2 }, 0, 0, 0));
      var holeEdge = pr.px(rawPoint({ x: mesh.width / 2 + mesh.holeRadius,
        y: mesh.height / 2 }, 0, 0, 0));
      svg.appendChild(svgEl('circle', { cx: fixed(holeCenter[0], 2),
        cy: fixed(holeCenter[1], 2), r: fixed(Math.abs(holeEdge[0] - holeCenter[0]), 2),
        'class': 'svg-plate-hole' }));
      var holeTag = svgEl('text', { 'class': 'svg-dim-text',
        x: fixed(holeCenter[0], 2), y: fixed(holeCenter[1] + 4, 2),
        'text-anchor': 'middle' });
      holeTag.textContent = 'HOLE';
      svg.appendChild(holeTag);
    }
    if (state.reactions && mesh.leftNodes && mesh.leftNodes.length) {
      var arrowG = svgEl('g');
      var left = mesh.leftNodes;
      var mid = left[Math.floor(left.length / 2)];
      if (mid !== undefined && p.nodes[mid]) {
        var n = p.nodes[mid];
        var base = pr.px(rawPoint(n, 0, 0, 0));
        drawArrow(arrowG, base, [-1, 0], pal().reaction, null, false, 'svg-reaction-text');
        var right = mesh.rightNodes[Math.floor(mesh.rightNodes.length / 2)];
        if (right !== undefined && p.nodes[right]) {
          var rb = pr.px(rawPoint(p.nodes[right], 0, 0, 0));
          drawArrow(arrowG, rb, [1, 0], pal().load, null, true, 'svg-label');
        }
      }
      svg.appendChild(arrowG);
    }
    svg.appendChild(hitG);
    var cap = svgEl('text', { 'class': 'svg-title', x: 14, y: 20 });
    cap.textContent = 'CST PLANE STRESS / ' + mesh.elements.length + ' TRIANGLES / ' +
      mesh.nodes.length + ' NODES / HOLE r = ' + fixed(mesh.holeRadius, 2) + ' m';
    svg.appendChild(cap);
  }

  /* ---------- canvas root ---------- */

  function renderCanvas() {
    var svg = $('structure-svg');
    if (!svg) return;
    var active = document.activeElement && document.activeElement.getAttribute && document.activeElement.getAttribute('data-inspect');
    if (state.kind === 'field') { renderFieldCanvas(); restoreFocus(active); return; }
    clear(svg);
    inspectionTargets = [];
    var pr = buildProjection();
    svg.setAttribute('viewBox', '0 0 ' + pr.W + ' ' + pr.H);
    svg.setAttribute('aria-label', (state.kind === 'plate' ? 'Plate' : 'Truss') +
      ' result view: ' + state.view);
    drawGrid(svg, pr);
    if (state.kind === 'plate') renderPlateCanvas(svg, pr);
    else renderTrussCanvas(svg, pr);
    restoreFocus(active);
  }
  function restoreFocus(index) {
    if (index !== null && index !== undefined && inspectionTargets[Number(index)]) {
      inspectionTargets.forEach(function (target, i) { target.element.setAttribute('tabindex', String(i) === index ? '0' : '-1'); });
      inspectionTargets[Number(index)].element.focus();
    }
  }

  /* ---------- header and stat tiles ---------- */

  function renderHeader() {
    var p = M.getPresets().filter(function (x) { return x.key === state.key; })[0];
    $('model-title').textContent = p.name;
  }

  function statTile(label, value, unit, note, hero) {
    var tile = hEl('div', 'stat-tile' + (hero ? ' hero' : ''));
    tile.appendChild(hEl('span', 'stat-label', label));
    var v = hEl('span', 'stat-value', value);
    if (unit) {
      var u = hEl('span', 'stat-unit', unit);
      v.appendChild(u);
    }
    tile.appendChild(v);
    if (note) tile.appendChild(hEl('span', 'stat-note', note));
    return tile;
  }

  function renderStats() {
    var g = $('stat-grid');
    clear(g);
    if (state.kind === 'plate') {
      var p = state.plate;
      var nominal = Math.abs(state.model.plate.loads.rightNodesForce) /
        (p.mesh.height * p.material.t);
      var ratio = nominal > 0 ? p.maxVonMises / nominal : null;
      g.appendChild(statTile('PEAK VON MISES', fixed(p.maxVonMises / 1e6, 2), 'MPa',
        'element maximum, sequential ramp', true));
      g.appendChild(statTile('CONCENTRATION RATIO', ratio === null ? 'unloaded' : fixed(ratio, 3), ratio === null ? '' : 'x',
        'finite width, stair step hole boundary'));
      g.appendChild(statTile('MAX DISPLACEMENT', fixed(p.maxDisplacement * 1000, 4), 'mm',
        'nodal translation magnitude'));
      g.appendChild(statTile('MAX SIGMA X', fixed(p.maxSigmaX / 1e6, 2), 'MPa',
        'axial component'));
      g.appendChild(statTile('CST ELEMENTS', String(p.elements.length), '',
        p.nodes.length + ' mesh nodes'));
      g.appendChild(statTile('EQUILIBRIUM', fixed(Math.max(
        Math.abs(p.equilibrium.residualX), Math.abs(p.equilibrium.residualY)) / 1000, 6), 'kN',
        'largest force residual'));
      return;
    }
    var sol = state.solution, sf = state.safety;
    var peakStress = 0;
    sol.members.forEach(function (m) { peakStress = Math.max(peakStress, Math.abs(m.stress)); });
    g.appendChild(statTile('MAX DISPLACEMENT', fixed(sol.maxDisplacement * 1000, 4), 'mm',
      state.view === 'deformed' ? 'exaggerated ' + state.scale + 'x in the canvas' : 'maximum nodal translation', true));
    g.appendChild(statTile('MIN FACTOR OF SAFETY',
      isFinite(sf.minFos) ? fixed(sf.minFos, 3) : 'unloaded', '',
      sf.criticalCount + ' critical members'));
    g.appendChild(statTile('STRAIN ENERGY', fixed(sol.strainEnergy, 4), 'J',
      'sum of 0.5 N delta'));
    g.appendChild(statTile('PEAK AXIAL STRESS', fixed(peakStress / 1e6, 2), 'MPa',
      'absolute member stress'));
    g.appendChild(statTile('NODES / MEMBERS', sol.nodes.length + ' / ' + sol.members.length, '',
      sol.dim + 'D, ' + sol.determinacy.classification));
    g.appendChild(statTile('FIRST MODE',
      state.modes.modes.length ? fixed(state.modes.modes[0].frequencyHz, 2) : 'n/a', 'Hz',
      state.modes.numModes + ' modes solved'));
  }

  /* ---------- side cards ---------- */

  function renderEquilibrium() {
    var t = $('equilibrium-table');
    var tb = t.getElementsByTagName ? t.getElementsByTagName('tbody')[0] : t.firstChild;
    clear(tb);
    var e = state.kind === 'plate' ? state.plate.equilibrium : state.solution.equilibrium;
    appendRow(tb, 'Sum applied Fx', fixed(e.sumFx / 1000, 4) + ' kN');
    appendRow(tb, 'Sum reaction Rx', fixed(e.sumRx / 1000, 4) + ' kN');
    appendRow(tb, 'Residual X', signed(e.residualX / 1000, 6) + ' kN');
    appendRow(tb, 'Sum applied Fy', fixed(e.sumFy / 1000, 4) + ' kN');
    appendRow(tb, 'Sum reaction Ry', fixed(e.sumRy / 1000, 4) + ' kN');
    appendRow(tb, 'Residual Y', signed(e.residualY / 1000, 6) + ' kN');
    if (state.kind !== 'plate') {
      appendRow(tb, 'Residual Z', signed((e.residualZ || 0) / 1000, 6) + ' kN');
      appendRow(tb, 'External work', fixed(state.solution.externalWork, 6) + ' J');
    }
    var res = Math.max(Math.abs(e.residualX || 0), Math.abs(e.residualY || 0), Math.abs(e.residualZ || 0));
    var badge = $('equilibrium-badge');
    badge.textContent = res < 1e-3 ? 'BALANCED' : 'CHECK';
    badge.setAttribute('data-state', res < 1e-3 ? 'good' : 'critical');
  }

  function renderDeterminacy() {
    var t = $('determinacy-table');
    var tb = t.getElementsByTagName ? t.getElementsByTagName('tbody')[0] : t.firstChild;
    clear(tb);
    if (state.kind === 'plate') {
      var p = state.plate;
      appendRow(tb, 'Mesh nodes', String(p.nodes.length));
      appendRow(tb, 'CST triangles', String(p.elements.length));
      appendRow(tb, 'Degrees of freedom', String(p.nodes.length * 2));
      appendRow(tb, 'Formulation', 'Plane stress');
      appendRow(tb, 'Thickness', fixed(p.material.t, 4) + ' m');
      appendRow(tb, 'Poisson ratio', fixed(p.material.nu, 3));
      return;
    }
    var d = state.solution.determinacy;
    var factor = state.solution.dim === 3 ? 3 : 2;
    appendRow(tb, 'Members m', String(d.m));
    appendRow(tb, 'Reactions r', String(d.r));
    appendRow(tb, 'Joints j', String(d.j));
    appendRow(tb, 'm + r - ' + factor + 'j', String(d.value));
    appendRow(tb, 'Classification', d.classification);
    appendRow(tb, 'Kinematically stable', d.stable ? 'YES' : 'NO');
  }

  function renderSafety() {
    var box = $('safety-meter');
    clear(box);
    var badge = $('safety-badge');
    if (state.kind === 'plate') {
      var p = state.plate;
      var ratio = p.maxVonMises / (Math.abs(state.model.plate.loads.rightNodesForce) /
        (p.mesh.height * p.material.t));
      var fos = state.model.material.sigmaYield / p.maxVonMises;
      var top = hEl('div', 'meter-top');
      top.appendChild(hEl('span', null, 'Yield FoS at peak stress'));
      top.appendChild(hEl('strong', null, isFinite(fos) ? fixed(fos, 3) : 'unloaded'));
      box.appendChild(top);
      var track = hEl('div', 'meter-track');
      var fill = hEl('div', 'meter-fill');
      fill.setAttribute('style', 'width:' + fixed(clamp(fos / 3 * 100, 2, 100), 2) +
        '%;background:' + (fos >= 1.5 ? STATUS.good : fos >= 1 ? STATUS.warning : STATUS.critical) + ';');
      track.appendChild(fill);
      var limit = hEl('div', 'meter-limit');
      limit.setAttribute('style', 'left:33.33%;');
      track.appendChild(limit);
      box.appendChild(track);
      var foot = hEl('div', 'meter-foot');
      foot.appendChild(hEl('span', null, 'FoS 0'));
      foot.appendChild(hEl('span', null, '1.0 LIMIT'));
      foot.appendChild(hEl('span', null, 'FoS 3'));
      box.appendChild(foot);
      var st = hEl('div', 'meter-status');
      var ic = hEl('span', 'icon', fos >= 1.5 ? 'OK' : fos >= 1 ? 'WATCH' : 'CRITICAL');
      ic.setAttribute('style', 'color:' + (fos >= 1.5 ? STATUS.good :
        fos >= 1 ? STATUS.warning : STATUS.critical) + ';font-family:var(--font-mono);font-size:10px;');
      st.appendChild(ic);
      st.appendChild(hEl('span', null, p.maxVonMises > 0 ? 'Peak / nominal ratio ' + fixed(ratio, 3) +
        '. Approximate finite width hole boundary.' : 'No applied load. Stress and displacement are zero.'));
      box.appendChild(st);
      badge.textContent = fos >= 1 ? 'OK' : 'CRITICAL';
      badge.setAttribute('data-state', fos >= 1 ? 'good' : 'critical');
      return;
    }
    var sf = state.safety;
    var fos = sf.minFos;
    var top = hEl('div', 'meter-top');
    top.appendChild(hEl('span', null, 'Governing minimum'));
    top.appendChild(hEl('strong', null, isFinite(fos) ? fixed(fos, 3) : 'unloaded'));
    box.appendChild(top);
    var track = hEl('div', 'meter-track');
    var fill = hEl('div', 'meter-fill');
    var col = fos >= 1.5 ? STATUS.good : fos >= 1 ? STATUS.warning : STATUS.critical;
    fill.setAttribute('style', 'width:' + fixed(clamp((isFinite(fos) ? fos : 3) / 3 * 100, 2, 100), 2) +
      '%;background:' + col + ';');
    track.appendChild(fill);
    var lim = hEl('div', 'meter-limit');
    lim.setAttribute('style', 'left:33.33%;');
    track.appendChild(lim);
    box.appendChild(track);
    var foot = hEl('div', 'meter-foot');
    foot.appendChild(hEl('span', null, 'FoS 0'));
    foot.appendChild(hEl('span', null, '1.0 LIMIT'));
    foot.appendChild(hEl('span', null, 'FoS 3'));
    box.appendChild(foot);
    var st = hEl('div', 'meter-status');
    var label = sf.safe ? (fos >= 1.5 ? 'OK' : 'WATCH') : 'CRITICAL';
    var ic = hEl('span', 'icon', label);
    ic.setAttribute('style', 'color:' + col + ';font-family:var(--font-mono);font-size:10px;');
    st.appendChild(ic);
    st.appendChild(hEl('span', null, label + '. ' + sf.criticalCount +
      ' members below FoS 1. Yield and Euler checks are combined.'));
    box.appendChild(st);
    badge.textContent = label;
    badge.setAttribute('data-state', sf.safe ? 'good' : 'critical');
  }

  function renderModeBars() {
    var box = $('mode-bars');
    clear(box);
    if (state.kind === 'plate' || !state.modes || !state.modes.modes.length) {
      box.appendChild(hEl('p', 'panel-sub',
        'Eigenmodes are computed for the truss presets. The plate view uses the CST contour.'));
      return;
    }
    var modes = state.modes.modes;
    var peak = Math.max.apply(null, modes.map(function (m) { return m.frequencyHz; })) || 1;
    modes.forEach(function (m, i) {
      var row = hEl('div', 'mode-bar-row');
      row.appendChild(hEl('span', 'mode-bar-idx', 'M' + m.index));
      var track = hEl('div', 'mode-bar-track');
      var fill = hEl('div', 'mode-bar-fill');
      fill.setAttribute('style', 'width:' + fixed(m.frequencyHz / peak * 100, 2) + '%;' +
        (i === state.modeIndex && state.view === 'modes' ? 'background:' + pal().reaction + ';' : ''));
      track.appendChild(fill);
      row.appendChild(track);
      row.appendChild(hEl('span', 'mode-bar-val', fixed(m.frequencyHz, 2) + ' Hz'));
      box.appendChild(row);
    });
  }

  /* ---------- legend ---------- */

  function legendSwatch(color, label, isLine) {
    var item = hEl('div', 'legend-item');
    var sw = hEl('span', isLine ? 'legend-line' : 'legend-swatch');
    sw.setAttribute('style', 'background:' + color + ';');
    item.appendChild(sw);
    item.appendChild(hEl('span', 'legend-label', label));
    return item;
  }

  function renderLegend() {
    var row = $('legend-row');
    clear(row);
    if (state.kind === 'plate' || state.view === 'stress') {
      var scale = hEl('div', 'legend-scale');
      scale.appendChild(hEl('span', 'legend-label', 'von Mises'));
      var lo, hi;
      if (state.kind === 'plate') {
        lo = state.plate.minVonMises / 1e6; hi = state.plate.maxVonMises / 1e6;
      } else { lo = 0; hi = maxAbsStress() / 1e6; }
      scale.appendChild(hEl('span', 'legend-tick', fixed(lo, 1) + ' MPa'));
      var r = hEl('div', 'legend-ramp');
      ramp().forEach(function (c) {
        var i = hEl('i');
        i.setAttribute('style', 'background:' + c + ';');
        r.appendChild(i);
      });
      scale.appendChild(r);
      scale.appendChild(hEl('span', 'legend-tick', fixed(hi, 1) + ' MPa'));
      row.appendChild(scale);
      if (state.kind === 'plate') {
        row.appendChild(legendSwatch(pal().load, 'APPLIED EDGE TRACTION', true));
        row.appendChild(legendSwatch(pal().reaction, 'SUPPORT REACTION', true));
      }
      if (state.ghost && (state.kind !== 'plate' || state.view === 'deformed')) row.appendChild(legendSwatch('var(--ghost)', 'UNDEFORMED GHOST', true));
      return;
    }
    if (state.view === 'force') {
      row.appendChild(legendSwatch(pal().tension, 'TENSION N > 0'));
      row.appendChild(legendSwatch(pal().compression, 'COMPRESSION N < 0'));
      row.appendChild(hEl('span', 'legend-label', 'LINE WEIGHT = |N|'));
    } else if (state.view === 'deformed') {
      row.appendChild(legendSwatch(pal().accent, 'DEFORMED SHAPE'));
      row.appendChild(legendSwatch('var(--ghost)', 'UNDEFORMED GHOST', true));
    } else {
      row.appendChild(legendSwatch(pal().accent, 'MODE SHAPE'));
      row.appendChild(legendSwatch('var(--ghost)', 'UNDEFORMED GHOST', true));
    }
    row.appendChild(legendSwatch(pal().load, 'APPLIED LOAD', true));
    row.appendChild(legendSwatch(pal().reaction, 'REACTION', true));
    if (state.safety && state.safety.criticalCount) row.appendChild(legendSwatch('var(--critical)', 'DASHED = BELOW SAFETY LIMIT', true));
  }

  /* ---------- mode scrubber ---------- */

  function renderModeScrub() {
    var scrub = $('mode-scrub');
    var show = state.view === 'modes' && state.kind !== 'plate';
    if (!show) {
      scrub.setAttribute('hidden', '');
      stopAnimation();
      var pb = $('mode-play');
      pb.setAttribute('aria-pressed', 'false');
      pb.textContent = 'PLAY';
      return;
    }
    scrub.removeAttribute('hidden');
    var sel = $('mode-select');
    if (sel.getAttribute('data-mode-key') !== state.key || sel.getElementsByTagName('option').length !== state.modes.modes.length) {
      clear(sel);
      state.modes.modes.forEach(function (m, i) {
        var o = hEl('option', null, 'Mode ' + m.index + ' / ' + fixed(m.frequencyHz, 2) + ' Hz');
        o.setAttribute('value', String(i)); sel.appendChild(o);
      });
      sel.setAttribute('data-mode-key', state.key);
    }
    sel.value = String(state.modeIndex);
    $('mode-phase').value = String(Math.round(state.phase));
    var cur = state.modes.modes[state.modeIndex];
    $('mode-readout').textContent = 'f = ' + fixed(cur.frequencyHz, 2) +
      ' Hz / T = ' + fixed(cur.period * 1000, 2) + ' ms / phase ' +
      Math.round(state.phase) + ' deg';
  }

  function stopAnimation() {
    state.playing = false; previousTime = null;
    if (animationFrame !== null && typeof cancelAnimationFrame === 'function') cancelAnimationFrame(animationFrame);
    animationFrame = null;
    if ($('mode-play')) { $('mode-play').setAttribute('aria-pressed', 'false'); $('mode-play').textContent = 'PLAY'; }
  }
  function tickAnimation(time) {
    animationFrame = null;
    if (!state.playing || document.hidden || reducedMotion) { stopAnimation(); return; }
    var elapsed = previousTime === null ? 0 : clamp(time - previousTime, 0, 100);
    previousTime = time;
    state.phase = (state.phase + elapsed * 0.09) % 360;
    $('mode-phase').value = String(Math.round(state.phase));
    var cur = state.modes.modes[state.modeIndex];
    $('mode-readout').textContent = 'f = ' + fixed(cur.frequencyHz, 2) + ' Hz / phase ' + Math.round(state.phase) + ' deg';
    renderCanvas();
    if (typeof requestAnimationFrame === 'function') animationFrame = requestAnimationFrame(tickAnimation);
    else stopAnimation();
  }

  /* ---------- member table ---------- */

  var TRUSS_COLS = [
    { key: 'id', label: 'MEMBER' }, { key: 'length', label: 'L (m)' },
    { key: 'axialForce', label: 'AXIAL N (kN)' }, { key: 'stress', label: 'STRESS (MPa)' },
    { key: 'strain', label: 'STRAIN' }, { key: 'euler', label: 'P_CR (kN)' },
    { key: 'slenderness', label: 'L / r' }, { key: 'fos', label: 'FoS GOV.' },
    { key: 'state', label: 'STATE' }
  ];
  var PLATE_COLS = [
    { key: 'id', label: 'ELEMENT' }, { key: 'area', label: 'AREA (m2)' },
    { key: 'sigmaX', label: 'SIGMA X (MPa)' }, { key: 'sigmaY', label: 'SIGMA Y (MPa)' },
    { key: 'tauMax', label: 'TAU MAX (MPa)' }, { key: 'vonMises', label: 'VON MISES (MPa)' },
    { key: 'sigma1', label: 'SIGMA 1 (MPa)' }, { key: 'centroid', label: 'CENTROID (m)' },
    { key: 'state', label: 'STATE' }
  ];

  function tableRows() {
    if (state.kind === 'plate') {
      var p = state.plate;
      var hi = p.maxVonMises;
      return p.elements.map(function (e) {
        return {
          id: e.id + 1, area: e.area, sigmaX: e.stress.sigmaX,
          sigmaY: e.stress.sigmaY, tauMax: e.tauMax, vonMises: e.vonMises,
          sigma1: e.sigma1, centroid: e.centroid.x, state: e.vonMises / hi,
          critical: e.vonMises > state.model.material.sigmaYield
        };
      });
    }
    var byId = {};
    state.safety.members.forEach(function (r) { byId[r.id] = r; });
    return state.solution.members.map(function (m) {
      var s = byId[m.id] || {};
      return {
        id: m.id + 1, label: m.label, length: m.length,
        axialForce: m.axialForce, stress: m.stress, strain: m.strain,
        euler: s.eulerCriticalLoad, slenderness: s.slendernessRatio,
        fos: s.fosGoverning, critical: !!s.critical,
        tension: m.axialForce > 1e-9, compression: m.axialForce < -1e-9, state: s.fosGoverning
      };
    });
  }

  function filterRows(rows) {
    if (state.kind === 'plate') {
      if (state.filter === 'critical') return rows.filter(function (r) { return r.critical; });
      if (state.filter === 'tension' || state.filter === 'compression') return rows;
      return rows;
    }
    if (state.filter === 'tension') return rows.filter(function (r) { return r.tension; });
    if (state.filter === 'compression') return rows.filter(function (r) { return r.compression; });
    if (state.filter === 'critical') return rows.filter(function (r) { return r.critical; });
    return rows;
  }

  function sortRows(rows) {
    var k = state.sortKey, dir = state.sortDir;
    return rows.slice().sort(function (a, b) {
      var x = a[k], y = b[k];
      if (typeof x === 'string' || typeof y === 'string') {
        return String(x).localeCompare(String(y)) * dir;
      }
      if (!isFinite(x)) x = Infinity;
      if (!isFinite(y)) y = Infinity;
      return (x === y ? 0 : x < y ? -1 : 1) * dir;
    });
  }

  function renderTableHead() {
    var thead = $('member-table');
    var head = thead.getElementsByTagName ? thead.getElementsByTagName('thead')[0] : null;
    if (!head) return;
    clear(head);
    var tr = hEl('tr');
    var cols = state.kind === 'plate' ? PLATE_COLS : TRUSS_COLS;
    cols.forEach(function (c) {
      var th = hEl('th');
      th.setAttribute('scope', 'col');
      th.setAttribute('aria-sort', c.key === state.sortKey ? (state.sortDir === 1 ? 'ascending' : 'descending') : 'none');
      var b = hEl('button', null, c.label);
      b.setAttribute('type', 'button');
      b.setAttribute('data-sort', c.key);
      if (c.key === state.sortKey) b.setAttribute('data-active', 'true');
      (function (key) {
        b.addEventListener('click', function () {
          if (state.sortKey === key) state.sortDir *= -1;
          else { state.sortKey = key; state.sortDir = 1; }
          renderTable();
        });
      })(c.key);
      th.appendChild(b);
      tr.appendChild(th);
    });
    head.appendChild(tr);
  }

  function renderTable() {
    renderTableHead();
    var tbody = $('member-tbody');
    clear(tbody);
    var rows = sortRows(filterRows(tableRows()));
    var cols = state.kind === 'plate' ? PLATE_COLS : TRUSS_COLS;
    rows.forEach(function (r) {
      var tr = hEl('tr');
      if (r.critical) tr.setAttribute('class', 'is-critical');
      cols.forEach(function (c) {
        var td = hEl('td');
        var k = c.key;
        var text = '';
        if (state.kind === 'plate') {
          if (k === 'id') text = 'E' + r.id;
          else if (k === 'area') text = fixed(r.area, 5);
          else if (k === 'centroid') text = fixed(r.centroid, 3);
          else if (k === 'state') text = '';
          else text = fixed(r[k] / 1e6, 3);
        } else {
          if (k === 'id') text = r.label;
          else if (k === 'length') text = fixed(r.length, 3);
          else if (k === 'axialForce') text = fixed(r.axialForce / 1000, 2);
          else if (k === 'stress') text = fixed(r.stress / 1e6, 2);
          else if (k === 'strain') text = fixed(r.strain, 7);
          else if (k === 'euler') text = fixed(r.euler / 1000, 1);
          else if (k === 'slenderness') text = fixed(r.slenderness, 1);
          else if (k === 'fos') text = isFinite(r.fos) ? fixed(r.fos, 2) : 'n/a';
          else text = '';
        }
        if (k === 'state') {
          if (state.kind === 'plate') {
            var hi = state.plate.maxVonMises;
            var color = hi > 0 ? rampColor(r.vonMises / hi) : rampColor(0);
            var mark = hEl('span', 'state-mark');
            mark.setAttribute('style', 'background:' + color + ';');
            var wrap = hEl('span', 'state-cell');
            wrap.appendChild(mark);
            wrap.appendChild(hEl('span', 'state-text',
              r.critical ? 'YIELD LIMIT' : 'BELOW YIELD'));
            td.appendChild(wrap);
          } else {
            var mark2 = hEl('span', 'state-mark');
            mark2.setAttribute('style', 'background:' +
              (r.tension ? pal().tension : r.compression ? pal().compression : pal().neutral) + ';');
            var wrap2 = hEl('span', 'state-cell');
            wrap2.appendChild(mark2);
            wrap2.appendChild(hEl('span', 'state-text', r.critical ? 'CRITICAL' :
              (r.tension ? 'TENSION' : r.compression ? 'COMPRESSION' : 'ZERO FORCE')));
            td.appendChild(wrap2);
          }
        } else {
          td.textContent = text;
          if (state.kind !== 'plate' && (k === 'axialForce' || k === 'stress')) {
            td.setAttribute('class', r.tension ? 'tension' : 'compression');
          }
        }
        tr.appendChild(td);
      });
      tbody.appendChild(tr);
    });
    $('table-foot').textContent = rows.length + ' of ' + tableRows().length +
      ' rows shown / click any column header to sort / all values are also in the .JSON report';
  }

  /* ---------- exports ---------- */

  function download(filename, text, mime) {
    var url;
    try {
      var blob = new Blob([text], { type: mime || 'text/plain' });
      url = URL.createObjectURL(blob);
    } catch (err) {
      url = 'data:' + (mime || 'text/plain') + ';charset=utf-8,' + encodeURIComponent(text);
    }
    var a = document.createElement('a');
    a.setAttribute('href', url);
    a.setAttribute('download', filename);
    if (document.body && document.body.appendChild) document.body.appendChild(a);
    if (typeof a.click === 'function') a.click();
    if (document.body && document.body.removeChild && a.parentNode) document.body.removeChild(a);
    if (url.indexOf('blob:') === 0 && typeof URL.revokeObjectURL === 'function') {
      setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
    }
  }

  function plateInp() {
    var p = state.plate, mesh = p.mesh, cfg = state.model.plate;
    var L = [];
    L.push('*HEADING');
    L.push('MeshLab CST plane stress plate with circular hole');
    L.push('*NODE, NSET=NALL');
    mesh.nodes.forEach(function (n) {
      L.push((n.id + 1) + ', ' + fixed(n.x, 8) + ', ' + fixed(n.y, 8) + ', 0.0');
    });
    L.push('*ELEMENT, TYPE=CPS3, ELSET=EALL');
    mesh.elements.forEach(function (e) {
      L.push((e.id + 1) + ', ' + e.nodes.map(function (n) { return n + 1; }).join(', '));
    });
    L.push('*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL');
    L.push(fixed(p.material.t, 8) + ',');
    L.push('*MATERIAL, NAME=STEEL');
    L.push('*ELASTIC');
    L.push(fixed(p.material.E, 4) + ', ' + fixed(p.material.nu, 4));
    L.push('*BOUNDARY');
    p.fixedMask.forEach(function (fixed, i) { if (fixed) L.push((Math.floor(i / 2) + 1) + ', ' + (i % 2 + 1) + ', ' + (i % 2 + 1) + ', 0.0'); });
    L.push('*STEP');
    L.push('*STATIC');
    L.push('*CLOAD');
    p.forceVector.forEach(function (force, i) { if (force) L.push((Math.floor(i / 2) + 1) + ', ' + (i % 2 + 1) + ', ' + fixed(force, 6)); });
    L.push('*NODE PRINT, NSET=NALL');
    L.push('U');
    L.push('*NODE FILE');
    L.push('U');
    L.push('*EL PRINT, ELSET=EALL');
    L.push('S');
    L.push('*END STEP');
    return L.join('\n') + '\n';
  }

  function plateSvgExport() {
    var p = state.plate, mesh = p.mesh;
    var W = 900, H = 620, pad = 60;
    var unit = Math.min((W - 2 * pad) / mesh.width, (H - 2 * pad) / mesh.height);
    var lo = p.minVonMises, hi = p.maxVonMises;
    var den = Math.max(hi - lo, 1e-9);
    function pr(n) { return [pad + n.x * unit, H - pad - n.y * unit]; }
    var out = [];
    out.push('<?xml version="1.0" encoding="UTF-8"?>');
    out.push('<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H +
      '" viewBox="0 0 ' + W + ' ' + H + '">');
    out.push('<title>MeshLab CST von Mises contour</title>');
    out.push('<desc>' + mesh.elements.length + ' constant strain triangles. Deformation scale ' +
      state.scale + 'x.</desc>');
    out.push('<rect width="100%" height="100%" fill="#09111E"/>');
    p.elements.forEach(function (e) {
      var pts = e.nodes.map(function (k) {
        var n = p.nodes[k];
        var q = pr({ x: n.x + n.ux * state.scale, y: n.y + n.uy * state.scale });
        return fixed(q[0], 2) + ',' + fixed(q[1], 2);
      }).join(' ');
      out.push('<polygon points="' + pts + '" fill="' +
        rampColor((e.vonMises - lo) / den) + '" stroke="#223550" stroke-width="0.4"/>');
    });
    out.push('<g fill="#cfe0f5" font-family="monospace" font-size="13">');
    out.push('<text x="20" y="26">MeshLab / CST plane stress contour</text>');
    out.push('<text x="20" y="46">Peak von Mises ' + fixed(hi / 1e6, 2) + ' MPa</text>');
    out.push('</g>');
    out.push('</svg>');
    return out.join('\n') + '\n';
  }

  function plateJsonExport() {
    var p = state.plate;
    var report = {
      format: 'MeshLab CST plate report v1.0.0',
      units: { length: 'm', force: 'N', stress: 'Pa' },
      mesh: { nx: p.mesh.nx, ny: p.mesh.ny, width: p.mesh.width,
        height: p.mesh.height, holeRadius: p.mesh.holeRadius,
        nodeCount: p.nodes.length, elementCount: p.elements.length },
      material: p.material,
      maxVonMises: p.maxVonMises, minVonMises: p.minVonMises,
      maxDisplacement: p.maxDisplacement, maxSigmaX: p.maxSigmaX,
      equilibrium: p.equilibrium,
      verification: p.verification,
      diagnostics: { method: p.diagnostics.method, iterations: p.diagnostics.iterations, relativeResidual: p.diagnostics.relativeResidual },
      forceVector: p.forceVector, fixedMask: p.fixedMask,
      nodes: p.nodes.map(function (n) {
        return { id: n.id, x: n.x, y: n.y, ux: n.ux, uy: n.uy,
          rx: n.rx, ry: n.ry, vonMises: n.vonMises };
      }),
      elements: p.elements.map(function (e) {
        return { id: e.id, nodes: e.nodes, area: e.area,
          centroid: e.centroid, strain: e.strain, stress: e.stress,
          sigma1: e.sigma1, sigma2: e.sigma2, tauMax: e.tauMax,
          vonMises: e.vonMises, strainEnergy: e.strainEnergy };
      })
    };
    return JSON.stringify(report, null, 2) + '\n';
  }

  function vtkData() {
    if (state.kind === 'field') {
      return { mesh: state.field.mesh, data: { title: 'MeshLab ' + state.model.name,
        pointScalars: [{ name: 'temperature', values: state.field.values }],
        pointVectors: [{ name: 'nodal_flux', values: state.field.nodeFlux }],
        cellVectors: [{ name: 'heat_flux', values: state.field.elements.map(function (e) { return [e.fluxX, e.fluxY, 0]; }) }] } };
    }
    if (state.kind === 'plate') {
      var p = state.plate;
      return { mesh: { nodes: p.mesh.nodes, elements: p.mesh.elements }, data: { title: 'MeshLab ' + state.model.name,
        pointScalars: [{ name: 'vonMises', values: p.nodes.map(function (n) { return n.vonMises; }) }],
        pointVectors: [{ name: 'displacement', values: p.nodes.map(function (n) { return [n.ux, n.uy, 0]; }) }],
        cellScalars: [{ name: 'vonMises_element', values: p.elements.map(function (e) { return e.vonMises; }) }] } };
    }
    var d = state.model.dim;
    return { mesh: { nodes: state.model.nodes, elements: state.solution.members.map(function (m) { return { nodes: [m.a, m.b] }; }) },
      data: { title: 'MeshLab ' + state.model.name,
        pointVectors: [{ name: 'displacement', values: state.solution.nodes.map(function (n) {
          return d === 3 ? n.displacement : n.displacement.concat([0]);
        }) }],
        cellScalars: [{ name: 'axial_force', values: state.solution.members.map(function (m) { return m.axialForce; }) },
          { name: 'axial_stress', values: state.solution.members.map(function (m) { return m.stress; }) }] } };
  }

  function doExport(kind) {
    var base = state.key;
    if (kind === 'vtk') {
      var v = vtkData();
      download(base + '.vtk', M.exportVtk(v.mesh, v.data), 'text/plain');
      return;
    }
    if (state.kind === 'field') {
      if (kind === 'json') download(base + '.json', JSON.stringify({ format: 'MeshLab Poisson report v2.0.0',
        units: { length: 'm', temperature: 'C', flux: 'W/m2' }, model: state.model, mesh: state.field.mesh, values: state.field.values,
        elements: state.field.elements, verification: state.field.verification,
        diagnostics: { method: state.field.diagnostics.method, iterations: state.field.diagnostics.iterations, relativeResidual: state.field.diagnostics.relativeResidual } }, null, 2) + '\n', 'application/json');
      else if (kind === 'svg') download(base + '.svg', fieldSvgExport(), 'image/svg+xml');
      return;
    }
    if (kind === 'inp') {
      var text = state.kind === 'plate' ? plateInp() : M.exportCalculiXInp(state.model);
      download(base + '.inp', text, 'text/plain');
    } else if (kind === 'svg') {
      var svg = state.kind === 'plate' ? plateSvgExport() :
        M.exportDeformedSvg(state.model, state.solution, state.scale);
      download(base + '.svg', svg, 'image/svg+xml');
    } else {
      var json = state.kind === 'plate' ? plateJsonExport() :
        M.exportReportJson(state.model, state.solution);
      download(base + '.json', json, 'application/json');
    }
  }

  /* ---------- full render ---------- */

  function renderAll() {
    $('determinacy-title').textContent = state.kind === 'plate' ? 'DISCRETIZATION' : 'DETERMINACY';
    $('safety-title').textContent = 'SAFETY ENVELOPE';
    $('modes-title').textContent = 'NATURAL FREQUENCIES';
    $('members-title').textContent = state.kind === 'plate' ? 'ELEMENT STRESSES & YIELD SAFETY' : 'MEMBER FORCES & BUCKLING SAFETY';
    $('preset-select').value = state.key;
    $('view-select').value = state.view;
    $('view-select').disabled = state.kind === 'field';
    $('member-filter').value = state.filter;
    $('load-range').value = String(state.load);
    $('load-output').textContent = state.load.toFixed(1) + 'x';
    var scaleR = $('scale-range');
    if (state.kind === 'field') {
      scaleR.setAttribute('min', '0.5'); scaleR.setAttribute('max', '3'); scaleR.setAttribute('step', '0.25');
      scaleR.value = String(state.resolution);
      scaleR.setAttribute('aria-label', 'Mesh resolution multiplier from 0.5x to 3x');
      $('scale-output').textContent = state.resolution + 'x';
      $('scale-label').textContent = 'MESH RESOLUTION';
      scaleR.disabled = false;
    } else {
      scaleR.setAttribute('min', '1'); scaleR.setAttribute('max', '500'); scaleR.setAttribute('step', '1');
      scaleR.value = String(state.scale);
      scaleR.setAttribute('aria-label', 'Deformation exaggeration from 1x to 500x');
      $('scale-output').textContent = state.scale + 'x';
      $('scale-label').textContent = state.view === 'modes' ? 'MODE DISPLAY AMPLITUDE' : 'DEFORMATION SCALE';
      scaleR.disabled = state.view !== 'deformed' && state.view !== 'modes';
    }
    $('load-label').textContent = state.kind === 'field'
      ? (state.model.field && state.model.field.heater !== undefined ? 'HEATER MULTIPLIER' : 'HOT EDGE MULTIPLIER')
      : 'LOAD MULTIPLIER';
    $('reactions-label').textContent = state.kind === 'field' ? 'Heat flux' : 'Reactions';
    $('view-select').getElementsByTagName('option')[1].textContent = state.kind === 'field' ? 'Temperature field' : 'von Mises stress';
    Array.prototype.forEach.call($('view-select').getElementsByTagName('option'), function (option) { option.disabled = state.kind === 'plate' && (option.value === 'force' || option.value === 'modes'); });
    $('toggle-ghost').disabled = state.kind === 'field' || (state.kind === 'plate' && state.view !== 'deformed');
    var filterOptions = $('member-filter').getElementsByTagName('option');
    var labels = state.kind === 'field' ? ['All triangles', 'At least 50 C', 'Below 50 C', 'Above 90% of peak'] :
      state.kind === 'plate' ? ['All triangles', 'All stress states', 'All stress states', 'Above yield limit'] : ['All members', 'Tension only', 'Compression only', 'Critical only'];
    Array.prototype.forEach.call(filterOptions, function (option, i) { option.textContent = labels[i]; option.disabled = state.kind === 'plate' && (i === 1 || i === 2); });
    Array.prototype.forEach.call(document.getElementsByClassName('dl-btn'), function (b) { b.disabled = state.kind === 'field' && b.getAttribute('data-export') === 'inp'; });
    if (state.kind === 'field') { renderFieldAll(); renderVerification(); renderModeScrub(); resetInspection(); return; }
    renderHeader();
    renderStats();
    renderEquilibrium();
    renderDeterminacy();
    renderSafety();
    renderModeBars();
    renderLegend();
    renderModeScrub();
    renderCanvas();
    renderTable();
    resetInspection(); renderVerification();
  }

  function syncViewSelect() {
    var sel = $('view-select');
    if (sel) sel.value = state.view;
  }

  /* ---------- scalar field experiment and visible verification ---------- */

  var HEAT_STOPS = ['#244970', '#287d9d', '#63c5cc', '#f0d68f', '#ed8c72'];
  var HEAT_RGB = null;
  function heatStops() {
    if (!HEAT_RGB) HEAT_RGB = HEAT_STOPS.map(function (h) {
      return [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
    });
    return HEAT_RGB;
  }
  function heatColor(value) {
    var t = clamp((value - state.field.min) / Math.max(state.field.max - state.field.min, 1e-12), 0, 1);
    var stops = heatStops(), x = t * (stops.length - 1), i0 = Math.floor(x), f = x - i0;
    var a = stops[i0], b = stops[Math.min(i0 + 1, stops.length - 1)];
    return 'rgb(' + Math.round(a[0] + (b[0] - a[0]) * f) + ',' +
      Math.round(a[1] + (b[1] - a[1]) * f) + ',' + Math.round(a[2] + (b[2] - a[2]) * f) + ')';
  }

  var DOMAIN_INFO = {
    rect: { title: 'HOT LEFT EDGE / COLD RIGHT EDGE / INTERIOR HEATER / ARROW LENGTH = |q|',
      boundary: 'Left/right prescribed', patch: 'rect' },
    l: { title: 'HOT LEFT EDGE / COLD RIGHT EDGE / INTERIOR HEATER / ARROW LENGTH = |q|',
      boundary: 'Left/right prescribed', patch: 'none' },
    annulus: { title: 'HOT INNER RING / COLD OUTER RING / PURE CONDUCTION',
      boundary: 'Inner/outer ring prescribed', patch: 'annulus' },
    sink: { title: 'HOT BASE / COLD FIN SURFACES / INTERIOR HEATER / ARROW LENGTH = |q|',
      boundary: 'Base + fin walls prescribed', patch: 'none' }
  };
  function fieldDomain() {
    return (state.field.mesh.domainSpec && DOMAIN_INFO[state.field.mesh.domainSpec.type]) || DOMAIN_INFO.rect;
  }

  function fieldProjection() {
    var mesh = state.field.mesh;
    var minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    mesh.nodes.forEach(function (n) {
      minX = Math.min(minX, n.x); maxX = Math.max(maxX, n.x);
      minY = Math.min(minY, n.y); maxY = Math.max(maxY, n.y);
    });
    var w = Math.max(maxX - minX, 1e-9), h = Math.max(maxY - minY, 1e-9);
    var unit = Math.min(880 / w, 500 / h);
    var ox = (1000 - w * unit) / 2, oy = (640 - h * unit) / 2;
    return function (n) { return [ox + (n.x - minX) * unit, 640 - oy - (n.y - minY) * unit]; };
  }

  function fieldNodeRows(node, i) {
    var fl = state.field.nodeFlux[i] || { x: 0, y: 0 };
    var q = Math.hypot(fl.x, fl.y);
    return [
      { name: 'x', value: fixed(node.x, 4) + ' m' },
      { name: 'y', value: fixed(node.y, 4) + ' m' },
      { name: 'Temperature', value: fixed(state.field.values[i], 3) + ' C', color: heatColor(state.field.values[i]) },
      { name: '|grad T|', value: fixed(q / Math.max(state.field.conductivity, 1e-12), 3) + ' K/m' },
      { name: '|flux q|', value: fixed(q, 3) + ' W/m2' }
    ];
  }

  function renderFieldCanvas() {
    var svg = $('structure-svg'); clear(svg); inspectionTargets = [];
    svg.setAttribute('viewBox', '0 0 1000 640');
    svg.setAttribute('aria-label', 'Temperature and heat flux on a Delaunay finite element mesh');
    var projection = fieldProjection(), field = state.field, mesh = field.mesh;
    var hitLayer = svgEl('g');
    var maxMag = 0;
    field.elements.forEach(function (e) { maxMag = Math.max(maxMag, Math.hypot(e.fluxX, e.fluxY)); });
    var stride = Math.max(1, Math.round(field.elements.length / 26));
    var midT = (field.min + field.max) / 2;
    field.elements.forEach(function (e, i) {
      var points = e.nodes.map(function (id) { var p = projection(mesh.nodes[id]); return fixed(p[0], 2) + ',' + fixed(p[1], 2); }).join(' ');
      var fill = heatColor(e.value);
      svg.appendChild(svgEl('polygon', { points: points, fill: fill, stroke: state.grid ? '#263b55' : fill, 'stroke-width': state.grid ? 0.7 : 0.2 }));
      if (state.reactions && i % stride === 0) {
        var m2 = Math.hypot(e.fluxX, e.fluxY);
        if (m2 > maxMag * 0.02) drawArrow(svg, projection(e.centroid), [e.fluxX, e.fluxY],
          e.value < midT ? '#e6f5ff' : '#182b3f', null, false, 'svg-field-flux', 10 + 26 * Math.pow(m2 / maxMag, 0.6));
      }
      var hit = svgEl('polygon', { points: points, 'class': 'svg-plate-hit', tabindex: '-1', role: 'button', 'aria-label': 'Triangle ' + (e.id + 1) + ', temperature ' + fixed(e.value, 2) + ' degrees C' });
      attachTip(hit, 'Field triangle ' + (e.id + 1), function () { return [
        { name: 'Temperature', value: fixed(e.value, 2) + ' C' },
        { name: 'Flux x', value: signed(e.fluxX, 2) + ' W/m2' },
        { name: 'Flux y', value: signed(e.fluxY, 2) + ' W/m2' },
        { name: 'Area', value: fixed(e.area, 5) + ' m2' }
      ]; });
      hitLayer.appendChild(hit);
    });
    svg.appendChild(hitLayer);
    if (mesh.nodes.length <= 800) {
      var nodeG = svgEl('g');
      mesh.nodes.forEach(function (n, i) {
        var p = projection(n);
        var dot = svgEl('circle', { cx: fixed(p[0], 2), cy: fixed(p[1], 2), r: 2.4, 'class': 'svg-field-node' });
        (function (node, idx, dotEl) {
          dotEl.addEventListener('pointerenter', function (evt) {
            var rows = fieldNodeRows(node, idx);
            showTooltip('Node ' + (idx + 1), rows, evt);
            if (state.pinned < 0) renderInspection('FIELD NODE ' + (idx + 1), rows);
          });
          dotEl.addEventListener('pointerleave', hideTooltip);
          dotEl.addEventListener('click', function () {
            var rows = fieldNodeRows(node, idx);
            renderInspection('FIELD NODE ' + (idx + 1), rows);
          });
        })(n, i, dot);
        nodeG.appendChild(dot);
      });
      svg.appendChild(nodeG);
    }
    var title = svgEl('text', { x: 20, y: 28, 'class': 'svg-title' });
    title.textContent = fieldDomain().title;
    svg.appendChild(title);
  }
  function fieldSvgExport() {
    var project = fieldProjection(), out = ['<?xml version="1.0" encoding="UTF-8"?>', '<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="640" viewBox="0 0 1000 640"><title>MeshLab Poisson temperature field</title><rect width="1000" height="640" fill="#0e1a2a"/>'];
    state.field.elements.forEach(function (e) {
      var points = e.nodes.map(function (id) { var p = project(state.field.mesh.nodes[id]); return fixed(p[0], 2) + ',' + fixed(p[1], 2); }).join(' ');
      out.push('<polygon points="' + points + '" fill="' + heatColor(e.value) + '" stroke="#263b55" stroke-width="0.5"/>');
    });
    var cfg = state.model.field;
    out.push('<text x="20" y="28" fill="#eef4fc" font-family="monospace" font-size="14">Poisson FEM / temperature C / hot edge ' +
      fixed(cfg.hotEdge, 0) + ' / cold edge ' + fixed(cfg.coldEdge || 0, 0) + '</text></svg>');
    return out.join('\n') + '\n';
  }
  function fieldTable() {
    var columns = [{ key: 'id', label: 'TRIANGLE' }, { key: 'value', label: 'TEMP (C)' }, { key: 'fluxX', label: 'FLUX X (W/m2)' }, { key: 'fluxY', label: 'FLUX Y (W/m2)' }, { key: 'area', label: 'AREA (m2)' }];
    var head = $('member-table').getElementsByTagName('thead')[0]; clear(head);
    var tr = hEl('tr');
    columns.forEach(function (c) {
      var th = hEl('th'); th.setAttribute('scope', 'col'); th.setAttribute('aria-sort', state.sortKey === c.key ? (state.sortDir > 0 ? 'ascending' : 'descending') : 'none');
      var button = hEl('button', null, c.label); button.setAttribute('type', 'button');
      button.addEventListener('click', function () { state.sortDir = state.sortKey === c.key ? -state.sortDir : 1; state.sortKey = c.key; fieldTable(); });
      th.appendChild(button); tr.appendChild(th);
    }); head.appendChild(tr);
    var rows = state.field.elements.slice();
    if (state.filter === 'tension') rows = rows.filter(function (e) { return e.value >= 50; });
    if (state.filter === 'compression') rows = rows.filter(function (e) { return e.value < 50; });
    if (state.filter === 'critical') rows = rows.filter(function (e) { return e.value > state.field.max * 0.9; });
    rows = sortRows(rows); var body = $('member-tbody'); clear(body);
    rows.forEach(function (row) { var tr = hEl('tr'); columns.forEach(function (c) { tr.appendChild(hEl('td', null, c.key === 'id' ? 'E' + (row.id + 1) : fixed(row[c.key], c.key === 'area' ? 5 : 3))); }); body.appendChild(tr); });
    $('table-foot').textContent = rows.length + ' of ' + state.field.elements.length + ' triangles / ' + fieldDomain().boundary + ' / conductivity ' + fixed(state.field.conductivity, 2) + ' W/(m K)';
  }
  function elementMinAngle(mesh, e) {
    var p = e.nodes.map(function (id) { return mesh.nodes[id]; }), worst = 180;
    for (var i = 0; i < 3; i++) {
      var o = p[i], u = p[(i + 1) % 3], v = p[(i + 2) % 3];
      var ax = u.x - o.x, ay = u.y - o.y, bx = v.x - o.x, by = v.y - o.y;
      var cs = (ax * bx + ay * by) / Math.max(Math.sqrt((ax * ax + ay * ay) * (bx * bx + by * by)), 1e-30);
      worst = Math.min(worst, Math.acos(clamp(cs, -1, 1)) * 180 / Math.PI);
    }
    return worst;
  }

  function renderFieldAll() {
    var field = state.field, mesh = field.mesh, diagnostics = field.diagnostics;
    var dom = fieldDomain(), cfg = state.model.field;
    $('model-title').textContent = state.model.name;
    $('determinacy-title').textContent = 'DISCRETIZATION'; $('safety-title').textContent = 'MESH QUALITY'; $('modes-title').textContent = 'PCG CONVERGENCE'; $('members-title').textContent = 'TEMPERATURE & HEAT FLUX';
    var grid = $('stat-grid'); clear(grid);
    grid.appendChild(statTile('PEAK TEMPERATURE', fixed(field.max, 2), 'C', fixed(cfg.hotEdge, 0) + ' C hot / ' + fixed(cfg.coldEdge, 0) + ' C cold', true));
    grid.appendChild(statTile('DELAUNAY TRIANGLES', String(field.elements.length), '', mesh.nodes.length + ' nodes'));
    grid.appendChild(statTile('PCG ITERATIONS', String(diagnostics.iterations), '', diagnostics.nonzeros + ' stored coefficients'));
    grid.appendChild(statTile('RELATIVE RESIDUAL', diagnostics.relativeResidual.toExponential(2), '', 'actual residual, checked at convergence'));
    grid.appendChild(statTile('MIN TRIANGLE ANGLE', fixed(mesh.quality.minimumAngle, 1), 'deg', 'max radius-edge ' + fixed(mesh.quality.maxRadiusEdge || 0, 2)));
    if (cfg.heater !== undefined) {
      grid.appendChild(statTile('HEATER INPUT', fixed(field.verification.sourceTotal, 2), 'W/m', 'heat input per unit depth'));
    } else {
      grid.appendChild(statTile('DOMAIN AREA', fixed(mesh.quality.area, 4), 'm2', 'Euler V-E+F ' + (mesh.chi === undefined ? 1 : mesh.chi)));
    }
    var equilibrium = $('equilibrium-table').getElementsByTagName('tbody')[0]; clear(equilibrium);
    appendRow(equilibrium, 'Integrated source', fixed(field.verification.sourceTotal, 6) + ' W/m');
    appendRow(equilibrium, 'Boundary flux', fixed(field.verification.boundaryFlux, 6) + ' W/m');
    appendRow(equilibrium, 'Balance residual', field.verification.balanceResidual.toExponential(2) + ' W/m');
    appendRow(equilibrium, 'Free node residual', field.verification.freeResidual.toExponential(2) + ' W/m');
    $('equilibrium-badge').textContent = diagnostics.converged ? 'BALANCED' : 'CHECK'; $('equilibrium-badge').setAttribute('data-state', diagnostics.converged ? 'good' : 'critical');
    var meshTable = $('determinacy-table').getElementsByTagName('tbody')[0]; clear(meshTable);
    appendRow(meshTable, 'Nodes / triangles', mesh.nodes.length + ' / ' + mesh.elements.length); appendRow(meshTable, 'Free scalar DOFs', String(field.freeDofs.length));
    appendRow(meshTable, 'Total mesh area', fixed(mesh.quality.area, 6) + ' m2');
    appendRow(meshTable, 'Euler V-E+F', String(mesh.chi === undefined ? 1 : mesh.chi) + (dom.patch === 'annulus' ? ' (one hole)' : ''));
    appendRow(meshTable, 'Minimum angle', fixed(mesh.quality.minimumAngle, 2) + ' deg'); appendRow(meshTable, 'Boundary', dom.boundary);
    var meter = $('safety-meter'); clear(meter);
    meter.appendChild(hEl('p', 'panel-sub', 'Field min ' + fixed(field.min, 2) + ' C / max ' + fixed(field.max, 2) + ' C. Multiplier ' + state.load.toFixed(1) + 'x. Element minimum-angle histogram below.'));
    var bins = [0, 0, 0, 0, 0, 0], peak = 1;
    field.elements.forEach(function (e) { bins[Math.min(5, Math.floor(elementMinAngle(mesh, e) / 15))]++; });
    bins.forEach(function (b) { peak = Math.max(peak, b); });
    bins.forEach(function (count, i) {
      var row = hEl('div', 'mode-bar-row');
      row.appendChild(hEl('span', 'mode-bar-idx', (i * 15) + '-' + ((i + 1) * 15)));
      var track = hEl('div', 'mode-bar-track'), fill = hEl('div', 'mode-bar-fill');
      fill.setAttribute('style', 'width:' + fixed(count / peak * 100, 2) + '%;' + (i < 2 ? 'background:var(--warning);' : ''));
      track.appendChild(fill); row.appendChild(track);
      row.appendChild(hEl('span', 'mode-bar-val', String(count)));
      meter.appendChild(row);
    });
    meter.appendChild(hEl('p', 'panel-sub', 'Degrees / ' + field.elements.length + ' triangles / min ' + fixed(mesh.quality.minimumAngle, 1) + ' deg'));
    $('safety-badge').textContent = 'THERMAL'; $('safety-badge').removeAttribute('data-state');
    var bars = $('mode-bars'); clear(bars);
    var hist = (diagnostics.history && diagnostics.history.length > 1) ? diagnostics.history : [1, diagnostics.relativeResidual];
    var W2 = 268, H2 = 64, pad = 5, lastX = 0, lastY = 0;
    var pts = hist.map(function (v, i) {
      var lg = Math.max(-16, Math.log10(Math.max(v, 1e-30)));
      lastX = pad + (W2 - 2 * pad) * (hist.length === 1 ? 1 : i / (hist.length - 1));
      lastY = pad + (H2 - 2 * pad) * (1 + lg / 16);
      return fixed(lastX, 1) + ',' + fixed(lastY, 1);
    }).join(' ');
    var spark = svgEl('svg', { 'class': 'conv-spark', viewBox: '0 0 ' + W2 + ' ' + H2, role: 'img', 'aria-label': 'PCG relative residual versus iteration, log scale' });
    spark.appendChild(svgEl('rect', { 'class': 'conv-frame', x: 0.5, y: 0.5, width: W2 - 1, height: H2 - 1 }));
    spark.appendChild(svgEl('polyline', { points: pts, 'class': 'conv-line' }));
    spark.appendChild(svgEl('circle', { 'class': 'conv-dot', cx: fixed(lastX, 1), cy: fixed(lastY, 1), r: 2.6 }));
    bars.appendChild(spark);
    bars.appendChild(hEl('div', 'conv-meta', 'rel residual ' + diagnostics.relativeResidual.toExponential(2) + ' in ' + diagnostics.iterations + ' iters / log10 axis 1 to 1e-16'));
    var cacheKey = state.key + '|' + mesh.nodes.length + '|' + mesh.elements.length;
    var checkText = 'No closed form on this domain. Flux balance and the PCG residual certify the solve.';
    if (dom.patch === 'rect' || dom.patch === 'annulus') {
      if (!(cacheKey in patchCheckCache)) {
        var hot2 = new Set(mesh.bcSets.hot), cold2 = new Set(mesh.bcSets.cold);
        var check = M.solvePoisson(mesh, { conductivity: 1, boundary: function (n, i) { if (hot2.has(i)) return 100; if (cold2.has(i)) return 0; } });
        var err = 0;
        mesh.nodes.forEach(function (n, i) {
          var exact;
          if (dom.patch === 'rect') exact = 100 * (1 - n.x / mesh.width);
          else {
            var rr = Math.hypot(n.x - mesh.meta.center.x, n.y - mesh.meta.center.y);
            exact = 100 * Math.log(rr / mesh.meta.rOut) / Math.log(mesh.meta.rIn / mesh.meta.rOut);
          }
          err = Math.max(err, Math.abs(check.values[i] - exact));
        });
        patchCheckCache[cacheKey] = err;
      }
      checkText = dom.patch === 'rect'
        ? 'With the heater off, T = 100(1 - x / width) is the exact solution. Maximum nodal error on this irregular mesh: ' + patchCheckCache[cacheKey].toExponential(2) + ' C.'
        : 'T(r) = 100 ln(r/Rout) / ln(Rin/Rout) is the exact radial solution. Maximum nodal error on this mesh: ' + patchCheckCache[cacheKey].toExponential(2) + ' C.';
    }
    bars.appendChild(hEl('p', 'panel-sub', checkText));
    var legend = $('legend-row'); clear(legend);
    var scale2 = hEl('div', 'legend-scale');
    scale2.appendChild(hEl('span', 'legend-tick', fixed(field.min, 0) + ' C'));
    var grad = hEl('div', 'legend-ramp legend-grad');
    for (var gi = 0; gi < 24; gi++) {
      var seg = hEl('i'); seg.setAttribute('style', 'background:' + heatColor(field.min + (field.max - field.min) * gi / 23) + ';');
      grad.appendChild(seg);
    }
    scale2.appendChild(grad);
    scale2.appendChild(hEl('span', 'legend-tick', fixed(field.max, 0) + ' C'));
    legend.appendChild(scale2);
    legend.appendChild(hEl('span', 'legend-label', 'ARROW LENGTH = |q| / DOTS = NODES'));
    fieldTable(); renderFieldCanvas();
  }
  function renderVerification() {
    var v = state.kind === 'field' ? state.field.verification : (state.kind === 'plate' ? state.plate.verification : state.solution.verification);
    var relative = state.kind === 'field' ? state.field.diagnostics.relativeResidual : v.relativeResidual;
    var passed = state.kind === 'field' ? state.field.diagnostics.converged : v.passed;
    $('verification-summary').textContent = (passed ? 'VERIFIED' : 'CHECK') + ' / residual ' + relative.toExponential(2);
    $('verification-summary').setAttribute('data-state', passed ? 'good' : 'critical');
    $('status-message').textContent = state.kind === 'field' ? 'PCG converged. Source and boundary heat flux balance. Run the round-trip check to reproduce the result.' :
      'Free DOF residual ' + v.relativeResidual.toExponential(2) + '. Energy identity error ' + v.energyError.toExponential(2) + '. Run the round-trip check to reproduce the result.';
  }
  function runRoundTrip() {
    try {
      var duplicate = JSON.parse(JSON.stringify(state.model)), expected, actual;
      if (state.kind === 'field') { expected = state.field.values; actual = solveField(duplicate.field.mesh, duplicate.field).values; }
      else if (state.kind === 'plate') {
        expected = state.plate.displacements;
        actual = M.solvePlateCST(duplicate.plate.mesh, duplicate.plate.material, duplicate.plate.bc, duplicate.plate.loads).displacements;
      } else { expected = state.solution.displacements; actual = M.solveTruss(duplicate).displacements; }
      var error = 0, peak = 0;
      expected.forEach(function (value, i) { error = Math.max(error, Math.abs(value - actual[i])); peak = Math.max(peak, Math.abs(value)); });
      var relative = error / Math.max(peak, 1e-20), passed = relative < 1e-8;
      $('status-message').textContent = (passed ? 'PASS' : 'CHECK') + ' / serialized model round-trip / relative error ' + relative.toExponential(3) +
        (state.kind === 'plate' ? ' / independent Gaussian solver comparison' : ' / reproduced numerical solution');
    } catch (error) { $('status-message').textContent = 'CHECK / ' + error.message; }
  }

  /* ---------- wiring ---------- */

  function init() {
    var presetSel = $('preset-select');
    if (!presetSel) return;
    M.getPresets().forEach(function (p) {
      var o = hEl('option', null, p.name);
      o.setAttribute('value', p.key);
      presetSel.appendChild(o);
    });
    presetSel.value = state.key;

    presetSel.addEventListener('change', function (e) {
      state.filter = 'all';
      state.resolution = 1;
      Array.prototype.forEach.call(document.getElementsByClassName('scenario-btn'), function (b) { b.setAttribute('aria-pressed', 'false'); });
      updateAnalysis((e.target && e.target.value) || state.key);
    });
    $('view-select').addEventListener('change', function (e) {
      state.view = (e.target && e.target.value) || 'force';
      stopAnimation();
      if ((state.view === 'modes' || state.view === 'force') && state.kind === 'plate') {
        state.view = 'stress'; syncViewSelect();
      }
      renderAll();
    });
    $('scale-range').addEventListener('input', function (e) {
      var v = Number((e.target && e.target.value) || 100);
      if (!isFinite(v)) return;
      if (state.kind === 'field') {
        state.resolution = clamp(v, 0.5, 3);
        $('scale-output').textContent = state.resolution + 'x';
        updateAnalysis(state.key);
        return;
      }
      state.scale = clamp(v, 1, 500);
      $('scale-output').textContent = state.scale + 'x';
      renderStats();
      renderCanvas();
    });
    $('load-range').addEventListener('input', function (e) {
      var load = Number(e.target.value);
      if (!isFinite(load)) return;
      state.load = clamp(load, 0, 5); updateAnalysis(state.key);
    });
    $('run-checks').addEventListener('click', runRoundTrip);
    Array.prototype.forEach.call(document.getElementsByClassName('scenario-btn'), function (button) {
      button.addEventListener('click', function () {
        var scenarios = { bridge: ['pratt', 'force', 1, 100], crane: ['crane', 'deformed', 5, 50], plate: ['plate', 'stress', 1, 100], heat: ['heat', 'stress', 1, 100] };
        var scenario = scenarios[button.getAttribute('data-scenario')];
        state.view = scenario[1]; state.load = scenario[2]; state.scale = scenario[3]; state.filter = 'all';
        state.resolution = 1;
        updateAnalysis(scenario[0]);
        Array.prototype.forEach.call(document.getElementsByClassName('scenario-btn'), function (b) { b.setAttribute('aria-pressed', b === button ? 'true' : 'false'); });
      });
    });
    $('toggle-ghost').addEventListener('change', function (e) {
      state.ghost = !!(e.target && e.target.checked); renderLegend(); renderCanvas();
    });
    $('toggle-reactions').addEventListener('change', function (e) {
      state.reactions = !!(e.target && e.target.checked); renderCanvas();
    });
    $('toggle-grid').addEventListener('change', function (e) {
      state.grid = !!(e.target && e.target.checked); renderCanvas();
    });
    $('member-filter').addEventListener('change', function (e) {
      state.filter = (e.target && e.target.value) || 'all'; if (state.kind === 'field') fieldTable(); else renderTable();
    });
    $('mode-select').addEventListener('change', function (e) {
      var index = Number((e.target && e.target.value) || 0);
      if (!isFinite(index) || !state.modes || !state.modes.modes.length) return;
      state.modeIndex = clamp(Math.floor(index), 0, state.modes.modes.length - 1);
      renderModeScrub(); renderModeBars(); renderCanvas();
    });
    $('mode-phase').addEventListener('input', function (e) {
      var phase = Number((e.target && e.target.value) || 0);
      if (!isFinite(phase)) return;
      state.phase = clamp(phase, 0, 360);
      renderModeScrub(); renderCanvas();
    });
    $('mode-play').addEventListener('click', function () {
      if (reducedMotion) { $('status-message').textContent = 'Reduced motion is enabled. Use the phase slider to inspect each mode.'; return; }
      if (state.playing) { stopAnimation(); return; }
      state.playing = true; previousTime = null;
      var b = $('mode-play');
      b.setAttribute('aria-pressed', state.playing ? 'true' : 'false');
      b.textContent = state.playing ? 'PAUSE' : 'PLAY';
      if (typeof requestAnimationFrame === 'function' && animationFrame === null) animationFrame = requestAnimationFrame(tickAnimation);
    });
    $('theme-toggle').addEventListener('click', function () {
      state.theme = state.theme === 'dark' ? 'light' : 'dark';
      var root = document.documentElement;
      if (root && root.setAttribute) root.setAttribute('data-theme', state.theme);
      var b = $('theme-toggle');
      b.setAttribute('aria-pressed', state.theme === 'dark' ? 'true' : 'false');
      $('theme-label').textContent = state.theme === 'dark' ? 'LIGHT MODE' : 'DARK MODE';
      b.setAttribute('title', state.theme === 'dark' ?
        'Switch to light blueprint mode' : 'Switch to dark slate mode');
      if (state.kind === 'field') renderFieldAll();
      else { renderLegend(); renderCanvas(); renderTable(); renderModeBars(); }
    });
    var dlButtons = document.getElementsByClassName
      ? document.getElementsByClassName('dl-btn') : [];
    Array.prototype.forEach.call(dlButtons, function (b) {
      b.addEventListener('click', function () {
        doExport(b.getAttribute('data-export'));
      });
    });

    document.addEventListener('visibilitychange', function () { if (document.hidden) stopAnimation(); });
    updateAnalysis(state.key);
  }

  init();
})();
