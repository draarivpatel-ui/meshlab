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
  var STATUS = { good: '#0ca30c', warning: '#fab219', serious: '#ec835a', critical: '#d03b3b' };

  var state = {
    key: 'pratt', kind: 'truss', model: null, solution: null,
    safety: null, modes: null, plate: null,
    view: 'force', scale: 100, ghost: true, reactions: true, grid: true,
    theme: 'dark', modeIndex: 0, phase: 0, playing: false,
    filter: 'all', sortKey: 'id', sortDir: 1, hovered: -1
  };

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
    state.model = model;
    state.plate = null; state.solution = null;
    state.safety = null; state.modes = null;
    if (state.kind === 'plate') {
      state.plate = M.solvePlateCST(model.plate.mesh, model.plate.material,
        model.plate.bc, model.plate.loads);
      if (state.view === 'modes') { state.view = 'stress'; syncViewSelect(); }
      clear($('mode-select'));
    } else {
      var sol = M.solveTruss(model);
      state.solution = sol;
      state.safety = M.analyzeSafety(sol.model, sol);
      state.modes = M.solveModes(model, 4);
    }
    state.modeIndex = 0;
    state.phase = 0;
    state.hovered = -1;
    state.sortKey = 'id';
    state.sortDir = 1;
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
    var s = state.scale;
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
          var amp = Math.max(state.solution.maxDisplacement, 1e-9) * 0.55;
          var mx = modeDisplacement(amp, d * i), my = modeDisplacement(amp, d * i + 1),
            mz = d === 3 ? modeDisplacement(amp, d * i + 2) : 0;
          pts.push(rawPoint(ns[i], mx, my, mz));
          pts.push(rawPoint(ns[i], -mx, -my, -mz));
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
      var amp = Math.max(state.solution.maxDisplacement, 1e-9) * 0.55;
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
    var x = (evt && evt.clientX ? evt.clientX - rect.left : 24) + 14;
    var y = (evt && evt.clientY ? evt.clientY - rect.top : 24) + 14;
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
    function open(evt) { showTooltip(title, rowsFn(), evt); }
    ['pointermove', 'pointerenter', 'focus'].forEach(function (ev) {
      el.addEventListener(ev, open);
    });
    ['pointerleave', 'blur'].forEach(function (ev) {
      el.addEventListener(ev, hideTooltip);
    });
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

  function drawArrow(g, at, vec, color, label, into, textClass) {
    var mag = Math.sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
    if (!(mag > 1e-12)) return;
    var len = 34;
    var ux = vec[0] / mag, uy = -vec[1] / mag;
    var tail = into ? [at[0] - ux * len, at[1] - uy * len] : [at[0], at[1]];
    var tip = into ? [at[0], at[1]] : [at[0] + ux * len, at[1] + uy * len];
    var line = svgEl('line', { x1: fixed(tail[0], 2), y1: fixed(tail[1], 2),
      x2: fixed(tip[0], 2), y2: fixed(tip[1], 2), stroke: color,
      'stroke-width': 2, 'stroke-linecap': 'round' });
    line.setAttribute('class', textClass === 'svg-reaction-text' ? 'svg-reaction' : 'svg-load');
    g.appendChild(line);
    var bx = tip[0] - ux * 9, by = tip[1] - uy * 9;
    var px = -uy * 4.5, py = ux * 4.5;
    var pts = fixed(tip[0], 2) + ',' + fixed(tip[1], 2) + ' ' +
      fixed(bx + px, 2) + ',' + fixed(by + py, 2) + ' ' +
      fixed(bx - px, 2) + ',' + fixed(by - py, 2);
    var head = svgEl('polygon', { points: pts, fill: color });
    head.setAttribute('class', textClass === 'svg-reaction-text' ? 'svg-reaction-head' : 'svg-load-head');
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
      line.setAttribute('class', 'svg-member' + (state.hovered === i ? ' hovered' : ''));
      memberG.appendChild(line);

      var hit = svgEl('line', { x1: fixed(q0[0], 2), y1: fixed(q0[1], 2),
        x2: fixed(q1[0], 2), y2: fixed(q1[1], 2),
        'class': 'svg-member-hit', tabindex: '0', role: 'img',
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
            { value: mm.axialForce >= 0 ? 'TENSION' : 'COMPRESSION', name: 'State' }
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
        var amp = Math.max(sol.maxDisplacement, 1e-9) * 0.55;
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
      var lx = n.load[0], ly = n.load[1];
      if (Math.abs(lx) + Math.abs(ly) > 1e-9) {
        drawArrow(arrowG, base, [lx / maxL, ly / maxL], pal().load,
          fixed(Math.sqrt(lx * lx + ly * ly) / 1000, 1) + ' kN', true, 'svg-label');
      }
      if (state.reactions && (Math.abs(n.reaction[0]) + Math.abs(n.reaction[1])) > 1e-6) {
        drawArrow(arrowG, base, [n.reaction[0] / maxR, n.reaction[1] / maxR],
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
        tabindex: '0', role: 'img',
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
    });
    svg.appendChild(triG);
    svg.appendChild(edgeG);
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
    clear(svg);
    var pr = buildProjection();
    svg.setAttribute('viewBox', '0 0 ' + pr.W + ' ' + pr.H);
    svg.setAttribute('aria-label', (state.kind === 'plate' ? 'Plate' : 'Truss') +
      ' result view: ' + state.view);
    drawGrid(svg, pr);
    if (state.kind === 'plate') renderPlateCanvas(svg, pr);
    else renderTrussCanvas(svg, pr);
  }

  /* ---------- header and stat tiles ---------- */

  function renderHeader() {
    var p = M.getPresets().filter(function (x) { return x.key === state.key; })[0];
    $('model-title').textContent = p.name;
    var caption;
    if (state.kind === 'plate') {
      caption = '2D plane stress CST mesh. Kirsch stress concentration around a circular hole. SI units.';
    } else if (state.model.dim === 3) {
      caption = '3D space frame with isometric projection. Direct stiffness solution, SI units.';
    } else {
      caption = '2D pin jointed truss. Direct stiffness solution, SI units.';
    }
    $('model-caption').textContent = caption;
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
      var ratio = nominal > 0 ? p.maxVonMises / nominal : 0;
      g.appendChild(statTile('PEAK VON MISES', fixed(p.maxVonMises / 1e6, 2), 'MPa',
        'element maximum, sequential ramp', true));
      g.appendChild(statTile('KIRSCH RATIO', fixed(ratio, 3), 'x',
        'peak / nominal membrane stress'));
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
      'exaggerated ' + state.scale + 'x in the canvas', true));
    g.appendChild(statTile('MIN FACTOR OF SAFETY',
      isFinite(sf.minFos) ? fixed(sf.minFos, 3) : 'n/a', '',
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
    var res = Math.max(Math.abs(e.residualX || 0), Math.abs(e.residualY || 0));
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
      var fos = 2.5e8 / p.maxVonMises;
      var top = hEl('div', 'meter-top');
      top.appendChild(hEl('span', null, 'Yield FoS at peak stress'));
      top.appendChild(hEl('strong', null, fixed(fos, 3)));
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
      st.appendChild(hEl('span', null, 'Peak / nominal ratio ' + fixed(ratio, 3) +
        '. Stress concentration at the hole boundary.'));
      box.appendChild(st);
      badge.textContent = fos >= 1 ? 'OK' : 'CRITICAL';
      badge.setAttribute('data-state', fos >= 1 ? 'good' : 'critical');
      return;
    }
    var sf = state.safety;
    var fos = sf.minFos;
    var top = hEl('div', 'meter-top');
    top.appendChild(hEl('span', null, 'Governing minimum'));
    top.appendChild(hEl('strong', null, isFinite(fos) ? fixed(fos, 3) : 'n/a'));
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
        'Eigenmodes are computed for the truss and frame presets. The plate view uses the CST contour.'));
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
      if (state.ghost) row.appendChild(legendSwatch('var(--ghost)', 'UNDEFORMED GHOST', true));
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
  }

  /* ---------- mode scrubber ---------- */

  function renderModeScrub() {
    var scrub = $('mode-scrub');
    var show = state.view === 'modes' && state.kind !== 'plate';
    if (!show) {
      scrub.setAttribute('hidden', '');
      state.playing = false;
      var pb = $('mode-play');
      pb.setAttribute('aria-pressed', 'false');
      pb.textContent = 'PLAY';
      return;
    }
    scrub.removeAttribute('hidden');
    var sel = $('mode-select');
    clear(sel);
    state.modes.modes.forEach(function (m, i) {
      var o = hEl('option', null, 'Mode ' + m.index + ' / ' + fixed(m.frequencyHz, 2) + ' Hz');
      o.setAttribute('value', String(i));
      if (i === state.modeIndex) o.setAttribute('selected', 'selected');
      sel.appendChild(o);
    });
    sel.setAttribute('value', String(state.modeIndex));
    $('mode-phase').setAttribute('value', String(Math.round(state.phase)));
    var cur = state.modes.modes[state.modeIndex];
    $('mode-readout').textContent = 'f = ' + fixed(cur.frequencyHz, 2) +
      ' Hz / T = ' + fixed(cur.period * 1000, 2) + ' ms / phase ' +
      Math.round(state.phase) + ' deg';
  }

  function tickAnimation() {
    if (!state.playing) return;
    state.phase = (state.phase + 5) % 360;
    renderModeScrub();
    renderCanvas();
    if (typeof requestAnimationFrame === 'function') requestAnimationFrame(tickAnimation);
    else state.playing = false;
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
          critical: e.vonMises > 0.7 * hi
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
        tension: m.axialForce > 0, state: s.fosGoverning
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
    if (state.filter === 'compression') return rows.filter(function (r) { return !r.tension; });
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
      return (x - y) * dir;
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
            var color = rampColor(r.vonMises / hi);
            var mark = hEl('span', 'state-mark');
            mark.setAttribute('style', 'background:' + color + ';');
            var wrap = hEl('span', 'state-cell');
            wrap.appendChild(mark);
            wrap.appendChild(hEl('span', 'state-text',
              r.critical ? 'HIGH' : 'NOMINAL'));
            td.appendChild(wrap);
          } else {
            var mark2 = hEl('span', 'state-mark');
            mark2.setAttribute('style', 'background:' +
              (r.tension ? pal().tension : pal().compression) + ';');
            var wrap2 = hEl('span', 'state-cell');
            wrap2.appendChild(mark2);
            wrap2.appendChild(hEl('span', 'state-text', r.critical ? 'CRITICAL' :
              (r.tension ? 'TENSION' : 'COMPRESSION')));
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
    $('table-caption').textContent = state.kind === 'plate'
      ? 'Per triangle CST stress state. Sequential ramp shows von Mises magnitude. Sortable.'
      : 'Per member axial force, stress and governing factor of safety. Sortable.';
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
    (cfg.bc.leftNodesFixed || mesh.leftNodes).forEach(function (n) {
      L.push((n + 1) + ', 1, 1, 0.0');
    });
    (cfg.bc.fixedNodes || []).forEach(function (n) {
      L.push((n + 1) + ', 1, 2, 0.0');
    });
    L.push('*STEP');
    L.push('*STATIC');
    L.push('*CLOAD');
    function edge(list, total) {
      if (!list || list.length < 2 || !total) return;
      var q = total / mesh.height;
      for (var i = 0; i < list.length - 1; i++) {
        var a = list[i], b = list[i + 1];
        if (a < 0 || b < 0) continue;
        var seg = Math.abs(mesh.nodes[b].y - mesh.nodes[a].y);
        var half = q * seg / 2;
        if (half) {
          L.push((a + 1) + ', 1, ' + fixed(half, 6));
          L.push((b + 1) + ', 1, ' + fixed(half, 6));
        }
      }
    }
    edge(mesh.rightNodes, cfg.loads.rightNodesForce);
    edge(mesh.leftNodes, cfg.loads.leftNodesForce);
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

  function doExport(kind) {
    var base = state.key;
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
  }

  function syncViewSelect() {
    var sel = $('view-select');
    if (sel) sel.setAttribute('value', state.view);
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
    presetSel.setAttribute('value', state.key);

    presetSel.addEventListener('change', function (e) {
      compute((e.target && e.target.value) || state.key);
      renderAll();
    });
    $('view-select').addEventListener('change', function (e) {
      state.view = (e.target && e.target.value) || 'force';
      if (state.view === 'modes' && state.kind === 'plate') {
        state.view = 'stress'; syncViewSelect();
      }
      renderAll();
    });
    $('scale-range').addEventListener('input', function (e) {
      state.scale = Number((e.target && e.target.value) || 100);
      $('scale-output').textContent = state.scale + 'x';
      renderStats();
      renderCanvas();
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
      state.filter = (e.target && e.target.value) || 'all'; renderTable();
    });
    $('mode-select').addEventListener('change', function (e) {
      state.modeIndex = Number((e.target && e.target.value) || 0);
      renderModeScrub(); renderModeBars(); renderCanvas();
    });
    $('mode-phase').addEventListener('input', function (e) {
      state.phase = Number((e.target && e.target.value) || 0);
      renderModeScrub(); renderCanvas();
    });
    $('mode-play').addEventListener('click', function () {
      state.playing = !state.playing;
      var b = $('mode-play');
      b.setAttribute('aria-pressed', state.playing ? 'true' : 'false');
      b.textContent = state.playing ? 'PAUSE' : 'PLAY';
      if (state.playing) tickAnimation();
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
      renderLegend(); renderCanvas(); renderTable(); renderModeBars();
    });
    var dlButtons = document.getElementsByClassName
      ? document.getElementsByClassName('dl-btn') : [];
    Array.prototype.forEach.call(dlButtons, function (b) {
      b.addEventListener('click', function () {
        doExport(b.getAttribute('data-export'));
      });
    });

    compute(state.key);
    renderAll();
  }

  init();
})();
