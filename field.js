/* MeshLab scalar FEM tools. SI units, deterministic meshes, no dependencies. */
(function (root, factory) {
  'use strict';
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.MeshField = factory();
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';

  function finite(v, label) {
    if (typeof v !== 'number' || !Number.isFinite(v)) throw new Error('MeshLab: ' + label + ' must be finite');
    return v;
  }
  function dot(a, b) { var s = 0; for (var i = 0; i < a.length; i++) s += a[i] * b[i]; return s; }

  function toCSR(rows) {
    var n = rows.length, ptr = [0], cols = [], vals = [], diagonal = new Float64Array(n);
    rows.forEach(function (row, i) {
      if (!(row instanceof Map) && (!row || row.length !== n)) throw new Error('MeshLab: matrix dimensions do not match');
      var entries = row instanceof Map ? Array.from(row.entries()) :
        Array.from(row).map(function (v, j) { return [j, v]; });
      entries.sort(function (a, b) { return a[0] - b[0]; });
      entries.forEach(function (entry) {
        var j = entry[0], value = finite(entry[1], 'matrix entry');
        if (!Number.isInteger(j) || j < 0 || j >= n) throw new Error('MeshLab: matrix column out of range');
        if (value !== 0) { cols.push(j); vals.push(value); }
        if (j === i) diagonal[i] = value;
      });
      ptr.push(cols.length);
    });
    return { size: n, rowPtr: Int32Array.from(ptr), columns: Int32Array.from(cols), values: Float64Array.from(vals), diagonal: diagonal };
  }
  function multiply(A, x) {
    var out = new Float64Array(A.size);
    for (var i = 0; i < A.size; i++) for (var k = A.rowPtr[i]; k < A.rowPtr[i + 1]; k++) out[i] += A.values[k] * x[A.columns[k]];
    return out;
  }
  function solvePCG(input, rhs, options) {
    var A = input.rowPtr ? input : toCSR(input), n = A.size, opt = options || {};
    if (!Number.isInteger(n) || n < 0 || A.rowPtr.length !== n + 1 || A.values.length !== A.columns.length ||
        A.rowPtr[0] !== 0 || A.rowPtr[n] !== A.values.length || !A.diagonal || A.diagonal.length !== n) throw new Error('MeshLab: invalid sparse matrix');
    var sparseRows = [];
    for (var row = 0; row < n; row++) {
      if (A.rowPtr[row] > A.rowPtr[row + 1]) throw new Error('MeshLab: invalid sparse row pointers');
      var entries = new Map();
      for (var entry = A.rowPtr[row]; entry < A.rowPtr[row + 1]; entry++) {
        var column = A.columns[entry], value = finite(A.values[entry], 'sparse coefficient');
        if (!Number.isInteger(column) || column < 0 || column >= n) throw new Error('MeshLab: sparse column out of range');
        entries.set(column, value);
      }
      sparseRows.push(entries);
    }
    sparseRows.forEach(function (entries, row) {
      entries.forEach(function (value, column) {
        var mirror = sparseRows[column].get(row) || 0;
        if (Math.abs(value - mirror) > 1e-12 * Math.max(Math.abs(value), Math.abs(mirror), 1e-30)) throw new Error('MeshLab: PCG matrix must be symmetric');
      });
    });
    if (rhs.length !== n) throw new Error('MeshLab: right hand side dimension mismatch');
    var b = Float64Array.from(rhs), x = new Float64Array(n), r = Float64Array.from(b), z = new Float64Array(n), p = new Float64Array(n);
    var tolerance = opt.tolerance === undefined ? 1e-11 : finite(opt.tolerance, 'PCG tolerance');
    if (!(tolerance > 0 && tolerance < 1)) throw new Error('MeshLab: PCG tolerance must be between zero and one');
    var limit = opt.maxIterations === undefined ? Math.max(32, n * 4) : opt.maxIterations;
    if (!Number.isInteger(limit) || limit < 1) throw new Error('MeshLab: PCG iteration limit must be a positive integer');
    for (var i = 0; i < n; i++) {
      finite(b[i], 'right hand side');
      if (!(A.diagonal[i] > 0)) throw new Error('MeshLab: PCG requires a positive diagonal');
      z[i] = r[i] / A.diagonal[i]; p[i] = z[i];
    }
    var normB = Math.sqrt(dot(b, b)), rz = dot(r, z), residual = normB, iteration = 0;
    var target = tolerance * (normB || 1);
    while (residual > target && iteration < limit) {
      var Ap = multiply(A, p), curvature = dot(p, Ap);
      if (!(curvature > 0) || !Number.isFinite(curvature)) throw new Error('MeshLab: PCG encountered a mechanism or non positive stiffness');
      var alpha = rz / curvature;
      for (i = 0; i < n; i++) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
      iteration++;
      /* Recompute the actual residual before accepting convergence. */
      residual = Math.sqrt(dot(r, r));
      if (residual <= target || iteration % 40 === 0) {
        var Ax = multiply(A, x);
        for (i = 0; i < n; i++) r[i] = b[i] - Ax[i];
        residual = Math.sqrt(dot(r, r));
      }
      if (residual <= target) break;
      for (i = 0; i < n; i++) z[i] = r[i] / A.diagonal[i];
      var next = dot(r, z), beta = next / rz;
      for (i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
      rz = next;
    }
    var finalProduct = multiply(A, x);
    for (i = 0; i < n; i++) r[i] = b[i] - finalProduct[i];
    residual = Math.sqrt(dot(r, r));
    return { solution: x, converged: residual <= target, iterations: iteration,
      residual: residual, relativeResidual: residual / (normB || 1), nonzeros: A.values.length, method: 'Jacobi PCG' };
  }

  function triangulate(points) {
    if (!Array.isArray(points)) throw new Error('MeshLab: triangulation needs points[]');
    var nodes = [], seen = new Set(), duplicates = 0;
    points.forEach(function (p) {
      if (!p || typeof p !== 'object') throw new Error('MeshLab: invalid triangulation point');
      var x = finite(Array.isArray(p) ? p[0] : p.x, 'point x');
      var y = finite(Array.isArray(p) ? p[1] : p.y, 'point y');
      var key = x + ',' + y;
      if (seen.has(key)) { duplicates++; return; }
      seen.add(key); nodes.push({ id: nodes.length, x: x, y: y });
    });
    if (nodes.length < 3) throw new Error('MeshLab: triangulation needs three distinct points');
    var xs = nodes.map(function (p) { return p.x; }), ys = nodes.map(function (p) { return p.y; });
    var xmin = Math.min.apply(null, xs), ymin = Math.min.apply(null, ys);
    var width = Math.max.apply(null, xs) - xmin, height = Math.max.apply(null, ys) - ymin;
    var span = Math.max(width, height);
    if (!(span > 0) || !Number.isFinite(span)) throw new Error('MeshLab: triangulation has invalid extent');
    var pts = nodes.map(function (p) { return { x: (p.x - xmin) / span, y: (p.y - ymin) / span }; });
    var count = pts.length;
    var bestCross = 0, p0 = pts[0];
    for (var bi = 1; bi < count - 1; bi++) {
      for (var bj = bi + 1; bj < count; bj++) {
        var cVal = Math.abs((pts[bi].x - p0.x) * (pts[bj].y - p0.y) - (pts[bj].x - p0.x) * (pts[bi].y - p0.y));
        if (cVal > bestCross) bestCross = cVal;
      }
    }
    if (bestCross < 1e-13) throw new Error('MeshLab: points are collinear');

    function orient(a, b, c) {
      var pa = pts[a], pb = pts[b], pc = pts[c];
      return (pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x);
    }

    function inCircle(a, b, c, d) {
      var pa = pts[a], pb = pts[b], pc = pts[c], pd = pts[d];
      var ax = pa.x - pd.x, ay = pa.y - pd.y;
      var bx = pb.x - pd.x, by = pb.y - pd.y;
      var cx = pc.x - pd.x, cy = pc.y - pd.y;
      var a2 = ax * ax + ay * ay;
      var b2 = bx * bx + by * by;
      var c2 = cx * cx + cy * cy;
      var det = ax * (by * c2 - b2 * cy) - ay * (bx * c2 - b2 * cx) + a2 * (bx * cy - by * cx);
      var s = a2 + b2 + c2;
      return det > 1e-12 * s * s;
    }

    var order = Array.from({ length: count }, function (_, i) { return i; });
    order.sort(function (a, b) {
      return pts[a].x !== pts[b].x ? pts[a].x - pts[b].x : pts[a].y - pts[b].y;
    });

    var seedEnd = 2;
    while (seedEnd < count) {
      var oBase = Math.abs(orient(order[0], order[1], order[seedEnd]));
      var d01 = Math.hypot(pts[order[1]].x - pts[order[0]].x, pts[order[1]].y - pts[order[0]].y);
      if (oBase > 1e-13 * d01) break;
      seedEnd++;
    }
    if (seedEnd >= count) {
      var maxAreaIdx = 2, maxAreaVal = 0;
      for (var si = 2; si < count; si++) {
        var av = Math.abs(orient(order[0], order[1], order[si]));
        if (av > maxAreaVal) { maxAreaVal = av; maxAreaIdx = si; }
      }
      if (maxAreaIdx !== 2) {
        var tmpIdx = order[maxAreaIdx];
        order.splice(maxAreaIdx, 1);
        order.splice(2, 0, tmpIdx);
      }
      seedEnd = 2;
    }

    var tris = [];
    var hull = [];
    var apex = order[seedEnd];
    var seedSign = orient(order[0], order[seedEnd - 1], apex);
    if (seedSign > 0) {
      for (var cj = 0; cj < seedEnd - 1; cj++) {
        tris.push([order[cj], order[cj + 1], apex]);
      }
      for (var hj = 0; hj <= seedEnd; hj++) hull.push(order[hj]);
    } else {
      for (var cj2 = 0; cj2 < seedEnd - 1; cj2++) {
        tris.push([order[cj2 + 1], order[cj2], apex]);
      }
      hull.push(apex);
      for (var hj2 = seedEnd - 1; hj2 >= 0; hj2--) hull.push(order[hj2]);
    }

    for (var oi = seedEnd + 1; oi < count; oi++) {
      var pid = order[oi], m = hull.length, vis = new Array(m), anyVis = false;
      for (var e = 0; e < m; e++) {
        var u = hull[e], v = hull[(e + 1) % m];
        var ov = orient(u, v, pid);
        var scaleEps = 1e-14 * (Math.abs(pts[v].x - pts[u].x) + Math.abs(pts[v].y - pts[u].y) +
          Math.abs(pts[pid].x - pts[v].x) + Math.abs(pts[pid].y - pts[v].y));
        vis[e] = ov < -scaleEps;
        if (vis[e]) anyVis = true;
      }
      if (!anyVis) {
        for (var e2 = 0; e2 < m; e2++) {
          vis[e2] = orient(hull[e2], hull[(e2 + 1) % m], pid) < 0;
          if (vis[e2]) anyVis = true;
        }
      }
      if (!anyVis) continue;
      var start = -1;
      for (var e3 = 0; e3 < m; e3++) {
        if (vis[e3] && !vis[(e3 - 1 + m) % m]) { start = e3; break; }
      }
      if (start === -1) start = 0;
      var visCount = 0;
      while (visCount < m && vis[(start + visCount) % m]) {
        var eu = hull[(start + visCount) % m], ev = hull[(start + visCount + 1) % m];
        tris.push([eu, pid, ev]);
        visCount++;
      }
      var nextHull = [pid];
      var endVertexIdx = (start + visCount) % m;
      for (var keepStep = 0; keepStep <= m - visCount; keepStep++) {
        nextHull.push(hull[(endVertexIdx + keepStep) % m]);
      }
      hull = nextHull;
    }

    var edgeOwner = new Map();
    function edgeKey(u, v) { return u * count + v; }
    function registerTri(ti) {
      var t = tris[ti];
      edgeOwner.set(edgeKey(t[0], t[1]), { t: ti, opp: t[2] });
      edgeOwner.set(edgeKey(t[1], t[2]), { t: ti, opp: t[0] });
      edgeOwner.set(edgeKey(t[2], t[0]), { t: ti, opp: t[1] });
    }
    for (var ti = 0; ti < tris.length; ti++) registerTri(ti);

    var queue = [];
    for (var qti = 0; qti < tris.length; qti++) {
      var qt = tris[qti];
      queue.push([qt[0], qt[1]], [qt[1], qt[2]], [qt[2], qt[0]]);
    }
    var maxFlips = Math.max(100, count * count * 4), flips = 0;
    while (queue.length && flips < maxFlips) {
      var pair = queue.pop(), eu = pair[0], ev = pair[1];
      var left = edgeOwner.get(edgeKey(eu, ev));
      var right = edgeOwner.get(edgeKey(ev, eu));
      if (!left || !right) continue;
      var cVert = left.opp, dVert = right.opp;
      if (inCircle(eu, ev, cVert, dVert) && orient(dVert, cVert, eu) > 1e-15 && orient(cVert, dVert, ev) > 1e-15) {
        var t1 = left.t, t2 = right.t;
        edgeOwner.delete(edgeKey(eu, ev));
        edgeOwner.delete(edgeKey(ev, eu));
        tris[t1] = [dVert, cVert, eu];
        tris[t2] = [cVert, dVert, ev];
        registerTri(t1);
        registerTri(t2);
        queue.push([eu, dVert], [dVert, ev], [ev, cVert], [cVert, eu]);
        flips++;
      }
    }

    var elements = tris.map(function (t, id) { return { id: id, nodes: t }; });
    if (!elements.length) throw new Error('MeshLab: points are collinear');
    var edgeCounts = new Map(), minimumAngle = 180, totalArea = 0;
    elements.forEach(function (el) {
      var p = el.nodes.map(function (id) { return pts[id]; });
      var area = Math.abs((p[1].x - p[0].x) * (p[2].y - p[0].y) - (p[2].x - p[0].x) * (p[1].y - p[0].y)) / 2;
      totalArea += area * span * span;
      for (var e = 0; e < 3; e++) {
        var a = el.nodes[e], b = el.nodes[(e + 1) % 3], key = Math.min(a, b) + ':' + Math.max(a, b);
        edgeCounts.set(key, (edgeCounts.get(key) || 0) + 1);
        var u = p[(e + 1) % 3], v = p[(e + 2) % 3], origin = p[e];
        var ux = u.x - origin.x, uy = u.y - origin.y, vx = v.x - origin.x, vy = v.y - origin.y;
        var cosine = (ux * vx + uy * vy) / Math.sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy));
        minimumAngle = Math.min(minimumAngle, Math.acos(Math.max(-1, Math.min(1, cosine))) * 180 / Math.PI);
      }
    });
    var boundarySet = new Set();
    edgeCounts.forEach(function (count, key) { if (count === 1) key.split(':').forEach(function (id) { boundarySet.add(Number(id)); }); });
    return { nodes: nodes, elements: elements, boundaryNodes: Array.from(boundarySet).sort(function (a, b) { return a - b; }),
      width: width, height: height, quality: { minimumAngle: minimumAngle, area: totalArea }, duplicatesRemoved: duplicates };
  }

  function generateFieldMesh(nx, ny, width, height) {
    nx = nx === undefined ? 14 : nx; ny = ny === undefined ? 10 : ny;
    width = width === undefined ? 2.4 : width; height = height === undefined ? 1.6 : height;
    if (!Number.isInteger(nx) || !Number.isInteger(ny) || nx < 2 || ny < 2 || nx > 60 || ny > 60 || !(width > 0 && height > 0)) {
      throw new Error('MeshLab: invalid field mesh dimensions');
    }
    var pts = [];
    for (var j = 0; j <= ny; j++) for (var i = 0; i <= nx; i++) {
      var boundary = i === 0 || j === 0 || i === nx || j === ny;
      pts.push({ x: width * (i + (boundary ? 0 : 0.18 * Math.sin(i * 17 + j * 31))) / nx,
        y: height * (j + (boundary ? 0 : 0.18 * Math.cos(i * 23 + j * 13))) / ny });
    }
    var mesh = triangulate(pts);
    mesh.nx = nx; mesh.ny = ny; mesh.leftNodes = []; mesh.rightNodes = [];
    mesh.nodes.forEach(function (n) { if (n.x === 0) mesh.leftNodes.push(n.id); if (Math.abs(n.x - width) < 1e-12) mesh.rightNodes.push(n.id); });
    return mesh;
  }

  function solvePoisson(mesh, config) {
    var cfg = config || {}, conductivity = cfg.conductivity === undefined ? 1 : finite(cfg.conductivity, 'conductivity');
    if (!(conductivity > 0)) throw new Error('MeshLab: conductivity must be positive');
    if (!mesh || !Array.isArray(mesh.nodes) || !Array.isArray(mesh.elements) || !mesh.nodes.length || !mesh.elements.length) throw new Error('MeshLab: scalar FEM needs a mesh');
    var n = mesh.nodes.length, rows = Array.from({ length: n }, function () { return new Map(); }), b = new Float64Array(n), data = [];
    var parent = Array.from({ length: n }, function (_, i) { return i; });
    function component(i) { while (parent[i] !== i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; }
    mesh.elements.forEach(function (el) {
      if (el.nodes.length !== 3 || el.nodes.some(function (id) { return !Number.isInteger(id) || !mesh.nodes[id]; })) throw new Error('MeshLab: invalid scalar triangle');
      var p = el.nodes.map(function (id) { return mesh.nodes[id]; });
      parent[component(el.nodes[1])] = component(el.nodes[0]); parent[component(el.nodes[2])] = component(el.nodes[0]);
      var cross = (p[1].x - p[0].x) * (p[2].y - p[0].y) - (p[2].x - p[0].x) * (p[1].y - p[0].y);
      if (!Number.isFinite(cross) || Math.abs(cross) <= 1e-20) throw new Error('MeshLab: degenerate scalar triangle');
      var area = Math.abs(cross) / 2;
      var gx = [(p[1].y - p[2].y) / cross, (p[2].y - p[0].y) / cross, (p[0].y - p[1].y) / cross];
      var gy = [(p[2].x - p[1].x) / cross, (p[0].x - p[2].x) / cross, (p[1].x - p[0].x) / cross];
      var centroid = { x: (p[0].x + p[1].x + p[2].x) / 3, y: (p[0].y + p[1].y + p[2].y) / 3 };
      var source = typeof cfg.source === 'function' ? cfg.source(centroid.x, centroid.y) : (cfg.source === undefined ? 0 : cfg.source);
      finite(source, 'heat source');
      el.nodes.forEach(function (a, i) {
        b[a] += source * area / 3;
        el.nodes.forEach(function (c, j) { rows[a].set(c, (rows[a].get(c) || 0) + conductivity * area * (gx[i] * gx[j] + gy[i] * gy[j])); });
      });
      data.push({ nodes: el.nodes, gx: gx, gy: gy, area: area, centroid: centroid });
    });
    var prescribed = new Map();
    mesh.nodes.forEach(function (node, i) {
      var value = typeof cfg.boundary === 'function' ? cfg.boundary(node, i) : cfg.boundary && cfg.boundary[i];
      if (value !== undefined && value !== null) prescribed.set(i, finite(value, 'boundary value'));
    });
    if (!prescribed.size) throw new Error('MeshLab: Poisson problem needs a prescribed boundary value');
    var anchored = new Set(); prescribed.forEach(function (_, id) { anchored.add(component(id)); });
    for (var node = 0; node < n; node++) if (!anchored.has(component(node))) throw new Error('MeshLab: every mesh component needs a prescribed boundary value');
    var free = [], index = new Map();
    for (var i = 0; i < n; i++) if (!prescribed.has(i)) { index.set(i, free.length); free.push(i); }
    var reduced = free.map(function (a) {
      var row = new Map();
      rows[a].forEach(function (v, c) { if (index.has(c)) row.set(index.get(c), v); });
      return row;
    });
    var rhs = free.map(function (a) {
      var value = b[a]; rows[a].forEach(function (v, c) { if (prescribed.has(c)) value -= v * prescribed.get(c); }); return value;
    });
    var result = solvePCG(reduced, rhs, cfg), values = new Float64Array(n);
    if (!result.converged) throw new Error('MeshLab: Poisson PCG did not converge');
    prescribed.forEach(function (v, id) { values[id] = v; }); free.forEach(function (id, i) { values[id] = result.solution[i]; });
    var elements = data.map(function (d, i) {
      var dx = 0, dy = 0, value = 0;
      d.nodes.forEach(function (id, j) { dx += d.gx[j] * values[id]; dy += d.gy[j] * values[id]; value += values[id] / 3; });
      return { id: i, nodes: d.nodes.slice(), centroid: d.centroid, area: d.area, value: value,
        fluxX: -conductivity * dx, fluxY: -conductivity * dy, energy: conductivity * d.area * (dx * dx + dy * dy) / 2 };
    });
    var residuals = multiply(toCSR(rows), values), reaction = 0, sourceTotal = 0, freeResidual = 0;
    for (i = 0; i < n; i++) {
      residuals[i] -= b[i]; sourceTotal += b[i];
      if (prescribed.has(i)) reaction += residuals[i]; else freeResidual = Math.max(freeResidual, Math.abs(residuals[i]));
    }
    return { mesh: mesh, values: Array.from(values), elements: elements,
      min: Math.min.apply(null, Array.from(values)), max: Math.max.apply(null, Array.from(values)),
      energy: elements.reduce(function (s, e) { return s + e.energy; }, 0), diagnostics: result,
      verification: { freeResidual: freeResidual, balanceResidual: sourceTotal + reaction, sourceTotal: sourceTotal, boundaryFlux: reaction },
      freeDofs: free, conductivity: conductivity };
  }

  return { toCSR: toCSR, multiply: multiply, solvePCG: solvePCG, triangulate: triangulate,
    generateFieldMesh: generateFieldMesh, solvePoisson: solvePoisson };
});
