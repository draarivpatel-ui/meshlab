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
    var normB = Math.sqrt(dot(b, b)), rz = dot(r, z), residual = normB, iteration = 0, history = [1];
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
      history.push(residual / (normB || 1));
      if (residual <= target) break;
      for (i = 0; i < n; i++) z[i] = r[i] / A.diagonal[i];
      var next = dot(r, z), beta = next / rz;
      for (i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
      rz = next;
    }
    var finalProduct = multiply(A, x);
    for (i = 0; i < n; i++) r[i] = b[i] - finalProduct[i];
    residual = Math.sqrt(dot(r, r));
    history.push(residual / (normB || 1));
    return { solution: x, converged: residual <= target, iterations: iteration,
      residual: residual, relativeResidual: residual / (normB || 1), nonzeros: A.values.length, method: 'Jacobi PCG', history: history };
  }

  function pointInPolygon(x, y, polygon) {
    var inside = false;
    for (var i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
      var xi = polygon[i][0], yi = polygon[i][1], xj = polygon[j][0], yj = polygon[j][1];
      if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) inside = !inside;
    }
    return inside;
  }

  function triangulate(points, opts) {
    opts = opts || {};
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

    /* Domain filtering: drop triangles whose centroid leaves the polygon or enters a hole. */
    var keepTris = tris;
    var outerPoly = null, holePolys = [];
    if (opts.polygon) {
      if (!Array.isArray(opts.polygon) || opts.polygon.length < 3) throw new Error('MeshLab: polygon domain needs three vertices');
      outerPoly = opts.polygon.map(function (v) { return [finite(Array.isArray(v) ? v[0] : v.x, 'polygon x'), finite(Array.isArray(v) ? v[1] : v.y, 'polygon y')]; });
    }
    if (opts.holes) {
      if (!Array.isArray(opts.holes)) throw new Error('MeshLab: holes must be polygons[]');
      opts.holes.forEach(function (hole) {
        if (!Array.isArray(hole) || hole.length < 3) throw new Error('MeshLab: each hole needs three vertices');
        holePolys.push(hole.map(function (v) { return [finite(Array.isArray(v) ? v[0] : v.x, 'hole x'), finite(Array.isArray(v) ? v[1] : v.y, 'hole y')]; }));
      });
    }
    if (outerPoly || holePolys.length) {
      /* A triangle survives only if its centroid and all three edge midpoints stay in the
         domain. The midpoint test stops Delaunay edges that bridge notches or gaps (an edge
         whose endpoints are both on the boundary can still cross void space). Midpoints that
         sit on a polygon edge are kept via a small distance tolerance. */
      var bx0 = Infinity, bx1 = -Infinity, by0 = Infinity, by1 = -Infinity;
      nodes.forEach(function (n) {
        if (n.x < bx0) bx0 = n.x; if (n.x > bx1) bx1 = n.x;
        if (n.y < by0) by0 = n.y; if (n.y > by1) by1 = n.y;
      });
      var epsClip = 1e-7 * Math.max(bx1 - bx0, by1 - by0, 1e-9);
      var distToPoly = function (px, py, poly) {
        var best = Infinity;
        for (var i = 0, j = poly.length - 1; i < poly.length; j = i++) {
          var ax = poly[j][0], ay = poly[j][1], dx = poly[i][0] - ax, dy = poly[i][1] - ay;
          var L = dx * dx + dy * dy;
          var tt = L ? Math.max(0, Math.min(1, ((px - ax) * dx + (py - ay) * dy) / L)) : 0;
          var ex = px - (ax + tt * dx), ey = py - (ay + tt * dy);
          if (ex * ex + ey * ey < best) best = ex * ex + ey * ey;
        }
        return Math.sqrt(best);
      };
      var inDomain = function (px, py) {
        if (outerPoly && !pointInPolygon(px, py, outerPoly) && distToPoly(px, py, outerPoly) > epsClip) return false;
        for (var h = 0; h < holePolys.length; h++)
          if (pointInPolygon(px, py, holePolys[h]) && distToPoly(px, py, holePolys[h]) > epsClip) return false;
        return true;
      };
      keepTris = tris.filter(function (t) {
        var cx = 0, cy = 0;
        t.forEach(function (id) { cx += nodes[id].x; cy += nodes[id].y; });
        cx /= 3; cy /= 3;
        if (!inDomain(cx, cy)) return false;
        for (var e = 0; e < 3; e++) {
          var mx = (nodes[t[e]].x + nodes[t[(e + 1) % 3]].x) / 2, my = (nodes[t[e]].y + nodes[t[(e + 1) % 3]].y) / 2;
          if (!inDomain(mx, my)) return false;
        }
        return true;
      });
      /* Stair-step cleanup: keep the largest edge-connected component, then drop degenerate
         boundary slivers until none remain. Removing a boundary sliver merges its pocket
         into the outside, so both connectivity and Euler characteristic converge. */
      var minAngleOf = function (t) {
        var p = t.map(function (id) { return nodes[id]; });
        var worst = 180;
        for (var e = 0; e < 3; e++) {
          var u = p[(e + 1) % 3], v = p[(e + 2) % 3], o = p[e];
          var ux = u.x - o.x, uy = u.y - o.y, vx = v.x - o.x, vy = v.y - o.y;
          var cs = (ux * vx + uy * vy) / Math.sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy) || 1);
          worst = Math.min(worst, Math.acos(Math.max(-1, Math.min(1, cs))) * 180 / Math.PI);
        }
        return worst;
      };
      for (var clean = 0; clean < 120; clean++) {
        var triByEdge = new Map();
        keepTris.forEach(function (t, i) {
          for (var ei = 0; ei < 3; ei++) {
            var ea = Math.min(t[ei], t[(ei + 1) % 3]), eb = Math.max(t[ei], t[(ei + 1) % 3]), ekey = ea + ':' + eb;
            if (!triByEdge.has(ekey)) triByEdge.set(ekey, []);
            triByEdge.get(ekey).push(i);
          }
        });
        var comp = new Int32Array(keepTris.length).fill(-1), stack = [], compCount = 0;
        keepTris.forEach(function (t, i) {
          if (comp[i] >= 0) return;
          stack.push(i); comp[i] = compCount;
          while (stack.length) {
            var ti = stack.pop();
            for (var e2 = 0; e2 < 3; e2++) {
              var a2 = Math.min(keepTris[ti][e2], keepTris[ti][(e2 + 1) % 3]), b2 = Math.max(keepTris[ti][e2], keepTris[ti][(e2 + 1) % 3]);
              (triByEdge.get(a2 + ':' + b2) || []).forEach(function (t2) { if (comp[t2] < 0) { comp[t2] = compCount; stack.push(t2); } });
            }
          }
          compCount++;
        });
        if (compCount > 1) {
          var sizes = new Map();
          for (var ci = 0; ci < keepTris.length; ci++) sizes.set(comp[ci], (sizes.get(comp[ci]) || 0) + 1);
          var bestComp = 0, bestSize = -1;
          sizes.forEach(function (sz, c) { if (sz > bestSize) { bestSize = sz; bestComp = c; } });
          keepTris = keepTris.filter(function (t, i) { return comp[i] === bestComp; });
          continue;
        }
        var drop = -1, dropAngle = 4.5;
        keepTris.forEach(function (t, i) {
          if (drop !== -1) return;
          var bound = 0;
          for (var e3 = 0; e3 < 3; e3++) {
            var a3 = Math.min(t[e3], t[(e3 + 1) % 3]), b3 = Math.max(t[e3], t[(e3 + 1) % 3]);
            if ((triByEdge.get(a3 + ':' + b3) || []).length === 1) bound++;
          }
          if (!bound) return;
          var ang = minAngleOf(t);
          if (ang < dropAngle) { drop = i; dropAngle = ang; }
        });
        if (drop !== -1) { keepTris.splice(drop, 1); continue; }
        /* Phantom-void check: walk directed boundary edges head to tail. A void that touches
           the outside at one vertex still counts as its own loop, matching the Euler count. */
        var boundDir = [];
        keepTris.forEach(function (t, i) {
          for (var e4 = 0; e4 < 3; e4++) {
            var u4 = t[e4], v4 = t[(e4 + 1) % 3];
            var k4 = Math.min(u4, v4) + ':' + Math.max(u4, v4);
            if ((triByEdge.get(k4) || []).length === 1) boundDir.push({ u: u4, v: v4, t: i });
          }
        });
        var outMap = new Map();
        boundDir.forEach(function (be, idx) {
          if (!outMap.has(be.u)) outMap.set(be.u, []);
          outMap.get(be.u).push(idx);
        });
        var visited = new Uint8Array(boundDir.length), loopId = new Int32Array(boundDir.length).fill(-1), loopSizes = [], loops = 0;
        boundDir.forEach(function (_, idx) {
          if (visited[idx]) return;
          var count = 0, cur = idx, guard = 0;
          while (!visited[cur] && guard++ <= boundDir.length) {
            visited[cur] = 1; loopId[cur] = loops; count++;
            /* Half-edge boundary rule: next edge is CCW-first after the reversed incoming
               direction, keeping pinched loops (vertex-touching voids) separate. */
            var be = boundDir[cur], pv = nodes[be.v], pu = nodes[be.u];
            var phi = Math.atan2(pu.y - pv.y, pu.x - pv.x);
            var outs = outMap.get(be.v) || [], nxtIdx = -1, bestTurn = Infinity;
            for (var q4 = 0; q4 < outs.length; q4++) {
              var oi = outs[q4];
              if (visited[oi]) continue;
              var pw = nodes[boundDir[oi].v];
              var turn = Math.atan2(pw.y - pv.y, pw.x - pv.x) - phi;
              while (turn <= 0) turn += 2 * Math.PI;
              while (turn > 2 * Math.PI) turn -= 2 * Math.PI;
              if (turn < bestTurn) { bestTurn = turn; nxtIdx = oi; }
            }
            if (nxtIdx === -1) break;
            cur = nxtIdx;
          }
          loopSizes.push(count); loops++;
        });
        if (loops <= (outerPoly ? 1 : 0) + holePolys.length) break;
        var smallestLoop = 0;
        loopSizes.forEach(function (sz, l) { if (sz < loopSizes[smallestLoop]) smallestLoop = l; });
        var drop2 = -1, dropArea = Infinity;
        boundDir.forEach(function (be, idx) {
          if (loopId[idx] !== smallestLoop) return;
          var p = keepTris[be.t].map(function (id) { return nodes[id]; });
          var ar = Math.abs((p[1].x - p[0].x) * (p[2].y - p[0].y) - (p[2].x - p[0].x) * (p[1].y - p[0].y)) / 2;
          if (ar < dropArea) { dropArea = ar; drop2 = be.t; }
        });
        if (drop2 === -1) break;
        keepTris.splice(drop2, 1);
      }
    }
    keepTris = keepTris.filter(function (t) { return Math.abs(orient(t[0], t[1], t[2])) > 1e-18; });
    var usedNodes = new Set();
    keepTris.forEach(function (t) { t.forEach(function (id) { usedNodes.add(id); }); });
    if (usedNodes.size < nodes.length) {
      var remap = new Int32Array(nodes.length).fill(-1), keptNodes = [];
      nodes.forEach(function (node) {
        if (usedNodes.has(node.id)) { remap[node.id] = keptNodes.length; keptNodes.push({ id: keptNodes.length, x: node.x, y: node.y }); }
      });
      keepTris = keepTris.map(function (t) { return [remap[t[0]], remap[t[1]], remap[t[2]]]; });
      nodes = keptNodes;
    }
    var elements = keepTris.map(function (t, id) { return { id: id, nodes: t }; });
    if (!elements.length) throw new Error('MeshLab: domain leaves no triangles inside');
    var edgeCounts = new Map(), minimumAngle = 180, totalArea = 0, minArea = Infinity, maxRadiusEdge = 0;
    elements.forEach(function (el) {
      var p = el.nodes.map(function (id) { return nodes[id]; });
      var cross = (p[1].x - p[0].x) * (p[2].y - p[0].y) - (p[2].x - p[0].x) * (p[1].y - p[0].y);
      var area = Math.abs(cross) / 2;
      totalArea += area;
      if (area < minArea) minArea = area;
      var l01 = Math.hypot(p[1].x - p[0].x, p[1].y - p[0].y), l12 = Math.hypot(p[2].x - p[1].x, p[2].y - p[1].y), l20 = Math.hypot(p[0].x - p[2].x, p[0].y - p[2].y);
      var shortest = Math.min(l01, l12, l20);
      if (area > 0 && shortest > 0) maxRadiusEdge = Math.max(maxRadiusEdge, (l01 * l12 * l20) / (4 * area * shortest));
      for (var e = 0; e < 3; e++) {
        var a = el.nodes[e], b = el.nodes[(e + 1) % 3], key = Math.min(a, b) + ':' + Math.max(a, b);
        edgeCounts.set(key, (edgeCounts.get(key) || 0) + 1);
        var u = p[(e + 1) % 3], v = p[(e + 2) % 3], origin = p[e];
        var ux = u.x - origin.x, uy = u.y - origin.y, vx = v.x - origin.x, vy = v.y - origin.y;
        var cosine = (ux * vx + uy * vy) / Math.sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy));
        minimumAngle = Math.min(minimumAngle, Math.acos(Math.max(-1, Math.min(1, cosine))) * 180 / Math.PI);
      }
    });
    var boundaryEdges = 0, boundarySet = new Set();
    edgeCounts.forEach(function (count, key) {
      if (count === 1) { boundaryEdges++; key.split(':').forEach(function (id) { boundarySet.add(Number(id)); }); }
    });
    return { nodes: nodes, elements: elements, boundaryNodes: Array.from(boundarySet).sort(function (a, b) { return a - b; }),
      width: width, height: height, edges: edgeCounts.size, boundaryEdges: boundaryEdges,
      chi: nodes.length - edgeCounts.size + elements.length,
      quality: { minimumAngle: minimumAngle, area: totalArea, minArea: minArea, maxRadiusEdge: maxRadiusEdge }, duplicatesRemoved: duplicates };
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
    mesh.bcSets = { hot: mesh.leftNodes.slice(), cold: mesh.rightNodes.slice() };
    mesh.domainSpec = { type: 'rect', nx: nx, ny: ny, width: width, height: height };
    return mesh;
  }

  /* Boundary-edge sampling at spacing <= ds so domain edges become mesh edges. */
  function samplePolygonEdge(pts, ax, ay, bx, by, ds) {
    var steps = Math.max(1, Math.ceil(Math.hypot(bx - ax, by - ay) / ds));
    for (var s = 0; s <= steps; s++) pts.push({ x: ax + (bx - ax) * s / steps, y: ay + (by - ay) * s / steps });
  }
  function samplePolygon(pts, poly, ds) {
    for (var i = 0; i < poly.length; i++) {
      var a = poly[i], b = poly[(i + 1) % poly.length];
      samplePolygonEdge(pts, a[0], a[1], b[0], b[1], ds);
    }
  }
  /* Minimum distance from (x, y) to the segments of every polygon in `polys`. */
  function distToEdges(x, y, polys) {
    var best = Infinity;
    polys.forEach(function (poly) {
      for (var i = 0; i < poly.length; i++) {
        var a = poly[i], b = poly[(i + 1) % poly.length];
        var vx = b[0] - a[0], vy = b[1] - a[1], wx = x - a[0], wy = y - a[1];
        var t = Math.max(0, Math.min(1, (wx * vx + wy * vy) / (vx * vx + vy * vy || 1)));
        var dx = wx - t * vx, dy = wy - t * vy, d = dx * dx + dy * dy;
        if (d < best) best = d;
      }
    });
    return Math.sqrt(best);
  }

  function gridInterior(pts, xmin, ymin, xmax, ymax, nx, ny, inside, clear) {
    for (var j = 1; j < ny; j++) for (var i = 1; i < nx; i++) {
      var bx = xmin + (xmax - xmin) * i / nx, by = ymin + (ymax - ymin) * j / ny;
      if (!inside(bx, by) || (clear && !clear(bx, by))) continue;
      var jx = 0.16 * Math.sin(i * 17 + j * 31) * (xmax - xmin) / nx;
      var jy = 0.16 * Math.cos(i * 23 + j * 13) * (ymax - ymin) / ny;
      var x = bx + jx, y = by + jy;
      if (inside(x, y) && (!clear || clear(x, y))) pts.push({ x: x, y: y });
      else if (!clear || clear(bx, by)) pts.push({ x: bx, y: by });
    }
  }

  function generateDomainMesh(spec) {
    spec = spec || {};
    var type = spec.type || 'rect';
    if (type === 'rect') {
      return generateFieldMesh(spec.nx === undefined ? 14 : spec.nx, spec.ny === undefined ? 10 : spec.ny,
        spec.width === undefined ? 2.4 : spec.width, spec.height === undefined ? 1.6 : spec.height);
    }
    var nx = spec.nx === undefined ? 20 : spec.nx, ny = spec.ny === undefined ? 16 : spec.ny;
    var width = spec.width === undefined ? 2.4 : spec.width, height = spec.height === undefined ? 1.6 : spec.height;
    if (!Number.isInteger(nx) || !Number.isInteger(ny) || nx < 2 || ny < 2 || nx > 60 || ny > 60 ||
        !(width > 0 && height > 0)) throw new Error('MeshLab: invalid domain mesh dimensions');
    var ds = Math.min(width / nx, height / ny), pts = [], polygon = null, holes = [], inside = null, bcSets = null, meta = {};

    if (type === 'l') {
      var x1 = width * 0.5, y1 = height * 0.5;
      polygon = [[0, 0], [width, 0], [width, y1], [x1, y1], [x1, height], [0, height]];
      inside = function (x, y) { return pointInPolygon(x, y, polygon); };
      var clearL = function (x, y) { return distToEdges(x, y, [polygon]) > 0.35 * ds; };
      samplePolygon(pts, polygon, ds);
      gridInterior(pts, 0, 0, width, height, nx, ny, inside, clearL);
      bcSets = 'lr';
      meta = { notchX: x1, notchY: y1 };
    } else if (type === 'annulus') {
      var cx = width / 2, cy = height / 2, rOut = Math.min(width, height) / 2;
      var rIn = rOut * (spec.ratio === undefined ? 0.35 : finite(spec.ratio, 'annulus ratio'));
      if (!(rIn > 0 && rIn < rOut * 0.95)) throw new Error('MeshLab: invalid annulus ratio');
      var nOut = Math.max(48, Math.ceil(Math.PI * 2 * rOut / ds)), nIn = Math.max(28, Math.ceil(Math.PI * 2 * rIn / ds));
      polygon = [];
      for (var oi = 0; oi < nOut; oi++) {
        var th = 2 * Math.PI * oi / nOut;
        polygon.push([cx + rOut * Math.cos(th), cy + rOut * Math.sin(th)]);
        pts.push({ x: polygon[oi][0], y: polygon[oi][1] });
      }
      var hole = [];
      for (var hi = 0; hi < nIn; hi++) {
        var th2 = 2 * Math.PI * hi / nIn;
        hole.push([cx + rIn * Math.cos(th2), cy + rIn * Math.sin(th2)]);
        pts.push({ x: hole[hi][0], y: hole[hi][1] });
      }
      holes.push(hole);
      inside = function (x, y) { var r = Math.hypot(x - cx, y - cy); return r < rOut - 1e-9 && r > rIn + 1e-9; };
      var clearA = function (x, y) { return distToEdges(x, y, [polygon, hole]) > 0.3 * ds; };
      var rings = Math.max(2, Math.round((rOut - rIn) / ds));
      for (var rk = 1; rk < rings; rk++) {
        var rr = rIn + (rOut - rIn) * rk / rings;
        var count = Math.max(6, Math.round(2 * Math.PI * rr / ds));
        for (var ri = 0; ri < count; ri++) {
          var base = 2 * Math.PI * (ri + 0.5 * (rk % 2)) / count;
          var jr = rr + 0.10 * ds * Math.sin(ri * 19 + rk * 29), jt = base + 0.10 * ds / rr * Math.cos(ri * 13 + rk * 7);
          var ax = cx + jr * Math.cos(jt), ay = cy + jr * Math.sin(jt);
          if (inside(ax, ay) && clearA(ax, ay)) pts.push({ x: ax, y: ay });
        }
      }
      bcSets = 'annulus';
      meta = { center: { x: cx, y: cy }, rIn: rIn, rOut: rOut };
    } else if (type === 'sink') {
      var baseH = height * 0.45, fins = spec.fins === undefined ? 5 : spec.fins;
      if (!Number.isInteger(fins) || fins < 1 || fins > 12) throw new Error('MeshLab: invalid fin count');
      var period = width / fins;
      polygon = [[0, 0], [width, 0], [width, baseH]];
      for (var fi = fins - 1; fi >= 0; fi--) {
        var fr = fi * period + period * 0.75, fl = fi * period + period * 0.25;
        polygon.push([fr, baseH], [fr, height], [fl, height], [fl, baseH]);
      }
      polygon.push([0, baseH]);
      inside = function (x, y) { return pointInPolygon(x, y, polygon); };
      var clearS = function (x, y) { return distToEdges(x, y, [polygon]) > 0.35 * ds; };
      samplePolygon(pts, polygon, ds);
      gridInterior(pts, 0, 0, width, height, nx, ny, inside, clearS);
      bcSets = 'sink';
      meta = { baseH: baseH, fins: fins };
    } else {
      throw new Error('MeshLab: unknown domain type');
    }

    var mesh = triangulate(pts, { polygon: polygon, holes: holes });
    mesh.domainSpec = { type: type, nx: nx, ny: ny, width: width, height: height, ratio: spec.ratio, fins: spec.fins };
    mesh.meta = meta;
    var boundary = new Set(mesh.boundaryNodes);
    var hot = [], cold = [];
    mesh.nodes.forEach(function (node) {
      if (!boundary.has(node.id)) return;
      if (bcSets === 'lr') { if (node.x < ds * 0.5) hot.push(node.id); else if (node.x > width - ds * 0.5) cold.push(node.id); }
      else if (bcSets === 'annulus') {
        var r = Math.hypot(node.x - meta.center.x, node.y - meta.center.y);
        if (Math.abs(r - meta.rIn) < Math.abs(r - meta.rOut)) hot.push(node.id); else cold.push(node.id);
      } else if (bcSets === 'sink') {
        if (node.y < ds * 0.5) hot.push(node.id); else cold.push(node.id);
      }
    });
    if (!hot.length || !cold.length) throw new Error('MeshLab: domain boundary sets are incomplete');
    mesh.bcSets = { hot: hot, cold: cold };
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
      var source;
      if (typeof cfg.source === 'function') {
        /* Mid-edge rule: integrates quadratics exactly, tighter than a centroid sample for the Gaussian heater. */
        var m01x = (p[0].x + p[1].x) / 2, m01y = (p[0].y + p[1].y) / 2;
        var m12x = (p[1].x + p[2].x) / 2, m12y = (p[1].y + p[2].y) / 2;
        var m20x = (p[2].x + p[0].x) / 2, m20y = (p[2].y + p[0].y) / 2;
        source = (finite(cfg.source(m01x, m01y), 'heat source') + finite(cfg.source(m12x, m12y), 'heat source') + finite(cfg.source(m20x, m20y), 'heat source')) / 3;
      } else {
        source = cfg.source === undefined ? 0 : cfg.source;
        finite(source, 'heat source');
      }
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
    /* Area-weighted nodal flux for exact per-node (x, y, u, |grad u|) inspection. */
    var nodeFlux = mesh.nodes.map(function () { return { x: 0, y: 0, w: 0 }; });
    elements.forEach(function (e) {
      e.nodes.forEach(function (id) {
        nodeFlux[id].x += e.fluxX * e.area; nodeFlux[id].y += e.fluxY * e.area; nodeFlux[id].w += e.area;
      });
    });
    nodeFlux = nodeFlux.map(function (nf) { return nf.w > 0 ? { x: nf.x / nf.w, y: nf.y / nf.w } : { x: 0, y: 0 }; });
    return { mesh: mesh, values: Array.from(values), elements: elements, nodeFlux: nodeFlux,
      min: Math.min.apply(null, Array.from(values)), max: Math.max.apply(null, Array.from(values)),
      energy: elements.reduce(function (s, e) { return s + e.energy; }, 0), diagnostics: result,
      verification: { freeResidual: freeResidual, balanceResidual: sourceTotal + reaction, sourceTotal: sourceTotal, boundaryFlux: reaction },
      freeDofs: free, conductivity: conductivity };
  }

  /* VTK legacy ASCII UNSTRUCTURED_GRID. Elements: 2 nodes -> VTK_LINE (3), 3 -> VTK_TRIANGLE (5), 4 -> VTK_QUAD (9). */
  function exportVtk(mesh, data) {
    if (!mesh || !Array.isArray(mesh.nodes) || !Array.isArray(mesh.elements) || !mesh.nodes.length || !mesh.elements.length) {
      throw new Error('MeshLab: VTK export needs a mesh');
    }
    data = data || {};
    var n = mesh.nodes.length, m = mesh.elements.length;
    var cellTypes = { 2: 3, 3: 5, 4: 9 };
    var cellInts = 0;
    mesh.elements.forEach(function (el) {
      if (!el.nodes || !cellTypes[el.nodes.length]) throw new Error('MeshLab: VTK export supports line/triangle/quad cells');
      el.nodes.forEach(function (id) { if (!Number.isInteger(id) || id < 0 || id >= n) throw new Error('MeshLab: VTK export found an invalid node reference'); });
      cellInts += el.nodes.length + 1;
    });
    var title = String(data.title === undefined ? 'MeshLab unstructured grid' : data.title).replace(/[\r\n]+/g, ' ').slice(0, 255);
    var num = function (v) {
      var x = finite(Number(v), 'VTK value');
      if (x === 0) x = 0;
      return String(x);
    };
    var lines = ['# vtk DataFile Version 3.0', title, 'ASCII', 'DATASET UNSTRUCTURED_GRID', 'POINTS ' + n + ' float'];
    mesh.nodes.forEach(function (node) { lines.push(num(node.x) + ' ' + num(node.y) + ' ' + num(node.z || 0)); });
    lines.push('CELLS ' + m + ' ' + cellInts);
    mesh.elements.forEach(function (el) { lines.push(el.nodes.length + ' ' + el.nodes.join(' ')); });
    lines.push('CELL_TYPES ' + m);
    mesh.elements.forEach(function (el) { lines.push(String(cellTypes[el.nodes.length])); });
    function block(setKey, count, scalars, vectors) {
      scalars = scalars || []; vectors = vectors || [];
      if (!scalars.length && !vectors.length) return;
      lines.push(setKey + ' ' + count);
      scalars.concat(vectors).forEach(function (field) {
        if (!field || typeof field.name !== 'string' || !field.name.length) throw new Error('MeshLab: VTK field needs a name');
        if (!/^[-\w.]+$/.test(field.name)) throw new Error('MeshLab: VTK field name must be a single word');
        var isVector = vectors.indexOf(field) !== -1;
        if (!Array.isArray(field.values) || field.values.length !== count) throw new Error('MeshLab: VTK field ' + field.name + ' length mismatch');
        if (isVector) {
          lines.push('VECTORS ' + field.name + ' float');
          field.values.forEach(function (v) {
            if (!v) throw new Error('MeshLab: VTK vector entry missing');
            var vx = Array.isArray(v) ? v[0] : v.x, vy = Array.isArray(v) ? v[1] : v.y, vz = Array.isArray(v) ? v[2] : (v.z === undefined ? 0 : v.z);
            lines.push(num(vx) + ' ' + num(vy) + ' ' + num(vz));
          });
        } else {
          lines.push('SCALARS ' + field.name + ' float 1');
          lines.push('LOOKUP_TABLE default');
          field.values.forEach(function (v) { lines.push(num(v)); });
        }
      });
    }
    block('POINT_DATA', n, data.pointScalars, data.pointVectors);
    block('CELL_DATA', m, data.cellScalars, data.cellVectors);
    return lines.join('\n') + '\n';
  }

  return { toCSR: toCSR, multiply: multiply, solvePCG: solvePCG, triangulate: triangulate,
    generateFieldMesh: generateFieldMesh, generateDomainMesh: generateDomainMesh, solvePoisson: solvePoisson,
    exportVtk: exportVtk, pointInPolygon: pointInPolygon };
});
