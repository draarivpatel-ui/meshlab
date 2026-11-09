/*
 * MeshLab - finite element structural analysis for trusses, frames and plates
 * MIT License. Copyright (c) 2026 draarivpatel-ui
 *
 * Units are SI throughout: metres, newtons, pascals, kilograms.
 * Works in Node via require('./fem.js') and in the browser via window.MeshLab.
 */
(function (root, factory) {
  'use strict';
  if (typeof module === 'object' && module.exports) {
    module.exports = factory();
  } else {
    root.MeshLab = factory();
  }
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';

  // ---------------------------------------------------------------------------
  // Small dense linear algebra toolbox
  // ---------------------------------------------------------------------------

  function zeros(n) {
    var a = new Float64Array(n);
    for (var i = 0; i < n; i++) a[i] = 0;
    return a;
  }

  function matrix(r, c) {
    var m = new Array(r), i;
    for (i = 0; i < r; i++) m[i] = zeros(c);
    return m;
  }

  function identity(n) {
    var m = matrix(n, n), i;
    for (i = 0; i < n; i++) m[i][i] = 1;
    return m;
  }

  function matVec(A, x) {
    var n = A.length, out = zeros(n), i, j, s;
    for (i = 0; i < n; i++) {
      s = 0;
      for (j = 0; j < n; j++) s += A[i][j] * x[j];
      out[i] = s;
    }
    return out;
  }

  function matMul(A, B) {
    var n = A.length, k = B.length, c = B[0].length, C = matrix(n, c), i, j, p, s;
    for (i = 0; i < n; i++) {
      for (j = 0; j < c; j++) {
        s = 0;
        for (p = 0; p < k; p++) s += A[i][p] * B[p][j];
        C[i][j] = s;
      }
    }
    return C;
  }

  function transpose(A) {
    var r = A.length, c = A[0].length, T = matrix(c, r), i, j;
    for (i = 0; i < r; i++) for (j = 0; j < c; j++) T[j][i] = A[i][j];
    return T;
  }

  /* Gaussian elimination with partial pivoting. Returns null on a singular matrix. */
  function solveLinear(Ain, bin) {
    var n = bin.length, A = new Array(n), b = Float64Array.from(bin), i, j, k;
    for (i = 0; i < n; i++) A[i] = Float64Array.from(Ain[i]);
    var scale = zeros(n);
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) scale[i] = Math.max(scale[i], Math.abs(A[i][j]));
      if (scale[i] === 0) return null;
    }
    for (k = 0; k < n; k++) {
      var piv = k, best = Math.abs(A[k][k]) / scale[k];
      for (i = k + 1; i < n; i++) {
        var q = Math.abs(A[i][k]) / scale[i];
        if (q > best) { best = q; piv = i; }
      }
      if (!(best > 1e-13) || !isFinite(best)) return null;
      if (piv !== k) {
        var tmp = A[k]; A[k] = A[piv]; A[piv] = tmp;
        var tb = b[k]; b[k] = b[piv]; b[piv] = tb;
        var ts = scale[k]; scale[k] = scale[piv]; scale[piv] = ts;
      }
      for (i = k + 1; i < n; i++) {
        var f = A[i][k] / A[k][k];
        if (f === 0) continue;
        A[i][k] = 0;
        for (j = k + 1; j < n; j++) A[i][j] -= f * A[k][j];
        b[i] -= f * b[k];
      }
    }
    var x = zeros(n);
    for (i = n - 1; i >= 0; i--) {
      var s = b[i];
      for (j = i + 1; j < n; j++) s -= A[i][j] * x[j];
      var d = A[i][i];
      if (Math.abs(d) < 1e-300) return null;
      x[i] = s / d;
    }
    for (i = 0; i < n; i++) if (!isFinite(x[i])) return null;
    return x;
  }

  /* Jacobi eigenvalue solver for a real symmetric matrix. */
  function jacobiEigen(Ain, maxSweeps) {
    var n = Ain.length, A = new Array(n), V = identity(n), i, j, p, q, sweep;
    for (i = 0; i < n; i++) A[i] = Float64Array.from(Ain[i]);
    var limit = maxSweeps || 60 * n * n;
    for (sweep = 0; sweep < limit; sweep++) {
      var off = 0, biggest = 0;
      for (i = 0; i < n; i++) for (j = i + 1; j < n; j++) {
        off += A[i][j] * A[i][j];
        biggest = Math.max(biggest, Math.abs(A[i][j]));
      }
      var norm = 0;
      for (i = 0; i < n; i++) norm += A[i][i] * A[i][i];
      if (biggest <= 1e-12 * Math.sqrt(norm + 1) || off === 0) break;
      for (p = 0; p < n - 1; p++) {
        for (q = p + 1; q < n; q++) {
          var apq = A[p][q];
          var threshold = 1e-14 * Math.sqrt(norm + 1);
          if (Math.abs(apq) <= threshold) continue;
          var theta = (A[q][q] - A[p][p]) / (2 * apq);
          var t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
          if (!isFinite(t) || t === 0) continue;
          var c = 1 / Math.sqrt(t * t + 1), s = t * c, tau = s / (1 + c);
          A[p][p] -= t * apq;
          A[q][q] += t * apq;
          A[p][q] = 0; A[q][p] = 0;
          for (j = 0; j < n; j++) {
            if (j === p || j === q) continue;
            var ajp = A[j][p], ajq = A[j][q];
            A[j][p] = A[p][j] = ajp - s * (ajq + tau * ajp);
            A[j][q] = A[q][j] = ajq + s * (ajp - tau * ajq);
          }
          for (j = 0; j < n; j++) {
            var vjp = V[j][p], vjq = V[j][q];
            V[j][p] = c * vjp - s * vjq;
            V[j][q] = s * vjp + c * vjq;
          }
        }
      }
    }
    var pairs = [];
    for (i = 0; i < n; i++) {
      var vec = zeros(n);
      for (j = 0; j < n; j++) vec[j] = V[j][i];
      pairs.push({ value: A[i][i], vector: vec });
    }
    pairs.sort(function (a, b) { return a.value - b.value; });
    return pairs;
  }

  function dot(a, b) {
    var s = 0, i;
    for (i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }

  function clamp(v, lo, hi) { return Math.min(hi, Math.max(lo, v)); }
  function num(v, fallback) { return (typeof v === 'number' && isFinite(v)) ? v : fallback; }

  // ---------------------------------------------------------------------------
  // Model normalisation
  // ---------------------------------------------------------------------------

  var BC_DOF = {
    pin: function (d) { return [true, true].concat(d === 3 ? [true] : []); },
    fixed: function (d) { return [true, true].concat(d === 3 ? [true] : []); },
    rollerX: function (d) { return [true, false].concat(d === 3 ? [false] : []); },
    rollerY: function (d) { return [false, true].concat(d === 3 ? [false] : []); }
  };

  function normalizeModel(input) {
    if (!input || !Array.isArray(input.nodes) || !Array.isArray(input.members)) {
      throw new Error('MeshLab: model needs nodes[] and members[] arrays');
    }
    var dim = (input.dim === 3 || (input.dim === 'auto' && input.nodes.some(function (n) {
      return n && !Array.isArray(n) && isFinite(n.z) && Math.abs(n.z) > 1e-12;
    }))) ? 3 : 2;
    var nodes = input.nodes.map(function (n, i) {
      var p;
      if (Array.isArray(n)) {
        p = n;
      } else if (n.coords || n.xyz || n.p) {
        p = n.coords || n.xyz || n.p;
      } else {
        p = [n.x, n.y, n.z];
      }
      var id = num(Array.isArray(n) ? i : n.id, i);
      var load = (!Array.isArray(n) && n.load) ? n.load.slice(0, dim) : [0, 0, 0].slice(0, dim);
      return {
        id: id,
        x: num(p[0], 0), y: num(p[1], 0), z: dim === 3 ? num(p[2], 0) : 0,
        support: (!Array.isArray(n) && n.support) || 'free',
        load: load
      };
    });
    var defaultMat = input.material || {};
    var E = num(defaultMat.E, 2.0e11);
    var rho = num(defaultMat.rho, 7850);
    var yieldStress = num(defaultMat.sigmaYield, num(defaultMat.sigma_yield, 2.5e8));
    var members = input.members.map(function (m, i) {
      var a = num(m.a !== undefined ? m.a : m.n1, 0);
      var b = num(m.b !== undefined ? m.b : m.n2, 0);
      var mat = m.material || {};
      return {
        id: i, a: a, b: b,
        label: m.label || ('M' + (i + 1)),
        E: num(m.E, num(mat.E, E)),
        A: num(m.A, num(mat.A, 0.01)),
        I: num(m.I, num(mat.I, 1e-6)),
        rho: num(m.rho, num(mat.rho, rho)),
        sigmaYield: num(m.sigmaYield, num(mat.sigmaYield, yieldStress)),
        Keff: clamp(num(m.Keff, num(mat.Keff, 1.0)), 0.25, 4),
        rGyr: num(m.rGyr, num(mat.rGyr, null))
      };
    });
    return { dim: dim, nodes: nodes, members: members, material: {
      E: E, rho: rho, sigmaYield: yieldStress
    }, name: input.name || 'Untitled structure' };
  }

  function nodeCoords(n, dim) { return dim === 3 ? [n.x, n.y, n.z] : [n.x, n.y]; }

  function memberGeometry(model, m) {
    var A = model.nodes[m.a], B = model.nodes[m.b];
    if (!A || !B) throw new Error('MeshLab: member ' + m.label + ' references a missing node');
    var p = nodeCoords(A, model.dim), q = nodeCoords(B, model.dim);
    var d = [], L = 0, i;
    for (i = 0; i < model.dim; i++) { d.push(q[i] - p[i]); L += d[i] * d[i]; }
    L = Math.sqrt(L);
    if (L <= 1e-12) throw new Error('MeshLab: member ' + m.label + ' has zero length');
    var c = d.map(function (v) { return v / L; });
    return { L: L, cx: c[0], cy: c[1], cz: model.dim === 3 ? c[2] : 0, c: c, p: p, q: q };
  }

  /* Supports may be declared on the node or in a top level supports map. */
  function applySupports(model, raw) {
    var list = raw.supports;
    if (Array.isArray(list)) {
      list.forEach(function (s) {
        var n = model.nodes[num(s.node !== undefined ? s.node : s.n, -1)];
        if (n) n.support = s.type || s.support || 'pin';
      });
    } else if (list && typeof list === 'object') {
      Object.keys(list).forEach(function (key) {
        var n = model.nodes[Number(key)];
        if (n) n.support = list[key];
      });
    }
    var loads = raw.loads;
    if (Array.isArray(loads)) {
      loads.forEach(function (L) {
        var n = model.nodes[num(L.node !== undefined ? L.node : L.n, -1)];
        if (!n) return;
        var f = L.f || L.load || L.force || [];
        for (var i = 0; i < model.dim; i++) n.load[i] += num(f[i], 0);
      });
    }
  }

  function buildSystem(model, opts) {
    var d = model.dim, N = model.nodes.length, ndof = d * N;
    var K = matrix(ndof, ndof), F = zeros(ndof), M = matrix(ndof, ndof);
    var fixedMask = new Array(ndof), i, j, k;
    for (i = 0; i < ndof; i++) fixedMask[i] = false;
    for (i = 0; i < N; i++) {
      var n = model.nodes[i];
      for (j = 0; j < d; j++) F[d * i + j] = num(n.load[j], 0);
      var rule = BC_DOF[n.support];
      if (rule) {
        var mask = rule(d);
        for (j = 0; j < d; j++) fixedMask[d * i + j] = !!mask[j];
      }
    }
    var geo = [];
    for (k = 0; k < model.members.length; k++) {
      var m = model.members[k], g = memberGeometry(model, m);
      geo.push(g);
      var stiffness = m.E * m.A / g.L;
      var idx = [], a;
      for (a = 0; a < d; a++) idx.push(d * m.a + a);
      for (a = 0; a < d; a++) idx.push(d * m.b + a);
      var n2 = 2 * d, ke = matrix(n2, n2);
      for (i = 0; i < d; i++) for (j = 0; j < d; j++) {
        var v = stiffness * g.c[i] * g.c[j];
        ke[i][j] += v; ke[i][d + j] -= v;
        ke[d + i][j] -= v; ke[d + i][d + j] += v;
      }
      for (i = 0; i < n2; i++) for (j = 0; j < n2; j++) K[idx[i]][idx[j]] += ke[i][j];
      if (!opts || opts.mass !== false) {
        var total = m.rho * m.A * g.L;
        var mode = (opts && opts.massType) || 'lumped';
        if (mode === 'consistent') {
          var me = matrix(n2, n2), w = total / 6;
          for (i = 0; i < d; i++) for (j = 0; j < d; j++) {
            var delta = i === j ? 1 : 0;
            var local = w * (2 * delta);
            me[i][j] += local; me[i][d + j] += delta * w;
            me[d + i][j] += delta * w; me[d + i][d + j] += local;
          }
          for (i = 0; i < n2; i++) for (j = 0; j < n2; j++) M[idx[i]][idx[j]] += me[i][j];
        } else {
          for (i = 0; i < d; i++) { M[idx[i]][idx[i]] += total / 2; M[idx[d + i]][idx[d + i]] += total / 2; }
        }
      }
    }
    return { K: K, F: F, M: M, fixedMask: fixedMask, geo: geo, ndof: ndof };
  }

  function reducedSystem(sys) {
    var free = [], i, j;
    for (i = 0; i < sys.ndof; i++) if (!sys.fixedMask[i]) free.push(i);
    var n = free.length, Kr = matrix(n, n), Fr = zeros(n), Mr = matrix(n, n);
    for (i = 0; i < n; i++) {
      Fr[i] = sys.F[free[i]];
      for (j = 0; j < n; j++) {
        Kr[i][j] = sys.K[free[i]][free[j]];
        Mr[i][j] = sys.M[free[i]][free[j]];
      }
    }
    return { free: free, K: Kr, F: Fr, M: Mr, n: n };
  }

  // ---------------------------------------------------------------------------
  // Direct stiffness truss / frame solver
  // ---------------------------------------------------------------------------

  function solveTruss(modelInput) {
    var model = normalizeModel(modelInput);
    applySupports(model, modelInput);
    var sys = buildSystem(model, { mass: true });
    var red = reducedSystem(sys);
    var ur = solveLinear(red.K, red.F);
    if (!ur) {
      throw new Error('MeshLab: stiffness matrix is singular (mechanism or insufficient supports)');
    }
    var u = zeros(sys.ndof), i, j;
    for (i = 0; i < red.n; i++) u[red.free[i]] = ur[i];
    var Ku = matVec(sys.K, u);
    var R = zeros(sys.ndof);
    for (i = 0; i < sys.ndof; i++) R[i] = Ku[i] - sys.F[i];
    var nodes = model.nodes.map(function (n, k) {
      var disp = [], react = [], f = [];
      for (j = 0; j < model.dim; j++) {
        disp.push(u[model.dim * k + j]);
        react.push(R[model.dim * k + j]);
        f.push(sys.F[model.dim * k + j]);
      }
      return { id: n.id, x: n.x, y: n.y, z: n.z, support: n.support,
        load: f, displacement: disp, reaction: react };
    });
    var members = model.members.map(function (m, k) {
      var g = sys.geo[k];
      var delta = 0;
      for (j = 0; j < model.dim; j++) {
        delta += g.c[j] * (u[model.dim * m.b + j] - u[model.dim * m.a + j]);
      }
      var N = m.E * m.A / g.L * delta;
      var sigma = N / m.A;
      var eps = delta / g.L;
      var U = 0.5 * N * delta;
      return {
        id: m.id, label: m.label, a: m.a, b: m.b, length: g.L,
        directionCosines: { cx: g.cx, cy: g.cy, cz: g.cz },
        E: m.E, A: m.A, I: m.I, Keff: m.Keff, sigmaYield: m.sigmaYield,
        delta: delta, axialForce: N, stress: sigma, strain: eps,
        strainEnergy: U, tension: N > 0
      };
    });
    var r = 0;
    for (i = 0; i < sys.ndof; i++) if (sys.fixedMask[i]) r++;
    var totalU = members.reduce(function (s, m) { return s + m.strainEnergy; }, 0);
    var dimFactor = model.dim === 3 ? 3 : 2;
    var determinacy = model.members.length + r - dimFactor * model.nodes.length;
    var extFx = 0, extFy = 0, extFz = 0, rx = 0, ry = 0, rz = 0, work = 0;
    for (i = 0; i < model.nodes.length; i++) {
      var nd = nodes[i];
      extFx += nd.load[0]; extFy += nd.load[1];
      rx += nd.reaction[0]; ry += nd.reaction[1];
      for (j = 0; j < model.dim; j++) work += nd.load[j] * nd.displacement[j];
      if (model.dim === 3) { extFz += nd.load[2]; rz += nd.reaction[2]; }
    }
    return {
      name: model.name, dim: model.dim, model: model,
      nodes: nodes, members: members, displacements: Array.from(u),
      forceVector: Array.from(sys.F), reactions: Array.from(R),
      stiffnessMatrix: sys.K, massMatrix: sys.M, fixedMask: sys.fixedMask,
      maxDisplacement: Math.max.apply(null, Array.from(u).map(Math.abs).concat([0])),
      strainEnergy: totalU, externalWork: work,
      determinacy: {
        m: model.members.length, r: r, j: model.nodes.length,
        value: determinacy,
        classification: determinacy < 0 ? 'unstable mechanism' :
          (determinacy === 0 ? 'statically determinate' : 'statically indeterminate'),
        stable: determinacy >= 0
      },
      equilibrium: {
        sumFx: extFx, sumFy: extFy, sumFz: extFz,
        sumRx: rx, sumRy: ry, sumRz: rz,
        residualX: extFx + rx, residualY: extFy + ry, residualZ: extFz + rz
      }
    };
  }

  // ---------------------------------------------------------------------------
  // Euler buckling, yield and factor of safety analyzer
  // ---------------------------------------------------------------------------

  function analyzeSafety(modelInput, solutionInput) {
    var model = solutionInput ? solutionInput.model : normalizeModel(modelInput);
    if (!model || !model.members) model = normalizeModel(modelInput);
    var sol = solutionInput || solveTruss(modelInput);
    var results = sol.members.map(function (mem) {
      var m = model.members[mem.id] || mem;
      var E = num(m.E, mem.E), I = num(m.I, mem.I), L = mem.length;
      var A = num(m.A, mem.A), sy = num(m.sigmaYield, mem.sigmaYield);
      var keff = num(m.Keff, mem.Keff || 1);
      var rGyr = num(m.rGyr, Math.sqrt(I / A));
      var N = mem.axialForce, sigma = mem.stress;
      var compression = N < 0;
      var Pcr = Math.PI * Math.PI * E * I / Math.pow(keff * L, 2);
      var fosYield = Math.abs(sigma) > 1e-12 ? sy / Math.abs(sigma) : Infinity;
      var fosBuckling = compression && Math.abs(N) > 1e-12 ? Pcr / Math.abs(N) : Infinity;
      var governing = Math.min(fosYield, fosBuckling);
      return {
        id: mem.id, label: mem.label, length: L, axialForce: N,
        stress: sigma, compression: compression,
        eulerCriticalLoad: Pcr, slendernessRatio: rGyr > 0 ? L / rGyr : Infinity,
        radiusOfGyration: rGyr, effectiveLength: keff * L,
        fosYield: fosYield, fosBuckling: fosBuckling, fosGoverning: governing,
        yieldCritical: fosYield < 1, bucklingCritical: fosBuckling < 1,
        critical: governing < 1,
        warning: governing < 1.5,
        governingMode: fosBuckling < fosYield ? 'buckling' : 'yield'
      };
    });
    var critical = results.filter(function (r) { return r.critical; });
    var finite = results.filter(function (r) { return isFinite(r.fosGoverning); })
      .map(function (r) { return r.fosGoverning; });
    return {
      members: results, criticalMembers: critical,
      criticalCount: critical.length,
      minFos: finite.length ? Math.min.apply(null, finite) : Infinity,
      safe: critical.length === 0
    };
  }

  // ---------------------------------------------------------------------------
  // Natural frequencies and mode shapes
  // ---------------------------------------------------------------------------

  function solveModes(modelInput, numModes) {
    var model = normalizeModel(modelInput);
    applySupports(model, modelInput);
    var sys = buildSystem(model, { mass: true, massType: modelInput.massType || 'lumped' });
    var red = reducedSystem(sys);
    var n = red.n;
    var count = clamp(Math.floor(num(numModes, 4)), 1, n);
    if (n === 0) return { modes: [], frequencies: [], numModes: 0, ndof: sys.ndof, massType: 'lumped' };
    /* Transform to the standard symmetric problem A = M^(-1/2) K M^(-1/2). */
    var invSqrt = zeros(n), i, j;
    for (i = 0; i < n; i++) {
      var mi = red.M[i][i];
      if (!(mi > 1e-14)) throw new Error('MeshLab: mass matrix has a zero diagonal entry');
      invSqrt[i] = 1 / Math.sqrt(mi);
    }
    var A = matrix(n, n);
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
      A[i][j] = invSqrt[i] * red.K[i][j] * invSqrt[j];
    }
    var eig = jacobiEigen(A);
    var modes = [];
    for (var k = 0; k < count; k++) {
      var lambda = Math.max(eig[k].value, 0);
      var omega = Math.sqrt(lambda);
      var w = eig[k].vector;
      var phiRed = zeros(n);
      for (i = 0; i < n; i++) phiRed[i] = invSqrt[i] * w[i];
      /* Mass normalise: phi^T M phi = 1. */
      var mNorm = 0;
      for (i = 0; i < n; i++) {
        var sm = 0;
        for (j = 0; j < n; j++) sm += red.M[i][j] * phiRed[j];
        mNorm += phiRed[i] * sm;
      }
      if (mNorm > 1e-14) {
        var f = 1 / Math.sqrt(mNorm);
        for (i = 0; i < n; i++) phiRed[i] *= f;
      }
      var phi = zeros(sys.ndof);
      for (i = 0; i < n; i++) phi[red.free[i]] = phiRed[i];
      modes.push({
        index: k + 1, eigenvalue: lambda, omega: omega,
        frequencyHz: omega / (2 * Math.PI), period: omega > 0 ? 2 * Math.PI / omega : Infinity,
        modeShape: Array.from(phi), shape: Array.from(phi)
      });
    }
    return {
      modes: modes, frequencies: modes.map(function (m) { return m.frequencyHz; }),
      numModes: modes.length, ndof: sys.ndof,
      massType: modelInput.massType || 'lumped',
      massMatrix: sys.M, reducedMassMatrix: red.M, freeDofs: red.free
    };
  }

  // ---------------------------------------------------------------------------
  // Constant strain triangle plane stress solver
  // ---------------------------------------------------------------------------

  function generatePlateMesh(nx, ny, width, height, holeRadius) {
    nx = Math.max(2, Math.floor(num(nx, 8)));
    ny = Math.max(2, Math.floor(num(ny, 6)));
    width = num(width, 2); height = num(height, 1);
    holeRadius = num(holeRadius, 0);
    var nodes = [], id = 0, i, j;
    var xc = width / 2, yc = height / 2;
    var grid = [];
    for (j = 0; j <= ny; j++) {
      grid[j] = [];
      for (i = 0; i <= nx; i++) {
        var x = width * i / nx, y = height * j / ny;
        var insideHole = holeRadius > 0 &&
          Math.pow(x - xc, 2) + Math.pow(y - yc, 2) < holeRadius * holeRadius;
        if (insideHole) { grid[j][i] = -1; continue; }
        grid[j][i] = id;
        nodes.push({ id: id, x: x, y: y });
        id++;
      }
    }
    var elements = [];
    function addTri(a, b, c) {
      if (a < 0 || b < 0 || c < 0) return;
      var na = nodes[a], nb = nodes[b], nc = nodes[c];
      var area2 = (nb.x - na.x) * (nc.y - na.y) - (nc.x - na.x) * (nb.y - na.y);
      if (Math.abs(area2) < 1e-14) return;
      if (area2 < 0) { var t = b; b = c; c = t; }
      elements.push({ id: elements.length, nodes: [a, b, c] });
    }
    for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
      var a = grid[j][i], b = grid[j][i + 1], c = grid[j + 1][i], d = grid[j + 1][i + 1];
      addTri(a, b, d);
      addTri(a, d, c);
    }
    var left = [], right = [];
    for (j = 0; j <= ny; j++) {
      left.push(grid[j][0]); right.push(grid[j][nx]);
    }
    return {
      nodes: nodes, elements: elements, nx: nx, ny: ny,
      width: width, height: height, holeRadius: holeRadius,
      leftNodes: left.filter(function (v) { return v >= 0; }),
      rightNodes: right.filter(function (v) { return v >= 0; }),
      grid: grid
    };
  }

  function cstB(mesh, el) {
    var p = el.nodes.map(function (k) { return mesh.nodes[k]; });
    var x1 = p[0].x, y1 = p[0].y, x2 = p[1].x, y2 = p[1].y, x3 = p[2].x, y3 = p[2].y;
    var area2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    var A = area2 / 2;
    var inv = 1 / area2;
    var b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
    var c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;
    var B = [
      [b1 * inv, 0, b2 * inv, 0, b3 * inv, 0],
      [0, c1 * inv, 0, c2 * inv, 0, c3 * inv],
      [c1 * inv, b1 * inv, c2 * inv, b2 * inv, c3 * inv, b3 * inv]
    ];
    return { B: B, area: Math.abs(A), areaSigned: A, coords: p };
  }

  function planeStressD(E, nu) {
    var f = E / (1 - nu * nu);
    return [[f, f * nu, 0], [f * nu, f, 0], [0, 0, f * (1 - nu) / 2]];
  }

  function stressMeasures(sx, sy, txy) {
    var mean = (sx + sy) / 2;
    var rad = Math.sqrt(Math.pow((sx - sy) / 2, 2) + txy * txy);
    return {
      sigma1: mean + rad, sigma2: mean - rad,
      tauMax: rad,
      vonMises: Math.sqrt(Math.max(0, sx * sx - sx * sy + sy * sy + 3 * txy * txy))
    };
  }

  function solvePlateCST(meshInput, material, bc, loads) {
    var mesh = meshInput.nodes ? meshInput : generatePlateMesh(8, 6, 2, 1, 0);
    var mat = material || {};
    var E = num(mat.E, 2.0e11), nu = clamp(num(mat.nu, 0.3), -0.49, 0.49);
    var t = num(mat.t, num(mat.thickness, 0.01));
    var D = planeStressD(E, nu);
    var N = mesh.nodes.length, ndof = 2 * N;
    var K = matrix(ndof, ndof), F = zeros(ndof), fixed = new Array(ndof);
    var i, j, k, l;
    for (i = 0; i < ndof; i++) fixed[i] = false;
    var elemData = [];
    for (k = 0; k < mesh.elements.length; k++) {
      var el = mesh.elements[k], c = cstB(mesh, el);
      var BT = transpose(c.B), DB = matMul(D, c.B);
      var ke = matMul(BT, DB);
      var idx = [];
      for (i = 0; i < 3; i++) idx.push(2 * el.nodes[i], 2 * el.nodes[i] + 1);
      var factor = t * c.area;
      for (i = 0; i < 6; i++) for (j = 0; j < 6; j++) K[idx[i]][idx[j]] += ke[i][j] * factor;
      elemData.push({ B: c.B, DB: DB, area: c.area, coords: c.coords, idx: idx });
    }
    var cfg = bc || {};
    function fixSet(list, dofMask) {
      if (!list) return;
      (Array.isArray(list) ? list : [list]).forEach(function (nodeId) {
        var n = num(nodeId, -1);
        if (n < 0 || n >= N) return;
        if (dofMask[0]) fixed[2 * n] = true;
        if (dofMask[1]) fixed[2 * n + 1] = true;
      });
    }
    if (cfg.left === 'fixed' || cfg.left === 'pin') fixSet(mesh.leftNodes, [true, true]);
    else if (cfg.left === 'rollerX') fixSet(mesh.leftNodes, [true, false]);
    else if (cfg.left === 'rollerY') fixSet(mesh.leftNodes, [false, true]);
    else if (Array.isArray(cfg.left)) fixSet(cfg.left, [true, true]);
    if (cfg.right === 'rollerY') fixSet(mesh.rightNodes, [false, true]);
    else if (cfg.right === 'rollerX') fixSet(mesh.rightNodes, [true, false]);
    else if (cfg.right === 'fixed') fixSet(mesh.rightNodes, [true, true]);
    else if (Array.isArray(cfg.right)) fixSet(cfg.right, [true, true]);
    (cfg.fixedNodes || []).forEach(function (n) { fixSet(n, [true, true]); });
    (cfg.constraints || []).forEach(function (q) {
      var node = num(q.node, -1);
      if (node >= 0 && node < N) {
        if (q.x || q.ux === 0) fixed[2 * node] = true;
        if (q.y || q.uy === 0) fixed[2 * node + 1] = true;
      }
    });
    /* Default: restrain the left edge in x and the bottom left corner in y. */
    var anyFixed = fixed.some(function (v) { return v; });
    if (!anyFixed) {
      fixSet(mesh.leftNodes, [true, false]);
      var bottomLeft = mesh.grid[0][0];
      if (bottomLeft >= 0) fixed[2 * bottomLeft + 1] = true;
    }
    var loadCfg = loads || {};
    function applyNodal(list, fx, fy) {
      if (!list) return;
      (Array.isArray(list) ? list : [list]).forEach(function (nodeId) {
        var n = num(nodeId, -1);
        if (n < 0 || n >= N) return;
        F[2 * n] += fx; F[2 * n + 1] += fy;
      });
    }
    /* Consistent uniform edge traction: each boundary segment contributes
       q * segmentLength / 2 to both of its end nodes. This reproduces exact
       uniform stress for a constant strain triangle patch. */
    function applyEdgeTraction(list, total, sign) {
      if (!list || list.length < 2 || !total) return;
      var q = total / mesh.height;
      for (var e = 0; e < list.length - 1; e++) {
        var na = list[e], nb = list[e + 1];
        if (na < 0 || nb < 0 || na >= N || nb >= N) continue;
        var seg = Math.abs(mesh.nodes[nb].y - mesh.nodes[na].y);
        var half = q * seg / 2;
        F[2 * na] += sign * half;
        F[2 * nb] += sign * half;
      }
    }
    if (loadCfg.rightNodesForce) {
      applyEdgeTraction(mesh.rightNodes, Math.abs(num(loadCfg.rightNodesForce, 0)),
        Math.sign(num(loadCfg.rightNodesForce, 0)));
    }
    if (loadCfg.leftNodesForce) {
      applyEdgeTraction(mesh.leftNodes, Math.abs(num(loadCfg.leftNodesForce, 0)),
        Math.sign(num(loadCfg.leftNodesForce, 0)));
    }
    (loadCfg.nodal || []).forEach(function (q) {
      var n = num(q.node, -1);
      if (n >= 0 && n < N) { F[2 * n] += num(q.fx, 0); F[2 * n + 1] += num(q.fy, 0); }
    });
    if (loadCfg.nodalArray && Array.isArray(loadCfg.nodalArray)) {
      applyNodal(loadCfg.nodalArray.nodes || [], num(loadCfg.nodalArray.fx, 0), num(loadCfg.nodalArray.fy, 0));
    }
    var free = [], r = 0;
    for (i = 0; i < ndof; i++) { if (!fixed[i]) free.push(i); else r++; }
    var n = free.length, Kr = matrix(n, n), Fr = zeros(n);
    for (i = 0; i < n; i++) {
      Fr[i] = F[free[i]];
      for (j = 0; j < n; j++) Kr[i][j] = K[free[i]][free[j]];
    }
    var ur = solveLinear(Kr, Fr);
    if (!ur) throw new Error('MeshLab: plate stiffness matrix is singular');
    var u = zeros(ndof);
    for (i = 0; i < n; i++) u[free[i]] = ur[i];
    var Ku = matVec(K, u), R = zeros(ndof);
    for (i = 0; i < ndof; i++) R[i] = Ku[i] - F[i];
    var elements = elemData.map(function (ed, k) {
      var ue = ed.idx.map(function (q) { return u[q]; });
      var strain = [0, 0, 0], stress = [0, 0, 0];
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 6; j++) strain[i] += ed.B[i][j] * ue[j];
      }
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) stress[i] += D[i][j] * strain[j];
      }
      var sm = stressMeasures(stress[0], stress[1], stress[2]);
      return {
        id: k, nodes: mesh.elements[k].nodes.slice(), area: ed.area,
        centroid: {
          x: (ed.coords[0].x + ed.coords[1].x + ed.coords[2].x) / 3,
          y: (ed.coords[0].y + ed.coords[1].y + ed.coords[2].y) / 3
        },
        strain: { exx: strain[0], eyy: strain[1], gxy: strain[2] },
        stress: { sigmaX: stress[0], sigmaY: stress[1], tauXY: stress[2] },
        sigma1: sm.sigma1, sigma2: sm.sigma2, tauMax: sm.tauMax, vonMises: sm.vonMises,
        strainEnergy: 0.5 * t * ed.area * (
          strain[0] * stress[0] + strain[1] * stress[1] + strain[2] * stress[2])
      };
    });
    var nodalStress = [];
    for (i = 0; i < N; i++) nodalStress.push({ count: 0, sx: 0, sy: 0, txy: 0, vm: 0 });
    elements.forEach(function (e) {
      e.nodes.forEach(function (q) {
        var s = nodalStress[q];
        s.count++; s.sx += e.stress.sigmaX; s.sy += e.stress.sigmaY;
        s.txy += e.stress.tauXY; s.vm += e.vonMises;
      });
    });
    var nodeResults = mesh.nodes.map(function (p, k) {
      var s = nodalStress[k], f = s.count || 1;
      return {
        id: p.id, x: p.x, y: p.y,
        ux: u[2 * k], uy: u[2 * k + 1],
        rx: R[2 * k], ry: R[2 * k + 1],
        vonMises: s.vm / f,
        stress: stressMeasures(s.sx / f, s.sy / f, s.txy / f)
      };
    });
    var vm = elements.map(function (e) { return e.vonMises; });
    var sumFx = 0, sumFy = 0, sumRx = 0, sumRy = 0;
    for (i = 0; i < ndof; i += 2) {
      sumFx += F[i]; sumFy += F[i + 1]; sumRx += R[i]; sumRy += R[i + 1];
    }
    return {
      mesh: mesh, material: { E: E, nu: nu, t: t }, D: D,
      nodes: nodeResults, elements: elements, displacements: Array.from(u),
      reactions: Array.from(R), stiffnessMatrix: K,
      maxVonMises: vm.length ? Math.max.apply(null, vm) : 0,
      minVonMises: vm.length ? Math.min.apply(null, vm) : 0,
      maxDisplacement: Math.max.apply(null, Array.from(u).map(Math.abs).concat([0])),
      maxSigmaX: Math.max.apply(null, elements.map(function (e) { return e.stress.sigmaX; }).concat([0])),
      equilibrium: { sumFx: sumFx, sumFy: sumFy, sumRx: sumRx, sumRy: sumRy,
        residualX: sumFx + sumRx, residualY: sumFy + sumRy }
    };
  }

  // ---------------------------------------------------------------------------
  // Built in structural presets
  // ---------------------------------------------------------------------------

  var STEEL = { E: 2.0e11, rho: 7850, sigmaYield: 2.5e8 };

  function presetMaterial(over) {
    return Object.assign({}, STEEL, over || {});
  }

  function makeModel(name, dim, nodes, members, material) {
    return { name: name, dim: dim, nodes: nodes, members: members,
      material: presetMaterial(material), supports: nodes.map(function (n, i) {
        return n.support ? { node: i, type: n.support } : null;
      }).filter(Boolean), loads: nodes.map(function (n, i) {
        return n.load ? { node: i, f: n.load } : null;
      }).filter(Boolean) };
  }

  function prattBridge() {
    var bays = 7, L = 3, H = 3, nodes = [], members = [], i;
    for (i = 0; i <= bays; i++) nodes.push({ x: i * L, y: 0,
      support: i === 0 ? 'pin' : (i === bays ? 'rollerY' : 'free') });
    for (i = 0; i <= bays; i++) nodes.push({ x: i * L, y: H });
    var top = bays + 1;
    for (i = 0; i < bays; i++) {
      members.push({ a: i, b: i + 1, A: 0.012, I: 6e-5, label: 'Bottom ' + (i + 1) });
      members.push({ a: top + i, b: top + i + 1, A: 0.010, I: 5e-5, label: 'Top ' + (i + 1) });
      members.push({ a: i, b: top + i, A: 0.008, I: 3e-5, label: 'Post ' + (i + 1) });
    }
    members.push({ a: bays, b: top + bays, A: 0.008, I: 3e-5, label: 'End post' });
    for (i = 1; i < bays; i++) {
      members.push({ a: i, b: top + i - 1, A: 0.007, I: 2e-5, label: 'Diagonal ' + i });
    }
    for (i = 1; i < bays; i++) {
      members.push({ a: i, b: top + i + 1, A: 0.007, I: 2e-5, label: 'Diagonal ' + (i + 7) });
    }
    for (i = 1; i < bays; i++) {
      nodes[i].load = [0, -60000];
    }
    return makeModel('Pratt Bridge Truss (7-Bay)', 2, nodes, members);
  }

  function warrenRoof() {
    var bays = 6, span = 12, rise = 3.2, eave = 0.45, nodes = [], members = [], i;
    function roofY(x) { return eave + rise * (1 - Math.abs(x - span / 2) / (span / 2)); }
    for (i = 0; i <= bays; i++) nodes.push({ x: span * i / bays, y: 0,
      support: i === 0 ? 'pin' : (i === bays ? 'rollerY' : 'free') });
    for (i = 0; i <= bays; i++) nodes.push({ x: span * i / bays, y: roofY(span * i / bays) });
    var top = bays + 1;
    for (i = 0; i < bays; i++) {
      members.push({ a: i, b: i + 1, A: 0.009, I: 4e-5, label: 'Tie ' + (i + 1) });
      members.push({ a: top + i, b: top + i + 1, A: 0.008, I: 3e-5, label: 'Rafter ' + (i + 1) });
      if (i % 2 === 0) members.push({ a: i, b: top + i + 1, A: 0.006, I: 2e-5, label: 'Web ' + i });
      else members.push({ a: i + 1, b: top + i, A: 0.006, I: 2e-5, label: 'Web ' + i });
    }
    for (i = 0; i <= bays; i++) {
      members.push({ a: i, b: top + i, A: 0.005, I: 1.5e-5, label: 'Vertical ' + i });
    }
    for (i = 1; i < bays; i++) nodes[top + i].load = [0, -24000];
    return makeModel('Warren Roof Truss', 2, nodes, members);
  }

  function cantileverCrane() {
    var bays = 5, L = 2.5, H = 2, nodes = [], members = [], i;
    for (i = 0; i <= bays; i++) nodes.push({ x: i * L, y: 0 });
    for (i = 0; i <= bays; i++) nodes.push({ x: i * L, y: H });
    var top = bays + 1;
    nodes[0].support = 'fixed';
    nodes[top].support = 'fixed';
    for (i = 0; i < bays; i++) {
      members.push({ a: i, b: i + 1, A: 0.014, I: 8e-5, label: 'Lower boom ' + (i + 1) });
      members.push({ a: top + i, b: top + i + 1, A: 0.014, I: 8e-5, label: 'Upper boom ' + (i + 1) });
      members.push({ a: i, b: top + i, A: 0.009, I: 4e-5, label: 'Post ' + (i + 1) });
      members.push({ a: i, b: top + i + 1, A: 0.009, I: 4e-5, label: 'Brace ' + (i + 1) });
    }
    members.push({ a: bays, b: top + bays, A: 0.009, I: 4e-5, label: 'Tip post' });
    nodes[bays].load = [0, -95000];
    nodes[top + bays].load = [0, -35000];
    return makeModel('Cantilever Crane Boom', 2, nodes, members);
  }

  function transmissionTower() {
    var levels = 6, H = 2.2, nodes = [], members = [], i, k;
    function halfWidth(y) { return 2.6 - 0.28 * y / H; }
    for (k = 0; k <= levels; k++) {
      var w = halfWidth(k * H);
      nodes.push({ x: -w, y: k * H });
      nodes.push({ x: w, y: k * H });
    }
    nodes[0].support = 'fixed'; nodes[1].support = 'fixed';
    for (k = 0; k < levels; k++) {
      var a = 2 * k, b = a + 1, c = a + 2, d = a + 3;
      members.push({ a: a, b: c, A: 0.011, I: 5e-5, label: 'Left leg ' + (k + 1) });
      members.push({ a: b, b: d, A: 0.011, I: 5e-5, label: 'Right leg ' + (k + 1) });
      members.push({ a: c, b: d, A: 0.007, I: 2e-5, label: 'Cross arm ' + (k + 1) });
      members.push({ a: a, b: d, A: 0.006, I: 1.8e-5, label: 'X brace A' + (k + 1) });
      members.push({ a: b, b: c, A: 0.006, I: 1.8e-5, label: 'X brace B' + (k + 1) });
      members.push({ a: a, b: b, A: 0.005, I: 1.5e-5, label: 'Base tie ' + (k + 1) });
    }
    var last = 2 * levels;
    nodes[last].load = [-22000, -12000];
    nodes[last + 1].load = [22000, -12000];
    return makeModel('Transmission Tower Lattice', 2, nodes, members);
  }

  function spaceTripodDome() {
    var nodes = [], members = [], i, j;
    var R = 5, ringCount = 8, crownY = 3.4;
    for (i = 0; i < 3; i++) {
      var ang = 2 * Math.PI * i / 3;
      nodes.push({ x: R * Math.cos(ang), y: 0, z: R * Math.sin(ang), support: 'pin' });
    }
    var level1 = 3;
    for (i = 0; i < ringCount; i++) {
      var a1 = 2 * Math.PI * i / ringCount;
      nodes.push({ x: 3.2 * Math.cos(a1), y: 1.8, z: 3.2 * Math.sin(a1) });
    }
    var level2 = level1 + ringCount;
    for (i = 0; i < ringCount; i++) {
      var a2 = 2 * Math.PI * (i + 0.5) / ringCount;
      nodes.push({ x: 1.7 * Math.cos(a2), y: 2.9, z: 1.7 * Math.sin(a2) });
    }
    var crown = nodes.length;
    nodes.push({ x: 0, y: crownY, z: 0, load: [0, -80000, 0] });
    for (i = 0; i < 3; i++) {
      for (j = 0; j < ringCount; j++) {
        if (j % 3 === i) members.push({ a: i, b: level1 + j, A: 0.010, I: 5e-5, label: 'Strut ' + i + '-' + j });
        if (j % 3 === (i + 1) % 3) members.push({ a: i, b: level1 + j, A: 0.008, I: 4e-5, label: 'Tie strut ' + i + '-' + j });
      }
    }
    for (i = 0; i < ringCount; i++) {
      members.push({ a: level1 + i, b: level1 + (i + 1) % ringCount, A: 0.008, I: 3e-5, label: 'Ring 1-' + i });
      members.push({ a: level1 + i, b: level2 + i, A: 0.007, I: 2e-5, label: 'Web A' + i });
      members.push({ a: level1 + i, b: level2 + (i + ringCount - 1) % ringCount, A: 0.007, I: 2e-5, label: 'Web B' + i });
      members.push({ a: level2 + i, b: level2 + (i + 1) % ringCount, A: 0.006, I: 2e-5, label: 'Ring 2-' + i });
      members.push({ a: level2 + i, b: crown, A: 0.009, I: 4e-5, label: 'Crown strut ' + i });
    }
    return makeModel('3D Space Tripod Dome', 3, nodes, members);
  }

  function perforatedPlate() {
    var mesh = generatePlateMesh(20, 16, 2.4, 1.6, 0.28);
    var model = makeModel('Perforated Plate in Tension (Kirsch Stress Concentration)', 2, [], [],
      { E: 2.0e11, sigmaYield: 2.5e8 });
    model.plate = {
      mesh: mesh,
      material: { E: 2.0e11, nu: 0.3, t: 0.012 },
      bc: { left: 'rollerX', fixedNodes: [mesh.grid[Math.floor(mesh.ny / 2)][0]] },
      loads: { rightNodesForce: 480000, leftNodesForce: -480000 }
    };
    return model;
  }

  var PRESET_BUILDERS = {
    pratt: prattBridge,
    warren: warrenRoof,
    crane: cantileverCrane,
    tower: transmissionTower,
    dome: spaceTripodDome,
    plate: perforatedPlate
  };

  var PRESET_LIST = [
    { key: 'pratt', name: 'Pratt Bridge Truss (7-Bay)', dim: 2, kind: 'truss' },
    { key: 'warren', name: 'Warren Roof Truss', dim: 2, kind: 'truss' },
    { key: 'crane', name: 'Cantilever Crane Boom', dim: 2, kind: 'truss' },
    { key: 'tower', name: 'Transmission Tower Lattice', dim: 2, kind: 'truss' },
    { key: 'dome', name: '3D Space Tripod Dome', dim: 3, kind: 'truss' },
    { key: 'plate', name: 'Perforated Plate in Tension (Kirsch Stress Concentration)', dim: 2, kind: 'plate' }
  ];

  function getPresets() { return PRESET_LIST.map(function (p) { return Object.assign({}, p); }); }

  function getPreset(key) {
    var builder = PRESET_BUILDERS[key];
    if (!builder) throw new Error('MeshLab: unknown preset "' + key + '"');
    return builder();
  }

  // ---------------------------------------------------------------------------
  // Exporters
  // ---------------------------------------------------------------------------

  function fmt(v) {
    if (!isFinite(v)) return '0';
    if (v !== 0 && (Math.abs(v) >= 1e6 || Math.abs(v) < 1e-4)) return v.toExponential(6);
    return v.toPrecision(10).replace(/\.?0+$/, '');
  }

  function exportCalculiXInp(modelInput) {
    var model = normalizeModel(modelInput);
    applySupports(model, modelInput);
    var d = model.dim, lines = [];
    lines.push('** MeshLab CalculiX / Abaqus input deck');
    lines.push('** Model: ' + model.name.replace(/[\r\n]/g, ' '));
    lines.push('*HEADING');
    lines.push('MeshLab linear static truss analysis, SI units');
    lines.push('*NODE, NSET=NALL');
    model.nodes.forEach(function (n, i) {
      var coords = d === 3 ? [n.x, n.y, n.z] : [n.x, n.y];
      lines.push((i + 1) + ', ' + coords.map(fmt).join(', '));
    });
    var groups = {};
    model.members.forEach(function (m, i) {
      var key = m.E + '|' + m.A + '|' + m.rho;
      if (!groups[key]) groups[key] = { E: m.E, A: m.A, rho: m.rho, ids: [] };
      groups[key].ids.push(i + 1);
    });
    var gi = 0;
    Object.keys(groups).forEach(function (key) {
      var g = groups[key];
      gi++;
      lines.push('*ELEMENT, TYPE=' + (d === 3 ? 'T3D2' : 'T2D2') + ', ELSET=EG' + gi);
      g.ids.forEach(function (id) {
        var m = model.members[id - 1];
        lines.push(id + ', ' + (m.a + 1) + ', ' + (m.b + 1));
      });
      lines.push('*SOLID SECTION, ELSET=EG' + gi + ', MATERIAL=STEEL' + gi);
      lines.push(fmt(g.A) + ',');
      lines.push('*MATERIAL, NAME=STEEL' + gi);
      lines.push('*ELASTIC');
      lines.push(fmt(g.E) + ', 0.3');
      lines.push('*DENSITY');
      lines.push(fmt(g.rho));
    });
    lines.push('*ELSET, ELSET=EALL');
    var allIds = [];
    model.members.forEach(function (m, i) { allIds.push(i + 1); });
    for (var q = 0; q < allIds.length; q += 8) {
      lines.push(allIds.slice(q, q + 8).join(', '));
    }
    lines.push('*BOUNDARY');
    var dofNames = d === 3 ? ['1', '2', '3'] : ['1', '2'];
    var masks = { pin: [1, 1, 1], fixed: [1, 1, 1],
      rollerX: [1, 0, 0], rollerY: [0, 1, 0] };
    model.nodes.forEach(function (n, i) {
      var mask = masks[n.support];
      if (!mask) return;
      for (var k = 0; k < d; k++) {
        if (mask[k]) lines.push((i + 1) + ', ' + dofNames[k]);
      }
    });
    lines.push('*STEP');
    lines.push('*STATIC');
    lines.push('*CLOAD');
    var loadNames = d === 3 ? [1, 2, 3] : [1, 2];
    model.nodes.forEach(function (n, i) {
      for (var k = 0; k < d; k++) {
        if (Math.abs(n.load[k]) > 0) lines.push((i + 1) + ', ' + loadNames[k] + ', ' + fmt(n.load[k]));
      }
    });
    lines.push('*NODE PRINT, NSET=NALL');
    lines.push('U');
    lines.push('*NODE FILE');
    lines.push('U');
    lines.push('*EL PRINT, ELSET=EALL');
    lines.push('S');
    lines.push('*END STEP');
    return lines.join('\n') + '\n';
  }

  function xmlEscape(s) {
    return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;')
      .replace(/>/g, '&gt;').replace(/"/g, '&quot;').replace(/'/g, '&apos;');
  }

  function exportDeformedSvg(modelInput, solutionInput, scale) {
    var sol = solutionInput || solveTruss(modelInput);
    var model = sol.model;
    var s = num(scale, 100);
    var d = model.dim;
    var xs = model.nodes.map(function (n) { return n.x; });
    var ys = model.nodes.map(function (n) { return n.y; });
    var minX = Math.min.apply(null, xs), maxX = Math.max.apply(null, xs);
    var minY = Math.min.apply(null, ys), maxY = Math.max.apply(null, ys);
    var w = Math.max(maxX - minX, 1), h = Math.max(maxY - minY, 1);
    var pad = 40, W = 900, H = 560;
    var unit = Math.min((W - 2 * pad) / w, (H - 2 * pad) / h);
    function project(x, y, dx, dy) {
      var px = pad + (x - minX + (dx || 0) * s) * unit;
      var py = H - pad - (y - minY + (dy || 0) * s) * unit;
      return [px, py];
    }
    var out = [];
    out.push('<?xml version="1.0" encoding="UTF-8"?>');
    out.push('<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '">');
    out.push('<title>' + xmlEscape(model.name) + ' deformed shape</title>');
    out.push('<desc>MeshLab direct stiffness result. Undeformed geometry in ghost lines, deformed shape at ' + fmt(s) + 'x scale.</desc>');
    out.push('<rect width="100%" height="100%" fill="#09111E"/>');
    out.push('<g stroke="#223550" stroke-width="1" opacity="0.5">');
    for (var gx = 0; gx <= W; gx += 40) out.push('<line x1="' + gx + '" y1="0" x2="' + gx + '" y2="' + H + '"/>');
    for (var gy = 0; gy <= H; gy += 40) out.push('<line x1="0" y1="' + gy + '" x2="' + W + '" y2="' + gy + '"/>');
    out.push('</g>');
    out.push('<g stroke="#5b7899" stroke-width="1.4" stroke-dasharray="5 4" fill="none" opacity="0.75">');
    model.members.forEach(function (m) {
      var p = project(model.nodes[m.a].x, model.nodes[m.a].y);
      var q = project(model.nodes[m.b].x, model.nodes[m.b].y);
      out.push('<line x1="' + fmt(p[0]) + '" y1="' + fmt(p[1]) + '" x2="' + fmt(q[0]) + '" y2="' + fmt(q[1]) + '"/>');
    });
    out.push('</g>');
    out.push('<g stroke-width="3" stroke-linecap="round">');
    sol.members.forEach(function (m) {
      var na = model.nodes[m.a], nb = model.nodes[m.b];
      var ua = sol.displacements.slice(d * m.a, d * m.a + d);
      var ub = sol.displacements.slice(d * m.b, d * m.b + d);
      var p = project(na.x, na.y, ua[0], ua[1]);
      var q = project(nb.x, nb.y, ub[0], ub[1]);
      var color = m.axialForce > 0 ? '#3987e5' : (m.axialForce < 0 ? '#e66767' : '#8898ac');
      out.push('<line x1="' + fmt(p[0]) + '" y1="' + fmt(p[1]) + '" x2="' + fmt(q[0]) + '" y2="' + fmt(q[1]) + '" stroke="' + color + '"/>');
    });
    out.push('</g>');
    out.push('<g fill="#cfe0f5" font-family="monospace" font-size="13">');
    out.push('<text x="20" y="26">MeshLab / ' + xmlEscape(model.name) + '</text>');
    out.push('<text x="20" y="46">Deformation scale ' + fmt(s) + 'x</text>');
    out.push('<text x="' + (W - 260) + '" y="26" fill="#3987e5">blue = tension</text>');
    out.push('<text x="' + (W - 260) + '" y="46" fill="#e66767">red = compression</text>');
    out.push('</g>');
    out.push('</svg>');
    return out.join('\n') + '\n';
  }

  function exportReportJson(modelInput, solutionInput) {
    var sol = solutionInput || solveTruss(modelInput);
    var model = sol.model;
    var safety = analyzeSafety(model, sol);
    var modes = null;
    try { modes = solveModes(model, 4); } catch (e) { modes = null; }
    var report = {
      format: 'MeshLab structural report v1.0.0',
      generated: 'static analysis',
      units: { length: 'm', force: 'N', stress: 'Pa', mass: 'kg', frequency: 'Hz' },
      model: {
        name: model.name, dim: model.dim, nodeCount: model.nodes.length,
        memberCount: model.members.length,
        material: { E: model.material.E, rho: model.material.rho, sigmaYield: model.material.sigmaYield }
      },
      determinacy: sol.determinacy,
      equilibrium: sol.equilibrium,
      maxDisplacement: sol.maxDisplacement,
      strainEnergy: sol.strainEnergy,
      externalWork: sol.externalWork,
      nodes: sol.nodes.map(function (n) {
        return { id: n.id, x: n.x, y: n.y, z: n.z, support: n.support,
          load: n.load, displacement: n.displacement, reaction: n.reaction };
      }),
      members: sol.members.map(function (m) {
        return { id: m.id, label: m.label, a: m.a, b: m.b, length: m.length,
          axialForce: m.axialForce, stress: m.stress, strain: m.strain,
          strainEnergy: m.strainEnergy, tension: m.tension };
      }),
      safety: {
        minFos: safety.minFos, safe: safety.safe, criticalCount: safety.criticalCount,
        members: safety.members.map(function (r) {
          return { id: r.id, label: r.label, eulerCriticalLoad: r.eulerCriticalLoad,
            slendernessRatio: r.slendernessRatio, fosYield: r.fosYield,
            fosBuckling: r.fosBuckling, fosGoverning: r.fosGoverning,
            governingMode: r.governingMode, critical: r.critical };
        })
      },
      modes: modes ? modes.modes.map(function (m) {
        return { index: m.index, frequencyHz: m.frequencyHz, omega: m.omega,
          eigenvalue: m.eigenvalue, modeShape: m.modeShape };
      }) : []
    };
    return JSON.stringify(report, null, 2) + '\n';
  }

  // ---------------------------------------------------------------------------
  // Public API
  // ---------------------------------------------------------------------------

  return {
    version: '1.0.0',
    solveTruss: solveTruss,
    analyzeSafety: analyzeSafety,
    solveModes: solveModes,
    solvePlateCST: solvePlateCST,
    generatePlateMesh: generatePlateMesh,
    planeStressD: planeStressD,
    stressMeasures: stressMeasures,
    getPresets: getPresets,
    getPreset: getPreset,
    presets: PRESET_LIST,
    exportCalculiXInp: exportCalculiXInp,
    exportDeformedSvg: exportDeformedSvg,
    exportReportJson: exportReportJson,
    normalizeModel: normalizeModel,
    solveLinear: solveLinear,
    jacobiEigen: jacobiEigen,
    matVec: matVec,
    matMul: matMul,
    dot: dot,
    STEEL: STEEL
  };
});
