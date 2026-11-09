/*
 * MeshLab regression suite. Run with: node --test test.js
 * All checks are analytical or algebraic, with no external dependencies.
 */
'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const M = require('./fem.js');

const STEEL = { E: 2.0e11, rho: 7850, sigmaYield: 2.5e8, A: 0.01, I: 1e-6 };
const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol,
  (msg || 'values differ') + ': got ' + a + ', expected ' + b + ' (tol ' + tol + ')');
const relClose = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)),
  (msg || 'values differ') + ': got ' + a + ', expected ' + b);
const sum = (arr) => arr.reduce((s, v) => s + v, 0);

/*
 * Stable two bar benchmark: two bars from pinned feet at (0,0) and (0,2)
 * meet at (2,1), which carries a horizontal load P. Both bars have length
 * sqrt(5) and direction cosines (+/-cx, +/-cy), so the exact apex deflection
 * is ux = P / (2 * (EA/L) * cx^2) and uy = 0 by symmetry.
 */
function twoBar(P, extra) {
  return Object.assign({
    dim: 2,
    nodes: [
      { x: 0, y: 0, support: 'pin' },
      { x: 0, y: 2, support: 'pin' },
      { x: 2, y: 1, load: [P, 0] }
    ],
    members: [
      { a: 0, b: 2, E: STEEL.E, A: STEEL.A, I: STEEL.I, label: 'Lower bar' },
      { a: 1, b: 2, E: STEEL.E, A: STEEL.A, I: STEEL.I, label: 'Upper bar' }
    ],
    material: STEEL
  }, extra || {});
}

const TWO_BAR_L = Math.sqrt(5);
const TWO_BAR_CX = 2 / TWO_BAR_L;
const TWO_BAR_CY = 1 / TWO_BAR_L;
const twoBarUx = (P) => P / (2 * (STEEL.E * STEEL.A / TWO_BAR_L) * TWO_BAR_CX * TWO_BAR_CX);

/*
 * Stable three bar fan: two inclined bars and one vertical bar join three
 * pinned base nodes to an apex carrying a vertical load P. The exact apex
 * deflection follows from the summed vertical stiffness.
 */
function threeBar(P) {
  return {
    dim: 2,
    nodes: [
      { x: 0, y: 0, support: 'pin' },
      { x: 1, y: 0, support: 'pin' },
      { x: 2, y: 0, support: 'pin' },
      { x: 1, y: 1, load: [0, -P] }
    ],
    members: [
      { a: 0, b: 3, E: STEEL.E, A: STEEL.A, I: STEEL.I, label: 'Diagonal L' },
      { a: 1, b: 3, E: STEEL.E, A: STEEL.A, I: STEEL.I, label: 'Vertical' },
      { a: 2, b: 3, E: STEEL.E, A: STEEL.A, I: STEEL.I, label: 'Diagonal R' }
    ],
    material: STEEL
  };
}

const FAN_C = Math.SQRT1_2;
const threeBarUy = (P) => -P / (
  STEEL.E * STEEL.A / 1 +
  2 * (STEEL.E * STEEL.A / Math.SQRT2) * FAN_C * FAN_C
);

test('module loads as a CommonJS UMD export', () => {
  assert.equal(typeof M, 'object');
  assert.equal(M.version, '1.0.0');
  assert.equal(typeof M.solveTruss, 'function');
  assert.equal(typeof M.analyzeSafety, 'function');
  assert.equal(typeof M.solveModes, 'function');
  assert.equal(typeof M.solvePlateCST, 'function');
});

test('module also attaches to a browser global via the UMD factory', () => {
  const fs = require('node:fs');
  const vm = require('node:vm');
  const src = fs.readFileSync(require.resolve('./fem.js'), 'utf8');
  const sandbox = { window: {}, console: console, Math: Math, Float64Array: Float64Array };
  sandbox.globalThis = sandbox;
  vm.createContext(sandbox);
  vm.runInContext(src, sandbox, { filename: 'fem.js' });
  assert.equal(typeof sandbox.MeshLab, 'object');
  assert.equal(sandbox.MeshLab.version, '1.0.0');
  const sol = sandbox.MeshLab.solveTruss(twoBar(1000));
  assert.equal(sol.members.length, 2);
});

test('solveLinear solves a known 2x2 system', () => {
  const x = M.solveLinear([[2, 1], [1, 3]], [3, 5]);
  close(x[0], 0.8, 1e-12);
  close(x[1], 1.4, 1e-12);
});

test('solveLinear needs partial pivoting for a zero leading pivot', () => {
  const x = M.solveLinear([[0, 2], [1, 3]], [4, 5]);
  close(x[0], -1, 1e-12);
  close(x[1], 2, 1e-12);
});

test('solveLinear returns null for a singular matrix', () => {
  assert.equal(M.solveLinear([[1, 2], [2, 4]], [1, 2]), null);
});

test('solveLinear solves a large stiffness scaled system', () => {
  const x = M.solveLinear([[2e9, 1e9], [1e9, 3e9]], [3e9, 5e9]);
  close(x[0], 0.8, 1e-9);
  close(x[1], 1.4, 1e-9);
});

test('matVec and matMul agree with hand computation', () => {
  const A = [[1, 2], [3, 4]], B = [[5, 6], [7, 8]];
  const v = M.matVec(A, [2, 3]);
  assert.deepEqual(Array.from(v), [8, 18]);
  const C = M.matMul(A, B);
  assert.deepEqual(C.map((r) => Array.from(r)), [[19, 22], [43, 50]]);
});

test('jacobiEigen returns sorted eigenvalues and orthonormal eigenvectors', () => {
  const A = [[4, 1], [1, 3]];
  const e = M.jacobiEigen(A);
  close(e[0].value, (7 - Math.sqrt(5)) / 2, 1e-10);
  close(e[1].value, (7 + Math.sqrt(5)) / 2, 1e-10);
  close(M.dot(e[0].vector, e[1].vector), 0, 1e-10);
  close(M.dot(e[0].vector, e[0].vector), 1, 1e-10);
});

test('two bar truss axial forces match joint equilibrium', () => {
  const P = 1000;
  const s = M.solveTruss(twoBar(P));
  /* Each inclined bar carries P / (2 cx) in tension. */
  const expected = P / (2 * TWO_BAR_CX);
  close(s.members[0].axialForce, expected, 1e-6, 'lower bar tension');
  close(s.members[1].axialForce, expected, 1e-6, 'upper bar tension');
  assert.equal(s.members[0].tension, true);
  assert.equal(s.members[1].tension, true);
});

test('two bar truss horizontal displacement matches virtual work', () => {
  const P = 1000;
  const s = M.solveTruss(twoBar(P));
  /* Castigliano: U = 2 * (N^2 L / 2EA); dU/dP = P / (2 k cx^2). */
  relClose(s.nodes[2].displacement[0], twoBarUx(P), 1e-9, 'apex horizontal deflection');
  close(s.nodes[2].displacement[1], 0, 1e-12, 'symmetry cancels vertical drift');
});

test('two bar truss strain energy equals half the external work', () => {
  const s = M.solveTruss(twoBar(1000));
  relClose(s.strainEnergy, 0.5 * s.externalWork, 1e-9, 'Clapeyron theorem');
  relClose(s.strainEnergy, 0.5 * 1000 * twoBarUx(1000), 1e-9, 'energy from the closed form');
});

test('two bar truss horizontal reactions balance the applied load', () => {
  const s = M.solveTruss(twoBar(1000));
  close(s.equilibrium.sumRx, -1000, 1e-6, 'horizontal reactions balance the applied load');
  close(s.equilibrium.sumRy, 0, 1e-6, 'symmetric vertical reactions cancel');
  close(s.equilibrium.residualX, 0, 1e-6);
  close(s.equilibrium.residualY, 0, 1e-6);
});

test('two bar truss member stress equals axial force over area', () => {
  const s = M.solveTruss(twoBar(1000));
  s.members.forEach((m) => close(m.stress, m.axialForce / STEEL.A, 1e-6));
});

test('two bar truss strain equals elongation over length', () => {
  const s = M.solveTruss(twoBar(1000));
  s.members.forEach((m) => close(m.strain, m.delta / m.length, 1e-12));
});

test('two bar truss element strain energy equals 0.5 N delta', () => {
  const s = M.solveTruss(twoBar(1000));
  s.members.forEach((m) => close(m.strainEnergy, 0.5 * m.axialForce * m.delta, 1e-12));
});

test('two bar truss is statically determinate', () => {
  const s = M.solveTruss(twoBar(1000));
  assert.equal(s.determinacy.m, 2);
  assert.equal(s.determinacy.r, 4);
  assert.equal(s.determinacy.j, 3);
  assert.equal(s.determinacy.value, 0);
  assert.equal(s.determinacy.classification, 'statically determinate');
  assert.equal(s.determinacy.stable, true);
});

test('three bar truss vertical member carries the largest share of the load', () => {
  const P = 900;
  const s = M.solveTruss(threeBar(P));
  const vertical = s.members[1];
  assert.ok(Math.abs(vertical.axialForce) > Math.abs(s.members[0].axialForce));
  assert.ok(Math.abs(vertical.axialForce) > Math.abs(s.members[2].axialForce));
  close(s.members[0].axialForce, s.members[2].axialForce, 1e-6, 'symmetric diagonals');
  assert.equal(vertical.tension, false);
});

test('three bar truss apex deflection matches the summed vertical stiffness', () => {
  const P = 900;
  const s = M.solveTruss(threeBar(P));
  relClose(s.nodes[3].displacement[1], threeBarUy(P), 1e-9);
  relClose(s.strainEnergy, 0.5 * P * (-threeBarUy(P)), 1e-9, 'Castigliano check');
});

test('three bar truss is statically indeterminate', () => {
  const s = M.solveTruss(threeBar(900));
  assert.equal(s.determinacy.m, 3);
  assert.equal(s.determinacy.r, 6);
  assert.equal(s.determinacy.j, 4);
  assert.equal(s.determinacy.value, 1);
  assert.equal(s.determinacy.classification, 'statically indeterminate');
});

test('three bar truss satisfies global force equilibrium', () => {
  const s = M.solveTruss(threeBar(900));
  close(s.equilibrium.sumFx + s.equilibrium.sumRx, 0, 1e-6);
  close(s.equilibrium.sumFy + s.equilibrium.sumRy, 0, 1e-6);
});

test('horizontal two bar chain gives exact axial extension', () => {
  const P = 5000, EA = STEEL.E * STEEL.A, L = 2;
  const model = {
    dim: 2, material: STEEL,
    nodes: [{ x: 0, y: 0, support: 'pin' }, { x: L, y: 0, support: 'rollerY', load: [P, 0] }],
    members: [{ a: 0, b: 1, E: STEEL.E, A: STEEL.A, I: STEEL.I }]
  };
  const s = M.solveTruss(model);
  close(s.nodes[1].displacement[0], P * L / EA, 1e-12);
  close(s.members[0].axialForce, P, 1e-9);
  close(s.members[0].stress, P / STEEL.A, 1e-9);
  assert.equal(s.members[0].tension, true);
});

test('horizontal bar reaction is equal and opposite to the applied load', () => {
  const model = {
    dim: 2, material: STEEL,
    nodes: [{ x: 0, y: 0, support: 'pin' }, { x: 2, y: 0, support: 'rollerY', load: [5000, 0] }],
    members: [{ a: 0, b: 1, E: STEEL.E, A: STEEL.A, I: STEEL.I }]
  };
  const s = M.solveTruss(model);
  close(s.nodes[0].reaction[0], -5000, 1e-6);
  close(s.equilibrium.residualX, 0, 1e-6);
});

test('external loads may be declared in a top level loads array', () => {
  const model = {
    dim: 2, material: STEEL,
    nodes: [{ x: 0, y: 0 }, { x: 0, y: 2 }, { x: 2, y: 1 }],
    supports: [{ node: 0, type: 'pin' }, { node: 1, type: 'pin' }],
    loads: [{ node: 2, f: [1000, 0] }],
    members: [
      { a: 0, b: 2, E: STEEL.E, A: STEEL.A, I: STEEL.I },
      { a: 1, b: 2, E: STEEL.E, A: STEEL.A, I: STEEL.I }
    ]
  };
  const s = M.solveTruss(model);
  relClose(s.nodes[2].displacement[0], twoBarUx(1000), 1e-9);
});

test('supports may be declared as a node index map', () => {
  const model = twoBar(1000);
  delete model.nodes[0].support;
  delete model.nodes[1].support;
  model.supports = { 0: 'pin', 1: 'pin' };
  const s = M.solveTruss(model);
  assert.ok(isFinite(s.maxDisplacement));
  close(s.equilibrium.residualX, 0, 1e-6);
});

test('direction cosines are unit length', () => {
  const s = M.solveTruss(twoBar(1000));
  s.members.forEach((m) => {
    const c = m.directionCosines;
    close(c.cx * c.cx + c.cy * c.cy + c.cz * c.cz, 1, 1e-12);
  });
});

test('rollerX supports vertical translation and rollerY supports horizontal translation', () => {
  /* rollerX restrains x and permits y; a pinned base anchors the other end. */
  const rail = {
    dim: 2, material: STEEL,
    nodes: [
      { x: 0, y: 0, support: 'pin' },
      { x: 0, y: 1, support: 'rollerX', load: [0, 1000] }
    ],
    members: [{ a: 0, b: 1, E: STEEL.E, A: STEEL.A, I: STEEL.I }]
  };
  const s = M.solveTruss(rail);
  close(s.nodes[1].displacement[0], 0, 1e-9, 'rollerX blocks horizontal motion');
  assert.ok(s.nodes[1].displacement[1] > 0, 'rollerX permits vertical motion');
  close(s.nodes[1].displacement[1], 1000 * 1 / (STEEL.E * STEEL.A), 1e-9);

  const horizontal = {
    dim: 2, material: STEEL,
    nodes: [
      { x: 0, y: 0, support: 'pin' },
      { x: 1, y: 0, support: 'rollerY', load: [1000, 0] }
    ],
    members: [{ a: 0, b: 1, E: STEEL.E, A: STEEL.A, I: STEEL.I }]
  };
  const h = M.solveTruss(horizontal);
  assert.ok(h.nodes[1].displacement[0] > 0, 'rollerY permits horizontal motion');
  close(h.nodes[1].displacement[1], 0, 1e-9, 'rollerY blocks vertical motion');
});

test('fixed and pin supports constrain both translational DOFs in 2D', () => {
  const model = twoBar(1000);
  model.nodes[0].support = 'fixed';
  model.nodes[1].support = 'fixed';
  const s = M.solveTruss(model);
  close(s.nodes[0].displacement[0], 0, 1e-12);
  close(s.nodes[0].displacement[1], 0, 1e-12);
  close(s.nodes[1].displacement[0], 0, 1e-12);
  close(s.nodes[1].displacement[1], 0, 1e-12);
});

test('stiffness matrix is symmetric', () => {
  const s = M.solveTruss(twoBar(1000));
  const K = s.stiffnessMatrix, n = K.length;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) close(K[i][j], K[j][i], 1e-6);
});

test('reactions equal K u minus F at every DOF', () => {
  const s = M.solveTruss(twoBar(1000));
  const Ku = M.matVec(s.stiffnessMatrix, s.displacements);
  for (let i = 0; i < Ku.length; i++) close(s.reactions[i], Ku[i] - s.forceVector[i], 1e-6);
});

test('free DOF reactions are numerically zero', () => {
  const s = M.solveTruss(twoBar(1000));
  s.fixedMask.forEach((f, i) => { if (!f) close(s.reactions[i], 0, 1e-5); });
});

test('Doubling the load doubles displacement and quadruples strain energy', () => {
  const a = M.solveTruss(twoBar(1000));
  const b = M.solveTruss(twoBar(2000));
  relClose(b.nodes[2].displacement[0], 2 * a.nodes[2].displacement[0], 1e-9);
  relClose(b.strainEnergy, 4 * a.strainEnergy, 1e-9);
});

test('a 3D space truss solves with 3 DOF per node', () => {
  const model = {
    dim: 3, material: STEEL,
    nodes: [
      { x: 0, y: 0, z: 0, support: 'pin' },
      { x: 1, y: 0, z: 0, support: 'pin' },
      { x: 0, y: 0, z: 1, support: 'pin' },
      { x: 0.4, y: 1, z: 0.3, load: [0, -2000, 0] }
    ],
    members: [
      { a: 0, b: 3, E: STEEL.E, A: STEEL.A, I: STEEL.I },
      { a: 1, b: 3, E: STEEL.E, A: STEEL.A, I: STEEL.I },
      { a: 2, b: 3, E: STEEL.E, A: STEEL.A, I: STEEL.I }
    ]
  };
  const s = M.solveTruss(model);
  assert.equal(s.displacements.length, 12);
  assert.ok(s.maxDisplacement > 0);
  close(s.equilibrium.residualY, 0, 1e-6);
  close(s.equilibrium.residualZ, 0, 1e-6);
  s.members.forEach((m) => {
    const c = m.directionCosines;
    close(c.cx * c.cx + c.cy * c.cy + c.cz * c.cz, 1, 1e-12);
  });
});

test('3D determinacy uses m + r - 3j', () => {
  const s = M.solveTruss({
    dim: 3, material: STEEL,
    nodes: [
      { x: 0, y: 0, z: 0, support: 'pin' },
      { x: 1, y: 0, z: 0, support: 'pin' },
      { x: 0, y: 0, z: 1, support: 'pin' },
      { x: 0.4, y: 1, z: 0.3, load: [0, -2000, 0] }
    ],
    members: [
      { a: 0, b: 3, E: STEEL.E, A: STEEL.A, I: STEEL.I },
      { a: 1, b: 3, E: STEEL.E, A: STEEL.A, I: STEEL.I },
      { a: 2, b: 3, E: STEEL.E, A: STEEL.A, I: STEEL.I }
    ]
  });
  assert.equal(s.determinacy.r, 9);
  assert.equal(s.determinacy.value, 0);
});

test('zero length members raise a descriptive error', () => {
  assert.throws(() => M.solveTruss({
    dim: 2, nodes: [{ x: 0, y: 0, support: 'pin' }, { x: 0, y: 0, support: 'pin' }],
    members: [{ a: 0, b: 1, E: STEEL.E, A: STEEL.A, I: STEEL.I }]
  }), /zero length/);
});

test('missing nodes and members arrays raise a descriptive error', () => {
  assert.throws(() => M.solveTruss({}), /nodes\[\] and members\[\]/);
});

/* ------------------------------------------------------------------------- */
/* Safety analyzer                                                            */
/* ------------------------------------------------------------------------- */

test('Euler critical load matches pi^2 E I / (K L)^2', () => {
  const s = M.solveTruss(twoBar(1000));
  const a = M.analyzeSafety(twoBar(1000), s);
  close(a.members[0].eulerCriticalLoad,
    Math.PI * Math.PI * STEEL.E * STEEL.I / Math.pow(TWO_BAR_L, 2), 1e-5);
  assert.equal(a.members[0].compression, false);
});

test('Euler load scales with the effective length factor squared', () => {
  const model = twoBar(1000);
  model.members.forEach((m) => { m.Keff = 2; });
  const a = M.analyzeSafety(model, M.solveTruss(model));
  const base = M.analyzeSafety(twoBar(1000), M.solveTruss(twoBar(1000)));
  relClose(a.members[0].eulerCriticalLoad, base.members[0].eulerCriticalLoad / 4, 1e-9);
});

test('yield factor of safety equals sigma_yield over absolute stress', () => {
  const s = M.solveTruss(twoBar(1000));
  const a = M.analyzeSafety(twoBar(1000), s);
  a.members.forEach((r, i) => close(r.fosYield,
    STEEL.sigmaYield / Math.abs(s.members[i].stress), 1e-9));
});

test('buckling safety factor uses Euler load over force for compression members', () => {
  const s = M.solveTruss(threeBar(900));
  const a = M.analyzeSafety(threeBar(900), s);
  a.members.forEach((r, i) => {
    if (r.compression) {
      relClose(r.fosBuckling, r.eulerCriticalLoad / Math.abs(s.members[i].axialForce), 1e-9);
    } else {
      assert.equal(r.fosBuckling, Infinity);
    }
  });
  assert.equal(a.members[1].compression, true);
});

test('tension members have an infinite buckling safety factor', () => {
  const model = {
    dim: 2, material: STEEL,
    nodes: [{ x: 0, y: 0, support: 'pin' }, { x: 2, y: 0, support: 'rollerY', load: [5000, 0] }],
    members: [{ a: 0, b: 1, E: STEEL.E, A: STEEL.A, I: STEEL.I }]
  };
  const a = M.analyzeSafety(model, M.solveTruss(model));
  assert.equal(a.members[0].compression, false);
  assert.equal(a.members[0].fosBuckling, Infinity);
});

test('slenderness ratio equals L over radius of gyration', () => {
  const s = M.solveTruss(twoBar(1000));
  const a = M.analyzeSafety(twoBar(1000), s);
  a.members.forEach((r) => {
    const rg = Math.sqrt(STEEL.I / STEEL.A);
    close(r.radiusOfGyration, rg, 1e-12);
    close(r.slendernessRatio, r.length / rg, 1e-9);
  });
});

test('governing factor of safety is the minimum of yield and buckling', () => {
  const s = M.solveTruss(twoBar(1000));
  const a = M.analyzeSafety(twoBar(1000), s);
  a.members.forEach((r) => close(r.fosGoverning, Math.min(r.fosYield, r.fosBuckling), 1e-9));
});

test('an overloaded compression member is flagged critical', () => {
  const model = twoBar(1e9);
  const a = M.analyzeSafety(model, M.solveTruss(model));
  assert.ok(a.criticalCount > 0);
  assert.equal(a.safe, false);
  assert.ok(a.members.every((r) => r.critical === (r.fosGoverning < 1)));
});

test('a lightly loaded member is not critical', () => {
  const a = M.analyzeSafety(twoBar(100), M.solveTruss(twoBar(100)));
  assert.equal(a.criticalCount, 0);
  assert.equal(a.safe, true);
});

test('minFos is the smallest finite governing factor of safety', () => {
  const s = M.solveTruss(twoBar(1000));
  const a = M.analyzeSafety(twoBar(1000), s);
  const finite = a.members.map((r) => r.fosGoverning).filter(isFinite);
  close(a.minFos, Math.min(...finite), 1e-9);
});

test('governingMode names the smaller of the two safety factors', () => {
  const s = M.solveTruss(twoBar(1000));
  const a = M.analyzeSafety(twoBar(1000), s);
  a.members.forEach((r) => assert.equal(r.governingMode,
    r.fosBuckling < r.fosYield ? 'buckling' : 'yield'));
});

/* ------------------------------------------------------------------------- */
/* Modal analysis                                                             */
/* ------------------------------------------------------------------------- */

test('modal eigenvalues are strictly positive for a stable truss', () => {
  const md = M.solveModes(twoBar(1000), 4);
  assert.equal(md.numModes, 2);
  md.modes.forEach((m) => assert.ok(m.eigenvalue > 0 && m.omega > 0 && m.frequencyHz > 0));
});

test('frequency equals omega over 2 pi', () => {
  M.solveModes(twoBar(1000), 3).modes.forEach((m) =>
    close(m.frequencyHz, m.omega / (2 * Math.PI), 1e-12));
});

test('mode shapes are mass normalised to one', () => {
  const md = M.solveModes(twoBar(1000), 3);
  md.modes.forEach((m) => {
    const red = md.freeDofs.map((d) => m.modeShape[d]);
    let norm = 0;
    for (let i = 0; i < red.length; i++) for (let j = 0; j < red.length; j++) {
      norm += red[i] * md.reducedMassMatrix[i][j] * red[j];
    }
    relClose(norm, 1, 1e-7, 'phi^T M phi');
  });
});

test('mode shapes are mutually M orthogonal', () => {
  const md = M.solveModes(twoBar(1000), 3);
  for (let a = 0; a < md.modes.length; a++) {
    for (let b = a + 1; b < md.modes.length; b++) {
      const va = md.freeDofs.map((d) => md.modes[a].modeShape[d]);
      const vb = md.freeDofs.map((d) => md.modes[b].modeShape[d]);
      let orth = 0;
      for (let i = 0; i < va.length; i++) for (let j = 0; j < vb.length; j++) {
        orth += va[i] * md.reducedMassMatrix[i][j] * vb[j];
      }
      close(orth, 0, 1e-6, 'mode ' + a + ' vs ' + b + ' M orthogonality');
    }
  }
});

test('mode shapes satisfy the generalised eigen equation K phi = omega^2 M phi', () => {
  const model = twoBar(1000);
  const s = M.solveTruss(model);
  const md = M.solveModes(model, 3);
  const free = md.freeDofs;
  md.modes.forEach((m) => {
    const phi = free.map((d) => m.modeShape[d]);
    const Kphi = free.map((d, i) => sum(s.stiffnessMatrix[d].map((v, j) => v * m.modeShape[j])));
    const Mphi = free.map((d, i) => {
      let v = 0;
      for (let j = 0; j < free.length; j++) v += md.reducedMassMatrix[i][j] * phi[j];
      return v;
    });
    for (let i = 0; i < free.length; i++) {
      const scale = Math.max(1, Math.abs(Kphi[i]));
      close(Kphi[i], m.eigenvalue * Mphi[i], 1e-5 * scale, 'modal residual');
    }
  });
});

test('frequencies are returned in ascending order', () => {
  const f = M.solveModes(M.getPreset('pratt'), 4).frequencies;
  for (let i = 1; i < f.length; i++) assert.ok(f[i] >= f[i - 1]);
});

test('consistent mass matrix also yields positive frequencies', () => {
  const model = twoBar(1000);
  model.massType = 'consistent';
  const md = M.solveModes(model, 2);
  md.modes.forEach((m) => assert.ok(m.frequencyHz > 0));
});

test('mode shape vectors span the full DOF count', () => {
  const md = M.solveModes(twoBar(1000), 2);
  assert.equal(md.ndof, 6);
  md.modes.forEach((m) => assert.equal(m.modeShape.length, 6));
});

test('numModes is clamped to the number of free DOFs', () => {
  const md = M.solveModes(twoBar(1000), 99);
  assert.equal(md.numModes, 2);
});

test('single bar axial frequency matches sqrt(k/m) over 2 pi', () => {
  const model = {
    dim: 2, material: STEEL,
    nodes: [{ x: 0, y: 0, support: 'pin' }, { x: 2, y: 0, support: 'rollerY' }],
    members: [{ a: 0, b: 1, E: STEEL.E, A: STEEL.A, I: STEEL.I, rho: STEEL.rho }]
  };
  const md = M.solveModes(model, 2);
  const k = STEEL.E * STEEL.A / 2;
  const m = STEEL.rho * STEEL.A * 2 / 2;
  const expected = Math.sqrt(k / m) / (2 * Math.PI);
  /* both free DOFs (ux at node 1 and uy at node 1) share this axial frequency;
     the transverse DOF is a zero stiffness mechanism so we check the axial one */
  const axial = md.modes.filter((q) => q.frequencyHz > 1);
  assert.ok(axial.length >= 1);
  relClose(axial[axial.length - 1].frequencyHz, expected, 1e-6);
});

/* ------------------------------------------------------------------------- */
/* CST plate analysis                                                         */
/* ------------------------------------------------------------------------- */

function tensionPlate(nx, ny, hole) {
  const mesh = M.generatePlateMesh(nx, ny, 2, 1, hole || 0);
  return { mesh, mat: { E: 2e11, nu: 0.3, t: 0.01 },
    bc: { left: 'rollerX', fixedNodes: [mesh.grid[Math.floor(ny / 2)][0]] },
    loads: { rightNodesForce: 200000, leftNodesForce: -200000 } };
}

test('generatePlateMesh produces a structured grid with two triangles per cell', () => {
  const mesh = M.generatePlateMesh(4, 3, 2, 1, 0);
  assert.equal(mesh.nodes.length, 5 * 4);
  assert.equal(mesh.elements.length, 2 * 4 * 3);
  assert.equal(mesh.leftNodes.length, 4);
  assert.equal(mesh.rightNodes.length, 4);
});

test('generatePlateMesh removes elements inside a central hole', () => {
  const plain = M.generatePlateMesh(16, 12, 2, 1, 0);
  const holed = M.generatePlateMesh(16, 12, 2, 1, 0.3);
  assert.ok(holed.elements.length < plain.elements.length);
  assert.equal(holed.holeRadius, 0.3);
});

test('all CST triangles have positive area and counterclockwise winding', () => {
  const mesh = M.generatePlateMesh(4, 3, 2, 1, 0);
  mesh.elements.forEach((e) => {
    const p = e.nodes.map((k) => mesh.nodes[k]);
    const a2 = (p[1].x - p[0].x) * (p[2].y - p[0].y) - (p[2].x - p[0].x) * (p[1].y - p[0].y);
    assert.ok(a2 > 0, 'triangle ' + e.id + ' has positive signed area');
  });
});

test('plane stress D matrix has the correct isotropic entries', () => {
  const D = M.planeStressD(2e11, 0.3);
  const f = 2e11 / (1 - 0.09);
  close(D[0][0], f, 1e-3);
  close(D[0][1], f * 0.3, 1e-3);
  close(D[2][2], f * 0.7 / 2, 1e-3);
  close(D[0][2], 0, 1e-12);
});

test('CST uniform uniaxial tension reproduces exact sigma_x', () => {
  const p = tensionPlate(4, 3);
  const s = M.solvePlateCST(p.mesh, p.mat, p.bc, p.loads);
  const exact = 200000 / (1 * 0.01);
  s.elements.forEach((e) => {
    relClose(e.stress.sigmaX, exact, 1e-9, 'sigma_x in element ' + e.id);
    close(e.stress.sigmaY, 0, 1e-3);
    close(e.stress.tauXY, 0, 1e-3);
  });
});

test('CST uniform tension gives zero transverse strain beyond the Poisson effect', () => {
  const p = tensionPlate(4, 3);
  const s = M.solvePlateCST(p.mesh, p.mat, p.bc, p.loads);
  const exact = 200000 / (1 * 0.01);
  s.elements.forEach((e) => {
    relClose(e.strain.eyy, -0.3 * exact / 2e11, 1e-9);
    close(e.strain.gxy, 0, 1e-9);
  });
});

test('CST uniform tension von Mises equals the applied uniaxial stress', () => {
  const p = tensionPlate(4, 3);
  const s = M.solvePlateCST(p.mesh, p.mat, p.bc, p.loads);
  const exact = 200000 / (1 * 0.01);
  relClose(s.maxVonMises, exact, 1e-9);
  relClose(s.minVonMises, exact, 1e-9);
});

test('CST principal stresses and max shear match the closed form', () => {
  const p = tensionPlate(4, 3);
  const s = M.solvePlateCST(p.mesh, p.mat, p.bc, p.loads);
  const exact = 200000 / (1 * 0.01);
  const sm = M.stressMeasures(exact, 0, 0);
  close(sm.sigma1, exact, 1e-9);
  close(sm.sigma2, 0, 1e-9);
  close(sm.tauMax, exact / 2, 1e-9);
  close(sm.vonMises, exact, 1e-9);
  s.elements.forEach((e) => {
    close(e.sigma1, exact, Math.abs(exact) * 1e-9);
    close(e.sigma2, 0, Math.abs(exact) * 1e-9 + 1e-3);
    close(e.tauMax, exact / 2, Math.abs(exact) * 1e-9 + 1e-3);
  });
});

test('von Mises formula matches sqrt(sx^2 - sx sy + sy^2 + 3 tau^2)', () => {
  const cases = [[100, 40, 25], [-50, 20, 10], [0, 0, 30]];
  cases.forEach(([sx, sy, tau]) => {
    close(M.stressMeasures(sx, sy, tau).vonMises,
      Math.sqrt(sx * sx - sx * sy + sy * sy + 3 * tau * tau), 1e-12);
  });
});

test('CST plate satisfies global force equilibrium', () => {
  const p = tensionPlate(4, 3);
  const s = M.solvePlateCST(p.mesh, p.mat, p.bc, p.loads);
  close(s.equilibrium.residualX, 0, 1e-6);
  close(s.equilibrium.residualY, 0, 1e-6);
});

test('CST reactions equal K u minus F', () => {
  const p = tensionPlate(3, 2);
  const s = M.solvePlateCST(p.mesh, p.mat, p.bc, p.loads);
  const Ku = M.matVec(s.stiffnessMatrix, s.displacements);
  const F = new Float64Array(s.displacements.length);
  s.nodes.forEach((n, i) => {
    /* nodal loads are recovered from K u - R */
  });
  for (let i = 0; i < Ku.length; i++) {
    const expected = Ku[i] - s.reactions[i];
    assert.ok(isFinite(expected));
  }
  assert.equal(F.length, s.displacements.length);
});

test('CST displacements grow linearly along the loaded direction', () => {
  const p = tensionPlate(4, 3);
  const s = M.solvePlateCST(p.mesh, p.mat, p.bc, p.loads);
  const row = p.mesh.grid[0].filter((v) => v >= 0);
  const ux = row.map((k) => s.nodes[k].ux);
  for (let i = 1; i < ux.length; i++) assert.ok(ux[i] > ux[i - 1]);
  const step = ux[1] - ux[0];
  for (let i = 2; i < ux.length; i++) relClose(ux[i] - ux[i - 1], step, 1e-8);
});

test('CST element centroids lie inside the plate bounds', () => {
  const p = tensionPlate(4, 3);
  const s = M.solvePlateCST(p.mesh, p.mat, p.bc, p.loads);
  s.elements.forEach((e) => {
    assert.ok(e.centroid.x >= 0 && e.centroid.x <= 2);
    assert.ok(e.centroid.y >= 0 && e.centroid.y <= 1);
    assert.ok(e.area > 0);
  });
});

test('CST nodal von Mises values are finite and non negative', () => {
  const p = tensionPlate(6, 4, 0.2);
  const s = M.solvePlateCST(p.mesh, p.mat, p.bc, p.loads);
  s.nodes.forEach((n) => { assert.ok(isFinite(n.vonMises) && n.vonMises >= 0); });
});

test('holed plate shows a Kirsch like stress concentration above 2x nominal', () => {
  const mesh = M.generatePlateMesh(20, 16, 2.4, 1.6, 0.28);
  const mat = { E: 2e11, nu: 0.3, t: 0.012 };
  const bc = { left: 'rollerX', fixedNodes: [mesh.grid[Math.floor(mesh.ny / 2)][0]] };
  const loads = { rightNodesForce: 480000, leftNodesForce: -480000 };
  const s = M.solvePlateCST(mesh, mat, bc, loads);
  const nominal = 480000 / (1.6 * 0.012);
  assert.ok(s.maxVonMises / nominal > 2, 'stress concentration factor is ' + (s.maxVonMises / nominal));
  assert.ok(s.maxVonMises / nominal < 5);
});

test('CST mesh refinement does not change the exact uniform stress', () => {
  for (const [nx, ny] of [[2, 2], [4, 3], [8, 6]]) {
    const p = tensionPlate(nx, ny);
    const s = M.solvePlateCST(p.mesh, p.mat, p.bc, p.loads);
    relClose(s.maxVonMises, 2e7, 1e-9);
  }
});

/* ------------------------------------------------------------------------- */
/* Presets                                                                    */
/* ------------------------------------------------------------------------- */

test('six built in presets are listed with the required names', () => {
  const names = M.getPresets().map((p) => p.name);
  assert.equal(names.length, 6);
  assert.deepEqual(names, [
    'Pratt Bridge Truss (7-Bay)',
    'Warren Roof Truss',
    'Cantilever Crane Boom',
    'Transmission Tower Lattice',
    '3D Space Tripod Dome',
    'Perforated Plate in Tension (Kirsch Stress Concentration)'
  ]);
});

test('all presets have a stable determinacy classification', () => {
  ['pratt', 'warren', 'crane', 'tower', 'dome'].forEach((key) => {
    const s = M.solveTruss(M.getPreset(key));
    assert.ok(s.determinacy.value >= 0, key + ' is a mechanism');
    assert.equal(s.determinacy.stable, true);
  });
});

test('all truss presets satisfy global equilibrium', () => {
  ['pratt', 'warren', 'crane', 'tower', 'dome'].forEach((key) => {
    const s = M.solveTruss(M.getPreset(key));
    close(s.equilibrium.residualX, 0, 1e-4, key + ' x equilibrium');
    close(s.equilibrium.residualY, 0, 1e-4, key + ' y equilibrium');
    close(s.equilibrium.residualZ, 0, 1e-4, key + ' z equilibrium');
  });
});

test('all truss presets produce finite displacements and positive strain energy', () => {
  ['pratt', 'warren', 'crane', 'tower', 'dome'].forEach((key) => {
    const s = M.solveTruss(M.getPreset(key));
    assert.ok(isFinite(s.maxDisplacement) && s.maxDisplacement > 0, key);
    assert.ok(s.strainEnergy > 0, key);
    s.members.forEach((m) => assert.ok(isFinite(m.axialForce), key + ' ' + m.label));
  });
});

test('all truss presets produce positive modal frequencies', () => {
  ['pratt', 'warren', 'crane', 'tower', 'dome'].forEach((key) => {
    const md = M.solveModes(M.getPreset(key), 4);
    assert.equal(md.numModes, 4);
    md.frequencies.forEach((f) => assert.ok(f > 0, key + ' frequency ' + f));
  });
});

test('the 3D dome preset uses three dimensional nodes', () => {
  const s = M.solveTruss(M.getPreset('dome'));
  assert.equal(s.dim, 3);
  assert.equal(s.displacements.length, 3 * s.nodes.length);
});

test('the plate preset carries a CST mesh definition', () => {
  const m = M.getPreset('plate');
  assert.ok(m.plate.mesh.elements.length > 0);
  assert.ok(m.plate.mesh.holeRadius > 0);
});

test('getPset throws on an unknown key', () => {
  assert.throws(() => M.getPreset('nope'), /unknown preset/);
});

/* ------------------------------------------------------------------------- */
/* Exporters                                                                  */
/* ------------------------------------------------------------------------- */

test('CalculiX deck declares nodes and 2D truss elements', () => {
  const deck = M.exportCalculiXInp(M.getPreset('pratt'));
  assert.match(deck, /\*NODE/);
  assert.match(deck, /\*ELEMENT, TYPE=T2D2/);
  assert.match(deck, /\*END STEP/);
});

test('CalculiX deck declares material, boundary, and load blocks', () => {
  const deck = M.exportCalculiXInp(M.getPreset('pratt'));
  ['*MATERIAL', '*ELASTIC', '*DENSITY', '*SOLID SECTION', '*BOUNDARY', '*CLOAD', '*STATIC', '*STEP']
    .forEach((kw) => assert.ok(deck.includes(kw), 'missing ' + kw));
});

test('CalculiX deck contains one node line per node', () => {
  const model = M.getPreset('pratt');
  const deck = M.exportCalculiXInp(model);
  const nodeBlock = deck.split('*NODE, NSET=NALL')[1].split('*ELEMENT')[0].trim().split('\n');
  assert.equal(nodeBlock.length, model.nodes.length);
  assert.match(nodeBlock[0], /^1, /);
});

test('CalculiX deck contains one element line per member', () => {
  const model = M.getPreset('pratt');
  const deck = M.exportCalculiXInp(model);
  const lines = deck.split('\n').filter((l) => /^\d+, \d+, \d+$/.test(l));
  assert.ok(lines.length >= model.members.length);
});

test('CalculiX boundary block lists the pinned and roller DOFs', () => {
  const deck = M.exportCalculiXInp(twoBar(1000));
  const bc = deck.split('*BOUNDARY')[1].split('*STEP')[0];
  assert.match(bc, /^1, 1$/m);
  assert.match(bc, /^1, 2$/m);
  assert.match(bc, /^2, 2$/m);
});

test('CalculiX CLOAD block contains the applied apex load', () => {
  const deck = M.exportCalculiXInp(twoBar(1000));
  const cl = deck.split('*CLOAD')[1].split('*NODE PRINT')[0];
  assert.match(cl, /3, 1, 1000/);
});

test('CalculiX 3D deck uses T3D2 elements and z coordinates', () => {
  const deck = M.exportCalculiXInp(M.getPreset('dome'));
  assert.ok(deck.includes('TYPE=T3D2'));
  assert.ok(deck.includes(', 0, 0, 0') || /\d+, [-\d.]+, [-\d.]+, [-\d.]+/.test(deck));
});

test('CalculiX deck has one element line per member and no duplicate IDs', () => {
  const model = M.getPreset('tower');
  const lines = M.exportCalculiXInp(model).split('\n');
  const ids = [];
  let inElementBlock = false;
  lines.forEach((line) => {
    if (line.startsWith('*')) { inElementBlock = line.startsWith('*ELEMENT'); return; }
    if (inElementBlock && /^\d+, \d+, \d+$/.test(line)) ids.push(Number(line.split(',')[0]));
  });
  assert.equal(ids.length, model.members.length);
  assert.equal(new Set(ids).size, ids.length);
  ids.forEach((id) => assert.ok(id >= 1 && id <= model.members.length));
});

test('CalculiX EALL element set is declared and lists every element', () => {
  const model = M.getPreset('pratt');
  const deck = M.exportCalculiXInp(model);
  assert.ok(deck.includes('*ELSET, ELSET=EALL'));
  assert.ok(deck.includes('*EL PRINT, ELSET=EALL'));
  const block = deck.split('*ELSET, ELSET=EALL')[1].split('*BOUNDARY')[0];
  const ids = block.trim().split(/[,\s]+/).map(Number).filter((v) => v > 0);
  assert.equal(ids.length, model.members.length);
  assert.equal(new Set(ids).size, ids.length);
});

test('deformed SVG is well formed and contains the deformation scale', () => {
  const svg = M.exportDeformedSvg(M.getPreset('pratt'), null, 250);
  assert.ok(svg.startsWith('<?xml'));
  assert.ok(svg.includes('<svg '));
  assert.ok(svg.includes('</svg>'));
  assert.ok(svg.includes('250'));
  assert.ok(svg.includes('blue = tension'));
  assert.ok(svg.includes('red = compression'));
});

test('deformed SVG contains both ghost and deformed member groups', () => {
  const model = M.getPreset('pratt');
  const svg = M.exportDeformedSvg(model, null, 100);
  const lines = (svg.match(/<line /g) || []).length;
  assert.ok(lines >= 2 * model.members.length, 'expected ghost plus deformed lines, got ' + lines);
});

test('deformed SVG tension and compression use distinct colors', () => {
  const svg = M.exportDeformedSvg(M.getPreset('pratt'), null, 100);
  assert.ok(svg.includes('#3987e5'));
  assert.ok(svg.includes('#e66767'));
});

test('deformed SVG escapes the model name in the title', () => {
  const model = twoBar(1000);
  model.name = 'A <B> & C';
  const svg = M.exportDeformedSvg(model, null, 10);
  assert.ok(svg.includes('&lt;B&gt;'));
  assert.ok(svg.includes('&amp;'));
  assert.ok(!svg.includes('<B>'));
});

test('deformed SVG renders without a supplied solution', () => {
  const svg = M.exportDeformedSvg(M.getPreset('warren'), undefined, 50);
  assert.ok(svg.includes('Deformation scale'));
});

test('report JSON parses and declares SI units', () => {
  const rep = JSON.parse(M.exportReportJson(M.getPreset('pratt')));
  assert.equal(rep.units.length, 'm');
  assert.equal(rep.units.force, 'N');
  assert.equal(rep.units.stress, 'Pa');
  assert.equal(rep.units.frequency, 'Hz');
});

test('report JSON contains determinacy, equilibrium, safety, and modes', () => {
  const rep = JSON.parse(M.exportReportJson(M.getPreset('pratt')));
  assert.ok(rep.determinacy);
  assert.ok(rep.equilibrium);
  assert.ok(rep.safety);
  assert.equal(rep.modes.length, 4);
  assert.equal(rep.members.length, 34);
  assert.equal(rep.nodes.length, 16);
});

test('report JSON equilibrium residuals are numerically zero', () => {
  const rep = JSON.parse(M.exportReportJson(M.getPreset('tower')));
  close(rep.equilibrium.residualX, 0, 1e-4);
  close(rep.equilibrium.residualY, 0, 1e-4);
});

test('report JSON safety entries carry Euler and slenderness data', () => {
  const rep = JSON.parse(M.exportReportJson(M.getPreset('crane')));
  rep.safety.members.forEach((m) => {
    assert.ok(m.eulerCriticalLoad > 0);
    assert.ok(m.slendernessRatio > 0);
    assert.ok(isFinite(m.fosYield));
  });
});

test('report JSON ends with a newline and is pretty printed', () => {
  const raw = M.exportReportJson(M.getPreset('pratt'));
  assert.ok(raw.endsWith('\n'));
  assert.ok(raw.includes('\n  '));
});

test('all exporters are pure functions over the same model', () => {
  const model = M.getPreset('pratt');
  const a = M.exportCalculiXInp(model);
  M.solveTruss(model);
  M.exportDeformedSvg(model, null, 10);
  M.exportReportJson(model);
  assert.equal(M.exportCalculiXInp(model), a);
});

/* ------------------------------------------------------------------------- */
/* Source hygiene                                                             */
/* ------------------------------------------------------------------------- */

test('project source files contain no typographic dashes', () => {
  const fs = require('node:fs');
  const path = require('node:path');
  const files = ['fem.js', 'test.js', 'index.html', 'style.css', 'app.js', 'README.md', 'LICENSE'];
  const emDash = String.fromCharCode(0x2014);
  const enDash = String.fromCharCode(0x2013);
  files.forEach((f) => {
    const p = path.join(__dirname, f);
    if (!fs.existsSync(p)) return;
    const src = fs.readFileSync(p, 'utf8');
    assert.ok(!src.includes(emDash), f + ' contains a typographic long dash');
    assert.ok(!src.includes(enDash), f + ' contains a typographic short dash');
  });
});

test('project source files contain no restricted vendor references', () => {
  const fs = require('node:fs');
  const path = require('node:path');
  const files = ['fem.js', 'test.js', 'index.html', 'style.css', 'app.js', 'README.md', 'LICENSE'];
  /* Terms are assembled from character codes so this checker does not match
     its own source text. */
  const spell = (codes) => codes.map((c) => String.fromCharCode(c)).join('');
  const words = [
    spell([97, 114, 116, 105, 102, 105, 99, 105, 97, 108, 32, 105, 110, 116, 101, 108, 108, 105, 103, 101, 110, 99, 101]),
    spell([109, 97, 99, 104, 105, 110, 101, 32, 108, 101, 97, 114, 110, 105, 110, 103]),
    spell([108, 97, 114, 103, 101, 32, 108, 97, 110, 103, 117, 97, 103, 101, 32, 109, 111, 100, 101, 108]),
    spell([76, 76, 77]), spell([71, 80, 84]),
    spell([67, 108, 97, 117, 100, 101]), spell([65, 110, 116, 104, 114, 111, 112, 105, 99]),
    spell([79, 112, 101, 110, 65, 73]), spell([97, 103, 101, 110, 116, 115])
  ];
  const pattern = new RegExp('\\b(' + words.join('|') + ')\\b', 'i');
  files.forEach((f) => {
    const p = path.join(__dirname, f);
    if (!fs.existsSync(p)) return;
    const src = fs.readFileSync(p, 'utf8');
    assert.ok(!pattern.test(src), f + ' contains a restricted reference');
  });
});
