'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const M = require('./fem.js');
const F = require('./field.js');
const close = (a, b, tolerance = 1e-9) => assert.ok(Math.abs(a - b) <= tolerance, `${a} differs from ${b}`);
const bar = () => ({ dim: 2, nodes: [{ x: 0, y: 0, support: 'pin' }, { x: 2, y: 0, support: 'rollerY', load: [1000, 0] }], members: [{ a: 0, b: 1, E: 2e11, A: 0.01, I: 1e-6 }] });
const patch = () => M.generatePlateMesh(4, 3, 2, 1, 0);

for (const [key, force] of [['pratt', -360000], ['warren', -120000], ['crane', -130000], ['tower', -24000], ['dome', -80000]]) {
  test(`${key} preset applies its declared load exactly once`, () => {
    const model = M.getPreset(key), solution = M.solveTruss(model);
    close(solution.equilibrium.sumFy, force, 1e-9);
    assert.equal(model.loads, undefined);
    assert.equal(solution.verification.passed, true);
    assert.ok(solution.verification.relativeResidual < 1e-10);
  });
}
test('short nodal load arrays retain additional top level components', () => {
  const model = bar(); model.nodes[1].load = [1000]; model.loads = [{ node: 1, f: [2000, 0] }];
  const solution = M.solveTruss(model);
  close(solution.forceVector[2], 3000); close(solution.forceVector[3], 0);
});
test('automatic dimension detection recognizes coordinate arrays and aliases', () => {
  for (const nodes of [[[0, 0, 0], [1, 1, 1]], [{ coords: [0, 0, 0] }, { xyz: [1, 1, 1] }]]) {
    assert.equal(M.normalizeModel({ dim: 'auto', nodes, members: [{ a: 0, b: 1 }] }).dim, 3);
  }
});
test('invalid model data produces descriptive errors', () => {
  for (const mutate of [m => { m.nodes[0].x = NaN; }, m => { m.members[0].A = 0; }, m => { m.members[0].E = -1; }, m => { m.members[0].rho = Infinity; }, m => { m.nodes[0].support = 'unknown'; }, m => { m.members[0].a = 0.5; }, m => { m.loads = [{ node: 12, f: [1, 0] }]; }]) {
    const model = bar(); mutate(model); assert.throws(() => M.solveTruss(model), /MeshLab:/);
  }
});
test('empty structures and incompatible linear dimensions are rejected', () => {
  assert.throws(() => M.solveTruss({ nodes: [], members: [] }), /needs two nodes/);
  assert.throws(() => M.solveLinear([[1]], [1, 2]), /dimensions/);
  assert.deepEqual(Array.from(M.solveLinear([], [])), []);
});
test('a fully constrained model returns loads as support reactions', () => {
  const model = bar(); model.nodes[1].support = 'pin';
  const s = M.solveTruss(model); close(s.maxDisplacement, 0); close(s.reactions[2], -1000); assert.equal(s.verification.passed, true);
});
test('maximum translation is the vector norm rather than a component', () => {
  const model = M.getPreset('dome'), s = M.solveTruss(model);
  const expected = Math.max(...s.nodes.map(n => Math.hypot(...n.displacement)));
  close(s.maxDisplacement, expected, 1e-15);
  assert.ok(s.maxDisplacement >= s.maxComponentDisplacement);
});
test('consistent mass modes satisfy the full generalized eigenproblem', () => {
  const model = { ...M.getPreset('pratt'), massType: 'consistent' };
  const s = M.solveTruss(model), modes = M.solveModes(model, 4);
  for (const mode of modes.modes) {
    assert.ok(mode.relativeResidual < 1e-8);
    const Kphi = M.matVec(s.stiffnessMatrix, mode.modeShape), Mphi = M.matVec(modes.massMatrix, mode.modeShape);
    const r = modes.freeDofs.map(i => Kphi[i] - mode.eigenvalue * Mphi[i]);
    assert.ok(Math.hypot(...r) / Math.hypot(...modes.freeDofs.map(i => Kphi[i])) < 1e-8);
  }
  for (const a of modes.modes) for (const b of modes.modes) {
    close(M.dot(a.modeShape, M.matVec(modes.massMatrix, b.modeShape)), a === b ? 1 : 0, 1e-8);
  }
});
test('consistent mass differs from lumped mass when two free nodes are coupled', () => {
  const model = M.getPreset('pratt');
  assert.ok(Math.abs(M.solveModes(model, 1).frequencies[0] - M.solveModes({ ...model, massType: 'consistent' }, 1).frequencies[0]) > 0.01);
  assert.throws(() => M.solveModes({ ...model, massType: 'unknown' }, 1), /mass formulation/);
});
test('small scale Jacobi systems preserve relative eigenvalue accuracy', () => {
  const eig = M.jacobiEigen([[4e-15, 1e-15], [1e-15, 3e-15]]);
  close(eig[0].value, (7 - Math.sqrt(5)) * 0.5e-15, 1e-26);
});
test('plate nodalArray applies the same load as explicit nodal entries', () => {
  const mesh = patch(), ids = mesh.rightNodes;
  const a = M.solvePlateCST(mesh, {}, {}, { nodalArray: { nodes: ids, fx: 1200, fy: 0 } });
  const b = M.solvePlateCST(mesh, {}, {}, { nodal: ids.map(node => ({ node, fx: 1200 })) });
  assert.deepEqual(a.forceVector, b.forceVector); a.displacements.forEach((v, i) => close(v, b.displacements[i], 1e-15));
});
test('plate solver returns an independently checkable force vector and energy', () => {
  const s = M.solvePlateCST(patch(), {}, {}, { rightNodesForce: 200000 });
  const Ku = M.matVec(s.stiffnessMatrix, s.displacements);
  Ku.forEach((v, i) => close(s.reactions[i], v - s.forceVector[i], 1e-6));
  close(s.strainEnergy * 2, s.externalWork, 1e-8); assert.equal(s.verification.passed, true);
});
test('plate PCG agrees with independent Gaussian elimination', () => {
  const cfg = M.getPreset('plate').plate;
  const a = M.solvePlateCST(cfg.mesh, cfg.material, cfg.bc, cfg.loads);
  const b = M.solvePlateCST(cfg.mesh, cfg.material, cfg.bc, cfg.loads, { method: 'pcg', tolerance: 1e-12 });
  assert.ok(b.diagnostics.converged); assert.ok(b.diagnostics.relativeResidual < 1e-12);
  a.displacements.forEach((v, i) => close(v, b.displacements[i], 1e-12));
  close(a.maxVonMises, b.maxVonMises, 0.01); assert.equal(b.verification.passed, true);
});
test('plate PCG reports a failed solve when its iteration budget is insufficient', () => {
  assert.throws(() => M.solvePlateCST(patch(), {}, {}, { rightNodesForce: 1000 }, { method: 'pcg', maxIterations: 1 }), /did not converge/);
});
test('zero load plates have zero stress and exact verification', () => {
  const s = M.solvePlateCST(patch(), {}, {}, {}, { method: 'pcg' });
  close(s.maxVonMises, 0); close(s.maxDisplacement, 0); assert.equal(s.verification.passed, true);
});
test('plate dimensions, material, and degenerate elements are validated', () => {
  for (const args of [[4, 3, 0, 1, 0], [4, 3, 2, 1, 0.5], [4, 3, 2, 1, -0.1], [61, 3, 2, 1, 0], [NaN, 3, 2, 1, 0]]) assert.throws(() => M.generatePlateMesh(...args), /MeshLab:/);
  for (const mat of [{ t: 0 }, { E: -1 }, { nu: 0.5 }, { nu: -1 }, { nu: NaN }]) assert.throws(() => M.solvePlateCST(patch(), mat), /MeshLab:/);
  const mesh = patch(); mesh.elements[0] = { nodes: [0, 0, 1] }; assert.throws(() => M.solvePlateCST(mesh), /zero area/);
});
test('custom CST meshes do not require structured grid metadata', () => {
  const mesh = { nodes: [{ id: 0, x: 0, y: 0 }, { id: 1, x: 1, y: 0 }, { id: 2, x: 0, y: 1 }], elements: [{ nodes: [0, 1, 2] }], leftNodes: [0, 2] };
  const s = M.solvePlateCST(mesh, {}, undefined, { nodal: [{ node: 1, fx: 1000 }] }); assert.ok(s.maxDisplacement > 0); assert.equal(s.verification.passed, true);
});
test('hole meshes contain no unconnected nodes or cavity crossing edges', () => {
  const mesh = M.generatePlateMesh(9, 7, 2, 1, 0.07), used = new Set(mesh.elements.flatMap(e => e.nodes));
  assert.equal(used.size, mesh.nodes.length); assert.ok(mesh.elements.length < 2 * 9 * 7);
  for (const e of mesh.elements) for (let i = 0; i < 3; i++) {
    const a = mesh.nodes[e.nodes[i]], b = mesh.nodes[e.nodes[(i + 1) % 3]], dx = b.x - a.x, dy = b.y - a.y;
    const u = Math.max(0, Math.min(1, ((1 - a.x) * dx + (0.5 - a.y) * dy) / (dx * dx + dy * dy)));
    assert.ok(Math.hypot(a.x + u * dx - 1, a.y + u * dy - 0.5) >= 0.07 * (1 - 1e-12));
  }
});
test('3D SVG export preserves a member that exists only along z', () => {
  const model = { dim: 3, nodes: [{ x: 0, y: 0, z: 0, support: 'pin' }, { x: 0, y: 0, z: 1, support: 'pin' }], members: [{ a: 0, b: 1 }] };
  const svg = M.exportDeformedSvg(model), ghost = svg.split('stroke="#5b7899"')[1].split('</g>')[0];
  const line = /x1="([^"]+)" y1="([^"]+)" x2="([^"]+)" y2="([^"]+)"/.exec(ghost);
  assert.notEqual(line[1], line[3]); assert.notEqual(line[2], line[4]);
});

test('PCG solves a known symmetric positive definite system', () => {
  const r = F.solvePCG([[2, 1], [1, 3]], [3, 5]); close(r.solution[0], 0.8); close(r.solution[1], 1.4); assert.ok(r.converged); assert.equal(r.nonzeros, 4);
});
test('PCG handles zero and empty right hand sides', () => {
  const r = F.solvePCG([[2, 0], [0, 3]], [0, 0]); assert.equal(r.iterations, 0); assert.equal(r.relativeResidual, 0);
  assert.equal(F.solvePCG([], []).converged, true);
});
test('PCG rejects nonsymmetric, indefinite and malformed systems', () => {
  assert.throws(() => F.solvePCG([[2, 1], [0, 3]], [1, 1]), /symmetric/);
  assert.throws(() => F.solvePCG([[1, 2], [2, 1]], [1, -1]), /non positive/);
  assert.throws(() => F.solvePCG([[1]], [NaN]), /finite/);
  assert.throws(() => F.solvePCG([[1]], [1, 2]), /dimension/);
  assert.throws(() => F.solvePCG([[1]], [1], { tolerance: 0 }), /tolerance/);
});
test('PCG reports actual residual and finite iteration exhaustion', () => {
  const A = [[4, -1, 0], [-1, 4, -1], [0, -1, 3]], b = [1, 2, 5];
  const r = F.solvePCG(A, b, { maxIterations: 1 }); assert.equal(r.converged, false); assert.equal(r.iterations, 1);
  close(r.relativeResidual, Math.hypot(...F.multiply(F.toCSR(A), r.solution).map((v, i) => b[i] - v)) / Math.hypot(...b), 1e-12);
});
test('Delaunay square with center has four triangles and correct hull', () => {
  const mesh = F.triangulate([[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]]);
  assert.equal(mesh.elements.length, 4); assert.equal(mesh.boundaryNodes.length, 4); close(mesh.quality.area, 1); close(mesh.quality.minimumAngle, 45);
});
test('Delaunay removes exact duplicate points and rejects collinear data', () => {
  const mesh = F.triangulate([[0, 0], [1, 0], [0, 1], [0, 1]]); assert.equal(mesh.duplicatesRemoved, 1); assert.equal(mesh.nodes.length, 3);
  assert.throws(() => F.triangulate([[0, 0], [1, 0], [2, 0]]), /collinear/);
  assert.throws(() => F.triangulate([[0, 0], [1, 0], [NaN, 1]]), /finite/);
});
test('Delaunay normalization supports microscopic geometry', () => {
  const mesh = F.triangulate([[0, 0], [1e-9, 0], [1e-9, 1e-9], [0, 1e-9], [0.5e-9, 0.5e-9]]);
  assert.equal(mesh.elements.length, 4); close(mesh.quality.area, 1e-18, 1e-30);
});
test('every irregular mesh triangle has counterclockwise winding and an empty circumcircle', () => {
  const mesh = F.generateFieldMesh(6, 4, 2, 1);
  for (const e of mesh.elements) {
    const [a, b, c] = e.nodes.map(id => mesh.nodes[id]);
    const cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); assert.ok(cross > 0);
    const d = 2 * cross, bx = b.x - a.x, by = b.y - a.y, cx = c.x - a.x, cy = c.y - a.y;
    const x = a.x + ((bx * bx + by * by) * cy - (cx * cx + cy * cy) * by) / d;
    const y = a.y + ((cx * cx + cy * cy) * bx - (bx * bx + by * by) * cx) / d;
    const radius2 = (x - a.x) ** 2 + (y - a.y) ** 2;
    mesh.nodes.forEach((p, id) => { if (!e.nodes.includes(id)) assert.ok((p.x - x) ** 2 + (p.y - y) ** 2 >= radius2 - 1e-10); });
  }
});
test('default field mesh has verified geometry and no isolated nodes', () => {
  const mesh = F.generateFieldMesh(); assert.equal(mesh.nodes.length, 165); assert.equal(mesh.elements.length, 280); close(mesh.quality.area, 3.84); assert.ok(mesh.quality.minimumAngle > 33);
  assert.equal(new Set(mesh.elements.flatMap(e => e.nodes)).size, mesh.nodes.length);
  assert.throws(() => F.generateFieldMesh(0, 4, 2, 1), /invalid/);
});
test('Poisson FEM reproduces an affine nonzero Dirichlet patch and constant flux', () => {
  const mesh = F.generateFieldMesh(8, 6, 2, 1), boundary = n => n.x === 0 ? 100 : n.x === 2 ? 0 : undefined;
  const s = F.solvePoisson(mesh, { conductivity: 2, boundary });
  mesh.nodes.forEach((n, i) => close(s.values[i], 100 * (1 - n.x / 2), 1e-7));
  s.elements.forEach(e => { close(e.fluxX, 100, 1e-6); close(e.fluxY, 0, 1e-6); });
  assert.ok(s.diagnostics.converged); close(s.verification.balanceResidual, 0, 1e-6);
});
test('Poisson constant boundary values produce a constant field with zero source', () => {
  const mesh = F.generateFieldMesh(5, 4, 2, 1), s = F.solvePoisson(mesh, { boundary: (n, id) => mesh.boundaryNodes.includes(id) ? 42 : undefined });
  s.values.forEach(v => close(v, 42, 1e-8)); assert.ok(s.energy < 1e-14);
});
test('Poisson positive source is positive under zero exterior boundary conditions', () => {
  const mesh = F.generateFieldMesh(8, 6, 2, 1), s = F.solvePoisson(mesh, { source: 10, boundary: (n, id) => mesh.boundaryNodes.includes(id) ? 0 : undefined });
  assert.ok(s.max > 0); assert.ok(s.min >= -1e-10); close(s.verification.sourceTotal, 20, 1e-10); close(s.verification.boundaryFlux, -20, 1e-8);
});
test('Poisson solution obeys conductivity and source scaling', () => {
  const mesh = F.generateFieldMesh(4, 3, 2, 1), boundary = (n, id) => mesh.boundaryNodes.includes(id) ? 0 : undefined;
  const a = F.solvePoisson(mesh, { source: 1, conductivity: 1, boundary }), b = F.solvePoisson(mesh, { source: 4, conductivity: 2, boundary });
  a.values.forEach((v, i) => close(b.values[i], 2 * v, 1e-10));
});
test('Poisson validates boundary conditions, conductivity and source values', () => {
  const mesh = F.generateFieldMesh(4, 3, 2, 1);
  assert.throws(() => F.solvePoisson(mesh), /prescribed boundary/);
  assert.throws(() => F.solvePoisson(mesh, { conductivity: 0 }), /positive/);
  assert.throws(() => F.solvePoisson(mesh, { source: () => NaN }), /finite/);
  assert.throws(() => F.solvePoisson(mesh, { source: NaN }), /finite/);
  assert.throws(() => F.solvePoisson(mesh, { boundary: () => 0, maxIterations: 0 }), /iteration/);
});
test('disconnected scalar meshes need a prescribed value in every component', () => {
  const mesh = { nodes: [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 0, y: 1 }, { x: 3, y: 0 }, { x: 4, y: 0 }, { x: 3, y: 1 }], elements: [{ nodes: [0, 1, 2] }, { nodes: [3, 4, 5] }] };
  assert.throws(() => F.solvePoisson(mesh, { boundary: (n, id) => id === 0 ? 0 : undefined }), /every mesh component/);
  const s = F.solvePoisson(mesh, { boundary: (n, id) => id === 0 ? 5 : id === 3 ? 10 : undefined });
  s.values.forEach((value, i) => close(value, i < 3 ? 5 : 10));
});
test('zero load PCG plates still reject insufficient supports', () => {
  assert.throws(() => M.solvePlateCST(patch(), {}, { constraints: [{ node: 0, x: true }] }, {}, { method: 'pcg' }), /singular/);
});
test('PCG rejects singular models even with nonzero loads (4x3 rollerX + right load)', () => {
  const plate = M.generatePlateMesh(1.2, 0.9, 4, 3);
  const rightLoad = { edgeLoads: [{ edge: 'right', fx: 1000, fy: 0 }] };
  assert.throws(() => M.solvePlateCST(plate, { E: 210e9, nu: 0.3, thickness: 0.01 }, { left: 'rollerX' }, rightLoad, { method: 'pcg' }), /singular/);
});
test('valid skinny triangles triangulate cleanly instead of failing as collinear', () => {
  const mesh = F.triangulate([[0, 0], [0.5, 0.001], [1, 0]]);
  assert.equal(mesh.elements.length, 1);
  assert.equal(mesh.boundaryNodes.length, 3);
  close(mesh.quality.area, 0.0005, 1e-10);
});
test('skinny triangles near a convex hull edge are preserved without dropping interior points', () => {
  const mesh = F.triangulate([[0, 0], [1, 0], [0, 1], [0.5, 0.00001]]);
  assert.equal(mesh.elements.length, 3);
  assert.deepEqual(mesh.boundaryNodes, [0, 1, 2]);
  close(mesh.quality.area, 0.5, 1e-12);
});

/* ------------------------------------------------------------------ */
/* Domain meshes, MMS convergence, PCG history, nodal flux, VTK        */
/* ------------------------------------------------------------------ */

test('generated domains are connected, oriented and carry the right Euler characteristic', () => {
  for (const [type, chi] of [['l', 1], ['annulus', 0], ['sink', 1]]) {
    const mesh = F.generateDomainMesh({ type, nx: 16, ny: 12, width: 2, height: 2, fins: 4 });
    assert.equal(mesh.chi, chi, type + ' Euler');
    assert.ok(mesh.quality.minimumAngle > 5, type + ' min angle ' + mesh.quality.minimumAngle);
    assert.ok(mesh.bcSets.hot.length && mesh.bcSets.cold.length, type + ' boundary sets');
    for (const e of mesh.elements) {
      const [a, b, c] = e.nodes.map(id => mesh.nodes[id]);
      assert.ok((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) > 0, type + ' CCW');
    }
    const edgeUse = new Map(), adj = mesh.nodes.map(() => []);
    mesh.elements.forEach(e => e.nodes.forEach((u, i) => {
      const v = e.nodes[(i + 1) % 3], k = Math.min(u, v) + ':' + Math.max(u, v);
      edgeUse.set(k, (edgeUse.get(k) || 0) + 1); adj[u].push(v); adj[v].push(u);
    }));
    for (const c of edgeUse.values()) assert.ok(c <= 2, type + ' non-manifold edge');
    const seen = new Set([0]), stack = [0];
    while (stack.length) for (const nxt of adj[stack.pop()]) if (!seen.has(nxt)) { seen.add(nxt); stack.push(nxt); }
    assert.equal(seen.size, mesh.nodes.length, type + ' connected');
    assert.equal(new Set(mesh.elements.flatMap(e => e.nodes)).size, mesh.nodes.length, type + ' no orphans');
  }
});

test('interior edges of generated domain meshes satisfy the local Delaunay criterion', () => {
  const mesh = F.generateDomainMesh({ type: 'l', nx: 12, ny: 12, width: 2, height: 2 });
  const inc = new Map();
  mesh.elements.forEach((e, k) => e.nodes.forEach((u, i) => {
    const v = e.nodes[(i + 1) % 3], w = e.nodes[(i + 2) % 3];
    const key = Math.min(u, v) + ':' + Math.max(u, v);
    if (!inc.has(key)) inc.set(key, []);
    inc.get(key).push(w);
  }));
  function circle(a, b, c) {
    const d = 2 * ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
    const bx = b.x - a.x, by = b.y - a.y, cx = c.x - a.x, cy = c.y - a.y;
    const x = a.x + ((bx * bx + by * by) * cy - (cx * cx + cy * cy) * by) / d;
    const y = a.y + ((cx * cx + cy * cy) * bx - (bx * bx + by * by) * cx) / d;
    return { x, y, r2: (x - a.x) ** 2 + (y - a.y) ** 2 };
  }
  let checked = 0;
  for (const [key, tri] of inc) {
    if (tri.length !== 2) continue;
    const [u, v] = key.split(':').map(Number);
    const A = mesh.nodes[tri[0]], B = mesh.nodes[tri[1]], U = mesh.nodes[u], V = mesh.nodes[v];
    const ca = circle(U, V, A), cb = circle(U, V, B);
    assert.ok((B.x - ca.x) ** 2 + (B.y - ca.y) ** 2 >= ca.r2 - 1e-9, 'Delaunay violated on edge ' + key);
    assert.ok((A.x - cb.x) ** 2 + (A.y - cb.y) ** 2 >= cb.r2 - 1e-9, 'Delaunay violated on edge ' + key);
    checked++;
  }
  assert.ok(checked > 50);
});

test('manufactured Poisson solution converges at second order in the mesh size', () => {
  const exact = (x, y) => Math.sin(Math.PI * x) * Math.sin(Math.PI * y);
  const source = (x, y) => 2 * Math.PI * Math.PI * Math.sin(Math.PI * x) * Math.sin(Math.PI * y);
  const boundary = n => (n.x === 0 || n.y === 0 || Math.abs(n.x - 1) < 1e-12 || Math.abs(n.y - 1) < 1e-12) ? 0 : undefined;
  const errors = [];
  for (const n of [8, 16, 32]) {
    const mesh = F.generateFieldMesh(n, n, 1, 1);
    const s = F.solvePoisson(mesh, { source, boundary });
    let num = 0, den = 0;
    s.values.forEach((v, i) => { const e = v - exact(mesh.nodes[i].x, mesh.nodes[i].y); num += e * e; den += exact(mesh.nodes[i].x, mesh.nodes[i].y) ** 2; });
    errors.push(Math.sqrt(num / den));
  }
  assert.ok(errors[0] > errors[1] && errors[1] > errors[2], 'errors decrease ' + errors.join(','));
  const order = Math.log2(errors[0] / errors[2]) / 2;
  assert.ok(order > 1.7 && order < 2.6, 'convergence order ' + order);
});

test('Dirichlet boundary values hold exactly and PCG reports its residual history', () => {
  const mesh = F.generateFieldMesh(10, 8, 2, 1);
  const s = F.solvePoisson(mesh, { conductivity: 1, source: 5, boundary: n => n.x === 0 ? 37.5 : Math.abs(n.x - 2) < 1e-12 ? -12.25 : undefined });
  mesh.nodes.forEach((n, i) => {
    if (n.x === 0) assert.equal(s.values[i], 37.5);
    if (Math.abs(n.x - 2) < 1e-12) assert.equal(s.values[i], -12.25);
  });
  const h = s.diagnostics.history;
  assert.ok(h.length >= 2); assert.equal(h[0], 1);
  assert.ok(Math.abs(h[h.length - 1] - s.diagnostics.relativeResidual) < 1e-14);
  assert.ok(s.diagnostics.relativeResidual < 1e-10);
  assert.ok(s.diagnostics.iterations <= h.length);
});

test('annulus mesh preserves the hole, boundary sets and the exact radial profile', () => {
  const mesh = F.generateDomainMesh({ type: 'annulus', nx: 24, ny: 20, width: 2, height: 2, ratio: 0.35 });
  assert.equal(mesh.chi, 0);
  const rIn = mesh.meta.rIn, rOut = mesh.meta.rOut, c = mesh.meta.center;
  close(mesh.quality.area, Math.PI * (rOut * rOut - rIn * rIn), 0.1);
  const hot = new Set(mesh.bcSets.hot), cold = new Set(mesh.bcSets.cold);
  mesh.nodes.forEach((n, i) => {
    const r = Math.hypot(n.x - c.x, n.y - c.y);
    if (hot.has(i)) assert.ok(r < (rIn + rOut) / 2);
    if (cold.has(i)) assert.ok(r > (rIn + rOut) / 2);
  });
  const s = F.solvePoisson(mesh, { boundary: (n, i) => hot.has(i) ? 100 : cold.has(i) ? 0 : undefined });
  let err = 0;
  s.values.forEach((v, i) => {
    const r = Math.hypot(mesh.nodes[i].x - c.x, mesh.nodes[i].y - c.y);
    err = Math.max(err, Math.abs(v - 100 * Math.log(r / rOut) / Math.log(rIn / rOut)));
  });
  assert.ok(err < 0.5, 'radial solution error ' + err);
});

test('nodal flux is the area weighted element flux', () => {
  const mesh = F.generateFieldMesh(8, 6, 2, 1);
  const s = F.solvePoisson(mesh, { conductivity: 1, boundary: n => n.x === 0 ? 100 : Math.abs(n.x - 2) < 1e-12 ? 0 : undefined });
  assert.equal(s.nodeFlux.length, mesh.nodes.length);
  const acc = mesh.nodes.map(() => ({ x: 0, y: 0, w: 0 }));
  s.elements.forEach(e => e.nodes.forEach(id => { acc[id].x += e.fluxX * e.area; acc[id].y += e.fluxY * e.area; acc[id].w += e.area; }));
  acc.forEach((a, i) => { close(s.nodeFlux[i].x, a.x / a.w, 1e-12); close(s.nodeFlux[i].y, a.y / a.w, 1e-12); });
  mesh.nodes.forEach((n, i) => {
    if (n.x > 0.3 && n.x < 1.7) { close(s.nodeFlux[i].x, 50, 1e-9); close(s.nodeFlux[i].y, 0, 1e-9); }
  });
});

test('VTK export writes a valid legacy unstructured grid and guards its inputs', () => {
  const mesh = F.generateFieldMesh(6, 4, 2, 1);
  const s = F.solvePoisson(mesh, { boundary: n => mesh.boundaryNodes.includes(n.id) ? n.x * 10 : undefined });
  const vtk = F.exportVtk(mesh, { title: 'check', pointScalars: [{ name: 'u', values: s.values }], cellScalars: [{ name: 'area', values: s.elements.map(e => e.area) }] });
  const lines = vtk.trim().split('\n');
  assert.equal(lines[0], '# vtk DataFile Version 3.0'); assert.equal(lines[2], 'ASCII'); assert.equal(lines[3], 'DATASET UNSTRUCTURED_GRID');
  assert.equal(lines[4], 'POINTS ' + mesh.nodes.length + ' float');
  const ci = lines.findIndex(l => l.startsWith('CELLS ')), ti = lines.findIndex(l => l.startsWith('CELL_TYPES '));
  assert.equal(lines[ci], 'CELLS ' + mesh.elements.length + ' ' + mesh.elements.length * 4);
  lines.slice(ci + 1, ci + 1 + mesh.elements.length).forEach(l => assert.equal(l.split(' ').length, 4));
  assert.equal(lines[ti], 'CELL_TYPES ' + mesh.elements.length);
  for (let k = ti + 1; k < ti + 1 + mesh.elements.length; k++) assert.equal(lines[k], '5');
  const pi = lines.findIndex(l => l.startsWith('POINT_DATA'));
  assert.equal(lines[pi], 'POINT_DATA ' + mesh.nodes.length);
  assert.equal(lines[pi + 1], 'SCALARS u float 1');
  assert.ok(!/NaN|Infinity/.test(vtk));
  assert.throws(() => F.exportVtk(mesh, { pointScalars: [{ name: 'bad name', values: s.values }] }), /single word/);
  assert.throws(() => F.exportVtk(mesh, { pointScalars: [{ name: 'u', values: [1] }] }), /length mismatch/);
  assert.throws(() => F.exportVtk({ nodes: mesh.nodes, elements: [{ nodes: [0, 1, 9999] }] }), /invalid node reference/);
  assert.throws(() => F.exportVtk({ nodes: [], elements: [] }), /needs a mesh/);
});

test('triangulation terminates on clustered and near-collinear point sets', () => {
  const pts = [];
  for (let i = 0; i < 60; i++) pts.push([i * 0.1, (i % 3) * 1e-6 + (i % 7) * 0.01]);
  pts.push([3, 2]); pts.push([0, 1.5]);
  const mesh = F.triangulate(pts);
  assert.ok(mesh.elements.length > 0);
  assert.ok(Number.isFinite(mesh.quality.minimumAngle));
  for (const e of mesh.elements) {
    const [a, b, c] = e.nodes.map(id => mesh.nodes[id]);
    assert.ok((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) > 0);
  }
});

test('domain generation validates its spec', () => {
  assert.throws(() => F.generateDomainMesh({ type: 'bogus' }), /unknown domain/);
  assert.throws(() => F.generateDomainMesh({ type: 'annulus', nx: 1, ny: 20, width: 2, height: 2 }), /MeshLab:/);
});

/* ------------------------------------------------------------------ */
/* Moment equilibrium, SHA-256 integrity, CSV exports, determinism     */
/* ------------------------------------------------------------------ */

test('every truss preset balances moments about the origin', () => {
  for (const key of ['pratt', 'warren', 'crane', 'tower', 'dome']) {
    const s = M.solveTruss(M.getPreset(key));
    assert.ok(Math.abs(s.equilibrium.residualMz) < 1e-4, key + ' Mz residual ' + s.equilibrium.residualMz);
    assert.ok(s.equilibrium.relativeMomentResidual < 1e-10, key + ' relative moment residual');
    if (s.dim === 3) {
      assert.ok(Math.abs(s.equilibrium.residualMx) < 1e-4, key + ' Mx');
      assert.ok(Math.abs(s.equilibrium.residualMy) < 1e-4, key + ' My');
    }
  }
});

test('load and reaction moment sums decompose the balance', () => {
  const s = M.solveTruss(M.getPreset('crane'));
  let loadM = 0, reactM = 0;
  s.nodes.forEach((n) => {
    loadM += n.x * n.load[1] - n.y * n.load[0];
    reactM += n.x * n.reaction[1] - n.y * n.reaction[0];
  });
  close(s.equilibrium.sumMz, loadM, 1e-9);
  close(s.equilibrium.sumRMz, reactM, 1e-9);
  close(s.equilibrium.residualMz, loadM + reactM, 1e-9);
});

test('the CST plate balances moments about the origin', () => {
  const s = M.solvePlateCST(patch(), {}, {}, { rightNodesForce: 120000 });
  assert.ok(Math.abs(s.equilibrium.residualM) < 1e-4);
  assert.ok(s.equilibrium.relativeMomentResidual < 1e-10);
  assert.ok(Math.abs(s.equilibrium.momentResidual) < 1e-4);
});

test('edge traction applies its declared total even on partial edge coverage', () => {
  const mesh = M.generatePlateMesh(4, 4, 2, 1, 0);
  /* Keep only the bottom half of the right edge so coverage != mesh.height. */
  const partial = mesh.rightNodes.filter((id) => mesh.nodes[id].y <= 0.5);
  mesh.rightNodes = partial;
  const s = M.solvePlateCST(mesh, {}, { left: 'rollerX', fixedNodes: [mesh.grid[0][0]] },
    { rightNodesForce: 7000 });
  close(s.equilibrium.sumFx, 7000, 1e-9);
});

test('fixedNodes honors node id 0 instead of silently dropping the pin', () => {
  const mesh = patch();
  /* Node 0 is the bottom left corner. Before the fix, a falsy id fell through
     the list guard and the rigid body translation stayed free. */
  const s = M.solvePlateCST(mesh, {}, { left: 'rollerX', fixedNodes: [0] },
    { rightNodesForce: 7000 });
  assert.equal(s.nodes[0].ux, 0);
  assert.equal(s.nodes[0].uy, 0);
  assert.ok(Math.abs(s.nodes[0].ry) > 0 || Math.abs(s.reactions[1]) > 0, 'node 0 carries a reaction');
  assert.ok(Math.abs(s.equilibrium.residualM) < 1e-4);
});

test('sha256 matches the Node crypto implementation on edge cases', () => {
  const crypto = require('crypto');
  const cases = ['', 'a', 'abc', 'x'.repeat(55), 'y'.repeat(56), 'z'.repeat(64),
    'block boundary '.repeat(20), 'uni é÷▓ snow ☃', 'line\nbreaks\ttabs'];
  for (const s of cases) {
    assert.equal(M.sha256(s), crypto.createHash('sha256').update(s, 'utf8').digest('hex'), 'sha256 ' + JSON.stringify(s.slice(0, 20)));
  }
});

test('report JSON integrity hash verifies after removing the integrity block', () => {
  const crypto = require('crypto');
  for (const key of ['pratt', 'crane', 'dome']) {
    const rep = JSON.parse(M.exportReportJson(M.getPreset(key)));
    assert.equal(rep.integrity.algorithm, 'SHA-256');
    assert.match(rep.integrity.hash, /^[0-9a-f]{64}$/);
    const hash = rep.integrity.hash;
    delete rep.integrity;
    assert.equal(crypto.createHash('sha256').update(JSON.stringify(rep), 'utf8').digest('hex'), hash, key);
  }
});

test('every exporter is byte deterministic across fresh solves', () => {
  const a = M.getPreset('tower'), b = M.getPreset('tower');
  assert.equal(M.exportCalculiXInp(a), M.exportCalculiXInp(b));
  assert.equal(M.exportReportJson(a), M.exportReportJson(b));
  assert.equal(M.exportDeformedSvg(a, null, 120), M.exportDeformedSvg(b, null, 120));
  assert.equal(M.exportMembersCsv(M.solveTruss(a)), M.exportMembersCsv(M.solveTruss(b)));
  const m1 = F.generateFieldMesh(8, 6, 2, 1), m2 = F.generateFieldMesh(8, 6, 2, 1);
  assert.deepEqual(m1.elements, m2.elements);
  const s1 = F.solvePoisson(m1, { boundary: (n, i) => m1.boundaryNodes.includes(i) ? n.x : undefined });
  const s2 = F.solvePoisson(m2, { boundary: (n, i) => m2.boundaryNodes.includes(i) ? n.x : undefined });
  assert.deepEqual(s1.values, s2.values);
});

test('member CSV carries headers, member rows and CRLF line endings', () => {
  const s = M.solveTruss(M.getPreset('pratt'));
  const csv = M.exportMembersCsv(s);
  const lines = csv.trimEnd().split('\r\n');
  assert.equal(lines[0].split(',')[0], 'member_id');
  assert.equal(lines.length, s.members.length + 1);
  assert.match(lines[1], /^1,Bottom 1,/);
  assert.ok(csv.endsWith('\r\n'));
  /* numeric fields are machine readable */
  const axial = Number(lines[1].split(',')[5]);
  assert.ok(Math.abs(axial - s.members[0].axialForce) < 1e-6 * Math.abs(s.members[0].axialForce || 1));
});

test('CSV escapes spreadsheet formula injection and separators', () => {
  const model = {
    dim: 2, nodes: [{ x: 0, y: 0, support: 'pin' }, { x: 2, y: 0, support: 'rollerY', load: [1000, 0] }],
    members: [{ a: 0, b: 1, label: '=1+cmd|evil', E: 2e11, A: 0.01, I: 1e-6 }]
  };
  const csv = M.exportMembersCsv(M.solveTruss(model));
  assert.ok(csv.includes("'=1+cmd|evil"), 'formula label escaped');
  assert.ok(!/\r\n\d,=1/.test(csv), 'no raw formula cell');
  for (const evil of ['+2*3', '-9-9', '@SUM(1)', '\t=x()']) {
    assert.equal(M.csvEscape(evil), "'" + evil);
  }
  /* embedded separators additionally get RFC 4180 quoting */
  assert.equal(M.csvEscape('\rcalc'), '"\'\rcalc"');
  assert.equal(M.csvEscape('a,b'), '"a,b"');
  assert.equal(M.csvEscape('say "hi"'), '"say ""hi"""');
  assert.equal(M.csvEscape(-42.5), '-42.5');
  assert.equal(M.csvEscape(Infinity), 'Infinity');
  assert.equal(M.csvEscape('-Infinity'), "'-Infinity");
  assert.throws(() => M.toCsv(['a'], [[1, 2]]), /header/);
});

test('under-constrained trusses fail closed instead of returning mechanisms', () => {
  const unsupported = {
    dim: 2, nodes: [{ x: 0, y: 0 }, { x: 1, y: 0, load: [10, 0] }],
    members: [{ a: 0, b: 1 }]
  };
  assert.throws(() => M.solveTruss(unsupported), /singular/);
  const floatNode = {
    dim: 2, nodes: [{ x: 0, y: 0, support: 'pin' }, { x: 1, y: 0, support: 'pin' }, { x: 0.5, y: 1 }],
    members: [{ a: 0, b: 1 }]
  };
  assert.throws(() => M.solveTruss(floatNode), /singular/);
  const rollerOnly = {
    dim: 2, nodes: [{ x: 0, y: 0, support: 'rollerY' }, { x: 1, y: 0, support: 'rollerY', load: [5, 0] }],
    members: [{ a: 0, b: 1 }]
  };
  assert.throws(() => M.solveTruss(rollerOnly), /singular/);
});

test('a member joining a node to itself raises the zero length error', () => {
  assert.throws(() => M.solveTruss({
    dim: 2, nodes: [{ x: 0, y: 0, support: 'pin' }, { x: 1, y: 0 }],
    members: [{ a: 0, b: 0 }]
  }), /zero length/);
});

test('loads applied at fixed DOFs are absorbed into the reactions', () => {
  const model = {
    dim: 2, material: { E: 2e11 },
    nodes: [{ x: 0, y: 0, support: 'pin', load: [4000, 0] }, { x: 2, y: 0, support: 'rollerY', load: [4000, 0] }],
    members: [{ a: 0, b: 1, E: 2e11, A: 0.01, I: 1e-6 }]
  };
  const s = M.solveTruss(model);
  /* both 4 kN loads flow into the single x restraint at node 0 */
  close(s.reactions[0], -8000, 1e-6);
  close(s.nodes[1].displacement[0], 4000 * 2 / (2e11 * 0.01), 1e-12);
  close(s.equilibrium.residualX, 0, 1e-6);
  assert.ok(s.equilibrium.relativeMomentResidual < 1e-10);
});

test('Poisson rejects an unanchored disconnected component', () => {
  const mesh = {
    nodes: [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 0, y: 1 }, { x: 5, y: 0 }, { x: 6, y: 0 }, { x: 5, y: 1 }],
    elements: [{ nodes: [0, 1, 2] }, { nodes: [3, 4, 5] }]
  };
  assert.throws(() => F.solvePoisson(mesh, { boundary: (n, i) => i < 3 ? 10 : undefined }), /every mesh component/);
  const ok = F.solvePoisson(mesh, { boundary: (n, i) => i < 3 ? 10 : (i === 3 ? 4 : undefined) });
  ok.values.forEach((v, i) => close(v, i < 3 ? 10 : 4, 1e-9));
});

test('the cantilever crane reaction moment balances the applied tip moment', () => {
  const model = M.getPreset('crane');
  const s = M.solveTruss(model);
  /* Tip load total is vertical at x = 12.5; the built-in base supplies the
     counter moment through its reaction pair. */
  assert.ok(Math.abs(s.equilibrium.sumMz) > 1e6, 'applied moment exists');
  close(s.equilibrium.residualMz, 0, 1e-6);
  assert.ok(s.equilibrium.relativeMomentResidual < 1e-10);
});

/* ------------------------------------------------------------------ */
/* Adversarial audit regression tests                                  */
/* ------------------------------------------------------------------ */

test('sha256 encodes lone surrogates as U+FFFD like the WHATWG encoder', () => {
  const crypto = require('crypto');
  for (const s of ['\uD800', '\uDC00', 'a\uD800b', '\uD800\uD800', 'end\uDC00', '\uD83D\uDE00']) {
    assert.equal(M.sha256(s), crypto.createHash('sha256').update(s, 'utf8').digest('hex'), JSON.stringify(s));
  }
});

test('CSV escapes formula triggers after leading whitespace', () => {
  for (const evil of [' =1+1', '  -2+3', ' @SUM(1)', ' \t=x()', '\v=2', '\f+9']) {
    assert.equal(M.csvEscape(evil), "'" + evil);
  }
  assert.equal(M.csvEscape('  plain'), '  plain');
  assert.equal(M.csvEscape('  42'), '  42');
});

test('VTK export defaults only the missing third vector component', () => {
  const mesh = F.generateFieldMesh(3, 2, 1, 1);
  const deck = F.exportVtk(mesh, { title: 't',
    pointVectors: [{ name: 'u', values: mesh.nodes.map(() => [1, 2]) }] });
  assert.ok(deck.includes('1 2 0'), 'two component vector exports z 0');
  const deck2 = F.exportVtk(mesh, { title: 't',
    pointVectors: [{ name: 'u', values: mesh.nodes.map(() => ({ x: 1, y: 2 })) }] });
  assert.ok(deck2.includes('1 2 0'), 'object vector exports z 0');
  assert.throws(() => F.exportVtk(mesh, { title: 't',
    pointVectors: [{ name: 'u', values: mesh.nodes.map(() => ({ x: 1, y: 2, z: NaN })) }] }), /finite/);
  assert.throws(() => F.exportVtk({ nodes: [{ x: 0, y: 0, z: NaN }, { x: 1, y: 0 }, { x: 0, y: 1 }],
    elements: [{ nodes: [0, 1, 2] }] }, {}), /finite/);
});

test('INP export coerces a non-string model name instead of crashing', () => {
  const model = M.getPreset('pratt');
  model.name = 42;
  const deck = M.exportCalculiXInp(model);
  assert.ok(deck.length > 100);
  const named = { nodes: model.nodes, members: model.members, dim: 2, name: { evil: true } };
  assert.ok(M.exportCalculiXInp(named).length > 0);
});

test('extreme INP values keep ten significant digits', () => {
  const model = {
    dim: 2, name: 'big',
    nodes: [{ x: 0, y: 0, support: 'pin' }, { x: 123456789.123, y: 0, support: 'rollerY', load: [1000, 0] }],
    members: [{ a: 0, b: 1, E: 2e11, A: 0.01, I: 1e-6 }]
  };
  const deck = M.exportCalculiXInp(model);
  assert.ok(deck.includes('1.234567891'), 'full precision coordinate, got: ' + deck.match(/1\.23456[0-9e+]*/));
});

test('fractional and string node ids in plate BCs and loads fail closed', () => {
  const mesh = patch();
  const loads = { nodal: [{ node: 2.5, fx: 100 }] };
  assert.throws(() => M.solvePlateCST(mesh, {}, { left: 'fixed' }, loads), /integer/);
  assert.throws(() => M.solvePlateCST(mesh, {}, { left: 'fixed', fixedNodes: [0, 2.5] }, {}), /integer/);
  assert.throws(() => M.solvePlateCST(mesh, {}, { left: 'fixed', fixedNodes: ['2'] }, {}), /integer/);
  assert.throws(() => M.solvePlateCST(mesh, {}, { constraints: [{ node: 1.5, x: true }] }, {}), /integer/);
  const frac = patch(); frac.rightNodes = [1.5, 4];
  assert.throws(() => M.solvePlateCST(frac, {}, { left: 'fixed' }, { rightNodesForce: 100 }), /integer/);
  const gone = patch(); gone.rightNodes = [4, 99];
  assert.throws(() => M.solvePlateCST(gone, {}, { left: 'fixed' }, { rightNodesForce: 100 }), /integer/);
});

test('non-array plate BC and load containers produce validation errors', () => {
  const mesh = patch();
  assert.throws(() => M.solvePlateCST(mesh, {}, { fixedNodes: 0 }, {}), /must be an array/);
  assert.throws(() => M.solvePlateCST(mesh, {}, { constraints: 5 }, {}), /must be an array/);
  assert.throws(() => M.solvePlateCST(mesh, {}, { left: 'bogus' }, {}), /unknown plate left/);
  assert.throws(() => M.solvePlateCST(mesh, {}, { right: 3 }, {}), /unknown plate right/);
  assert.throws(() => M.solvePlateCST(mesh, {}, { left: 'fixed' }, { nodal: 7 }), /must be an array/);
});

test('a single boundary node carries its declared traction as a point load', () => {
  const mesh = patch();
  mesh.rightNodes = [mesh.rightNodes[0]];
  const s = M.solvePlateCST(mesh, {}, { left: 'fixed' }, { rightNodesForce: 5000 });
  close(s.equilibrium.sumFx, 5000, 1e-9);
  const empty = patch(); empty.rightNodes = [];
  assert.throws(() => M.solvePlateCST(empty, {}, { left: 'fixed' }, { rightNodesForce: 5000 }), /at least one boundary node/);
  const missing = { nodes: patch().nodes, elements: patch().elements };
  assert.throws(() => M.solvePlateCST(missing, {}, { left: 'fixed' }, { rightNodesForce: 5000 }), /boundary node list/);
});

test('constraints with x 0 pin to zero and nonzero prescriptions throw', () => {
  const mesh = patch();
  const s = M.solvePlateCST(mesh, {},
    { constraints: [{ node: 0, x: 0, y: 0 }].concat(mesh.leftNodes.slice(1).map(n => ({ node: n, x: true }))) },
    { rightNodesForce: 3000 });
  assert.equal(s.nodes[0].ux, 0);
  assert.equal(s.nodes[0].uy, 0);
  assert.throws(() => M.solvePlateCST(patch(), {}, { constraints: [{ node: 0, ux: 5 }] }, {}), /nonzero prescribed/);
  assert.throws(() => M.solvePlateCST(patch(), {}, { constraints: [{ node: 0, uy: -1 }] }, {}), /nonzero prescribed/);
});

test('sub micron plates keep their elements under a relative sliver tolerance', () => {
  const mesh = M.generatePlateMesh(4, 4, 5e-7, 5e-7, 0);
  assert.equal(mesh.nodes.length, 25);
  assert.ok(mesh.elements.length >= 32, 'elements survived, got ' + mesh.elements.length);
});

test('stress measures reject overflow instead of serializing null', () => {
  assert.throws(() => M.stressMeasures(1e308, 1e308, 1e308), /not finite/);
  const fine = M.stressMeasures(1e8, -1e8, 5e7);
  assert.ok(isFinite(fine.vonMises) && fine.tauMax >= 0);
});

test('scalar load vectors in the top level loads array fail closed', () => {
  const base = M.getPreset('pratt');
  assert.throws(() => M.solveTruss({ ...base, loads: [{ node: 2, f: 500 }] }), /load vector must be an array/);
  const ok = M.solveTruss({ ...base, loads: [{ node: 2, f: [0, -500] }] });
  close(ok.equilibrium.sumFy, ok.nodes.reduce((s, n) => s + n.load[1], 0), 1e-9);
});

