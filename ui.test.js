'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const { boot } = require('./dom-fixture.js');
const fs = require('node:fs');
const path = require('node:path');

test('actual static shell boots all three browser modules without reference errors', () => {
  const app = boot();
  assert.equal(app.$('solve-status').textContent, 'SOLVED');
  assert.equal(app.$('preset-select').value, 'pratt');
  assert.equal(app.$('member-tbody').children.length, 34);
  assert.match(app.$('verification-summary').textContent, /VERIFIED/);
  assert.equal(app.$('preset-select').getElementsByTagName('option').length, 7);
});
test('every structural preset and every supported result view runs against actual DOM bindings', () => {
  const app = boot();
  for (const key of ['pratt', 'warren', 'crane', 'tower', 'dome', 'plate']) {
    app.input('preset-select', key, 'change');
    for (const view of ['force', 'stress', 'deformed', 'modes']) {
      app.input('view-select', view, 'change');
      assert.ok(app.$('structure-svg').children.length > 0, `${key} ${view}`);
      assert.ok(app.$('member-tbody').children.length > 0);
      if (key === 'plate' && view === 'modes') assert.equal(app.$('view-select').value, 'stress');
    }
    assert.equal(app.$('solve-status').textContent, 'SOLVED');
  }
});
test('showcase buttons synchronize preset, live select values and scenario loads', () => {
  const app = boot(), buttons = app.document.getElementsByClassName('scenario-btn');
  for (const [i, key, view, load] of [[0, 'pratt', 'force', '1'], [1, 'crane', 'deformed', '5'], [2, 'plate', 'stress', '1'], [3, 'heat', 'stress', '1']]) {
    buttons[i].click();
    assert.equal(app.$('preset-select').value, key); assert.equal(app.$('view-select').value, view); assert.equal(app.$('load-range').value, load);
    assert.equal(buttons[i].getAttribute('aria-pressed'), 'true'); assert.equal(app.$('solve-status').textContent, 'SOLVED');
    if (key === 'crane') assert.equal(app.$('safety-badge').textContent, 'CRITICAL');
  }
});
test('load sweeps change solved displacement and zero loads remain fully inspectable', () => {
  const app = boot(), first = app.$('stat-grid').textContent;
  app.input('load-range', 2); assert.notEqual(app.$('stat-grid').textContent, first); assert.equal(app.$('load-output').textContent, '2.0x');
  app.input('load-range', 0); assert.match(app.$('stat-grid').textContent, /unloaded/); assert.match(app.$('verification-summary').textContent, /VERIFIED/);
  app.input('member-filter', 'compression', 'change'); assert.equal(app.$('member-tbody').children.length, 0);
});
test('plate zero load, filters and PCG verification render without invalid values', () => {
  const app = boot(); app.input('preset-select', 'plate', 'change'); app.input('load-range', 0);
  assert.match(app.$('stat-grid').textContent, /unloaded/); assert.match(app.$('safety-meter').textContent, /No applied load/);
  app.input('member-filter', 'critical', 'change'); assert.equal(app.$('member-tbody').children.length, 0);
  assert.ok(!app.$('stat-grid').textContent.includes('NaN'));
});
test('temperature experiment renders actual Delaunay nodes, PCG diagnostics and flux', () => {
  const app = boot(); app.input('preset-select', 'heat', 'change');
  assert.match(app.$('stat-grid').textContent, /280/); assert.match(app.$('stat-grid').textContent, /165 nodes/);
  assert.match(app.$('model-title').textContent, /Delaunay/); assert.equal(app.$('member-tbody').children.length, 280);
  assert.equal(app.$('view-select').disabled, true); assert.match(app.$('verification-summary').textContent, /VERIFIED/);
  app.input('load-range', 0); assert.match(app.$('mode-bars').textContent, /exact solution/);
  app.input('member-filter', 'compression', 'change'); assert.ok(app.$('member-tbody').children.length < 280);
});
test('round-trip checks reproduce every experiment and independently compare plate solvers', () => {
  const app = boot();
  for (const key of ['pratt', 'crane', 'dome', 'plate', 'heat']) {
    app.input('preset-select', key, 'change'); app.$('run-checks').click(); assert.match(app.$('status-message').textContent, /^PASS/);
    if (key === 'plate') assert.match(app.$('status-message').textContent, /independent Gaussian/);
  }
});
test('keyboard inspection uses one diagram tab stop and arrow navigation', () => {
  const app = boot(), marks = app.elements().filter(e => e.getAttribute('data-inspect') !== null);
  assert.equal(marks.filter(e => e.getAttribute('tabindex') === '0').length, 1);
  marks[0].focus(); marks[0].dispatch('keydown', { key: 'ArrowRight' }); assert.equal(app.document.activeElement, marks[1]);
  marks[1].dispatch('keydown', { key: 'Enter' }); assert.equal(marks[1].getAttribute('aria-pressed'), 'true'); assert.match(app.$('inspector').textContent, /PINNED RESULT/);
  marks[1].dispatch('keydown', { key: 'Escape' }); assert.equal(marks[1].getAttribute('aria-pressed'), 'false'); assert.match(app.$('inspector').textContent, /LIVE INSPECTOR/);
});
test('pointer inspection is bounded at the canvas edge and supports pinning', () => {
  const app = boot(), mark = app.elements().find(e => e.getAttribute('data-inspect') === '0');
  mark.dispatch('pointermove', { clientX: 1000, clientY: 1000 }); assert.equal(app.$('tooltip').getAttribute('hidden'), null);
  assert.equal(app.$('tooltip').getAttribute('style'), 'left:552px;top:232px;');
  mark.click(); assert.equal(app.$('tooltip').getAttribute('hidden'), ''); assert.match(app.$('inspector').textContent, /PINNED/);
});
test('deformation scale changes geometry and preserves focused inspection after render', () => {
  const app = boot(); app.input('view-select', 'deformed', 'change');
  const mark = app.elements().find(e => e.getAttribute('data-inspect') === '1'); mark.focus(); const x = mark.getAttribute('x1');
  app.input('scale-range', 350); const updated = app.elements().find(e => e.getAttribute('data-inspect') === '1');
  assert.notEqual(updated.getAttribute('x1'), x); assert.equal(app.document.activeElement, updated); assert.equal(app.$('scale-output').textContent, '350x');
});
test('modal selection and phase update live values without rebuilding the selector', () => {
  const app = boot(); app.input('view-select', 'modes', 'change');
  const option = app.$('mode-select').children[0]; app.input('mode-select', 2, 'change'); app.input('mode-phase', 90);
  assert.equal(app.$('mode-select').value, '2'); assert.equal(app.$('mode-phase').value, '90'); assert.equal(app.$('mode-select').children[0], option);
  assert.match(app.$('mode-readout').textContent, /phase 90 deg/);
  app.input('preset-select', 'plate', 'change'); app.input('preset-select', 'pratt', 'change'); app.input('view-select', 'modes', 'change');
  assert.equal(app.$('mode-select').getElementsByTagName('option').length, 4); assert.equal(app.$('mode-select').value, '0');
});
test('modal animation has one frame loop, time based phase and visibility cleanup', () => {
  const app = boot(); app.input('view-select', 'modes', 'change'); app.$('mode-play').click(); assert.equal(app.frames.size, 1);
  app.frame(0); app.frame(100); assert.equal(app.$('mode-phase').value, '9'); assert.equal(app.frames.size, 1);
  app.$('mode-play').click(); assert.equal(app.frames.size, 0); app.$('mode-play').click(); assert.equal(app.frames.size, 1);
  app.document.hidden = true; app.document.dispatch('visibilitychange'); assert.equal(app.frames.size, 0); assert.equal(app.$('mode-play').getAttribute('aria-pressed'), 'false');
});
test('modal geometry stays inside a stable viewBox through the phase cycle', () => {
  const app = boot(); app.input('view-select', 'modes', 'change'); const bounds = app.$('structure-svg').getAttribute('viewBox');
  for (const phase of [45, 90, 180, 270]) { app.input('mode-phase', phase); assert.equal(app.$('structure-svg').getAttribute('viewBox'), bounds); }
});
test('reduced motion uses manual mode inspection', () => {
  const app = boot({ reducedMotion: true }); app.input('view-select', 'modes', 'change'); app.$('mode-play').click(); assert.equal(app.frames.size, 0); assert.match(app.$('status-message').textContent, /Reduced motion/);
});
test('overlay controls and theme work in structural and scalar views', () => {
  const app = boot();
  for (const key of ['pratt', 'plate', 'heat']) {
    app.input('preset-select', key, 'change');
    for (const id of ['toggle-ghost', 'toggle-reactions', 'toggle-grid']) { app.$(id).checked = false; app.$(id).dispatch('change'); app.$(id).checked = true; app.$(id).dispatch('change'); }
    app.$('theme-toggle').click(); assert.equal(app.$('solve-status').textContent, 'SOLVED');
  }
});
test('table sort reports direction and filters agree with signed force', () => {
  const app = boot(); app.input('member-filter', 'tension', 'change'); const count = app.$('member-tbody').children.length; assert.ok(count > 0 && count < 34);
  const header = app.$('member-table').getElementsByTagName('thead')[0]; header.getElementsByTagName('button')[2].click();
  assert.equal(header.children[0].children[2].getAttribute('aria-sort'), 'ascending'); header.getElementsByTagName('button')[2].click(); assert.equal(header.children[0].children[2].getAttribute('aria-sort'), 'descending');
});
test('exports deliver valid structural reports, unique plate loads and release Blob URLs', async () => {
  const app = boot(), buttons = app.document.getElementsByClassName('dl-btn');
  buttons[2].click(); assert.equal(app.downloads[0].filename, 'pratt.json'); assert.equal(JSON.parse(await app.downloads[0].blob.text()).members.length, 34);
  app.input('preset-select', 'plate', 'change'); buttons[0].click();
  const deck = await app.downloads[1].blob.text(), block = deck.split('*CLOAD')[1].split('*NODE PRINT')[0];
  const lines = block.trim().split('\n'), ids = lines.map(line => line.split(',').slice(0, 2).join(',')); assert.equal(new Set(ids).size, ids.length);
  const force = lines.reduce((sum, line) => sum + Number(line.split(',')[2]), 0); assert.ok(Math.abs(force) < 1e-5);
  assert.ok(app.blobs.size > 0); app.runTimers(); assert.equal(app.blobs.size, 0);
});
test('scalar exports contain temperatures, verification and standalone SVG geometry', async () => {
  const app = boot(); app.input('preset-select', 'heat', 'change'); const buttons = app.document.getElementsByClassName('dl-btn');
  assert.equal(buttons[0].disabled, true); buttons[0].click(); assert.equal(app.downloads.length, 0);
  buttons[1].click(); buttons[2].click(); const svg = await app.downloads[0].blob.text(), json = JSON.parse(await app.downloads[1].blob.text());
  assert.equal((svg.match(/<polygon /g) || []).length, 280); assert.equal(json.values.length, 165); assert.ok(json.diagnostics.relativeResidual < 1e-10);
  assert.equal(json.model.field.heater, 1600); assert.equal(json.model.field.hotEdge, 100);
});
test('invalid range data and solver failures are visible without stale data corruption', () => {
  const app = boot(), geometry = app.$('structure-svg').textContent; app.input('scale-range', 'invalid'); assert.equal(app.$('structure-svg').textContent, geometry);
  const original = app.sandbox.MeshLab.solveTruss; app.sandbox.MeshLab.solveTruss = () => { throw new Error('Example mechanism'); };
  app.input('load-range', 2); assert.equal(app.$('solve-status').textContent, 'CHECK INPUT'); assert.equal(app.$('load-range').value, '1'); assert.match(app.$('status-message').textContent, /Example mechanism/);
  app.sandbox.MeshLab.solveTruss = original; app.input('load-range', 1); assert.equal(app.$('solve-status').textContent, 'SOLVED');
});
test('generated numeric SVG attributes are finite in every experiment', () => {
  const app = boot();
  for (const key of ['pratt', 'crane', 'dome', 'plate', 'heat']) {
    app.input('preset-select', key, 'change');
    for (const element of app.$('structure-svg').getElementsByTagName('*')) {
      for (const [name, value] of element.attributes) if (['x', 'y', 'x1', 'y1', 'x2', 'y2', 'cx', 'cy', 'r', 'points', 'viewBox'].includes(name)) assert.ok(!/NaN|Infinity|n\/a/.test(value), `${key}: ${name} ${value}`);
    }
  }
});
test('static assets are offline, ordered and force line widths remain data driven', () => {
  const html = fs.readFileSync(path.join(__dirname, 'index.html'), 'utf8'), css = fs.readFileSync(path.join(__dirname, 'style.css'), 'utf8');
  assert.ok(!/<(?:script|link)[^>]+(?:src|href)="https?:/.test(html));
  assert.ok(html.indexOf('src="field.js"') < html.indexOf('src="fem.js"')); assert.ok(html.indexOf('src="fem.js"') < html.indexOf('src="app.js"'));
  const rule = css.match(/\.svg-member \{([^}]+)\}/)[1]; assert.ok(!rule.includes('stroke-width'));
});
test('plate ghost boundary toggles draw real undeformed geometry', () => {
  const app = boot(); app.input('preset-select', 'plate', 'change'); app.input('view-select', 'deformed', 'change');
  const ghosts = () => app.$('structure-svg').getElementsByClassName('svg-ghost'); assert.ok(ghosts().length > 0);
  app.$('toggle-ghost').checked = false; app.$('toggle-ghost').dispatch('change'); assert.equal(ghosts().length, 0);
});
test('both themes meet 4.5 to 1 contrast for primary, secondary and status text', () => {
  const css = fs.readFileSync(path.join(__dirname, 'style.css'), 'utf8');
  const luminance = hex => { const rgb = [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16) / 255).map(v => v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4); return rgb[0] * 0.2126 + rgb[1] * 0.7152 + rgb[2] * 0.0722; };
  for (const theme of ['dark', 'light']) {
    const rule = css.match(new RegExp('html\\[data-theme="' + theme + '"\\] \\{([^}]+)\\}'))[1];
    const colors = Object.fromEntries([...rule.matchAll(/--([a-z0-9-]+): (#[a-f0-9]{6});/g)].map(m => [m[1], m[2]]));
    for (const key of ['ink', 'ink-2', 'ink-3', 'accent', 'good', 'warning', 'critical']) {
      const a = luminance(colors[key]), b = luminance(colors.surface), ratio = (Math.max(a, b) + 0.05) / (Math.min(a, b) + 0.05);
      assert.ok(ratio >= 4.5, `${theme} ${key} contrast ${ratio}`);
    }
  }
});
test('zero-load plate rows generate valid finite background colors', () => {
  const app = boot(); app.input('preset-select', 'plate', 'change'); app.input('load-range', 0);
  const marks = app.$('member-tbody').getElementsByClassName('state-mark');
  assert.equal(marks.length, 580);
  for (const mark of marks) {
    const style = mark.getAttribute('style') || '';
    assert.ok(!style.includes('undefined') && !style.includes('NaN'), style);
    assert.match(style, /background:/);
  }
});

