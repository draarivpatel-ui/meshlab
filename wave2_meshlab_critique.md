# MeshLab Wave 2 polish critique

Date: 2026-11-20. Auditor: Wave 2 polish judge (implementation sub-agent).
Scope: app.js, fem.js, field.js, style.css, index.html, README.md, DEVPOST_UPGRADE.md, plus live headless-Chrome verification at 1440x900 and 390x844.

## Baseline

- `node --test test.js` : 102 pass / 0 fail
- `node --test numerics.test.js` : 51 pass / 0 fail
- `node --test ui.test.js` : 29 pass / 0 fail
- Total 182/182, matching README.md:65 and DEVPOST_UPGRADE.md:56.
- README verification command `node --test --test-reporter=spec test.js numerics.test.js ui.test.js` reproduces 182/182.
- `node --check` passes on all seven scripts.

## Verified claims (all from real runs, no reliance on printed numbers)

- pratt sumFy -360000 N, peak translation 7.1582247 mm (README:15)
- crane at 5x: governing yield FoS 0.8615384615, member "Upper boom 1", 1 critical member (README:16)
- plate: 336 nodes / 580 CST triangles, peak von Mises 72.9961 MPa, peak/nominal 2.919844677, PCG 285 iterations, relative residual 7.03e-13 (README:17, :51-52, DEVPOST:63-65)
- heat preset: 165 nodes / 280 triangles, peak temperature 144.7337 C (README:18)
- lbeam 285/488, annulus 461/818 chi=0, sink 206 nodes / 304 triangles (README:53-55)
- Delaunay area 3.84 m2, minimum angle 33.6772 deg
- affine harmonic patch max nodal error 3.149e-10 C; annulus radial error 9.907e-2 C; MMS order 2.014
- warren/crane/tower/dome applied totals -120000 / -130000 / -24000 / -80000 N, each applied exactly once
- Contrast, offline-asset, script-order and fixture checks all pass.

## Independent verification beyond the suite

- 200 random trusses: every solve passes the residual + energy verification identity (worst relative residual 3.2e-9).
- 60 random point sets triangulated: CCW winding, manifold edges, no orphans in all; worst min angle 0.06 deg on pathological clouds, no crashes.
- Domain sweep (l/annulus/sink at 3 densities x 3 random sizes): correct Euler characteristic, PCG converged everywhere.
- Plate direct-vs-PCG differential on 6 mesh densities with holes: max displacement delta < 1e-11.
- Full interaction cycle in headless Chrome (all 10 presets x 4 views, load sweep incl. 0x, resolution re-mesh, round-trip checks, theme toggles, showcase buttons): zero JS exceptions; only console error is the favicon 404.
- Mobile (390x844) overflow probe on plate/dome/annulus/crane: no horizontal overflow; the wide table scrolls inside .table-scroll by design.
- annulus re-mesh at 3x resolution (2671 nodes / 5087 triangles): 46 ms mesh + 25 ms solve, no perf cliff.

## Findings

### [must]

- **M1 - Missing favicon produces a console error on every load.** index.html declares no icon, so Chrome requests /favicon.ico and logs `Failed to load resource: 404` on every page load (reproduced at both 1440x900 and 390x844; the required "consoleErrors": [] bar fails on this alone). Fix: an inline data-URI SVG icon keeps the zero-request, fully-offline property.

### [should]

- **S1 - Field SVG export misstates the boundary values.** app.js:1507 hardcodes `hot edge 100 / cold edge 0` in the exported SVG caption. The sink preset hot edge is 80 C (fem.js:1069), and for the annulus the hot edge scales with the multiplier (app.js:116), so exported artifacts can claim 100 while the solve used 80 or 200. Fix: print the actual `state.model.field.hotEdge/coldEdge`.
- **S2 - Zero-force members get the compression color in the result table.** app.js:1141-1143 colors the STATE mark `tension ? tension : compression`, so a zero-force member shows compression red while the canvas (app.js:283-291 uses neutral) and the tooltip ("ZERO FORCE") say neutral. pratt has 1 such member (Bottom 7), warren 7. Fix: three-way color.
- **S3 - README preset table miscounts the sink mesh.** README.md:55 says "206 nodes, 303 triangles"; deterministic runs give 304 triangles (nodes correct). Fix the count.

### [nit]

- N1 - style.css:227 hardcodes the live-dot glow `rgba(12,163,12,0.16)` instead of deriving from `--good` (dark #77dfb8 / light #146b45), so the glow hue mismatches the dot in light theme.
- N2 - app.js:1361 labels every field preset's multiplier "HEATER MULTIPLIER"; the annulus has no heater, the slider scales its hot boundary temperature.
- N3 - app.js:499 and :604 give hit targets `role="img"` plus `aria-pressed`; aria-pressed is not supported on role="img" (valid on role="button"). Field hits (app.js:1466) carry aria-pressed with no role at all.
- N4 - app.js:742 formats strain energy with `fixed(v,4)` which would print exponential notation for |v| >= 1e7; all current presets stay in fixed range (max 1237.1585 J), so no live impact.
- N5 - fem.js:1223 exportDeformedSvg hardcodes the dark background #09111E for all themes. Standalone artifact; acceptable.
- N6 - app.js:42 reads prefers-reduced-motion once at boot and never re-evaluates.

## Score

Pre-fix: **9.2 / 10**.
Engine and numerics: 9.5 (clean fuzz/differential verification, honest verification identities, exact closed-form benchmarks). UI: 9 (clean layouts at both viewports, one console error, minor label/ARIA nits). Honesty: 8.8 (README miscount, export caption). Accessibility: 9.

## Post-fix status (same day, applied and verified)

- M1 fixed: inline SVG data-URI favicon in index.html; both viewport checks now report `"consoleErrors": []` and exit 0.
- S1 fixed: fieldSvgExport prints the solved `hotEdge`/`coldEdge` (sink exports "hot edge 80", heat "hot edge 100").
- S2 fixed: zero-force members show the neutral color `#8898ac` in the table (regression test uses pratt Bottom 7; Bottom 4 as the tension control).
- S3 fixed: README sink count corrected to 304 triangles.
- N1 fixed: live-dot glow derives from `--good` via color-mix.
- N2 fixed: multiplier label reads "HOT EDGE MULTIPLIER" for heater-less field presets (annulus), "HEATER MULTIPLIER" otherwise.
- N3 fixed: member/plate hit targets use role="button" (aria-pressed valid); field triangle hits get role="button" + tabindex="-1".
- N6 fixed: prefers-reduced-motion is tracked live via matchMedia change listeners and stops a running animation.
- Regression coverage: 3 new UI tests (favicon data URI, neutral zero-force mark, solved boundary values in field SVG export). Suite: 185/185 (102 + 51 + 32).
- Post-fix browser runs: full interaction cycle (all presets, views, sliders, exports, theme) reports zero errors and no overflow; syntax checks pass on all seven scripts.
- README.md:65 and DEVPOST_UPGRADE.md:56 updated to 185 tests / 32 DOM checks.

Deferred (accepted, with reason): N4 (no live impact at any preset load), N5 (standalone SVG artifact, dark canvas is its documented style). N6-era one-shot read was upgraded rather than deferred.

Post-fix score: **9.5 / 10**.

## Changed files

- index.html (+1: data-URI favicon link)
- app.js (favicon-independent fixes S1, S2, N2, N3, N6)
- style.css (N1)
- ui.test.js (+3 regression tests)
- README.md, DEVPOST_UPGRADE.md (test counts, sink triangle count)
- wave2_meshlab_critique.md (this file)
