# MeshLab

## See the physics. Check the math.

MeshLab is an engineering workbench you can open in a browser and inspect immediately. Trace the load path through a bridge, push a crane through its yield limit, reveal the stress concentration around a hole, or follow heat across an irregular triangular mesh. Every result comes from a numerical solve with inspectable measurements and a reproducibility check.

**[Open the workbench](https://aariv.me/meshlab/)**

No package installation, external assets, build step, service, or runtime dependency is required. The complete application also works by opening `index.html` directly, including offline.

## A five second demo

| One click experiment | What to inspect | Verified result at the showcase load |
| --- | --- | --- |
| Trace a bridge | Blue tension, red compression, support reactions | 360 kN total vertical load; 7.1582 mm peak translation |
| Overload a crane | Dashed critical member, deformed outline, yield safety | At 5x load, Upper boom 1 governs with yield FoS 0.861538; one critical member |
| Stress around a hole | Individual CST stress tensors and the cavity boundary | 72.9961 MPa peak von Mises; peak/nominal ratio 2.919845 |
| Follow the heat | Interior heater, temperature, constant triangle heat flux | 165 nodes and 280 Delaunay triangles; 144.7337 C peak temperature |

Move the load or heater slider to compare scenarios. Hover any member or triangle for live measurements, select it to pin those measurements, then use **Run round-trip check** to serialize and reproduce the model. The plate check compares sparse PCG against independent Gaussian elimination.

## What is implemented

### Structural analysis

- **2D and 3D pin jointed trusses:** direct stiffness assembly, scaled partial pivot Gaussian elimination, nodal translations and support reactions.
- **Member mechanics:** axial force, stress, strain, elongation and strain energy.
- **Safety:** material yield and Euler buckling `Pcr = pi^2 E I / (Keff L)^2`. The lower factor of safety governs. Zero force is distinguished from compression.
- **Natural frequencies:** lumped and consistent mass formulations, full Cholesky reduction of the generalized eigenproblem, Jacobi eigenvectors and mass normalization. Every mode includes its relative eigen equation residual.
- **Plane stress plates:** constant strain triangles, isotropic material law, consistent edge traction normalized by the actually covered boundary span, principal stresses, maximum in-plane shear and von Mises stress. The browser solves the plate with sparse Jacobi PCG while the public API also exposes the direct solver.
- **Global equilibrium:** every solve reports applied forces, support reactions and their moments about the origin. The interface shows force and moment residuals for both 2D and 3D structures and the plate.

### Scalar field analysis

- **Delaunay meshing:** deterministic convex-hull sweep and Lawson edge-flip triangulation with normalized geometry, duplicate point removal, hull extraction, positive triangle winding, sliver cleanup and minimum angle / radius-edge / Euler characteristic reporting.
- **Non-convex and holed domains:** rectangle, L-bracket, annulus with a real hole (Euler V-E+F = 0) and a finned heat sink. Polygon and hole clipping keeps the boundary Delaunay valid; boundary sets mark hot and cold edges for each domain.
- **Poisson FEM:** linear scalar triangles solve `-div(k grad T) = q`. Function sources integrate on a three point edge midpoint rule (exact for quadratics, tighter than a centroid sample for the Gaussian heater). Nonzero prescribed boundary values enter the reduced right hand side correctly; unconstrained edges have the natural zero flux condition.
- **Sparse PCG:** compressed sparse row multiplication, Jacobi preconditioning, symmetry and positive curvature checks, finite iteration limits, a recomputed actual residual before reporting convergence, and a residual history the app plots as a log scale sparkline.
- **Heat flow:** prescribed hot and cold boundary sets with insulated free edges and an adjustable interior Gaussian heater. Flux is `-k grad T`, rendered as arrows whose length encodes local magnitude; heat balance is reported per unit depth and every node exposes exact `(x, y, T, |grad T|)` on hover.

All structural engine quantities use SI units: m, N, Pa, kg, Hz and J. Display values use mm, kN and MPa. The heat experiment uses C, W/m2 and W/m per unit depth.

## Presets

| Key | Model | Formulation |
| --- | --- | --- |
| `pratt` | Pratt Bridge Truss (7-Bay), 16 nodes, 34 members | 2D truss |
| `warren` | Warren Roof Truss | 2D truss |
| `crane` | Cantilever Crane Boom | 2D truss |
| `tower` | Transmission Tower Lattice | 2D truss |
| `dome` | 3D Space Tripod Dome with isometric projection | 3D truss |
| `plate` | Perforated plate, 336 nodes, 580 triangles | CST plane stress |
| `heat` | Irregular mesh, 165 nodes, 280 triangles | Delaunay Poisson FEM |
| `lbeam` | L-bracket non-convex domain, 285 nodes, 488 triangles | Delaunay Poisson FEM |
| `annulus` | Annulus with a real hole, 461 nodes, 818 triangles | Delaunay Poisson FEM |
| `sink` | Finned heat sink domain, 206 nodes, 304 triangles | Delaunay Poisson FEM |

The force count in `m + r - dj` is reported alongside the solved stiffness system. A nonnegative count alone is not a proof of stability. Singular direct systems raise a descriptive error. Under-constrained mechanisms, zero area triangles, collinear nodes, out of range node references and non-finite inputs are rejected before solving.

## Verification

```sh
node --test --test-reporter=spec test.js numerics.test.js ui.test.js
```

**218 tests passed, 0 failed, 0 skipped.** The suites contain 102 structural/export checks, 78 numerical and edge case checks, and 38 DOM/interface checks. All test code uses built in Node modules.

| Reproduced check | Verified value |
| --- | --- |
| Holed plate peak/nominal stress | 2.919844677 |
| Holed plate PCG convergence | 285 iterations; relative residual below 1e-12 |
| Delaunay mesh area | 3.84 m2 |
| Delaunay minimum angle | 33.6772 degrees |
| Harmonic field `T = 100(1 - x / width)` | Maximum nodal error 3.15e-10 C on the default irregular mesh |
| Annulus radial field `T(r)` | Maximum nodal error 9.91e-2 C against `100 ln(r/Rout)/ln(Rin/Rout)` |
| Manufactured solution `sin(pi x) sin(pi y)` | Second order L2 convergence, measured order near 2.0 |
| Domain Euler characteristic | Rectangle / L / sink V-E+F = 1; annulus V-E+F = 0 (one hole) |
| Consistent mass bridge modes | Full generalized eigen residual below 1e-8; mutual mass orthogonality checked |
| Global moment residual | Below 1e-10 relative on every truss preset and the CST plate |
| JSON report integrity | A receiver rehashes the compact report to the embedded SHA-256 value |
| Primary, secondary and status text in both themes | Automated contrast checks require at least 4.5:1 against the panel background |

DOM tests load the actual HTML and all browser scripts into a browser contract fixture. They exercise all ten presets, controls, numeric SVG attributes, mode animation, keyboard inspection, download generation, error recovery and round-trip checks. They verify bindings and behavior; browser layout and screen rendering are separate from this fixture.

Syntax can be checked without running the app:

```sh
node --check field.js
node --check fem.js
node --check app.js
node --check test.js
node --check numerics.test.js
node --check dom-fixture.js
node --check ui.test.js
```

All seven scripts pass syntax checks.

## Interface and accessibility

- Four instant experiments lead directly to real solved models.
- The load sweep covers 0x to 5x, including a fully unloaded state.
- The SVG scales without raster DPI loss. Force strokes and member hit targets keep their screen width as the viewport changes.
- One diagram tab stop avoids hundreds of sequential stops. Arrow keys move between results, Home/End jump to the extremes, Enter/Space pin a result, and Escape clears inspection.
- There are no single character keyboard shortcuts anywhere (WCAG 2.1.4); nothing on the page claims `accesskey` or `aria-keyshortcuts`, and character keys act on nothing.
- Text measurements duplicate every color encoding. Critical members have a dashed mark and a text state.
- The result table is sortable, reports sort direction, and can be scrolled with keyboard focus.
- Tooltips stay inside the plot; the pinned inspector remains readable on small screens.
- Modal animation uses elapsed time, one animation loop and a stable view box. Hiding the page stops it; reduced motion uses the phase slider.
- Light and dark themes have explicit palettes, visible focus states and automated text contrast checks.

## Exports

| Export | Trusses | Plate | Heat field |
| --- | --- | --- | --- |
| `.INP` | CalculiX / Abaqus style truss input deck | `CPS3` plate deck with the actual constrained DOFs and one consolidated load per DOF | Unavailable |
| `.SVG` | Standalone deformed result, including projected z coordinates | Standalone deformed von Mises contour | Standalone temperature contour |
| `.JSON` | Reactions, member results, safety, modes and verification | Nodal results, stress tensors, forces, constraints and PCG diagnostics | Mesh, temperatures, flux, balance and PCG diagnostics |
| `.VTK` | Legacy ASCII `UNSTRUCTURED_GRID` lines with displacement vectors and axial force/stress per member | Unstructured grid with nodal displacement and von Mises scalars | Unstructured grid with temperature, nodal and element heat flux |
| `.CSV` | Member forces, stresses, energy and safety factors | Element stress tensors and energy | Triangle temperatures, flux and energy |

Downloads are generated locally. Object URLs are released after downloading. Every JSON report embeds a SHA-256 integrity hash of its compact serialization so a received file can be rehashed for verification. CSV files are written with CRLF line endings and spreadsheet formula-injection hardening: text cells beginning with `=`, `+`, `-`, `@`, a tab or a carriage return are apostrophe prefixed, and separator-bearing cells are quoted.

## Security

A restrictive `Content-Security-Policy` meta confines the page (`default-src 'none'`, `script-src 'self'`, `style-src 'self' 'unsafe-inline'`, `style-src-elem 'self'`, `style-src-attr 'unsafe-inline'`, `object-src 'none'`, `connect-src 'none'`, `form-action 'none'`, `base-uri 'none'`). Scripts stay strictly same-directory. In CSP3 browsers the element/attribute split means an injected `<style>` element is blocked while data-driven inline style attributes still work; the `style-src` entry remains only as a fallback for older engines. Inline style attributes are allowed because member colors, force widths and the tooltip position are data-driven inline styles; `connect-src 'none'` prevents any exfiltration a style payload could attempt. The app performs no network requests, has no `eval`, no inline event handlers, and its only images are an inline SVG favicon and generated SVG downloads. All user reachable text flows through DOM text nodes, never HTML injection.

## Running and embedding

Open `index.html`, or serve this directory with any static file server. GitHub Pages needs no build command; `.nojekyll` is included.

Browser script order is `field.js`, `fem.js`, then `app.js`. In Node:

```js
const MeshLab = require('./fem.js');
const model = MeshLab.getPreset('pratt');
const result = MeshLab.solveTruss(model);
console.log(result.maxDisplacement, result.verification);
```

## Modeling scope

The structural models are linear elastic axial trusses and plane stress CST plates. They do not include beam bending, plastic redistribution, contact or large deformation. Euler buckling and yield ratios are screening calculations within those assumptions. Mode display amplitude is an inspection setting, not a prediction of forced vibration amplitude.

The perforated plate has finite width and a stair step approximation of the hole. Its measured concentration is compared with nominal membrane stress; it is not claimed as an exact solution for an infinite circular-hole plate. Delaunay field meshes cover rectangular, L-shaped, annular and finned domains; holed domains are discretized by polygon clipping, so the annulus inner radius is a polygonal approximation of the circle. Scalar function sources integrate on a three point edge midpoint rule per triangle.

## Files

| File | Purpose |
| --- | --- |
| `index.html`, `style.css` | Static semantic shell and responsive themes |
| `app.js` | Scenario state, SVG rendering, inspection, verification and exports |
| `fem.js` | Structural solver, presets and public UMD API |
| `field.js` | Delaunay mesh, sparse PCG and scalar Poisson FEM |
| `test.js` | Structural and export regression tests |
| `numerics.test.js` | Solver edge cases and independent math checks |
| `ui.test.js`, `dom-fixture.js` | Actual shell/browser-script contract tests |
| `DEVPOST_UPGRADE.md` | Submission narrative with reproduced numbers |

## License

MIT. Copyright (c) 2026 draarivpatel-ui. See [LICENSE](LICENSE).
