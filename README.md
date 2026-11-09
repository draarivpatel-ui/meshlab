# MeshLab

An interactive structural analysis workbench that runs entirely in the browser. MeshLab solves 2D and 3D pin jointed trusses by the direct stiffness method, checks yield and Euler buckling safety, computes natural frequencies and mode shapes, and solves a perforated 2D plate with constant strain triangles (CST) under plane stress. The interface is a blueprint engineering dashboard with a dark slate theme and an independently selected light blueprint theme.

## Run locally

No build step, no package install, and no server framework are required. The solver is a UMD module, so every file is a static asset.

```sh
git clone https://github.com/draarivpatel-ui/meshlab.git
cd meshlab
python3 -m http.server 8000
```

Open `http://localhost:8000/` in a browser. Opening `index.html` directly also works because the scripts are plain classic scripts. An internet connection is only needed for the three web fonts; the UI falls back to system fonts when offline.

## Test

```sh
node --test test.js
```

The suite has 102 tests and uses only the Node built in modules. It covers linear algebra, joint equilibrium, closed form deflections, strain energy, determinacy, support behaviour, 3D systems, Euler buckling, modal orthogonality, CST patch tests, global equilibrium, Kirsch like stress concentration, all presets, and the three exporters.

## How the solver works

`fem.js` exposes `MeshLab` on `window` in the browser and through `require('./fem.js')` in Node.

- **Direct stiffness:** each bar contributes a `2 x 2` (or `3 x 3` block) axial stiffness matrix assembled into the global system. Supports are applied by reduction, then the reduced system is solved with partial pivot Gaussian elimination. Reactions are recovered from `K u - F`.
- **Member results:** elongation, axial force, stress, strain, and element strain energy are derived from the nodal displacements. The report also includes determinacy (`m + r - 2j` or `m + r - 3j`) and global force residuals.
- **Safety:** every member gets a yield factor of safety and an Euler buckling factor of safety `pi^2 E I / (Keff L)^2`. The smaller of the two governs.
- **Modes:** a lumped or consistent mass matrix is reduced to the free degrees of freedom, transformed to a symmetric standard eigenproblem, and diagonalised with the Jacobi rotation method. Mode shapes are mass normalised and returned with frequency, angular frequency, eigenvalue, and period.
- **Plate:** `generatePlateMesh` builds a structured grid of two triangles per cell and removes elements inside a circular hole. Each CST element uses a plane stress constitutive matrix `D`, a constant strain displacement matrix `B`, and `K_e = B^T D B t A`. Element stress, principal stresses, maximum shear, and von Mises stress are post processed. Applied edge traction is distributed consistently along boundary segments, which reproduces exact uniform uniaxial stress and gives a stress concentration above 2x nominal at the hole.

All quantities are SI: metres, newtons, pascals, kilograms, hertz, and joules. The UI converts to kN, MPa, and mm for display.

## Presets

| Key | Structure | Solver |
| --- | --- | --- |
| `pratt` | Pratt Bridge Truss (7-Bay), 16 nodes, 34 members | 2D truss |
| `warren` | Warren Roof Truss | 2D truss |
| `crane` | Cantilever Crane Boom | 2D truss |
| `tower` | Transmission Tower Lattice | 2D truss |
| `dome` | 3D Space Tripod Dome, isometric projection | 3D truss |
| `plate` | Perforated Plate in Tension (Kirsch Stress Concentration), 336 nodes, 580 triangles | CST plane stress |

## Interface

- **Views:** axial force (diverging blue tension and red compression, line weight by magnitude), von Mises stress (single hue sequential ramp), deformed shape (undeformed ghost plus exaggerated shape), and modal shape (animated mode with a phase scrubber).
- **Controls:** structure preset, view, deformation scale from 1x to 500x, ghost geometry, reaction arrows, grid, and table filter.
- **Interaction:** hover or keyboard focus on any member or triangle opens a tooltip with the full numeric result. The sortable result table is an equivalent text view, legends always carry text labels, and the safety meter ships an icon plus label. The light and dark themes use independently validated palette steps.
- **Results:** stat tiles, global equilibrium, determinacy, a factor of safety meter with a FoS 1.0 limit, and natural frequency bars.

## Exports

| Button | Truss presets | Plate preset |
| --- | --- | --- |
| `.INP` | CalculiX and Abaqus style deck with `*NODE`, `*ELEMENT`, `*MATERIAL`, `*BOUNDARY`, `*CLOAD`, and `*END STEP` | `CPS3` plane stress deck with `*SOLID SECTION`, consistent edge traction, and output requests |
| `.SVG` | Standalone deformed shape SVG with ghost geometry, colour coded members, and the deformation scale | Standalone von Mises contour SVG with one filled polygon per triangle |
| `.JSON` | Full structural report with nodes, members, safety, determinacy, equilibrium, and modes | Full CST report with mesh metadata, nodal displacements, reactions, element stress tensors, and equilibrium |

Exports are generated in the browser with the Blob API and never leave the machine.

## Files

```
index.html    Static shell and semantic controls
style.css     Blueprint theme, layout, SVG mark specs, and accessibility states
app.js        UI state, SVG rendering, tooltips, table, theme, and exports
fem.js        UMD finite element solver, presets, and truss exporters
test.js       Node regression suite
.nojekyll     Disables Jekyll processing for static hosting
LICENSE       MIT license
```

## Accessibility

Colour never carries meaning alone: force and stress states have text labels, the table view repeats the numeric data, and status colours come with an icon and label. Focusable SVG hit targets are larger than the visible marks and show the same tooltip as pointer hover. Reduced motion preferences are honoured, and the light theme is a separate validated palette rather than an inversion of the dark theme.

## License

MIT. Copyright (c) 2026 draarivpatel-ui. See [LICENSE](LICENSE).
