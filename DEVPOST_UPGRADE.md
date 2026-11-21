# MeshLab: See the physics. Check the math.

## Inspiration

Engineering diagrams often show the final answer while hiding how that answer was obtained. We wanted a small laboratory where a visitor could change a physical input, see the load path or temperature field respond, inspect the numbers, and ask the solver to reproduce its own result. The first interaction should take seconds; the mathematical explanation should remain available as far as the visitor wants to explore.

## What it does

MeshLab opens as a complete browser workbench with four immediate experiments. Trace a bridge to distinguish tension and compression. Overload a crane to identify its governing yield limit. Inspect stress tensors around a perforated plate. Follow the flux from a heater across a Delaunay mesh.

The structural laboratory solves 2D and 3D pin jointed trusses, checks yield and Euler buckling, displays natural frequencies and mode shapes, and solves a plane stress plate with constant strain triangles. The field laboratory solves the Poisson equation with linear triangular finite elements and sparse preconditioned conjugate gradients.

Every member or triangle has live pointer and keyboard inspection. Selecting a result pins its measurements beside the plot. Load and heater sliders recompute the model, including the zero input state. Numerical badges report residuals and energy or heat balance. A round-trip button serializes and reproduces the model; for the plate, it compares the iterative result with an independent direct solver.

## How we built it

The complete application is classic JavaScript, HTML, CSS and SVG. It has zero runtime dependencies, external assets or build steps and works offline, including when opened directly from disk. GitHub Pages hosts the same files that run locally.

Structural bars contribute axial stiffness blocks to the assembled system. Constraints reduce the system to free degrees of freedom, scaled partial pivot elimination solves it, and `K u - F` recovers reactions. Member elongation determines axial force, stress and strain energy. Yield and Euler critical loads give separate factors of safety.

Modal analysis uses the full generalized eigenproblem. A Cholesky factor of the mass matrix reduces both lumped and consistent mass models to a symmetric problem; Jacobi rotations recover the modes, which are then mass normalized and checked against `K phi = omega^2 M phi`.

The plate assembles `B^T D B t A` for each constant strain triangle and distributes edge traction consistently along boundary segments. Its browser result uses sparse Jacobi PCG. The direct solver remains available for independent comparisons.

The thermal experiment creates a deterministic irregular mesh with normalized sweep-hull and Lawson edge-flip Delaunay triangulation. Scalar FEM assembles a sparse conductivity system, moves nonzero prescribed boundary values into the reduced right hand side, and solves with Jacobi PCG. The implementation checks symmetry, curvature, convergence and component boundary constraints. Flux comes from the element temperature gradient.

## Challenges

An audit found that preset loads were present both on nodes and in a top level load array, doubling the physical load. That error could survive equilibrium checks because the doubled model still balanced. Presets now have one load declaration and dedicated tests compare solved totals with the intended inputs.

Consistent mass was another subtle issue: treating its diagonal as the complete mass matrix gave positive frequencies while failing the coupled generalized eigenproblem. Full Cholesky reduction and independent residual and mass orthogonality checks now cover the coupled case.

The interface also needed honest data binding. Live select values had been updated through attributes, the force width encoding was overridden by CSS, modal amplitudes were almost invisible, and repeated selector rebuilding could disrupt focus. Controls now use live properties, force widths retain their data encoding, modal amplitudes are visible and independent of static load, and the animation has a single time based loop with stable bounds.

Additional repairs cover short load arrays, invalid materials and node references, plate load arrays, degenerate and orphaned mesh geometry, 3D SVG projection, consolidated export loads, zero load states, tooltip placement and download URL cleanup.

## Accomplishments

- Four instant experiments connect visual output to a physical question.
- Ten solved presets span axial structures, plate stress and scalar heat flow on four domain shapes.
- Numerical verification is visible while inspecting the model.
- The plate can compare two independent linear solvers in one interaction.
- The thermal mesh reports geometric quality and reproduces an exact affine field when the heater is off.
- Keyboard inspection uses a single diagram entry point with arrow navigation and persistent pinned measurements.
- Responsive SVG, screen-width member hit targets, explicit text states, reduced motion behavior and both theme palettes are implemented together.
- Input decks, standalone SVG and complete numerical JSON reports are generated locally.

## Verification numbers

Reproduced with:

```sh
node --test --test-reporter=spec test.js numerics.test.js ui.test.js
```

The run passes **182 tests with 0 failures, 0 skipped and 0 cancelled**: 102 structural/export checks, 51 numerical edge case checks and 29 DOM/interface checks. All seven JavaScript source and test files also pass `node --check`.

| Measurement | Reproduced value |
| --- | --- |
| Bridge applied vertical force | 360,000 N, applied once |
| Bridge peak translation | 7.1582247 mm |
| Crane showcase at 5x load | Yield FoS 0.86153846; Upper boom 1; one critical member |
| Holed plate mesh | 336 nodes, 580 CST triangles |
| Holed plate peak/nominal stress | 2.919844677 |
| Holed plate PCG solve | 285 iterations; relative residual below 1e-12 |
| Default Delaunay field mesh | 165 nodes, 280 triangles, area 3.84 m2 |
| Delaunay minimum triangle angle | 33.6772 degrees |
| Heater showcase peak temperature | 144.7337 C |
| Affine harmonic patch | Maximum nodal temperature error 3.15e-10 C |
| Annulus radial patch | Maximum nodal temperature error 9.91e-2 C against the exact `ln r` profile |
| Domain Euler characteristic | Rectangle / L / sink V-E+F = 1; annulus V-E+F = 0 (one hole) |
| Manufactured solution `sin(pi x) sin(pi y)` | Second order L2 convergence, measured order near 2.0 |
| Consistent mass bridge mode checks | Relative generalized eigen residual below 1e-8 and mutual mass orthogonality |
| Theme text contrast | Primary, secondary and status text meet automated 4.5:1 panel-background checks in both themes |

DOM verification loads the actual HTML and all browser scripts into a browser contract fixture. It exercises every preset and supported view, all sliders and overlays, keyboard inspection, mode animation, downloads, errors and round-trip checks. The fixture verifies behavior and bindings; it does not measure browser layout or screen rendering.

The plate's concentration is a finite width, stair step mesh result, not an exact infinite plate benchmark. Structural safety and modes use linear elastic axial and plane stress assumptions. Scalar function sources integrate on a three point edge midpoint rule; holed domains clip polygons, so the annulus inner boundary is a polygonal circle approximation.

## What's next

Add boundary-conforming curved-hole meshes and a posteriori adaptive refinement guided by field error, then expose model editing with the same input validation and verification checks. Extend the structural formulation with beam bending before introducing frame claims. Keep every addition inspectable, reproducible and compatible with static browser hosting.
