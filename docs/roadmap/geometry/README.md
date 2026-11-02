# Geometry, mesh, spatial

## Where it stands

Primitives and their intersections and distances (`geometry/`), meshes as the materialiser at the edge (`mesh/`: subdivision at the space's own midpoint, geodesics, transport, Voronoi, tessellation, marching cubes), and spatial indices (`spatial/`: a BVH over any bound, geodesic ball trees).

## Now

- nothing is taken up here now

## Later

- nothing is planned here yet; what was left open is in `--area`

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area geometry`.

## Reasoning and earlier history

### Geodesic algorithms

- MeshTopology: edges, vertex neighbors, face-edge map, boundary detection
- Geodesic distance field: Dijkstra with space.distance() — works on any Surface
- Shortest path: Dijkstra + predecessor trace
- Parallel transport: Schild's ladder (exp/log only)
- Geodesic Voronoi: multi-source Dijkstra + label propagation + face_labels

### Audit fixes, multi-mesh viewer, atom demo

- Marching cubes: Paul Bourke 256-entry tri_table (replaces centroid-fan)
- BVH stack: sized by the tree -- a fixed array of 64 on the call's stack while the tree is at most 62 deep, a vector past that. The line here used to say it had been a vector since long ago; the code had gone back to a bare `array<64>`, and a tree of centroids spread over many orders of magnitude overran it (2026-09-25, `tests/test_bvh.cpp`, caught by AddressSanitizer)
- MeshTopology: shared_ptr<const Mesh> ownership (replaces dangling reference)
- Matrix::inverse() → Result<Matrix> (was silent identity on singular)
- Quaternion::inverse() → Result<Quaternion> (was division by zero)
- identity() morphism: self-inverse (was nullopt)
- OBJ loader: std::from_chars (was std::stoi with exceptions)
- Viewer: multi-mesh rendering, per-mesh color push constants, alpha blending
- Physics: AtomModel, orbital wavefunctions, element database (118 elements), atom SVG
- atom_demo: CLI + console + Vulkan viewer + SVG export + shell toggle
- 437 tests, ~15K LOC

### Eigen interop, heat method, operator fix

- Eigen3 interop (optional, SPATIUM_EIGEN=ON): to_eigen/from_eigen, eigen_view (zero-copy Map)
- Heat method geodesics (Crane 2013): O(h^2) via Eigen::SimplicialLDLT. HeatSolver pre-factored.
- Discrete differential operators: cotangent Laplacian, lumped mass matrix, face gradients, integrated divergence
- GeodesicMethod enum: Dijkstra | Heat, backward-compatible default
- operator| disambiguation: Polygon union → operator+/union_of(), | reserved for intersection/pipe

### Geodesic procedural generation, render/ engine consolidation (2026-08-28)

- `examples/geodesic_procgen_demo.cpp` — the first two of three manifold-backlog ideas raised this session (the third, hyperbolic ray-marched rendering, is in Backlog → Manifold applications): `make_bumpy_sphere()` (new factory, `spaces/implicit.hpp`) → `marching_cubes()` → farthest-point sampling via repeated `geodesic_voronoi()` calls (each new site is the vertex the last call's own distance field says is farthest, so no separate single-source Dijkstra needed) → `face_labels()` → `BVH<Triangle3>::ray_cast()` against the real mesh, rendered through the same `render::Camera`/`parallel_for_rows`/`supersample_pixel`/`write_png_rgb` pipeline `tumbling_body_demo.cpp` established as the pattern. Zero new geometry algorithm — genuinely pure composition, as scoped.
- `render/color.hpp` (`hsv_to_rgb255`) and `render/sky.hpp` (`Sky`/`make_starfield`/`sample_sky_color`/`sample_sky`/`random_sky_dir`) promoted out of `examples/io_helpers.hpp` into the library proper — the same "found needed by another caller → promote to the engine" path `camera.hpp`/`parallel_for_rows.hpp`/`supersample.hpp` already went through, triggered this time by `geodesic_procgen_demo.cpp` needing the color/star-scatter pieces but not the whole-sky gradient. `examples/io_helpers.hpp` now holds only `confirm_overwrite()` (an example-CLI convention, not a rendering primitive).
- Real bug found and fixed, not just moved: `sample_sky_color()`'s background brightness gradient ("brighter near mid-latitude, dark at the poles") is correct for the wide-FOV whole-sky views the 3 GR raytracers use, but at `geodesic_procgen_demo`'s close, narrow-FOV single-object framing only a small patch of that gradient is ever visible — reading as an unexplained dark smear rather than a sky. Fixed with an additive `Sky::wide_sky` bool, default `true` (unchanged branch taken by every existing caller — spot-checked by re-rendering `blackhole_demo`, one frame of `blackhole_gr_demo`, and one frame of `tumbling_body_demo` after the move, all visually consistent with their established look); `false` gives a flat background at `tint` instead, used by `geodesic_procgen_demo`. New tests in `tests/test_render.cpp` cover both branches plus `hsv_to_rgb255` and per-star findability. 822/822 tests green.

### Mesh processing

- **[course]** Mesh simplification (edge collapse + QEM)
- **[course]** Isotropic remeshing on Surface
- **[want]** PLY import/export

### Tessellation density, and the rule that is not the obvious one

Measured 2026-09-22, `rsc/tools/tessellation_error.cpp`. The second
dispatch point after render level, and it is decided today by numbers
placed by hand: 240x120 for the donut's dough, 8x2 for its sprinkles,
48x24 and 12x4 as library defaults.

**For a surface with no displacement the rule is straightforward and
measurable.** A tessellation is wrong when its deviation from the true
surface is visible, so the quantity is the largest gap between the flat
mesh and the real surface, in pixels rather than in triangles. Against the
donut camera's pixel of 5.4e-3 world units, a plain torus 2.0/1.0 runs
123.5 / 37.2 / 9.7 / 2.5 / 0.62 / 0.099 pixels at 4x2 through 240x120, so a
half-pixel target lands near 120x60 and everything finer is paid for and
not delivered. The deviation is exact rather than estimated, since the
surface can project a point onto itself.

**For a displaced surface that rule does not apply, and the reason is
worth keeping.** The donut's dough is a torus under a noise offset, which
is a perfectly well-defined continuous surface -- `offset_surface`
evaluates its thickness field at (u, v), so it is approximated by a mesh
rather than defined by one. But it converges only once the tessellation
resolves the *displacement field's* frequency, not the base surface's
curvature, and the demo does not:

| dough tessellation | against 240x120 |
|---|---|
| 120x60 | 44.45 dB |
| 360x180 | 44.93 dB |

Going up changes the picture as much as going down, which is what
non-convergence looks like. The cause is in the field: the pore noise is
`pore_noise(u * 11.0, v * 4.5)` over a domain of 0..2π, so the argument
runs to 69 and there are about 69 pores across u. At 240 steps that is 3.5
samples per pore against a Nyquist minimum of two, so every density in the
range draws a different bread and none is near a limit.

**So the rule for a displaced surface is set by the displacement field,
and the frequency is already written down in the field itself** -- it does
not have to be guessed, which is exactly what hand-tuning a number by eye
is doing in its absence.

With one catch, and it is the same catch as everywhere else on this page.
The frequency is written down in the *source*, not in the field object.
The dough's `dough_surface` in `examples/donut_demo.cpp` is a single
`ScalarField` holding one opaque leaf, so `u * 11.0` lives inside a lambda
and nothing can read it back out: `field_report()` sees a closure, not a
frequency. Deriving this rule mechanically therefore needs the
displacement written structurally, which needs a noise operation in the
vocabulary -- the "Name the 32 opaque leaves against a vocabulary" row of
the index above. Until then the rule is derivable by a person reading the
code and not by the library, which is a weaker claim than the sentence
above makes on its own.

Not acted on for the demo: its numbers were tuned by eye and the picture
is what ships, so lowering the dough to 120x60 -- 78 104 triangles to
34 904, and the storage ratio from 36.4x to 79.7x -- changes the bread for
a 2.2x saving, which is a look decision rather than a correctness one.

### Geometry

- **[course]** Boolean ops on concave mesh (BSP tree)
