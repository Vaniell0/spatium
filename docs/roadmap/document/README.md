# Document

## Where it stands

What a person writes before anything is computed, a statement a line (`io/document.hpp`): a value, a parameter on an interval, a root with or without a bracket, an integral, a shape. Each reads only the names above it; it is printed back as written and as JSON read through the same reader, and its shapes are built as the scene its calls would build (`io/document_trace.hpp`). A question is a document and has one form, the rows of a LaTeX editor, each a formula as MathJSON read as a line, the last the question (`mathjson_row_text`, `rsc::respond_rows`); every answer gives it back as rows the library writes (`to_latex` of a statement, `question_latex`), which the page's parser reads as the same question. The scene DSL (`io/build.hpp`, a `Trace`) is the same layer for scenes.

## Now

- the steps of an answer written from what was computed, under the rows it was asked as

## Later

- a TypeScript builder that writes the document's JSON
- matrices and systems of equations as statements
- a workspace for documents and scenes, edited as a canvas

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area document`.

## Reasoning and earlier history

### DSL, concepts, Complex<T>, ray-quadric

- Vector UDLs: `3.0_x + 2.0_y + 1.0_z` via expression templates
- Uniform distance(): 7 new overloads (Box-Box, Circle-Circle, Triangle-Triangle, Polygon-Polygon, etc.)
- Polygon boolean operators: `operator&` (intersection), `operator-` (difference), `operator+` (union), `union_of()`
- Generalized pipe-unwrap: `Result<Vec> | Transform`, `Result<Point> | Morphism`, chains
- Lazy transform chains: `lazy(translate) * lazy(rotate)` with `.apply()` / `.collapse()`
- Concept constraints: Shape/Measurable requires in boolean.hpp
- Generic algebra functions: `power()` (Group), `commutator()` (Group), `adjoint()` (LieGroup), `poly_eval()` (Ring)
- Complex<T>: full arithmetic, conjugate, magnitude, phase, from_polar, sqrt, cbrt, constexpr
- Polynomial solvers: solve_quadratic, solve_cubic (Cardano), solve_quartic (Ferrari) → Complex roots
- Analytical ray-quadric: Quadric<T> (sphere/cylinder/cone/ellipsoid), ray_quadric (hits), ray_quadric_proximity (miss via imaginary part)
- 530 tests, ~18K LOC

### Declarative scene DSL: `io::build`, analytic-first, trace-based (2026-09-09)

`spatium::io::build::Trace<T>` — a getting-started demo (a donut, in the tradition of every 3D beginner's first Blender project) built out a real, general primitive along the way, not just example glue. Core decisions, all deliberate:

- **Analytic, not mesh-first.** `spaces/offset.hpp`'s `offset_surface(base, thickness)` composes a *new* `ParametricSurface` from `base(u,v) + thickness * base.normal_at(u,v)` — a function, not a second shape draped/projected onto the first (an earlier attempt built icing as a near-duplicate torus formula just to project it back onto the real one; redundant, and the mistake `offset_surface` exists to make impossible). `spaces/sample.hpp`'s `sample_surface_uniform` places points evenly *by actual surface area* via rejection sampling weighted by the first fundamental form's area element (`ParametricSurface::area_element()`, new) — no mesh, no geodesic-Voronoi graph, for the common case where the target is a `ParametricSurface`. Mesh only gets built at `materialize()`, for display, never as the scene's own representation.
- **Trace, not closures.** `Trace<T>` is a flat `vector<TraceNode<T>>`, each node a tagged operation (`Space`/`Offset`/`Scatter`/`Compose`/`Literal`) with explicit fields, addressed by index — not a tree of `std::function` closures (the first draft's mistake: opaque, un-inspectable, can't be walked or re-interpreted by a second pass later). `Handle<T>` wraps `{trace*, index}`; every factory method (`torus()`, `offset()`, `scatter()`, `cube()`, ...) appends one node and returns a handle. `Literal` is the deliberate escape hatch for shapes with no natural (u,v)→R³ formula (a cube has none) — holds a precomputed mesh directly.
- **One motion slot, not three.** Every node has `.moving(fn(point, t) -> point)` — a constant shift, a `Morphism` pipe, or genuine time-dependence all fit the same signature. The donut's opening act (delete the default cube, Blender-tutorial-style) uses it directly: displacement is `p*(1-smoothstep(t)) + turbulence(p,t)*smoothstep(t)`, a pure function evaluated wherever `materialize()` asks for it — no per-frame simulation state. (First attempt pushed vertices *outward* along their own direction for "explode" — on a mesh with fixed face connectivity that just inflates the solid into a blob, doesn't shatter it; collapsing inward reads as "deleted" correctly with zero per-fragment mesh work.)
- **Lazy, deliberately not reactive.** Nodes are immutable values; building one does no geometry work, only `materialize()` walking the trace does. No dependency-tracked incremental recompute (SolidJS/Vue-signals style) — that engine cost has no consumer yet in this scope. If a reactive layer gets built at all near-term, `gpu/derive_christoffel.py`'s declarative-derivation → auto-codegen pipeline (GPU rendering section below) is the more natural home for it than the scene graph.
- **New library primitive, reused nowhere else yet on purpose:** `algebra/noise.hpp`'s `PerlinNoise` (Ken Perlin's 2002 "improved noise", seeded — no shared mutable state, so each use gets an independent field) replaces the ad-hoc per-particle RNG velocity every hand-rolled particle effect in this tree would otherwise reach for.

**Real bug found and fixed while rendering it, not just tuning:** `examples/donut_demo.cpp`'s `--photo` path mixed two color scales — `render::sample_sky_color()` already returns [0,255] (matching every other offline demo's direct pixel write) but the new recursive raytracer's outer wrapper multiplied *everything* by 255 again, including sky contributions blended in through glossy reflection rays. Glossy icing reflecting the sky compounded the double-scale on top of itself, blowing the whole icing surface out to solid white. Fixed by keeping `trace_ray()` in a single [0,1] convention throughout (converting `sample_sky_color`'s result down once, at the point of use) and doing the one [0,255] conversion exactly once, at the outermost call.

`--photo` is now a real small raytracer, not a flat raycast: `Material<T>` gained a `roughness` field (1 = matte/Lambertian, the old implicit default; 0 = glossy) driving a Blinn-Phong specular lobe plus an actual recursive reflection ray for glossy surfaces (depth-limited to 3) — dough is rough, icing is glossy, and the difference is a real light-transport effect, not two different `base_color`s. Per-vertex smooth normals (area-weighted face-normal average, interpolated at the hit point via the BVH's own barycentric weights) replaced flat per-triangle shading, which was reading as faceted plastic on the demo's coarse-ish tessellation.

**Not done, explicitly deferred:** a live Vulkan window (CPU raytraces, Vulkan only presents the frame) — real new plumbing (`viewer::App` has no sampled-texture descriptor/pipeline path today, only mesh/point-cloud rasterization), too open-ended to land blind against a same-day deadline; console output + `--photo` PNG are the two supported modes for now, see Backlog below. Scattered items still sit slightly proud of the surface rather than lying fully flush (a small tangent-plane fix, not attempted this round).

1001 tests (was 989), 12 new (`tests/test_build_dsl.cpp`) covering `PerlinNoise`, `offset_surface`, `sample_surface_uniform`, and the `Trace` DSL's node-kind recording, materialization, and `.moving()` determinism. One pre-existing unrelated failure (`svd cross-validated against Eigen::JacobiSVD` — confirmed via `git stash` to already fail on a clean `main`, an Eigen `JacobiSVD` assertion about thin U/V needing a dynamic-column matrix, nothing to do with this work).

Docs: [`docs/getting-started-dsl.md`](getting-started-dsl.md) (new, zero-barrier-to-entry walkthrough), `docs/api-reference.md` and `CLAUDE.md` updated for every new primitive above — and, caught in the same pass, both were already stale against several *prior* PRs (`terminal_donut_demo`, `scene_demo`, `sound_synthesis_demo`, `spd_relationship_demo`, `cloth_sphere_probe`, and all of `discrete/combinatorics.hpp` were shipped but undocumented) — fixed together, not just today's additions.

### Object model as manifold substrate

Unlocked by the object→exact-`Surface` bridge (extending `geometry/surface_adapter.hpp`'s `ShapeSurface`, in flight as `feat/scene-object-model`): any renderable scene object exposing a real Surface, not an approximate tangent-plane one, means every manifold operation in the library (geodesics, Voronoi, parallel transport, DEC, contact, Riemannian optimization) becomes available on it for free instead of needing per-shape hand integration. Brainstormed 2026-09-06, honestly triaged (not everything below is equally close):

- **[course]** Recursive procedural world generation — geodesic Voronoi/Delaunay cell on a surface, one object per cell oriented by the normal, recurse per-object for a fractal scene. Mostly composition of what already exists (`mesh/voronoi.hpp`'s geodesic Voronoi, `Point`/`Morphism` pipes) once the object bridge lands; not new math.
- **[course]** Manifold-native RL as a new RSC domain — agent state = point on a `Surface`, action space = tangent space, transitions via `exp_map`. Fits RSC's existing REINFORCE-dispatcher methodology directly; also the concrete way to eventually fill the "geometric/manifold-constrained control theory" and "real-time control" candidates raised in the 2026-08-31 industry-directions session.
- **[course]** `GeodesicPath<S>` as a first-class 1D `MetricSpace` — a geodesic curve carries its own arc-length parametrization, so points/distance/Voronoi work along a route. Real routing use case: waypoints along a path, optimal stop placement.
- **[want]** Geodesic-guided growth simulation (cracks/roots/vessels) via existing parallel transport + distance fields — doesn't need the object bridge at all, already buildable on `Sphere`/`Hyperbolic`/`ParametricSurface` directly. Gallery-tier demo, not new library code.
- **[want]** Manifold-native game movement (walk on an asteroid, "up" = surface normal, Mario-Galaxy-style) — same status as above, a composition/demo of existing primitives (`Surface`, `exp_map`, `normal`, `render::Camera`), not a new subsystem. Sharpest concrete version, 2026-09-09: a `Hyperbolic<3>` dungeon crawler — exponential volume growth with distance from center makes rooms feel unbounded and corridors branch fractally, and since collision/physics/traversal already go through `exp_map()`/`log_map()`, enemies walking geodesics and light bending through the Poincaré-disk geometry falls out for free rather than needing hand-written non-Euclidean-movement logic. Worth building specifically because it's the cheapest of the "moonshot" ideas raised this session and it substantiates a claim the HN launch copy already makes ("build a game where up is whatever surface you're standing on") rather than leaving it as an unproven line.
- **[want]** Multiscale geometry pipeline (scan → mesh → `Quadric`/parametric fit → geodesic segmentation → recurse per segment) — not a separate idea, a continuation of the already-noted point-cloud → mesh reconstruction industry candidate (2026-08-31 session): reuses `ImplicitSurface`/`marching_cubes`/geodesic Voronoi, the new piece is only the fit-then-recurse loop.
- **[dare]** Cloth/shell dynamics intrinsic to a curved surface — solving the material's own equations of motion in tangent space accounting for curvature, distinct from the already-shipped "drape cloth onto a fixed obstacle" work (`examples/cloth_sphere_probe.cpp`, Contact physics section above). Real research-level lift, not a quick extension of XPBD.
- **[dare]** Topology optimization on an evolving surface (solve a stress/heat PDE on the current surface, remove material, iterate on the new one) — needs real elasticity FEM, well beyond the current DEC/heat-method infrastructure (which only covers scalar diffusion, not stress). Exciting, not scoped.
- Explicitly NOT a new subsystem, just a consequence of the type system once the object bridge lands: composing nested local frames (a point on a robot → a point on one of its sensors → a point on that sensor's own chip surface) through the same `Point<Space>` vocabulary at every level. Worth a demo to show it falls out for free; not separate work.
- **[want]** Reconstruction pipeline for real-world objects — photos/video → geometry → `ImplicitSurface` → existing Spatium post-processing (smoothing, simplification, geodesics, self-intersection checks, OBJ export). Corrected 2026-09-15: this used to end "Not researched yet" and carry `[dare]`. It has been researched, and the research changed the shape of the item. The learned-latents-to-SDF-decoder route (a JEPA-style encoder feeding a small SDF decoder) was tried and abandoned — the decoder produced a blob, not a surface. A frozen VGGT-1B was then validated locally with no training at all and produced a recognizable mesh from a DTU scan. **None of that work is in this repository** — no directory, no reference, nothing to build on — so what is actually open here is small and specific rather than a research programme: a load path in `io/`, a conversion utility, and somewhere for it to live. The reconstruction itself is a solved external step.

### Declarative scene DSL

#### Direction, settled 2026-09-14

Where this is going, so the individual items below read as one line of work rather than a pile of improvements.

**The authoring layer and the execution layer are different layers, and only the first one may be open.** Spatium's whole argument is that a space is an extension point, and the same should hold for scene operations — a user with a better geodesic solver or a new kind of node should not have to edit our headers. But openness is bought with indirection, and indirection is exactly what cannot be lowered to a GPU kernel or fused across operations. The resolution is not to pick one: type erasure, registries and anything else ergonomic belong to *describing* a scene, which happens once; *evaluating* it must be monomorphic, either at compile time or compiled from an IR. This is why the trace was built as a flat record of tagged operations rather than a tree of closures in the first place — closures cannot be lowered, data can. The trace is the IR.

This is not a novel architecture and should not be designed as though it were. Dr.Jit (Mitsuba 3's backend) traces C++ operations into an IR and JIT-compiles to CUDA/OptiX/LLVM for exactly this class of problem; Taichi captures an AST and compiles to CPU/GPU kernels; Halide separates the algorithm from its schedule; Kokkos gets one source onto both CPU and GPU through lambdas and execution spaces, and forbids `std::function` inside kernels for the reason above. What would be ours is the particular IR over Spatium's own operations, not the shape of the solution.

**Measured, not assumed** (`benchmarks/bench_trace.cpp`, on the shape `examples/donut_demo.cpp` actually builds: 19 800 single-mesh nodes, 12 vertices each, one motion hook per node). Three plausible performance hypotheses were tested and all three came back smaller than expected, which is the reason this section exists in this form:

| suspected cost | measured |
|---|---|
| the per-vertex indirect call | ~3.3 ns/vertex — about 4% of a realistic motion callee, invisible unless the callee is trivial |
| 9.7 MB of duplicated `PerlinNoise` tables captured into 19 800 closures | ~3% of frame time |
| materializing the nodes with no motion hook at all | **~20% of frame time** |

So the dominant structural cost is none of the things that looked like the problem: it is that the same 12-vertex mesh is stored and copied per instance, because the trace has no notion of instancing. `std::move_only_function` was still worth adopting, but for what it enables (hooks owning move-only or shared state, const-correct invocation through the `const Trace&` that `materialize()` holds) rather than for speed, and the commit says so.

The order that follows from this: instancing and the structural field representation together, then the full set of edge rules that fields make expressible, then opening the node set, then time. Revised 2026-09-15 — the field work moved ahead of the edge rules, because `EdgeRule`'s remaining values carry data (a cap radius; a reference to another node for "extend to this surface"), so they cannot be enum values and the rule set cannot be completed until fields exist. That argument is worth more than the original one: caching and GPU lowering will be true later, whereas a value with nowhere to put it is true now.

#### Chart versus manifold, recorded 2026-09-15

Splitting `offset()` from `offset_shell()` needed a closedness test, and writing it surfaced something larger than the function.

`periodic_u()` / `periodic_v()` describe the **chart**, not the surface. A sphere's chart is the rectangle [0,2π]×[0,π], and as a subset of R² that rectangle has an edge — which is what `periodic_v() == false` honestly reports. But the chart's `v` edges map to the poles: two *points*, not two curves. The surface has no boundary; the chart does. So a closedness test cannot read the flags, because the flags answer a question about the parametrization while the question being asked is topological. `is_closed()` therefore looks at the geometry — a direction closes by being periodic, or by both its edge curves collapsing to a point.

The general form of this, which is what `SurfaceWithBoundary` will actually need, is not "does the chart have an edge" but **what the chart's edge maps to**, and there are three answers:

| the chart's edge maps to | meaning | example |
|---|---|---|
| nothing (the map is periodic) | closed, seam only | torus, in both directions |
| a point | closed, pole — the surface passes through, it does not stop | sphere, in `v` |
| a curve | a genuine boundary | a band cut from a torus; a tube's open ends |

Only the third case needs a rule at all, and `EdgeRule::ZeroThickness` is one answer to it, not the general one. The first two need no rule and must not be made to ask for one — which is exactly the trap the pre-split `offset()` fell into by having a single operation cover all three.

`is_closed()` today collapses this to a yes/no because that is all `offset()` needs. The three-way classifier is what the boundary concept should be built on, and it should be built before the concept, not after.

#### Render levels, settled 2026-09-15 — and what the number actually says

A node now carries a `RenderLevel`: `Exact` (closed form), `Tessellated`
(triangles in a BVH), `Newton` (Newton's method on the `(u,v)` map).
Inferred in one place — an exact form when the node has one, tessellation
otherwise, and **Newton never inferred**, because at three orders of
magnitude nobody should pay it by accident. `.rendered_as()` overrides
and refuses on the spot a level the node cannot serve.

**The argument for it is build cost, not per-ray cost.** This was got
wrong twice on the way here and is worth stating flatly, because the
per-ray numbers point in whichever direction the leaf ratio happens to
push them:

| scene | exact leaf vs triangles | per ray | build |
|---|---|---|---|
| 19 800 specks | quadric vs 396k tris | exact **20% slower** | 11.8 ms vs 291 ms |
| one dough | torus vs 25 600 tris | exact **2× faster** | ~0 vs 10.8 ms |
| 64 tori, through the DSL | torus vs 147 456 tris | **a wash** (44 vs 49 ns) | 0.011 ms vs 78.8 ms |

The dough's 2× was read at the time as the speck result reversing. It
was not: one leaf against 25 600 is an extreme ratio, and adding 63 more
tori brings the per-ray cost back level. Nor is it true, as a related
guess had it, that the exact path's advantage grows with scene size in
*traversal* — the third row is the test of that and it came back flat.

What survives all three is that **the tessellation never exists**. Seven
thousand times the build cost in the DSL row, ~8× the memory in the
speck row. An animated scene pays the build on every frame while the
per-ray cost is paid once per pixel, and that is the whole case.

**Measured 2026-09-22: those three rows are three points, and there is a crossover between them.** The benchmark's tori count was a constant, which is why this entry had anecdotes from three different scenes and no boundary. Parameterised, per ray:

| tori | exact | tessellated |
|---|---|---|
| 1 | 22.3 ns | 29.8 ns |
| 16 | 40.3 | 42.1 |
| 64 | 57.4 | 57.5 |
| 256 | 76.5 | **61.1** |
| 1024 | 95.1 | **67.8** |

The crossover sits near a hundred objects, and the reason is in the growth rates: exact rises 4.3× across the sweep against tessellation's 2.3×. A ray-torus test is a quartic solve where a ray-triangle test is not, and a BVH saves a logarithm of the *number* of leaves, never the cost of one.

Build across the same sweep: 0.011 / 0.052 / 0.338 ms exact against 80.6 / 391 / 1822 ms tessellated, at 64 / 256 / 1024. So tessellation's per-ray saving pays its build back in **2.3 frames at 256 tori and 6.0 at 1024**, against the demo's own frame of 960×720 at 16 samples per pixel — 11.1 M rays.

**Which makes the rule in the code measurably wrong for large static scenes.** `Placed::render_level()` returns `exact.has_value() ? Exact : Tessellated` — a decision by availability. For a thousand static objects rendered for longer than six frames, availability picks the slower path.

Memory stays on exact's side and does not go away: 1024 tori is 2.36 M triangles, on the order of 170 MB against roughly a hundred kilobytes of closed forms.

So the decision has **three** features — object count, frame count, memory budget — and a rule reading only one of them cannot be right by construction. That is a dispatch problem with numbers attached, which is what `rsc/` is for and has never been pointed at.

A finding from the demo, kept because it is the mechanism working:
`donut_demo` reports `0 exact, 19 809 tessellated, 0 newton`. The dough
is an `Offset` carrying a noise bump, so it stopped being a torus the
moment it got bread texture. An exact form is a promise about the shape,
and a bumped torus cannot keep it.

#### `Chart`, settled 2026-09-15 — the requirement the operations already had

The DSL was open to new *operations* and closed to new *spaces*, and the
reason turned out not to be the one written down. It was not that `Kind`
is a closed enum, and not that the node type is fixed. It was that the
requirement every DSL operation shares had no name, and travelled through
the signatures disguised as the concrete type `ParametricSurface<T>`.

The claim is checkable, which is why it settled the design. Reading the
three operations for what they actually call:

| operation | what it calls |
|---|---|
| `offset_surface()` | `evaluate`, `normal_at`, `domain`, `periodic_u/v` |
| `sample_surface_uniform()` | `area_element`, `domain` |
| `parametric_mesh()` | `evaluate`, `domain` |

Not one of them calls `project()`, `normal()` or `exp_map()` — so not one
of them uses the `Surface` concept. Generalizing `resolve_surface` to
`Surface`, which is how this item was written for a week, would have
handed the operations methods nobody calls while losing every method they
need. The right name is `Chart`: a (u,v) parametrization, which is a
choice of coordinates *on* a space, not a kind of space. It is
deliberately not part of the `Set → Manifold → Surface` hierarchy — a
space can carry many charts or none.

`chart_of(space)` is the extension point, found by ADL, the same shape as
`point_to(p, surf)` in the contact concepts and for the same reason:
compile-time overload resolution, no runtime registry, no indirect call.
`ParametricSurface<T>` gets the identity overload and remains the erased
form every chart converts into, because a trace holds nodes of one type
and cannot be templated per node.

`Sphere<2, T>` gets its **first** chart here, not a generalization of an
existing one — as a `RiemannianManifold` it has no `evaluate(u, v)`
anywhere, which is the real reason it could never enter the DSL. Note
which sphere: `Sphere<N, T>` is the N-sphere in R^{N+1}, so the ordinary
sphere in R³ is `Sphere<2, T>`. `Sphere<3, T>` lives in R⁴ and gets no
overload, which is not an omission — there is no chart of that shape to
write, and `Chartable<Sphere<3, double>, double>` is `false` on purpose.

The chart's `periodic_v` is `false` and its `v` edges are the two poles:
the middle row of the chart-edge table above, and the reason `is_closed()`
reads geometry instead of these flags. The poles are also where the chart
degenerates — `area_element` is `r² sin v`, zero there, and the
finite-difference `normal_at` with it. That is a property of the chart,
not a defect of the sphere.

A sphere node deliberately leaves the `exact` slot empty. That slot means
"a renderer can hit this without triangles", nothing in the library can do
that for a `Sphere` today, and filling it would make `render_level()`
claim `Exact` for a shape no renderer can hit exactly — a picture right
about the shape and wrong about the scene. A `sphere()` factory recording
`BoundedQuadric` is available later and is a different change.

**`Kind` stays closed, and adding a space does not go near it.** Its five
values — Space, Offset, Scatter, Compose, Literal — are *operations*, not
shapes; a sphere is a new chart on the existing `Space` node.
`Kind::Literal` already covers shapes with no (u,v) map at all, so the
shape set was open through the chart or through `Literal` before any of
this. The single thing that would open `Kind` is **CSG boolean on
analytic surfaces** — `union(torus, sphere)` as an operation rather than a
mesh merge, whose result is a new analytic surface not expressible as any
of the five. It is not planned: `geometry/boolean.hpp` does this for
polygons, not surfaces, and there is no consumer. Recorded here so the
question is not reopened from scratch in a month.

#### A factor as an expression, and the two ops nobody costed (2026-09-23)

`scaled()` and `rotated()` now take a `Field<T>` as well as a callable,
and a node stores that expression beside the closure rather than instead
of it -- the closure stays the evaluation path, so nothing in `eval`, in
`is_placement()` or in `placement_at()` learned a second way to read a
factor. What the expression buys is not speed. It is that a factor in a
field can be read, exported and rewritten, and a factor in a closure can
be none of the three.

The convention is a choice and not a consequence, so it is written down
here as well as in the header: **a factor reads time as the field's first
parameter, and its second is zero.** `Field` names its inputs "first" and
"second" rather than u and v precisely so they can be bound to something
that is not a surface domain.

**The part that was missed in the estimate, which is the part worth
keeping.** The plan said these factors were "arithmetic in pure form" and
predicted 23 of 32 leaves would convert. The actual factors are
`smoothstep01`, and a smoothstep is `std::clamp` followed by a cubic: the
cubic was always expressible and the clamp was not, because the
vocabulary had no `Min` and no `Max`. With the overload alone, exactly one
leaf converts -- the table's constant yaw. The claim had been checked at
the level of "what kind of function is this" and not at the level of
reading it.

So `Min` and `Max` are now ops. They are the right kind of addition --
total on every input, no branch for a consumer to model, one native
instruction on every target an export would aim at -- and `clamp` and
`smoothstep` are built on them rather than being ops of their own.

Measured before and after, same scene, same configuration:

| | before | after |
|---|---|---|
| structural fields | 77 / 107 | **98 / 107** |
| opaque fields | 30 | **9** |
| opaque leaves | 32 | **10** |
| distinct closure types | 11 | 8 |
| placements / deformations | 33 / 0 | 33 / 0 |

The rendered frame is byte-identical to the baseline built from the
pre-change sources, which is the check this edit had to pass: it was
supposed to change what the library can *read*, and nothing else.

Two factors stay closures and neither is a shortfall to be fixed later.
The cube's `time < 0.12 ? 1.0 : 0.0` is a hard step, which needs a
comparison and therefore a decision about booleans in the IR. The dust's
tumble picks its axis through an integer hash and a modulo, which no
widening of an arithmetic vocabulary will ever express -- it is a leaf on
purpose, and `opaque_of_time` exists to say so.

One thing removed rather than translated: `dust_shrink`'s
`if (t < T_HOLD) return 1.0` guard. Below `T_HOLD` the argument is
negative, the clamp takes it to zero and the smoothstep with it, so the
branch and the expression agree -- checked over [-1, 5] at two million
points, difference exactly zero, rather than argued.

#### Open items

- **[course]** Classify what a chart's edge maps to — point, curve, or nothing — per direction, generalizing `is_closed()`'s yes/no. See "Chart versus manifold" above. This is the thing `SurfaceWithBoundary` (`boundary_distance`, `is_on_boundary`, what `exp_map` does past the edge) should be designed on top of, and it is cheap enough to have before there is a second consumer for the concept itself.
- ~~**[course]** A Debug job in CI.~~ — **done**, PR #34, 2026-09-15, and the entry stayed listed as open for the rest of that day. The reason it was worth doing, kept because it is the argument and not just the outcome: every job that ran tests built Release, so `NDEBUG` removed `assert` and a whole class of defect — misusing another library in a way that library catches with an assertion — passed green. Not hypothetical; it is how the `JacobiSVD` thin-U/V misuse survived, and how "the test passes now" got written down as "the old failure note went stale". Struck through rather than deleted, same as the `scatter()` orientation entry, because an item that outlives its own completion is the failure worth seeing.
- **[want]** `EdgeRule::ExtendTo(node)` as a structural node holding a trace index, exactly like `Offset::base` already does — independent of the field work, since a reference to another node is topology rather than a per-instance parameter. `RoundCap(radius)` is the one that genuinely needs a parameter and should not drag `ExtendTo` along with it.

- **Two boundaries on the scene work, written down 2026-09-16 because both were about to be treated as "while we are in there".**

  **Light, shadows and textures are edits; time is not.** A light is a field on `Material`, a shadow is a second ray in a recursion that already exists, a texture is a field plus a sampling call — three changes inside files that are already there. **Time is not in that list.** Animation already exists through `.moving()` and is enough for a cube that crumbles; "time as a property of space" is the series item about the metric — gravitational dilation, an object's own clock — and it needs a scene where the difference is *visible*, or there is nothing to check it against. Folded into a render pass it would become another multiplier on `t`, which is precisely the wrong shape.

  **Sound from geometry is gated on sparse linear algebra, not on the scene.** `io/wav.hpp` writes and `wave_string.hpp`/`wave_membrane.hpp` integrate, so the ends exist — but resonant modes of a surface are the generalised eigenproblem `L φ = λ M φ` on sparse matrices, and sparse is still `[course]` above with no native default at all. There is a narrower route that does not need a new sparse eigensolver: shift-invert Lanczos over the cotangent Laplacian and mass matrix in `mesh/differential.hpp`, reusing the `SimplicialLDLT` factorisation `HeatSolver` already performs. That is still its own numerical task, and it is `SPATIUM_EIGEN=ON`-gated, so a default build would not have it. Until then the honest options are a recorded sample — which defeats the point — or the wave equation on a 1D/2D domain, which is real physics but is not the surface, and calling it "the cube's sound" would be the same overclaim this file keeps catching.

- **[course]** **An IR evaluable in `Dual<T>` gives exact normals in `normal_at` and `offset_surface`.** A goal in its own right, not a by-product of structural fields, and recorded separately so it is not lost when someone later asks what fields were for.

  `ParametricSurface::normal_at` computes its normal by **finite differences** (`spaces/parametric.hpp:150-155`, step `1e-6 × extent`), and `offset_surface` is built directly on that normal — so every offset surface in the library, the donut's icing included, carries a differencing error that nothing currently removes. `Dual<T>` satisfies `Scalar`, so a field tree generic over its scalar evaluates in `Dual<T>` and yields an exact derivative with **no differentiation rule written per node**: the derivative comes from the arithmetic, not from a second implementation of the tree.

  What makes this worth its own line: it is the library used on itself. The IR exists to make scenes lowerable, and the same IR makes the library *more accurate* than it is without one. That is a different kind of argument from "20% of a frame", and a stronger one — an accuracy result does not depend on the scene's size.

- ~~**[course]** Instancing in the trace~~ and ~~**[course]** Structural fields~~ — **the library half landed 2026-09-16**, merged as #42 (fields), #45 (`cook()`), #44 (the instance leaf). What that means concretely, since "instancing" turned out to name three separable things:

  `Field`/`VecField` are flat pools of tagged ops with an opaque callable as a *leaf*; `cook()` turns a `Trace` into objects with deduplicated shapes; `Instanced<S>` is a reference to a shape plus a placement, so a BVH leaf is one test rather than N. Measured: 2 shapes for 19 801 objects and 138× fewer vertices held; 16.1× BVH build, 21.6× memory, 1.08× traversal.

  **What decides instanceability is not opacity but whether a motion reads the point.** `A(t) + p·s(t)` with both terms opaque is still a placement, because cost is per vertex and those terms are per object — so `PerlinNoise` never had to become an expression. The donut now reports `35209 placements, 21 deformations`.

  ~~`19804 placements, 9 deformations`, with the nine being the exploding cube's fragments, which genuinely deform.~~ — **both halves wrong, corrected 2026-09-17.** The counts rotted, which is routine. The explanation was wrong when it was written, which is not: there are no cube fragments. The cube is a single node whose motion is `p * (time < 0.12 ? 1 : 0)` — a uniform scale to nothing — and the explosion is the *dust*, which is separate `flake()` nodes with their own motions. Nothing about the cube deforms.

  What the 21 actually are, printed by the demo rather than asserted here, because a claim about this has already been wrong once:

  ```
  refused nodes: 1xLiteral 2xOffset 18xScatter
  ```

  The `Literal` is the cube, the two `Offset`s are the dough and the icing, the eighteen `Scatter`s are the sprinkle groups. The table is a `Literal` too and is *not* refused, because its motion is written structurally. Every refusal is over spelling, not over deformation.

  That also reconstructs the historical nine exactly: six `Scatter` (one per colour, before the bands were split three ways) plus two `Offset` plus one `Literal`. **The number was right and the explanation was wrong from the day it was written** — there were never any cube fragments; the cube is one node and the explosion is the dust.

  A number going stale is caught by re-running. A *reason* going stale is caught by nothing.

  ~~**Still open, and it is the part that moves a number:** no renderer consumes any of it.~~ — **wired 2026-09-16.** `donut_demo` builds a `BVH<Instanced<BoundedQuadric>>` beside its triangle tree and takes the nearer hit; `.moving()` stopped dropping the exact form for a *placement* (a deformation still drops it), which is what made the dust eligible at all. Measured on the same scene with only the render path toggled, at `t=1.5` where there is dust to see: **1.95 s instanced against 2.82 s tessellated, 1.44×, 31% of the frame.** The demo's own line now reads `render levels: 19800 exact, 9 tessellated`, where it had read `0 exact, 19809 tessellated` since it was written.

  Two measurement mistakes are recorded with it, because both were made here and both are easy to repeat. The first: `t=3.9` is a frame in which the dust has already dissolved, so a "56% of the frame" taken there measured traversal of geometry that is not in the picture. The second: comparing 2.02 s against 2.47 s changed the speck's *shape* and the render path at once, which measures neither.

  **Rotation in a placement, 2026-09-17.** A placement was a translation and a uniform scale, so every instance of a shared geometry sat not merely in the same pose but in the same *orientation*. Invisible while a speck was a ball; impossible to miss once it is a flake. `VecOp::Rotate` joins `Scale` as a second unary op — affine in the point, so a node that turns keeps its exact form and stays instanceable — and `Placement`/`Instanced` carry an orthonormal `R` beside the uniform `s`. The pair stays a pair rather than collapsing into one 3×3 because that is the contract: `Rᵀ` and a uniform scale are exactly what let a ray test enter local space without rescaling `t` or repairing a normal.

  **It also closed a disagreement nothing had noticed.** `materialize_mesh()` has oriented scattered items since the sprinkles complaint — each site gets its tangent frame plus a deterministic in-plane turn — while `cook()` kept only `site.position` and dropped the frame. Two renderings of one scene, differing on every scattered object. The reason it survived is plain: **`cook()` and `Instanced` had no tests at all**, so nothing ever compared the two paths. There is now a vertex-for-vertex agreement test, and it fails on 432 of its 443 assertions against the old code. The same fix corrects a second latent error in the same expression — the site position used to be added *outside* the node's scale.

  Raising the particle count is the next thing this unblocks; the frame is still linear in count (19 800 → 1.94 s, 59 400 → 4.27 s), but a flake is now one leaf rather than a dozen triangles.

- **The reasoning under those, kept because the decisions are still load-bearing.** The sub-sections below were written while the work was open; they record why each call was made and stay valid now that it has landed. A `Field` is an expression tree in which an opaque callable is a *leaf* rather than an alternative to the tree: under "tree or callable" a user picks a representation at the declaration site and silently forfeits caching, lowering and hoisting, whereas under "tree with opaque leaves" composing a structural field with an opaque one keeps the outer structure and loses only the one leaf.

  **Identity without names, measured 2026-09-15** (`benchmarks/bench_field_identity.cpp`). Sharing anything between fields requires deciding when two are the same function. The demo installs **39 611 field instances from 6 distinct closure types** — 19 800 dust motions and 19 800 dust colors each come from a single lambda written once inside a loop. By address they are 39 611 different functions; the discriminator that sees them correctly is `typeid` of the closure type, which is compiler-generated and so costs the user nothing and opens no door to registering or versioning function *names*.

  **The deduplication rule, and why it is narrower than identity.** Type identity alone is sound only for a *stateless* closure, where every instance really is interchangeable — `std::is_empty_v` decides this without asking the user, and the demo's 8 `grow_scale` instances are exactly that case. For a capturing closure the state that makes instances differ is hidden inside the object, so deduplicating on type would silently fuse 19 800 particles onto one trajectory: a correct answer under the assumption that the capture is empty, wrong otherwise, which is the same shape of defect as a NaN walking through a filter phrased to skip rather than to keep. So: deduplication is promised for stateless leaves, and for parameterised leaves on `(type, params)`; **a capturing leaf with empty params is opaque and is never deduplicated.** Wanting deduplication means putting the data in `params`.

  `params` is a flat POD map — numbers, strings, node references, arrays — and deliberately not nestable. Nesting would be a second IR living inside `params`, which is precisely how an IR becomes a language that has to be learned.

  **Where names stop being optional.** `type_index` does not survive serialization: once a scene reaches a browser as JSON there is no closure type anywhere in it, and none can be reconstructed. `io/build_json.hpp` is therefore the first real consumer of function *names*, not a hypothetical second contributor, and the date is already on the calendar even though it is not today. Recorded here so it does not get re-litigated from scratch: names are refused now for want of a consumer, not on principle.

  **Representation: a flat pool in topological order, not a linked tree. Settled 2026-09-15.** A `Field` is logically a tree; its representation is a `vector` of tagged ops addressed by index, children before parents, with a `Field` being a *range* in that pool. The reason is not the cheap structural hash and not the cheap deduplication — those are consequences. The reason is that a topologically ordered array evaluates in **one linear pass**: a loop over indices in which each node reads children that are already computed. No recursion, no explicit stack, no branching on depth. That is the shape a GPU can run, and it is the shape a linked tree cannot provide without a flattening walk first — so building the tree first would mean writing the flattening anyway.

  This buys the same question `Trace` already answered one level up: what is substitution, `f(g)`? In a linked tree it is a pointer swap, local and free. In a pool it is either a copy or shared structure. The answer is the one `compose()` already uses for nodes — **append, never mutate**: substitution pushes new ops onto the end of the pool and returns a new range, leaving every existing range valid. One principle applied at two levels, the same way type erasure and `chart_of` were.

  **The report must separate recognized from unknown.** A tally that says "39 600 fields, 1 distinct type" is indistinguishable between "they really are all one closure type" and "they all fell into one unknown bucket", and that is the class already named in `conventions.md` — a value wearing the costume of an answer. `field_report()` reports three numbers, not one: recognized, unknown, total. This is not hypothetical; a probe written while measuring the report's cost produced exactly that ambiguous line, counting `move_only_function` fields it could not classify alongside `std::function` fields it could.

  **Cost of the report, measured 2026-09-15.** A full walk classifying every field costs 2.8–5.1 ns/field — 112–204 µs at the demo's ~39 600 fields, against an 80.7 ms frame, i.e. 0.14–0.25%. So `field_report()` is called per frame and needs no incremental accumulation; the concern that an O(all fields) report would be too expensive to call, and would therefore die unused like the orphan mesh headers, does not survive the measurement.

  **Why `typeid` has to be captured by `Field` itself.** `PointField` is today a bare alias for `std::move_only_function`, which unlike `std::function` exposes no `target_type()` — so a motion field in the current tree cannot be counted by type, let alone deduplicated. That is a fact about an alias we do not own, not a property that blocks anything: the fix is `typeid(F)` captured in `Field`'s constructor, one line, and it works because a closure type is one type for every instance of one lambda (see the identity measurement above). It is recorded here rather than done because there is no `Field` constructor yet to put it in — this item is what creates one.

  **Ownership, settled 2026-09-15 and worth stating because the obvious phrase is wrong.** `Field` owns its pool: it is a value type, copies are deep, and `+=` appends into its own vector where nothing else can observe it. "Append, never mutate" is the right phrase for a *shared* pool — a `Field` that is a `(shared_ptr<Pool>, range)` and must keep other fields' ranges valid — and carrying it over to the owning version describes a guarantee the code neither needs nor gives. Here append-only protects exactly one thing, topological order, because a node can only name indices that already exist. If substitution across many fields later wants the shared form, the aliasing guarantee has to be added and proved, not inherited from the phrase.

  **Equal ops are merged on append, and lowering drops what the root does not read. Settled 2026-10-22.** Appending an op first looks for one with the same code, children and immediates and names that instead, so a repeated subexpression is one op; a field appended to another (`+=`, the binary operators, `noise`) is renumbered through a table, since a merged op lands at no fixed shift. Order and append-only are unchanged (a merged op is an earlier index) and the root stays the last op. What is compared is the whole meaning or nothing: constants by bit pattern (`-0.0` and `0.0` stay two ops, and so do NaNs of different payloads), noise by its 512-byte table and gather by every point of its table, and an `Opaque` leaf is never merged, since a callable cannot be compared and two closures of one type can capture different state. The lookup is a flat open-addressing table of op indices inside the `Field`, built on the first append and caught up by scanning only the ops added since, so a copy costs one `vector<uint32_t>` more (a node map would cost an allocation per op) and building by `+=` stays linear. `Field::reachable_from_root()` returns the ops the root reads, in order, renumbered; `lower()` calls it, so the POD interpreter and the GLSL emitter never see a dead op (the builders leave none; a pool assembled elsewhere can). `VecField` is untouched: its ops are a few per node and its factors are `Field`s, which merge. Measured on the donut scene (`tests/test_field_consing.cpp`, 42 scalar fields): 3 557 ops append-only against 1 019 merged (-71%), the 130 `VecField` ops unchanged; on `benchmarks/bench_trace.cpp`'s 19 800-term accumulation, 178 201 ops against 99 067, 48.2 ms before and 52.4 ms after (mean of three), and linear in the term count (2 000 terms 5.3 ms, 19 800 terms 52.4 ms); `acc = acc + term` for 1 000 terms 107 ms against 54 ms, copying a 99 067-op field 4.3 ms against 4.0 ms. Checked by a randomised oracle: 600 random trees with repeated subexpressions, opaque leaves of one closure type with different captures, signed zeros and NaN payloads evaluate bit for bit like the same expression built without merging, culled and through the POD interpreter. The oracle found a defect in the lowering it exercises: a gather table that differed from an earlier one only in a `-0.0` was joined with it, because tables were matched with `==`; they are matched by bit pattern now, in `field_pod.hpp` and in the GLSL data block.

  **Two cracks, known and not being solved now.**

  `typeid` gives one type per lambda *expression*, not per lambda body. Three `[noise]` lambdas written inside one loop are one type; three written in three places with identical bodies are three. That is correct for deduplication — different sites are different closures — but a report saying `types=3` cannot say *why* there are three, and the person asking "why did deduplication not fire on my scene" needs the report to answer rather than a person. Same seam as `unknown`.

  **Instancing: what is actually on the table, measured 2026-09-16 rather than estimated.** `cook()` on the donut's dust shape finds **2 shapes for 19 801 objects** — 159 552 vertices' worth of geometry if each object carried its own, against 1 160 actually stored, **138× fewer**. And it reports that **19 800 of those 19 801 are opaque-refused**: the sharing opportunity is total, and every instance of it is blocked by one thing, the motion being written as a deformation when it means a placement.

  The upper bound on what that is worth was then measured directly, by rendering the same frame with the dust removed: **2.02 s with, 0.895 s without — the dust is 56% of the frame.** That bound is what says the work is worth doing, and it is also what says *which* work: the ~45 ms of per-frame vertex building is the small part, and the large part is that 396 000 triangle leaves would become 19 800 instance leaves over one 12-triangle shape, a 20× reduction in what the BVH holds and traverses.

  Which means the three layers are not equal. Reading placements out of structural motions needs a scale op and a `t`-dependent leaf (`VecOp::Const` is a literal, so the cheap route is not available); extracting them in `cook()` is straightforward after that; but only the third — a per-leaf transform in the BVH — touches traversal, and `ray_hit(ray, shape)` takes no transform today, so that one is a change to load-bearing code rather than an addition. It is the one carrying most of the 56%.

  Layout. Once the threshold is payload against L3, the next question is arrangement: `FieldOp` carries a `std::function`, a `type_index` and a payload size inline, so walking a pool to evaluate it drags the opaque-leaf metadata through cache along with the arithmetic. Splitting the pool into arithmetic and leaf-metadata arrays (SoA) is the natural continuation of the same threshold, and it is not worth doing before there is a lowering pass that walks the arithmetic alone.

  **Scrubbing time backwards already works, and is a property rather than a feature.** `materialize(trace, idx, t)` is a pure function of `t`: a motion field takes `(p, t)` and the trace never stores a previous frame. So rewinding is passing a smaller `t`, with no history buffer and nothing to record — the thing games with time control build machinery for falls out of the trace being a description rather than a state. Worth writing down because it is easy to lose: it holds *exactly as long as motion stays a pure function of time*. The moment a solver writes a pose that depends on the previous one, rewind needs history again. That is the strongest practical argument for keeping fields pure and putting any statefulness somewhere it can be named, and it is an argument that costs nothing today.

  **BVH under mutation: built by SAH, then followed by insert, remove and a cost guard.** `spatial/bvh.hpp` builds by SAH and, until 2026-10-22, had no refit or incremental update (it now has both; the settled paragraph below gives the figures), which was exactly the weak spot a destructible scene would hit — a voxel grid updates one cell, a BVH has to rebuild the tree. This is not a defect to fix speculatively; it is the shape of the thing, and the interesting angle is ours rather than generic: an exact analytic leaf *replaces* leaves, and rebuild cost scales with leaf count, so the same measurement that showed a 7 000× build-time difference in the DSL row is also the argument that exact leaves make a rebuilding scene cheaper, not just a static one. If a mutating scene is ever a real consumer, refit comes first and exact leaves make refit smaller.

  **Raised again 2026-09-17 as "optimise the BVH across updates", and the answer is still "not yet", now with the reason being a number rather than a preference.** The donut rebuilds both trees every frame, so an animation pays the build once per frame and `--video` pays it 750 times. Refit — walk the tree bottom-up, widen each box, leave the topology alone — is the standard answer and is `O(n)` with no allocation. What makes it the wrong thing to reach for *here* is that **instancing already took the build down to 19.8 ms against a ~2 s frame**, so refit is optimising about 1%. The 16× that instancing won on build was the large factor, and it has been collected.

  **The part that does not change is reused rather than refit.** Regenerating a scene with the same torus and a different scatter seed or count paid the triangle tree of the non-instanceable geometry (a 240x120 torus, 57 600 triangles) every time, though it depends on nothing the scatter touches. `render::lay_out` takes an optional `StaticCache`, keyed by a hash of the triangle bytes and confirmed by comparing them, so any input that changes a triangle misses and a collision cannot return a wrong tree. Measured with `lay_out` alone, medians of 5 on E-cores, `-O3`: 55.7 -> 22.2 ms at 1e3 instances, 53.1 -> 19.5 ms at 1e4, 153.9 -> 115.1 ms at 1e5, where the instance tree is what is left. The remainder at small N is the triangles themselves, their normals and the copy of the cached tree, which are not cached.

  Two things would change that, and they are worth naming so the item is not re-raised from scratch a third time. **Traversal getting cheap** — shadows, textures and refraction all add rays, and every one of them amortises the build further, so build share *falls* rather than rises as the renderer grows. **Mutation within a frame** — a physics step that moves bodies many times per frame pays the build per step, and `ball_pit_demo`'s `step_physics` is exactly that consumer, already written and already not using a BVH.

  There is also a correctness-shaped caveat specific to this scene, and it is the interesting half. Refit never makes a tree *wrong* — the boxes still contain their geometry — it makes it *bad*, because the topology was chosen for where things were. The donut's dust is close to the worst case imaginable for that: particles burst outward from a single point, so a tree built at `t=0` has every leaf in one cluster and by `t=3` describes a cloud it no longer separates. So refit here cannot be unconditional; it needs a rebuild trigger, and the honest one is measured rather than a constant — track the SAH cost of the refitted tree against the cost it had when built, and rebuild when the ratio crosses a threshold. That is a number `BVH` can report about itself, which is a better thing to add than a refit that silently degrades.

  **Following a changing instance set, settled 2026-10-22.** `BVH` gained `insert(shape)`, `remove(index)`, `rebuild()` and a cost guard, all written against the `Bound` operations alone (`merge`, `sah_measure`, `bound_of`), so the same code serves boxes, balls and the geodesic-ball bound; the tests run it over all three and over a bound that counts its merges (O(depth) per operation). An insert descends by the smaller growth of the child's surface measure and splits a leaf past four shapes at the median of the widest centroid axis; a remove drops the shape from its leaf and, when the leaf empties, replaces the parent with the sibling; both then refit the bounds of the path bottom-up from the children (growing bounds on the way down is not enough for a ball, whose merge is not monotone). Indices: `insert` appends, `remove(i)` moves the last shape into slot i and returns its old index (the caller does the same swap on parallel data, as `CookedScene::looks` does), `rebuild()` keeps every index. Cost: the tree tracks its SAH cost incrementally (`cost()`, checked against `measured_cost()`); it is built again when the cost passes 1.25 times what the last build had, scaled by log2 of the size. That number is from measurement, not preference: on 20 000 scattered instances, 4N random replacements by instances from the same distribution raised the cost 3% over a fresh build and the ray time not at all beyond noise; replacements from a shifted distribution raised it 10% (rays +12-25%) and from a distribution shifted off the torus 40% (rays +27-34%), and a 1.25 limit lets the first two alone and rebuilds the third. The node grew from 56 to 64 bytes (explicit `left` and `parent` links, since a split cannot keep the left child at `index + 1`), with no change in ray time on the 1e5 scene (4.4 vs 4.4 us per instance ray, 1.5 vs 1.4 per triangle ray, E-cores, noise about 5%). An updated tree is no longer in depth-first order; code that copies `nodes()` wholesale, such as the device copy of the triangle tree, takes a built one.

  The consumer is `render::lay_out(previous, cooked, ...)`: it keeps the previous scene's quadric table if unchanged, matches the new instances against the old by shape and the bits of placement and look, and removes and inserts the difference when under a quarter of the instances differ, else builds. The result gives the ray hits of a fresh `lay_out` (tested on add one, remove one, new seed, remove half and a run of ten small changes); only the instance order differs, with `looks` following. `examples/scatter_regen.cpp` moves a scatter's count step by step and checks the two scenes against each other on random rays; `donut_live` cannot host it, because its scatters are moved on the device and it has no static instances, so no host instance tree is in its path. Measured with `benchmarks/bench_instance_update.cpp`, medians of 5, E-cores (a second session was running at load 5-8, so absolute values carry that), -O3, instances from a torus scatter; `k` is the number replaced (k removed plus k inserted):

  | N | k | rebuild | edit only | diff + edit |
  |---|---|---|---|---|
  | 1e3 | 1 / 10 / 100 | 519 us | 0.6 / 6.2 / 62 us | 110 / 107 / 148 us |
  | 1e4 | 1 / 10 / 100 | 7.2 ms | 0.8 / 8.3 / 80 us | 1.26 / 1.27 / 1.39 ms |
  | 1e5 | 1 / 10 / 100 | 104 ms | 1.7 / 17 / 190 us | 19.3 / 19.8 / 19.2 ms |

  The tree operations cost 1.7-1.9 us per replaced instance at 1e5 and are not the limit: finding which instances changed (one hash of each, one probe of each) is, at about 19 ms for 1e5, a fifth of the build it avoids. For the 475 ms the regenerating scatter costs at 1e5 (materialize 136, cook 100, instance tree 140 in the earlier single run), the instance tree goes from about a quarter to a twentieth of the total when one instance is added; materialize and cook are untouched and are now the larger share. Updating wins on time until about 45% of the instances are replaced (104 ms against 19 ms plus 1.9 us a pair); the 25% limit on removed plus inserted is set lower for the quality reason above.

  **Many sources pulling each other on a torus, with a Barnes-Hut tree on the BVH. Settled 2026-10-22.** `physics/mechanics/gravity_field.hpp` is gravity of many softened sources in N dimensions: the acceleration of a source of mass m is G m r / (r^2 + eps^2)^((p+1)/2) for an integer exponent p (N-1 is a point source's own field, 1/r in 2D and 1/r^2 in 3D; p = 1 in 3D is the softer 1/r law), with an optional cap on the length of the total acceleration (which breaks conservation, and says so) and a per-body ratio of gravitational charge to inertia, so two kinds of body are pulled differently without touching the integrator (a source's mass is its charge; entering a body as a source of charge x inertia keeps the third law exact in inertia). On a closed space Poisson's equation has a solution only for zero net source, so the periodic field is that of the sources plus a uniform background of the opposite total mass. It is summed by Ewald's method, not by the minimum image with a correcting term: the kernel is split by a Gaussian at 1/alpha into a short part over the nearest images (it is zero past 0.73 of the box side at the default alpha L = 6) and a long part over reciprocal vectors with k = 0 (the background) left out, both run to eta = 4.4 widths (tail about exp(-19)). The split is exact for every integer exponent from 1 to N-1, softened or not, because the Plummer factor exp(-t eps^2) goes into the Gaussian representation of the kernel; the reciprocal coefficient is then a one-dimensional integral per vector, computed once (Gauss-Legendre, 6 panels of 12 points), and the real-space factor is a regularized incomplete gamma function of half-integer order. **Checked against a lattice sum that shares no code with it:** the images of every source inside a ball around the point, weighted by a Gaussian-smoothed step 0.5 erfc((r - 6 delta) / delta), the background cancelling exactly because the ball is centred on the point; by Poisson summation its truncation error is the transform of the smoothed tail at the reciprocal lattice, about exp(-(pi delta / L)^2), and the test shows it converged by comparing delta = 1.5 L with 1.95 L (agreement 1e-9 relative). The Ewald field matches it to 1e-8 relative for 2D 1/r, 3D 1/r^2 and 3D 1/r, each with and without softening (`tests/test_gravity.cpp`), and two splits (alpha L 4.6 and 8) agree to 1e-7. The field is periodic and continuous across the box face, odd under reflection for a symmetric configuration, the total force of an isolated periodic set is zero to 1e-9 of the sum of |m a|, the force is the gradient of the periodic potential (centred differences, 1e-6), and 24 bodies advanced 400 leapfrog steps keep the energy within 2e-3 and the momentum within 1e-10.

  `physics/mechanics/gravity_tree.hpp` builds `spatial::BVH<Box>` over the sources (a source is a degenerate box) and keeps the monopole of every node (mass, centre of mass, longest side) in a side array indexed by node, computed bottom-up by a post-order over the links, so it holds for a tree that was edited and is not in depth-first order; the BVH knows nothing of it. A node is replaced by its monopole when the point is outside its box and size / distance < theta. Only the short-range part of the Ewald split needs the tree: the reciprocal sum goes through the structure factor of all the sources, exact, O(vectors) per point, so the tree error is the opening error of the short-range part alone. The short-range kernel is zero past the cutoff, so the walk is made once for each of the 3^N periodic images of the tree whose box comes within the cutoff of the point (the point and the sources are wrapped into the box first, so the nearest images are enough). With theta = 0 nothing is accepted and the tree equals the direct sum to 1e-10. Error of the force on every body against the direct sum, 4 000 bodies uniform in the box, softening 0.3 of the mean spacing, per-body relative error (mean, and worst body) and the same against the rms force (a uniform box has bodies whose forces nearly cancel, which inflates the per-body worst case):

  | theta | 3D torus mean / max | 3D torus mean / max, vs rms | 3D free mean / max | 2D torus mean / max |
  |---|---|---|---|---|
  | 0.2 | 7.0e-4 / 2.3e-2 | 4.8e-4 / 4.6e-3 | 2.8e-4 / 7.1e-3 | 5.5e-3 / 2.8e-1 |
  | 0.3 | 5.2e-3 / 1.6e-1 | 3.5e-3 / 1.6e-2 | 1.1e-3 / 1.3e-2 | 1.5e-2 / 2.3 |
  | 0.4 | 1.6e-2 / 4.9e-1 | 1.1e-2 / 4.8e-2 | 2.5e-3 / 5.4e-2 | 2.9e-2 / 2.9 |
  | 0.5 | 3.3e-2 / 6.0e-1 | 2.3e-2 / 9.7e-2 | 4.7e-3 / 1.7e-1 | 5.1e-2 / 5.1 |
  | 0.7 | 9.0e-2 / 1.6 | 6.2e-2 / 2.3e-1 | 1.3e-2 / 1.9e-1 | 1.2e-1 / 12.5 |
  | 1.0 | 3.0e-1 / 9.1 | 2.0e-1 / 8.3e-1 | 3.6e-2 / 7.3e-1 | 2.9e-1 / 26.5 |

  **The default is theta = 0.4**: about 1% of the force on a body on a torus (0.25% in free space), worst body 5% of the rms force. On a torus the error is about ten times that of free space at the same theta, and steeper than theta^2, because the short-range kernel falls off like a Gaussian on the scale 1/alpha, so a node that is small against its distance can still be large against the kernel's own variation. Holding the node size to theta d / (1 + c (alpha d)^2) and accepting any node whose nearest point is past alpha d = sqrt(7) was tried (c = 0.03, 0.1, 0.3, 1): at equal cost it was less accurate than the plain rule, so the plain rule stays. Time for the force on every body, one evaluation including the structure factor and the tree build, 3D, theta 0.4, medians of 5, E-cores (`taskset -c 4-11`; two long jobs of another session were running at load 6-7, so the absolute times are inflated and the ratios are what to read); the direct sum is timed on 256 targets and scaled where N is large:

  | N | torus: tree | torus: direct | speedup | free space: tree | free space: direct | speedup |
  |---|---|---|---|---|---|---|
  | 1e3 | 105 ms | 317 ms | 3.0x | 6.7 ms | 7.7 ms | 1.2x |
  | 3e3 | 445 ms | 2.85 s | 6.4x | 34 ms | 67 ms | 2.0x |
  | 1e4 | 1.81 s | 30.4 s | 16.8x | 167 ms | 762 ms | 4.6x |
  | 3e4 | 7.8 s | 277 s | 35x | 975 ms | 6.9 s | 7.1x |
  | 1e5 | 34.2 s | 3 026 s | 88x | 6.8 s | 87 s | 12.9x |

  The crossover is N of about 100 on a torus (each direct pair visits up to 2^N images and an incomplete gamma; at N = 100 and 200 the ratio is 1.1x and 1.2x) and about 700 in free space (0.5x at 100, 0.6x at 200, 0.7x at 400, 1.0x at 700). The alpha that balances the two sums is flat: alpha L = 4.6, 5.5, 7 and 9 took 1.38, 1.27, 1.21 and 1.45 s per step (two evaluations) at 1e4; the default is 6. The consumer is `examples/nbody_demo.cpp`: 2 000 bodies (2D or 3D, two kinds with charge-to-inertia 1 and 0.5) on a torus, kick-drift-kick leapfrog with the tree and with the direct sum, printing the drift of energy and momentum and ms per step of each, and an SVG of the final positions.

  **Gravity on the closed 3-sphere, flight under it, and the two stereographic charts. Settled 2026-10-23.** Poisson's equation on S^3 needs a source of zero total mass, so, as on the torus, the field is that of the sources plus a uniform background: the potential of a unit mass is the Green's function of Delta G = delta - 1/V, V = 2 pi^2. The closed form `G(psi) = -(pi - psi) cot(psi) / (4 pi^2)` at geodesic angle psi (unit radius) was taken as a claim and checked two ways that share no code with it, and it holds as written; only the constant needed stating: the mean of G over the sphere is -1/(8 pi^2), which is the constant between the closed form and the expansion. (a) The zonal eigenfunction expansion, eigenvalues n(n+2), degeneracy (n+1)^2, zonal functions sin((n+1) psi) / sin(psi), gives G = -(1/V) sum (n+1) / (n (n+2)) sin((n+1) psi) / sin(psi), conditionally convergent, summed with Fejer weights over 4e5 terms: agreement to 4e-6 at psi = 0.3, 1, 2, 3, 3.1 once the constant is added (the Fejer tail is the limit, not the formula). (b) The equation away from the source, G'' + 2 cot(psi) G' = -1/V, by differences: residual below 1e-6 on (0.5, 2.9), and a candidate without the (pi - psi) factor fails it at 1e-2; the flux of grad G through the geodesic sphere psi0 is 1 minus the share of the volume inside, 1 at a small sphere, to 1e-6 for psi0 = 0.01, 0.5, 1 and 2.5. Near the source G is -1/(4 pi psi), the flat law; at the antipode it is smooth: G is finite (1/(4 pi^2)), its slope vanishes, and the force of a source at the antipode is exactly zero and linear beside it (delta / (6 pi^2) at offset delta, 1e-3 relative at delta = 1e-2, 1e-4, 1e-6). On a sphere of radius R the potential is G / R and the acceleration scales as 1 / R^2. `physics/mechanics/sphere_gravity.hpp`: `SphereGravity` takes sources (unit vectors of R^4 times R, with masses), a Plummer softening in the angle written as tan^2(psi_s / 2) = tan^2(psi / 2) + eps^2 / 4 so that it stays smooth through the antipode (a softening of psi alone has a cone there), an optional clamp and a per-body coupling as `GravityField` does; the force is `G m H(psi) (s - cos(psi) x)`, H = G' / sin = [(pi - psi) + sin(psi) cos(psi)] / (4 pi^2 sin^3(psi)), tangent at the receiving point, the numerator by series where the two terms cancel. It is the gradient of the potential to 1e-6 (centred differences, softened or not), and at small angles it is Plummer's flat law to 2e-3. A source on the receiving point gives no force (no direction) and a finite softened potential; an empty set gives zero. **Flight.** `flight_step` in `manifold_body.hpp` is kick - exact great circle - kick: second order, symplectic for a gradient force and symmetric, and with no force it is `geodesic_step` itself. `geodesic_step` re-reads the speed from the velocity each step, so over a long free run the speed random-walks: 5e-14 after 2000 steps of 0.1 (200 rad), 1.2e-12 after 2000 of 0.01 and 4.6e-11 after 2e4 of 0.01 (a drift of the speed by 2e-14); with a force, the point is renormalised and the velocity projected each step. Six bodies of masses 0.2-0.7 with softening 0.15 over time 2: the largest energy deviation is 6.4e-8, 1.6e-8, 4.0e-9 at dt = 4e-3, 2e-3, 1e-3, a factor 4 each halving; 2000 steps forward, the velocities flipped, 2000 back returns to the start to 1.1e-14. The system is carried to itself by every rotation of R^4, so the six components L_ab = sum m (x_a v_b - x_b v_a) are conserved, not the sum of m v (a sphere has no translations): the torques of a pair cancel because the pair law is symmetric in x and s, and a drift kick adds none, so they hold to rounding (1.8e-15 over 20 000 steps, 8 bodies, 0.0105 ms a step), while a clamp, which scales each body separately, moves them by more than 1e-5 (the check must fail there, and does). **Antipodal focusing.** Every geodesic from P with the same speed reaches -P after an angle pi, whatever its direction (no masses, 32 shots, 2000 steps: 3.6e-14). With three fixed masses and softening 0.1 the arrivals spread in proportion to the mass, to first order: total mass 1e-3, 1e-2, 1e-1, 1 give a worst distance from -P of 1.96e-4, 1.96e-3, 1.96e-2, 0.19. **Charts.** `spaces/sphere_charts.hpp`: the north chart x = y' / (1 + y_N / R) is the hemisphere y_N > 0 around the north pole (projection from the south pole), the south chart x = y' / (1 - y_N / R) the hemisphere y_N < 0 around the south pole; both project along the same axes, so a point of the equator has the same coordinates p (|p| = R) in both charts, and the transition is the inversion x -> x / (kappa |x|^2) (an involution; the south chart is the mirror image of the north one seen from outside). Embedding to chart and back is identity to 1e-12 R in both hemispheres (kappa 1, 2.25, 0.3, in S^2 and S^3), the velocity pushforward is the exact derivative and is conformal (lambda |u| = |w| to 1e-10), the derivative of the inversion carries the north vector to the south one, a great circle of the embedding is the `ConstantCurvature::exp_map` geodesic of either chart (1e-10), and a ball crossing the equator at p is on the rim of the north disc at p and enters the south disc at p, where the south chart's own geodesic continues it to 1e-10. The consumer is `examples/sphere_nbody_demo.cpp` (the n-body run, the focusing table, an SVG of the two discs with a fan converging at the antipode). Not done: a Barnes-Hut tree on S^3 (the direct sum is O(sources) per point), and a periodic or Ewald-like sum is not needed there since the space is closed.

  **A scene version has an identity, so regeneration is told what changed. Settled 2026-10-22.** The update above found the changed instances by hashing every one of them (19 ms at 1e5, a fifth of the build it avoids), and `cook()` placed every copy again (about 100 ms) though a scatter's sampling is prefix-stable: a count of N+1 is the N sites already there and one more, and only a new seed moves them. Three pieces make that known instead of discovered. (1) `Field::content_id(i)` is a Merkle hash of op `i` (code, immediates, child ids; noise and gather tables whole), cached per field in a flat vector caught up like the lookup index and cleared with it where ops are edited in place; it does not depend on layout, so separately built fields, fields built with and without merging, and a field with dead ops all agree, and changing one leaf changes the ids on the path to the root and no others. A hash match only nominates: `same_value(a, b)` answers yes after comparing the canonical descriptions (`canonical_key`, equal ops numbered once), so a collision is a "no" (tested with every id forced equal), and anything without an exact description is equal to nothing, a false "changed" being a recomputation and a false "unchanged" a wrong picture. An `Opaque` leaf has a serial taken at construction and copied with the op, so a copy of a field is equal to it and a callable written again is not; a stateless closure is its type (interned, not hashed). `VecField` has the same ids over its ops and its factor fields. (2) `node_version(trace, node, t)` extends `content_key` to what the node produces at `t`: shape, material, motion, colour and emission fields, and for a Scatter the *generator* (item, target, seed, seat, all of those) is separate from the *count*; `Object::site` and `source_node` say which copy of which scatter an object is, and `Cooked::generators()` lists each scatter's run of objects with its version. A node with no exact description (a chart handed to `space()`, a closure the library cannot compare) has no version and keeps the old per-instance comparison. (3) `lay_out(previous, ...)` compares generators: an unchanged one keeps its first min(old, new) instances and adds or removes the tail, a changed one replaces all of its own (built again past a quarter, as before), and only instances without a version are hashed (`InstanceUpdate::hashed`). `cook(trace, root, t, &previous)` takes the leading objects of an unchanged generator from the previous cook and places only the rest; the sampler became the resumable `SurfaceSampler`, with its state kept every 1024 sites, so the draws for sites already known are skipped rather than repeated (on a 1e5 scatter the sampler is 75 of the 100 ms). Everything is checked against a fresh cook and lay-out: count +1, -1, +100, a new seed, a different torus radius, a different shape, a changed motion, the scatter removed and a second one added give the same objects and instances (as sets; removal reorders) and the same ray hits. Measured with `benchmarks/bench_instance_update.cpp`, medians of 5, E-cores, -O3, a 240x120 torus under the scatter, on a machine at load 8 from another session; k, the change in count, made no difference (1, 10 and 100 agree within noise), so one value per N:

  | N | cook, from scratch | cook, from the previous | `lay_out` update, before | after | fresh `lay_out` |
  |---|---|---|---|---|---|
  | 1e3 | 3.0 ms | 2.5 ms | 17.5-19.1 ms | 17.2-18.3 ms | 25-28 ms |
  | 1e4 | 12.2 ms | 3.0 ms | 19.1-20.4 ms | 17.8-18.6 ms | 34-38 ms |
  | 1e5 | 103 ms | 7.7-8.2 ms | 66.9-68.3 ms | 24.5-25.0 ms | 161-177 ms |

  A step of the slider at 1e5 is then 8 + 25 ms where it was 103 + 67 ms (a fresh cook and lay-out is 280 ms). The remaining 17 ms in the update was not the instances: it was the 57 600 triangles of the torus, emitted again from the cooked objects only to be looked up in the `StaticCache`, and it was the same at N = 1e3.

  **The static part is found before it is made, and a time that nothing reads replaces nothing. Settled 2026-10-22.** Two things closed the rest. (1) `lay_out` describes the static part before emitting a triangle: the words are, for each object that becomes triangles, its shape's `content_key` (kept on `Shape::key`), its placement and its material, as bits. `StaticCache::find(words)` answers from them, a hit being confirmed by comparing the words in full (a hash collision is a miss); on a hit the triangle emit, the corner normals and the BVH build are skipped and the tree and shaded triangles are copied out of the cache. A static object whose shape has no key (a placed mesh, a chart handed to `space()`) makes the whole part keyless, and then the triangles are made and compared as before, so a false miss costs time and a wrong hit cannot happen. (2) `Field::reads_time()` and `VecField::reads_time()` walk the ops the root reaches and answer whether the first input (the scene time) can matter: `U`, a metric's time coordinate and any callable leaf do (a callable cannot be asked), a constant, the point itself or an expression of those do not. `node_version` adds `t` to a node's version only if its motion, colour or emission reads it, so a scatter whose fields ignore time keeps its generator across times, and `cook(..., &previous)` reuses its placements at a new `t`; a scatter that reads time replaces its own instances and no one else's. For a move from `t` = 0 to `t` = 2.5 the cook and lay-out equal a fresh cook and lay-out on every instance and on random rays for a time-independent scatter (2000 of 2000 reused), a time-dependent one (0 reused) and a scene with one of each (only the time-independent 2000 reused). Update step of `examples/scatter_regen.cpp` (240x120 torus, E-cores, medians of 20 steps over 5 runs, -O3): 

  | N | before | after |
  |---|---|---|
  | 1e3 | 23.4 ms | 8.1 ms |
  | 1e4 | 24.6 ms | 9.2 ms |
  | 1e5 | 26.3 ms | 7.1 ms |

  What is left of the update is the copy of the cached tree and shaded triangles (a few megabytes) and, for growth steps, the inserts. Not done: the cooked sites of a scatter whose objects include a zero scale are compared by hash, not by generator; and `donut_live` still cannot use any of it, its scatters being moved on the device.

  **`cook()` and `Cooked<T>`, settled 2026-09-15.** The phase boundary the
  DSL already had in practice and never expressed: a `Trace` is mutated
  while it is built, then read-only while it is rendered, and nothing in
  the type system said so. `trace.cook()` names it.

  **A type, not a flag.** A `frozen_` bool would be a value that looks
  like a guarantee and does not hold one — nothing stops `trace.torus()`
  after it is set, and no compiler notices. That is the class named in
  `conventions.md`, applied to our own design. `Cooked<T>` holds the
  guarantee structurally: it simply has no `torus()`, `offset()`,
  `scatter()`. Mutation after cooking is not discouraged, it does not
  compile. Same shape as `Trace` describing and `Placed` evaluating,
  one level up.

  By value, and no `thaw()`. `Cooked` owns what it holds rather than
  pointing into a `Trace` that has to outlive it — the lifetime coupling
  is exactly what the owning-pool decision in `Field` rejected one level
  down. Going back to an editable trace means editing the source and
  cooking again, and if that is ever made a method it costs an O(n) copy
  and should say so, because the alternative is a flag under a new name.
  Note that "rebuild means a new trace" is *today* forced rather than
  chosen: `TraceNode` holds a `move_only_function`, so a `Trace` cannot
  be copied at all. Once motion is a `Field` that owns its pool, the
  constraint lifts and this becomes a decision on its merits.

  **What cook() does is one thing, and that is the argument for it.**
  Five loose ends that were being tracked separately turn out to be the
  same work, all of which can only happen once the whole scene is known:
  expanding operations into objects; deduplicating shapes (which *is*
  instancing — "one shape, N transforms" cannot be formed before you know
  there are N); compacting the pool; the field report, computed once
  instead of per frame; and AoS→SoA. One well-placed boundary closing
  several problems is the signal that they were one problem.

  Deduplication is part of the expansion, not an optimisation on top of
  it — without it there is no instancing at all. SoA and compaction are
  optimisations; those two are not the same kind of thing and the
  distinction is worth keeping.

  **Objects need their own index space, and `Scatter` proves it.** A
  trace node is an *operation*; a scene object is a *result*. One
  `Scatter` node materialises 600 sprinkles, so no node index names the
  37th. `Compose` is grouping and not an object at all — it has no
  `Placed` of its own — which is why `.moving()` on one is refused rather
  than implemented.

  **One layout, three sources of value.** Instancing and a runtime share
  a shape: a table of shapes plus N per-instance slots. They are not one
  mechanism, because the slots are filled three different ways, and
  collapsing that difference would hide the reason authored motion cannot
  serve physics:

  | source | when the value is produced | storage |
  |---|---|---|
  | instancing | once, in `cook()` | the transform, computed and kept |
  | authored motion (`.moving()`) | freshly, every frame | none — a pure function of `(p, t)` |
  | physics | accumulated across steps | the integrator's state |

  A pure function of time and an integrator are different mathematical
  objects and neither reduces to the other. `.moving()` cannot be
  extended into physics; physics is a second branch.

  **The consumer for all of this already exists and is not connected.**
  `examples/ball_pit_demo.cpp`, `tumbling_body_demo.cpp` and
  `native_collision_demo.cpp` include no `io/build.hpp` at all — the
  physics in this repository lives entirely outside the DSL. And
  `ball_pit_demo.cpp:349` is `step_physics(std::vector<SphereBody>&
  bodies, ...)`: a flat array of objects carrying mass, velocity and
  position, mutated every substep. That is a runtime, already written,
  already working, answering questions (1), (2) and (5) above before they
  were asked.

  So the next piece of work is not "design a runtime" but **express
  ball_pit's scene through the DSL and record where it breaks**. The
  breakages are the requirements, found rather than invented, and the
  demo is a deliverable either way. Questions (3) and (4) will not be
  answered by it — ball_pit uses a fixed `dt`, no event timers and no
  metric — and that is correct: a consumer answers the questions it
  actually has.

  **Five architecture questions, open and deliberately unanswered.** `PointField`'s eventual shape depends on them, and none has a consumer yet, so answering them now would be designing blind — the same call as `SurfaceWithBoundary`. Recorded so the questions are not rediscovered as surprises: (1) where mutable state lives, if `Trace` and `Cooked` are both immutable; (2) what the unit of mutation is — rebuilding a trace, or writing into a parallel structure; (3) that "time" today names two different things, an object's local clock and the metric's dilation, which are items 5 and 7 and not the same mechanism; (4) whether a timer is a function of `t` or state, which decides whether an object can wait, pause or start on an event; (5) whether physics bodies live in the trace or beside it.

  What follows from them *now* is only this, and it is cheap: a field's input should be a **named environment**, not a bare `t`. Today it holds one thing and behaves exactly as today. Tomorrow a local clock or a solver-written pose is a field added to a struct rather than a change to every signature and every call site. Item 3 taught the same lesson from the other side — the slot was right and what filled it changed.

  **The next crack, unsolved on purpose.** Two different `particle_motion`-shaped functions in one trace with different parameter counts — one takes three, the other five. Either the shorter pads with empties, or two genuinely different functions end up sharing a name. No answer yet, and it does not block the work, but it is the same seam and worth expecting.

  Justification is **not** call overhead, and not hoisting either: the indirect call measures ~3.3 ns/vertex (~4% of a realistic callee), so even eliminating it entirely is within noise. Hoisting `t`-only subexpressions is worth well under 1% on today's demo — `grow_scale`'s one `smoothstep` per vertex — and only becomes material when the hoisted subexpression is expensive, which is the local-time item below. What fields actually buy is the parameter table and deduplication above. Selling them on hoisting would be a fourth refuted hypothesis in this series. The one thing a closure genuinely cannot do at all, rather than merely more slowly, is lower to branchless GPU code.

  Composition becomes substitution of one tree into another's point slot — which is why `.moving()` was made to compose first, so the callable and structural forms agree. Degradation has to be **visible**, not felt: `is_structural()` plus a per-trace report ("8 opaque leaves, deduplication covers N of M fields"), for the same reason the polynomial solvers' lost root is pinned as a `[!shouldfail]` test rather than a comment. Silent degradation that the user discovers through a frame-rate drop is the NaN-through-the-filter failure one storey up. The risk to watch is that an expression IR quietly becomes a language users must learn; a plain lambda has to stay a normal form rather than a permitted exception, and the donut demo is the place that failure would show.
- **[want]** Time as a property of the space rather than a global scalar. A scene's `t` is currently one number shared by every node. Deriving each object's local time from where it sits degenerates to exactly today's behaviour in flat space at no cost, and in a curved one gravitational time dilation falls out of a metric the library already has — `physics/relativity/` ships Schwarzschild and Kerr as substitutable callables. What is missing is small and specific: no proper-time or local-clock-rate helper exists there yet. Deliberately preferred over a conventional keyframe timeline, which was considered and dropped.
- **[want]** Open the node set — `Kind` is a closed `enum class`, so a user cannot add an operation without editing `build.hpp`. `io/scene.hpp` already solved the same problem in the same namespace with an open registry (`register_shape_kind()`), so this is a consistency gap, not an unknown. Sequenced after the field work because whatever opens `Kind` has to survive lowering.
- **[course]** Vulkan live display — CPU raytraces (`donut_demo.cpp`'s `--photo` engine, extended), a small new Vulkan path just presents the frame (swapchain + one texture, uploaded and blitted every frame), not `viewer::App`'s mesh/point-cloud rasterizer. Real new plumbing (no sampled-texture descriptor/pipeline exists today), deliberately deferred rather than landed blind against a deadline.
- **[want]** Texture/UV mapping — **the sampling half landed 2026-09-25, the field half is open**: `render/texture.hpp` loads and samples an image, and `ball_pit_demo` lays one on the vat's UV net; `io::Material` still has no texture, so nothing in the DSL can reference one. `io::Material` currently has only `base_color`/`roughness`, no texture at all. Flagged back on 2026-09-07 alongside conform-to-surface (now shipped as `offset_surface`) as one of two primitives needed before the DSL; still open.

- **[course]** **A motion written as a plain lambda cannot be seen into, so almost every node in the donut demo is refused as a deformation when it is nothing of the kind.** Found 2026-09-17 by asking what the demo's 22 tessellated objects actually were.

  ```
  scene: 35222 objects -> 35200 instances + 429848 triangles;
         10995 refused (motion deforms), 35200 shared
  ```

  The 22 are eighteen sprinkle `Scatter` nodes (six colours times three bands, one merged object each), plus the cube, the dough and the icing, plus the table — which *is* a placement and goes to triangles only because it carries no closed form. The 10 995 counts refused *instances* rather than objects: 10 992 sprinkles and three whole nodes.

  The cause is not that those motions deform. It is that they are written as plain `(p, t)` lambdas:

  ```cpp
  auto grow_scale = [](const Vec<double, 3>& p, double time) {
      return Vec<double, 3>{p * e};        // a textbook placement
  };
  ```

  An opaque leaf that touches `p` sets `reads_point`, and `is_placement()` then answers "deformation" — correctly, since it cannot ask a closure what it does. Written structurally as `scaled(point(), e)` the very same motion is a placement. The donut's dust was rewritten that way when instancing landed; **nothing else was**, so everything but the dust is refused for a reason that is a matter of spelling.

  So sprinkle instancing has **two** preconditions, not one: the axis binding above, and structural motions here. Both are fixes rather than features, and they land together.

  **And rewriting the demo's motions does not close this**, which is the part worth separating. The demo is one user; the trap is in the API. Someone who writes `p * e(t)` as a lambda — the shape `.moving()` itself advertises, `(point, t) -> point` — gets a silent refusal, no instancing, and no way to find out why. `docs/conventions.md` records the general form under "the variant where the value is honest and the *name* is not": a predicate that can only see a syntactic property must not carry a name promising a semantic one. Three fixes, not exclusive:

  - **The refusal explains itself.** `VecField` can distinguish "refused because an opaque leaf reads the point" from "refused because the expression genuinely is not affine", and the first should say so with the remedy attached. `Cooked::refused_nodes()` already exists precisely so a report can point at nodes rather than state a number nobody can act on; it should also point at *why*.
  - **`is_placement()` is renamed to what it computes.** `is_recognisably_a_placement()` is ugly, and the ugliness is the point: it makes a caller ask "recognisably by whom?".
  - **The structural spelling becomes the obvious one.** `scaled(point(), e)` is already shorter than the lambda; it is simply undiscoverable, because the slot advertises a callable. That is an API problem and not a documentation one.

- ~~**[want]** Does collapsing the dust into one node remove `content_hash`'s O(vertices) cost?~~ — **measured 2026-09-17, and the answer is that the cost was never there.** The premise was wrong twice over. The dust is not a `Literal`: `flake()` builds a `Space` node, hashed by its exact form's type plus sixteen chart samples, so the O(vertices) path touches exactly two nodes in the whole scene (the cube and the table). The "158 000 vertices hashed across the dust" figure was inherited from when a speck was `literal(dust_speck(...))` and never re-derived.

  And the measurement settles it regardless of the premise:

  ```
  trace nodes   35 231 -> 34
  cook()         114.6 -> 105.6 ms
  ```

  A thousandfold fewer nodes and **8% less time**. So neither "proportional to the node count" nor "stronger than that" — hashing was not the cost at all. What dominates is per-*object* work, which the collapse does not reduce: 46 196 objects still each evaluate a placement and resolve a material through two opaque fields. The first line to look at, if cooking ever needs to be faster, is that — not the hash.

  The real win is the one the entry was not looking for: the trace went from 35 231 × 456 B ≈ **16 MB** to 34 × 456 B ≈ **15 KB**, which is what makes two million reachable at all.

- **[course]** **An exact form that does not agree with its own chart, created that way by the factory.** `docs/api-reference.md` already states the invariant — *"the exact form must agree with the map, so any operation that can move them apart clears it"* — and `.moving()` honours it. A factory can violate it at the moment of creation, and one does.

  | factory | exact form | chart it is recorded beside | agree |
  |---|---|---|---|
  | `torus(R, r)` | `Torus{R, r}` | `make_torus` | yes |
  | `cylinder(r, h)` | `BoundedQuadric::cylinder_z(r, 0, h)` | `make_cylinder` — lateral surface, no caps | yes |
  | `sphere(r)` | `BoundedQuadric::sphere(r)` | `chart_of(Sphere<2,T>)` | yes |
  | `flake(half, bulge)` | a sphere clipped to `Box{-half, half}` | a spherical cap out to `v_max`, derived from `rim = min(hx, hy)` | **no** |

  `flake`'s chart is the disc *inscribed* in the clip box; the exact form keeps the box's corners. For the donut's `flake({0.010, 0.010, 0.003})` the square's corner reaches 0.0141 against the disc's 0.010 — so the exact form is a rounded square and the tessellation is a circle. The exact form covers **more**, which is the opposite direction from the guess that prompted this entry.

  **Fixed 2026-09-17, and the measurement corrected the diagnosis twice on the way.** The disagreement is not in the corners — sampling both shows x and y agreeing to the bit — it is in **z**: the chart spans `[-0.001431, 0.003]` and the old clip claimed `[-0.003, 0.003]`. Below where the cap ends the sphere keeps widening past the slab's half-width and is cut by its sides, so the exact form carried a skirt the tessellation had never heard of.

  Making the chart match was not available: `ParametricSurface`'s domain is a rectangle in `(u, v)` and "sphere ∩ box" is not one. Making the clip match was, and is exactly right — clip to `rim` and to the height where the cap reaches it, and then at any `z` above that the sphere's radius is at most `rim`, the sides never cut, and the bottom cuts precisely where the cap ends. Sphere ∩ box is then *equal* to the cap. The clip stops describing the size the caller asked for and starts describing the region the surface occupies, which is what a clip is for.

  **The visual consequence, decided rather than absorbed.** Removing the skirt makes the dust finer: on a dust-heavy frame 10.5% of pixels change and 95% of those get brighter, which is less speck covering more background. Kept, for three reasons in order of weight. The skirt was surface the tessellated form never had, so keeping it would keep a disagreement between two renderings of one node — the thing this fix exists to remove. Finer specks read as dust rather than as gravel, which is the direction the scene has been moving anyway. And if the cloud wants more visual mass, the honest lever is the flake's `half`, not a loose clip. A frame hash cannot referee this: the representation changed honestly, so the pixels *should* differ, and only an eye can say whether the new shape is the wanted one.

  Guarded by a test over every factory that records both parts, comparing the sampled bounding box of the chart against the exact form's own. Bounding boxes rather than surfaces, deliberately: cheap, dependent on nothing that can itself be wrong, and it catches the failure that matters — one of the pair covering ground the other does not. A divergence is a bug until someone writes a documented exception, and there are none. Checked against the old clip before being trusted: four assertions fail there.

  The reason this is worth an entry rather than a fix-in-place: the invariant exists, is documented, and is enforced on the one operation that was thought to threaten it. Nothing checks it where the two are first written down together — which is also where it is cheapest to check, since both are right there.

  **Withdrawn along the way:** the claim that this class was demonstrated by `cylinder()` losing its end caps. It is not. `make_cylinder` is the lateral surface with no caps and `BoundedQuadric::cylinder_z` matches it exactly. The donut's sprinkles became open tubes because the demo swapped its own capped `solid_cylinder` helper for `scene.cylinder()`, which is a different *shape*, not a different *representation* of the same one.

- **[course]** **`scatter_frame` builds its frame from one normal, and orienting *along* a surface needs a full one.** Found 2026-09-17, by checking a claim rather than by hitting a bug: the plan said the donut's sprinkles would instance for free once the renderer read a cooked scene, since `cylinder()` records a `BoundedQuadric`. They did not, because the demo does not use `cylinder()` — the sprinkle is a `literal()` of a hand-built mesh whose axes have been permuted so its length points sideways.

  That permutation is the finding. A frame derived from a normal is enough for anything that points *along* the normal, or that has no orientation worth speaking of. A sprinkle lies *across* the surface, and there is no way to say so: `scatter()` binds the item's local z to the normal and offers no alternative, so the only remaining move is to rotate the mesh by hand — which a closed form cannot follow, since the exact `BoundedQuadric` is a cylinder about z and stays one. Hence `Literal`, hence no instancing.

  Two distinct gaps sit under that, worth separating because they block different things:

  - **The in-plane directions are arbitrary.** `basis_from_normal` (Duff et al. 2017) picks *some* pair perpendicular to the normal, and `scatter_frame` then turns it by a hash of the site index. Perfect for sprinkles, which really do lie every which way. Useless for anything that should follow the surface itself — hair along a curvature direction, scales along a flow, tiles along a parametrization.
  - **The axis binding is fixed.** Local z goes to the normal, always. This is the one the sprinkle is blocked by, and it is the cheaper of the two.

  Until a full frame exists — normal, tangent, bitangent, with the tangent meaning something — scattered items that must lie along a surface stay `Literal` meshes and cannot be instanced. The ceiling that puts on the donut demo is concrete and measured: 11 000 sprinkles are 429 848 triangles that a closed form would have made 11 000 instances of one shape.

- ~~**[want]** Scattered items lying fully flush to the target surface, not slightly proud of it~~ — **done 2026-09-17**, and the cause was not the one this entry assumed. They were not proud because of a tangent-plane approximation; they were proud because `offset()` silently rendered at 48×24 instead of the tessellation it was asked for, and a coarse mesh of a convex surface sits *below* the analytic surface that `scatter()` places against. The gap was lifting every item into view by accident. `scatter_lift()` now raises an item by the depth of its own geometry below its local origin, and `scatter()`'s new `seat` parameter says how deep it sits — 1 rests it on the surface, 0 puts its origin there, between is a press into something soft. Note the trap that made a named parameter necessary: since the lift is derived from the item's own lowest point, shifting the item's mesh down by *d* lowers `min_z` by *d* and raises the lift by *d*, so the two cancel exactly and the nudge does nothing. The donut demo had exactly such a nudge, with a comment explaining what it was for, and it had silently stopped working.

- **[course]** **Per-instance parameters: `MotionEnv` gains the instance's own origin.** The blocking step for two million particles, analysed 2026-09-17, and the analysis is the valuable part because the obvious answer is the wrong one.

  The obvious answer is "carry the exact form through `cook()` so a `Scatter` can be instanced". That is true and necessary and it is **not** what stops the dust, because the dust is not a `Scatter` at all: `Scatter` places an item on a *surface*, while dust flies through space with a per-particle burst direction, target letterform and swirl seed. The recorded measurement "one node, N instances → trace = 3 nodes" was about sprinkles.

  Restructure the dust as a scatter over an invisible sphere — so that `site.position` *is* the burst direction — and the real wall appears one line further in:

  ```cpp
  pl = n.transform.placement_at(MotionEnv<T>{Vec<T, 3>{}, t});   // cook(), once per node
  ```

  **A `Scatter` node has one motion field for all its instances.** For particles to fly differently the field must tell them apart, and the only thing it can tell them apart by today is `p`, the vertex position — which sets `reads_point`, makes the motion a deformation, and makes instancing refuse. Correctly: that test is doing its job. Per-instance variation smuggled through the vertex position is exactly what it exists to catch.

  The fix is the one `io/field.hpp` already wrote down the reason for — *"tomorrow a field is added to this struct and nothing else moves"*. `MotionEnv` gains the instance's **origin**, its rest position. A field reading `env.origin` stays a placement, because an origin is one value per object rather than one per vertex, so the affine test still holds and the instance path still applies. The dust then becomes a single node with N instances, each on its own trajectory.

  **A second reason, found 2026-09-17 while trying to avoid needing this at all.** The plan had a cheaper intermediate step: a purely *radial* burst needs no per-instance parameters, because `A(t) + p·s(t)` is affine in the point and each particle's trajectory falls out of its own site position. The instancing half of that is true — `is_placement()` returns true and the objects collapse onto one node. The animation half is not.

  A placement is `p ↦ b + R·s·p` applied to every vertex, and a scattered instance's vertex is `seat + F·local`. So:

  ```
  particle centre = b + R·s·seat     grows with s
  particle radius = s·|local|        grows with s, the same s
  ```

  One scalar drives both. Expanding the cloud tenfold inflates every fleck tenfold: not a burst, an inflating balloon. They cannot be separated, because `seat` and `local` reach the field already summed into one point — and separating them is exactly what an instance origin is.

  So `MotionEnv::origin` is not only what swirl and letterform convergence need. It is what a *plain radial burst at constant particle size* needs, which removes the cheaper intermediate step from the table: building it would mean building something thrown away one step later.

  **A prediction, registered before the measurement rather than after.** Stage 1 gave 36.4× on stored vertices and 11% off the frame, and the gap between those two is itself explained: traversal is logarithmic, and `log(46 000)` against `log(429 000)` is almost nothing. So at two million the expectation is **7–10× on storage and 10–15% on the frame** — storage *falls* because one node's worth of geometry is amortised over vastly more instances while the instance array itself grows, and the frame moves little for the same logarithmic reason. Written down here so that a frame growing by more than that is a signal to look for a cause rather than a number to accept.

  **Scored 2026-09-18, and the prediction was wrong in the direction nobody checks.** Measured at 2 021 984 instances: **1646×** on stored vertices, not the 7–10× predicted — and the prediction was that the ratio would *fall* from stage 1's 36.4×, where it rose by forty-five times instead.

  The reasoning recorded alongside it is where the error lives, and it is worth keeping visible rather than deleting: *"storage falls because one node's worth of geometry is amortised over vastly more instances while the instance array itself grows"*. The instance array does not appear in that ratio at all. Its numerator is vertices-without-instancing, which grows with N, and its denominator is the shape table, which does not — so the ratio must rise. Two different quantities were being conflated, the vertex ratio and total memory, and only the second behaves the way the sentence describes.

  Being wrong by forty-five times in the flattering direction went unnoticed for a day because 1646× is a pleasant number and nobody re-derives a pleasant number. The frame half of the prediction, 10–15%, has still not been scored at all; peak RSS at two million is 0.93 GB after the diagnostic loop was removed, and `cook()` is 2642.7 ms over 34 trace nodes.

  Scope, so this is not mistaken for a parameter change: `MotionEnv`, a new leaf factory in `VecField` alongside `opaque_of_time`, `is_placement`/`affine_in_point`, `cook()`, `materialize_mesh()`, `Shape::exact`, and the renderer moving from `materialize()` to `cook()`. That last one uncovers a further gap — `cook()` hands back *rest* geometry for a refused (deforming) object with no way to recover its motion, so a renderer driven by `Cooked` alone would draw the exploding cube unexploded. ~~Cost at two million, measured by extrapolation rather than guessed: roughly 350 MB for instances plus tree~~ — **wrong by a factor of three, corrected 2026-09-17 by measuring `sizeof` instead of recalling it**, and the correction is the point of checking before rather than after:

  ```
  sizeof(Object<double>)     192 B  ->  366 MB at 2M
  sizeof(Instanced<BQ>)      112 B  ->  214 MB
  the BVH's own copy of them         ->  214 MB
  nodes_, reserving 2N               ->  244 MB
                                         ~1.0 GB
  ```

  plus about 160 MB of build temporaries (`boxes`, `centroids`, `indices`). The 350 MB estimate forgot `Object<T>` entirely, and forgot that `BVH::build` takes its shapes **by value** — so the instance array exists twice.

  **The full breakdown, measured 2026-09-17, since a figure without its layout is not a figure.** `Object<double>` is 192 B: `shape` 8 at offset 0, `source_node` 8, `translation` 24, `scale` 8, `rotation` **72**, `instanceable` 1 at offset 120, `material` **64** at 128 — members summing to 185 with 7 bytes of padding. That padding is *not* reclaimable by reordering: 185 rounds to 192 under 8-byte alignment in any permutation. A `BVH` node is `{Box<3,double> bounds, uint32 first, uint32 count}` = 56 B, of which the bounds are 48.

  What that makes available, in order of return, with the projected `Object` alongside:

  | change | saves at 2M | `Object` after |
  |---|---|---|
  | ~~`rotation` 72 → 32 (quaternion)~~ — **done**, `Object` only; measured 192 → 152 B | 80 MB | 152 B ✓ |
  | `material` 64 → a `uint16` index into a palette — **measured first, and it is not what the estimate assumed**, see below | ~97 MB | 88 B |
  | `shape`, `source_node` → `uint32` | 16 MB | 80 B |
  | BVH bounds in `float`, rounded outward — a bound only has to be conservative | 96 MB | — |
  | `prim_indices_` `size_t` → `uint32` | 8 MB | — |
  | the BVH owning its shapes: hold indices into `Cooked` instead | 214 MB | — |

  **How the rotation change is staged, decided before it starts.** The rotation is a matrix the whole way down — `scatter_frame` assembles columns from a basis, `rotated()` produces one through `SO3::exp`, their product is one — so a quaternion can only come from `from_matrix`, and the round trip is unavoidable. Measured over four thousand real objects: **zero bit-exact**, worst element 1.72e-15, about 8 ulp. A byte-identical frame after the matrix goes is therefore impossible rather than unlikely.

  What that does *not* explain is any visible change, and the first two attempts to explain one both failed the same way — by sounding like a mechanism without ever being given a number.

  The donut's worst object radius is 10.3 world units, so the worst vertex displacement is 1.8e-14 against a pixel of 5.4e-3: eleven orders below one. *The picture cannot shift.* The second story — that exact ties at a silhouette flip instead — is arithmetic too, and it does not survive it either. The band in which a tie could fall the other way is 3.3e-12 of a pixel wide; a 960×720 frame at four samples a pixel is 2.8 million rays, of which perhaps 5% land on a silhouette, so the expected number of flipped ties is about **5e-7**. Not "occasionally" — not once.

  So the honest prediction is the stronger one: **step two should be byte-identical as well.** If it is not, neither of those stories covers it and the cause is something else — so the response is not a third explanation but a count: how many pixels differ, and where. One on a silhouette is a coincidence to accept and close. A dozen scattered across the frame is systematic, and no displacement of 1e-14 accounts for it under any reading.

  **Both steps came out byte-identical, including the second one** — the strengthened prediction held, and the picture did not move by a pixel even though every rotation now arrives through a quaternion.

  Four preconditions were checked rather than assumed, and each could have made the result meaningless:

  - **The renderer is deterministic.** Three runs of one binary, one hash. Without that, "one pixel differs" could not have been attributed to anything.
  - **The comparison can fail.** Substituting an identity quaternion for one object in forty changes the hash. A byte-identical result from a comparison that cannot detect a difference is not a result.
  - **The expansion stayed cold.** `to_world` runs per vertex and now takes the matrix as an argument; the quaternion is unpacked once per object, beside the instance that keeps it as a matrix. A caller writing `o.rotation_q.to_matrix()` inside the vertex loop would have spent the memory and bought nothing.
  - **The object count is asserted separately from any picture.** A lost rotation makes instances coalesce, and a slightly denser cluster is not a visible defect — neither a hash nor a pixel diff distinguishes forty thousand instances from forty thousand with two on top of each other. The test checks counts, pairwise distinct positions, pairwise distinct orientations, and that none of them is the identity.

  So, two steps and two rules:

  - **Step one: store both, and the matrix stays authoritative.** The renderer reads the matrix; the quaternion is carried and compared. A byte-identical hash then means something, because the thing being drawn has not changed — it tests the plumbing and nothing else. Two sources of truth for one rotation is otherwise exactly how they drift apart in silence.
  - **Step two: remove the matrix and introduce the pixel comparison in the same change.** Splitting them leaves a window where the matrix still exists but the renderer already reads the quaternion and nothing checks the difference. The comparison is what makes the removal safe; without it the removal is blind. It reports a count and locations rather than a pass/fail tolerance, because the expected count is zero and anything else is a lead rather than a threshold to clear.

  **The palette was costed before being built, and the costing changed it.** `Material` carries values rather than functions — `base_color`, `roughness`, `emissive`, `opacity` — so a palette is possible in principle; the question is whether per-object materials *collapse*, and that is a property of the scene rather than of the idea. The demo now reports it:

  ```
  t = 1.5    9 762 distinct materials over 46 196 objects   (21%)
  t = 3.9   13 188                                          (29%)
  ```

  Neither the two-to-five the estimate imagined nor the N that would kill it. Three things follow, and the mechanism matters more than the ratio:

  - **The collapse is saturation, not sharing.** `dust_color` clamps at `d/1.4` and the emissive at `d/0.75`, so every particle past those distances lands on exactly the same grey and exactly zero emission. The particles are not sharing a material because they mean the same thing; they are sharing one because two `clamp` calls flattened them.
  - **The saving is smaller than the 124 MB claimed**, which was computed as though the palette were tiny. At 21% on two million: 128 MB of inline materials becomes 4 MB of indices plus 27 MB of palette, so about **97 MB**.
  - **It is fragile, and cheaply broken.** Replace either falloff with something that does not saturate — an exponential, say — and the collapse disappears entirely, the palette becomes size N, and the change *costs* memory rather than saving it. The 97 MB currently rests on two `clamp`s in a demo.

  The ratio also moves with time (9 762 against 13 188 on two frames of one scene), so a palette built in `cook()` is a different size every frame. Not a defect, but not a constant to plan against either.

  If it is built: `Cooked`'s materials are immutable like everything else `cook()` produces, and that has to be said out loud rather than discovered — sharing one entry across N objects means recolouring one object recolours all of them, and someone will eventually want to flash a single instance. And the test that matters is not "the palette works" but "the palette *collapses*": two objects with equal materials must land on one index, asserted directly.

  **`Instanced<S>` is excluded from that list on purpose**, and `docs/conventions.md` carries the principle: compaction belongs to cold storage, and the hot path keeps whatever computes fastest. `Object`'s rotation is read once a frame; `Instanced`'s is applied at every leaf test the traversal reaches, and a quaternion costs more arithmetic to apply than a matrix. Same saving, traversal untouched.

  The trigger for revisiting it, written down instead of guessed: two million `Instanced` at 112 bytes is 214 MB that traversal walks and no cache holds. *If* traversal at that scale measures as memory-bound, the options are splitting the array — positions apart from rotations, so a ray touches only what it needs — or shrinking the element after all. Against a measurement showing the stall, not before it.

  Landing at roughly `2M × 80 B` for objects plus a much thinner tree. The next boundary after that is SoA — positions in one array, rotations in another — which is a different change and should not be started until the cheap ones are measured.

  **One of them was free and is already taken:** `BVH::build` takes its shapes by value, and the demo handed it an lvalue, so the instance array existed twice for the whole frame. `std::move` at the call site removes one copy — 214 MB at two million — and changes nothing else. Worth recording that checking *what still reads the vector afterwards* caught a real break: the report line counted `insts.size()` after the build and would have silently printed zeros.

  **Two million is a memory problem before it is a speed problem**, and that reframes the work: the reductions available are layout ones the entry below already separates from deduplication. `rotation` as a full 3x3 is 72 of `Object`'s 192 bytes and a quaternion would be 32; `shape` and `source_node` are `size_t` where `uint32_t` would do; `Material` is 64 bytes repeated per instance where most instances of a node share everything but one field. And the BVH taking its shapes by value is a duplication nothing needs. None of that is hard, and none of it should be guessed at either — the number above is what makes it worth doing, and a measured number after each change is what says it worked.

- **[course]** **One object deforming another.** Raised 2026-09-17, from wanting sprinkles to press visibly into liquid icing rather than merely sitting in it. Today a field is a function of `(u, v)` or of a `MotionEnv`, and it cannot see anything else in the scene — so "the glaze closes around each sprinkle" is not expressible, and neither is a footprint, a dent, a contact weld, or a drip running off a placed object. The concrete first case is narrow enough to build: `scatter()` knows its sites and does not expose them, so the target's own thickness field could read them. The cost is what makes it a real item rather than a tuning change — the naive form is a sum over every site per surface sample, 600 × 10 240 on today's donut and around 50 million at the sprinkle counts this scene now uses, so it needs a grid or similar index over the `(u, v)` domain before it is usable at all. Worth naming as the general capability rather than as "dimples": the same mechanism is what a scene needs before objects can be said to interact at all, and it is the one the demo keeps reaching for. Deliberately deferred past the first release.

  **What it breaks, which is not what it looks like it breaks.** Raised 2026-09-17 as a worry that this item ends the "`materialize(trace, t)` is a pure function of `t`" property — and therefore time scrubbing, which falls out of that purity for free. Worth being exact, because the worry is right about there being a cost and wrong about where:

  *Purity survives the dimple.* A scatter's sites are a pure function of (target surface, count, seed) with no history in them, so a thickness field reading those sites is still a pure function of `t`. Rewinding keeps working.

  *What the dimple actually breaks is acyclicity.* The icing's shape would depend on where the sprinkles are, and where the sprinkles are depends on the icing's shape — a cycle in a structure that is a DAG **by construction**, since a node can only name indices that already exist. The resolution is to scatter against the *undented* surface and apply the dent afterwards: one truncated iteration of a fixed point, which is what a solver would do anyway. That has to be stated here rather than discovered, because the alternative reading — make it self-consistent — has no answer, and someone will spend a day finding that out.

  *Purity dies one step later, at the thing anyone asks for next:* a mark that stays. A footprint, a dent that does not spring back, glaze that remembers. That is accumulation over history rather than a function of the current configuration, and it is where rewind genuinely stops — not here. Two different failure modes, one step apart, and only the second one is fatal to the property.

- **[want]** Removing the self-intersection an offset surface produces when its thickness exceeds the local radius of curvature. Seen in the donut's icing as pink fins around the hole, where the drip field pushed glaze down the inner wall of the torus and the offset overran its own centre of curvature. The demo's fix is to stop drips at the inner rim, which is correct for glaze — it runs off the outside, it does not climb into the middle — and is a scene decision, not a library one. The library answer is either a validity check (`offset_shell` refusing, or clamping, a thickness above the base's minimum radius of curvature, which `ParametricSurface` has the derivatives to compute) or self-intersection removal on the resulting mesh. The first is cheap and catches the case at the call that makes it; the second is a general mesh operation this tree does not have.
- ~~Per-instance orientation variance in `scatter()`~~ — **done**, and this entry was stale from the day it was written: `build.hpp`'s scatter branch already gives each instance its own spin about the normal, hashed from `(seed, site index)` so the same trace at the same `t` reproduces exactly. Kept here struck through rather than deleted, because it was listed as open for a week while the code was already in the tree, which is the failure worth remembering.
- **[want]** Generalize `Trace::space()`/`resolve_surface()` beyond `ParametricSurface<T>` to any `Surface` (`ImplicitSurface`, `Sphere`, `Hyperbolic`, ...) — today's `offset()`/`scatter()` chain is hand-specialized to `ParametricSurface`'s `evaluate`/`normal_at`/`area_element`, not the general `Surface` concept.
- ~~Exercise `mesh/conform.hpp` and `mesh/scatter.hpp` in a real demo~~ — **settled 2026-09-15, and the old entry was wrong in a way worth keeping visible.** It said they "compile and pass their own tests"; they had no tests at all, no caller anywhere, and are not reached by the umbrella header either — so their bodies are templates that had never been instantiated, only parsed. Both now have a test in `tests/test_mesh_ops.cpp` that instantiates them and checks the result, and both headers say in their own first lines that nothing calls them and why. They are kept rather than deleted because they serve the case the analytic path cannot — a target that is only ever a mesh — and they are explicitly **not** getting a demo. That is the whole of the maintenance they are owed.
- Explicitly not planned for the scene DSL: dependency-tracked reactive recompute (SolidJS/Vue-signals style). If that gets built at all, `gpu/derive_christoffel.py`'s declarative-derivation → auto-codegen pipeline is the intended home for it (GPU rendering section below), not this graph.
