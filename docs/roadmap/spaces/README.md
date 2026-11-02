# Spaces

## Where it stands

The concept hierarchy (`core/`) and the spaces on it (`spaces/`): Euclidean, sphere, hyperbolic as one constant-curvature family, Minkowski as a Lorentzian space with no distance, L^2 of an interval, Lie groups as spaces, products, parametric and implicit surfaces with exact partials, a space from its metric alone. `distance` is the geodesic distance or absent. The connectivity matrix grades every scalar against every space and algorithm.

## Now

- nothing is taken up here now

## Later

- nothing is planned here yet; what was left open is in `--area`

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area spaces`.

## Reasoning and earlier history

### Concept hierarchy, core geometry, mesh, viewer foundations

- Concept hierarchy: Set → TopologicalSpace → MetricSpace → ... → Surface (12 concepts)
- Concrete spaces: Euclidean<N>, Sphere<N>, Hyperbolic<N>, ProductSpace
- Geometry: 10 primitives, 11 intersect pairs, distance, contains, BoundedRegion + clip
- Boolean ops: intersection_region (Sutherland-Hodgman), difference/symdiff area
- Algebra: Vec, Matrix, AffineTransform, SO(3), SE(3), Group/Ring/Field concepts
- Mesh: Mesh<Surface>, midpoint subdivision + project, LodChain, icosahedron/tetrahedron
- Viewer: Vulkan 1.3, orbit camera, per-face coloring, edge overlay, dynamic resize, nix run
- Discrete: FiniteSet<T> (∪ ∩ ∖ △ ⊆, power_set, cartesian), GeometricSet, SVG Venn
- IO: Table, Svg (polygon fill), std::format for all types
- Extras: Point<Space>, Morphism pipes, Real50/Real100, axiom verification
- 295 tests, 2 audits (0 CRIT/HIGH), ~10K LOC

### Surfaces & usability

- ParametricSurface<T>: f(u,v)→R³ with auto metric. Torus, cylinder, cone, Mobius factories
- ImplicitSurface<T>: F(x,y,z)=0 with auto gradient. Sphere, torus, gyroid factories
- Marching cubes: uniform grid extraction → Mesh
- Mesh operations: merge, transform, flip_normals, compute_face/vertex_normals, centered
- Mesh primitives: grid_mesh, uv_sphere_mesh, box_mesh
- OBJ I/O: load_obj (polygon fan triangulation, negative indices), save_obj
- STL I/O: load_stl (auto-detect binary/ASCII), save_stl (binary)
- Quaternion<T>: from_axis_angle, from_matrix, to_matrix, rotate, slerp, inverse

### Hyperbolic ray marching (2026-08-29)

Third and smallest of the three manifold-backlog ideas raised 2026-08-28 (the other two, Riemannian optimization and geodesic procgen, are above). `examples/hyperbolic_tessellation_demo.cpp` — sphere-tracing through `Hyperbolic<3>`'s own hyperboloid-model metric, not the mesh+BVH stack (Euclidean ray-primitive intersection, no way to tessellate hyperbolic space into flat triangles) and not RK4 geodesic integration (unlike Schwarzschild/Kerr, `Hyperbolic<N>::exp_map()` is exact, no ODE to integrate). Camera fixed at `Hyperbolic<3>::origin()`; because its tangent space there is genuinely Euclidean (`metric_at` reduces to the ordinary dot product once the normal-orthogonality constraint zeroes the ambient time coordinate), the existing Euclidean `render::Camera`/`camera_ray_dir` pinhole formula gives the initial ray direction directly, embedded into the hyperboloid's `Vec<T,4>` with a leading zero — no new camera-basis math. Markers sit at three shells of hyperbolic distance along a once-subdivided icosahedron's 42 directions (`mesh::subdivide_once`, reused as unit directions, not as a mesh), sphere-traced via `space.distance()` to the nearest marker as the safe step size.

Two real, non-obvious things found while tuning it, not guessed: (1) a hyperbolic ball's apparent angular size follows `sin(alpha) = sinh(radius)/sinh(distance)`, not the Euclidean `radius/distance` — an initial `radius=0.25` at `distance=1.0` filled a third of the frame, corrected to `0.08`; (2) with only 12 base directions spread across the *entire* sphere of view (not clustered toward the camera's forward hemisphere the way objects in ordinary scenes are), a normal ~40-55° FOV only ever catches 1-3 of them by chance — needed a wide (130°) FOV plus the subdivision (12→42 directions) to show a genuinely scattered field in one still frame; the wide-FOV rectilinear projection itself stretches off-center markers into radiating streaks, an honest projection-formula side effect (not a hyperbolic-metric effect) left in because it reads as a striking, correctly-motivated part of the image rather than a defect. Palette also needed golden-ratio-conjugate hue decorrelation (`hue = fmod(i * 0.618034, 1)`) since consecutive icosahedron-vertex indices are often spatially adjacent, which otherwise clustered same-hued neighbors together on screen. Wired into `examples/CMakeLists.txt` and `flake.nix` (`nix run .#hyperbolic-tessellation`), `nix flake check --no-build` passes, 825/825 tests still green (no library code changed, only the new example + existing `subdivide_once`/`icosahedron`/render-engine reuse).

### SPD(n): log-Euclidean manifold on symmetric positive-definite matrices (2026-09-01)

`spaces/spd.hpp`'s `SPD<N, T>` (`N=2,3`) — the natural sibling of SO(3)/SE(3) rotation averaging: covariance-matrix descriptors in computer vision, diffusion tensors in DTI/medical imaging, spatial-covariance features in BCI/EEG. Chosen over the affine-invariant metric deliberately, as the cheap first step (a fuller/affine-invariant treatment is a real, explicitly deferred follow-up, not dropped).

The defining move of log-Euclidean (Arsigny, Fillard, Pennec, Ayache 2006): matrix log is a global diffeomorphism SPD(n) ↔ Sym(n) (ordinary symmetric matrices, a flat vector space), so a "point" is represented directly as `vech(log(S))` rather than `vech(S)` itself. Under that parametrization every `Manifold`/`RiemannianManifold` operation (`distance`/`exp_map`/`log_map`/`metric_at`) is the ordinary flat Euclidean formula — the curvature of SPD(n) is entirely absorbed into the log/exp change of coordinates, not left for `exp_map`/`log_map` to handle the way `Sphere`/`Hyperbolic` must. `vech()` weights off-diagonal entries by `√2` so the plain dot product on `PointType` equals the Frobenius inner product on the matrix itself.

Sharp, checkable consequence of that flatness: the Fréchet mean of SPD matrices under log-Euclidean has a **closed form** (average the logs, exponentiate back) — `frechet_mean()` uses it directly, no `riemannian_minimize()` call, unlike `Sphere`/`Hyperbolic` where retraction genuinely has curvature to contend with. This is the textbook reason log-Euclidean is called "cheap," not just a performance claim.

Deliberately **not** a `Surface` (no `project()`/`normal()`): under log-Euclidean, SPD(n) is flat and full-dimensional in its own log-space coordinates, with no meaningful embedding as a hypersurface of some larger ambient space the way `Sphere`/`Hyperbolic` have — so there's no normal direction to define, and forcing one on (the way `Euclidean<N>`'s own placeholder `normal()` does, purely to satisfy the concept for `Mesh<Euclidean<N>>` compatibility) would silently break `riemannian_minimize()`'s `project_tangent()` for any caller who did try to use it that way, by removing one real degree of freedom from every gradient step. `Manifold + RiemannianManifold` is the honest concept-hierarchy fit here, not `Surface`.

`N` restricted to 2 and 3: log-Euclidean's real numerical engine is the matrix logarithm/exponential of the (small) SPD matrix, computed via closed-form eigendecomposition — for symmetric matrices, eigenvalues are the real roots of the characteristic polynomial, so `solve_quadratic()`/`solve_cubic()` (already in `algebra/polynomial.hpp`) give them directly with no iterative eigensolver. This reuses existing infrastructure and deliberately stays out of the separately-tracked, much bigger general SVD/eigendecomposition backlog item (Backlog → Native math) — general `N` needs that to land first. 3×3 eigenvectors come from the null-space-via-cross-product of the largest-norm row pair of `S - λI`, with a Gram-Schmidt pass against already-placed columns to stay robust for near-repeated eigenvalues (a real, if secondary, numerical edge case, not one this scope tries to handle exactly).

11 new tests (`tests/test_spd.cpp`): eigendecomposition reconstruction + orthonormality for hand-verifiable 2×2/3×3 cases (a known 45°-rotated eigenbasis, not just the trivial diagonal case), `to_spd(from_spd(S))` roundtrips, and `frechet_mean()` checked against an independently-derivable ground truth (the per-eigenvalue geometric mean, both axis-aligned and under a shared rotated eigenbasis) rather than only self-consistency. 825→836 tests, same one pre-existing unrelated failure noted throughout this file.

### SPD(n): affine-invariant metric (2026-09-01)

`spaces/spd.hpp`'s `SPDAffineInvariant<N, T>` (`N=2,3`) — the fuller treatment explicitly deferred when log-Euclidean landed above ("do it seriously later, not dropped"). Same matrices, different metric: `<U,V>_S = trace(S^-1 U S^-1 V)`, invariant under every congruence `S -> A S A^T`. Log-Euclidean's `<U,V> = trace(UV)` via `vech(log(S))` is *not* invariant under that action — the real mathematical gap this closes, not just an alternate style.

Where `SPD<N,T>` keeps a point as `vech(log(S))` specifically so every `Manifold` operation collapses to flat arithmetic, `SPDAffineInvariant` keeps a point as the literal SPD matrix and lets the curvature show: `exp_map`/`log_map` are the standard closed forms `Exp_S(V) = S^(1/2) exp(S^(-1/2) V S^(-1/2)) S^(1/2)` and `Log_S(Q) = S^(1/2) log(S^(-1/2) Q S^(-1/2)) S^(1/2)`. Both reduce to the same eigendecomposition-based symmetric matrix log/exp `SPD<N,T>` uses (`S^(-1/2) V S^(-1/2)` is symmetric whenever `V` is), so this reuses `detail::apply_eigen_sym` directly — no new numerical machinery, same `N=2,3` restriction for the same reason. `distance()` is defined as `sqrt(metric_at(p, log_map(p,q), log_map(p,q)))` rather than the textbook direct generalized-eigenvalue formula — the general Riemannian identity (distance = norm of the initial geodesic velocity), keeping the three consistent by construction instead of risking independent hand-derivations drifting apart (see the `SE3::exp()` translation-Jacobian bug above for what that risk costs in this codebase). No `from_spd()`/`to_spd()` needed here: `PointType` already *is* the matrix.

SPD(n) under this metric is a Hadamard manifold (complete, non-positively curved, uniquely geodesic) — geodesics genuinely bend away from the SPD cone's boundary (`Log_S(Q)` diverges as `Q` approaches singular) instead of being straight lines in disguise. That extra fidelity costs the closed form: `frechet_mean_affine_invariant()` is a real fixed-point (Karcher-mean) iteration — average `log_map()` to every sample, retract, repeat until the average tangent vector is ~zero — genuinely exercising the exp/log pair for the first time in this file, unlike log-Euclidean's `frechet_mean()`. No ambient projection step needed (unlike `algebra::riemannian_minimize()` on `Sphere`/`Hyperbolic`), since this space's tangent space is already the unconstrained space of symmetric matrices; hand-rolled rather than calling `riemannian_minimize()`, which requires `HasNormal`, deliberately absent here for the same reason `SPD<N,T>` has none.

**A real, pre-existing bug in `detail::eigen_sym` found and fixed along the way** (not introduced by this feature — latent since the log-Euclidean commit, just never exercised): the 3×3 null-space-via-cross-product construction silently produced a **zero eigenvector** whenever an eigenvalue had multiplicity ≥ 2 (e.g. `diag(3,3,10)`, an axisymmetric-inertia-tensor shape) — not just the exact `S = cI` case. `A = S - λI` then has rank ≤ 1, so the cross product of *any* two rows vanishes regardless of which pair gets picked, and worse, a genuinely *repeated* root is only located by `solve_cubic()` to about `sqrt(machine epsilon)` (a double root's location is a square-root-sensitive function of the polynomial's coefficients — standard numerical-analysis fact, not solver imprecision), so the "vanishing" cross product isn't even close to exactly zero — it's small-but-numerically-unstable noise (~2e-7 in the `diag(3,3,10)` case), large enough to sail past a naive near-zero check while still being garbage. First exposed by `SPDAffineInvariant`'s exp/log-at-the-identity test (`S=I`, the fully degenerate `rank(A)=0` case) — chasing it down surfaced the more general `rank(A)=1` case too. Fixed by detecting multiplicity directly on the eigenvalues themselves (compare `values[i]` against every other root, with a tolerance scaled to `sqrt(epsilon<T>())` — matching the solver's actual achievable precision for a repeated root, not machine epsilon) rather than inferring it from the cross product's shakiness, then falling back to a construction that only needs the single largest-norm row of `A` (crossed with a slot-index-varied seed so two occurrences of the same repeated eigenvalue don't land on the same raw candidate before Gram-Schmidt gets a chance to separate them). The 2×2 path had the analogous exact-tie bug (`S = cI` making the "closest diagonal entry" tie-break pick the same eigenvector for both slots) — no solver-precision subtlety there since 2×2 eigenvalues come from `S`'s own entries directly, not a root-finder, so a plain equality check on `a == d` is enough. 3 new regression tests target `eigen_sym` directly (2×2 and 3×3 scalar-multiple-of-identity, 3×3 axisymmetric).

12 new tests (`tests/test_spd.cpp`): the 3 `eigen_sym` regressions above, concept checks, `contains()` (symmetric + positive-eigenvalue membership, unlike log-Euclidean's unconstrained log-space), exp/log roundtrip at the identity cross-checked directly against `SPD<N,T>::matrix_log_sym` (not just internal self-consistency), exp/log roundtrip at a non-identity base point (2×2 and 3×3), `distance()` symmetry, `distance()` against the independent generalized-eigenvalue closed form for a shared eigenbasis, and `frechet_mean_affine_invariant()` checked both against the log-Euclidean closed form (commuting samples) and against its own first-order optimality condition (non-commuting samples, where no closed form exists to check against directly). 842→854 tests, same one pre-existing unrelated failure noted throughout this file.

### SPD(n) outreach demo: relationship-matrix training, naive vs. geodesic (2026-09-01)

`examples/spd_relationship_demo.cpp` — first concrete piece of a planned outreach/contributor-recruitment demo (see project vision doc), scoped deliberately small: not the full game/event-generation idea, just the one structural claim underneath it, isolated and made legible. Same argument as "Mario Plays on a Manifold" (arXiv:2206.00106), which uses Riemannian geometry in a VAE latent space to guarantee generated game content stays valid against a Euclidean baseline that visibly breaks — here the "content" is a relationship matrix between two entities instead of a game level.

Setup: a 2×2 relationship matrix starts at the identity (uncorrelated). A sequence of similar events keeps nudging it in a fixed "strengthening" direction (an off-diagonal-only symmetric perturbation) — modeling repeated causally-similar shocks, not independent random rolls, applied two ways in parallel: naive plain matrix addition, and `SPDAffineInvariant::exp_map`'s geodesic retraction. Naive addition breaks positive-definiteness at step 13 (nothing about plain addition prevents an eigenvalue from crossing zero); the geodesic version is still comfortably valid at the same step. Enforced as a real regression (`tests/test_spd.cpp`, "Naive matrix addition can break SPD validity; SPDAffineInvariant::exp_map cannot"), not just a printed narrative — 14 steps, matching what's verified numerically healthy. Also emits an SVG line chart of `lambda_min` vs. step for both trajectories, meant to become the visual for whatever writeup eventually surfaces this to the geometry-processing/geometric-ML audience identified as the actual target (see project vision doc) — not r/cpp/HN as the primary channel.

**Real numerical-precision finding, not swept under the rug**: pushed past ~18 steps with the same fixed raw nudge, the geodesic trajectory's `lambda_min` asymptotically approaches the SPD cone's boundary and ordinary double-precision underflow there eventually returns NaN from the eigendecomposition. This is a genuine floating-point-precision phenomenon near a manifold's boundary, not a counterexample to the structural guarantee (in exact arithmetic a geodesic never reaches the boundary of a Hadamard manifold at all) — but it's real and worth remembering if this demo is ever pushed further: repeatedly applying the *same raw ambient tangent vector* at each new point isn't the same as following one smooth geodesic ray, and doesn't inherit that ray's asymptotic slowdown. Both the demo and the regression test stay within the verified-healthy range (14 steps) rather than papering over this.

Deliberately not yet built: the game/event-generation layer itself (player actions, entity hierarchy, causally-coherent generation) — that's the next, separate piece, once this one lands. "По немногу" (little by little), not one large PR.

---

## Backlog

Grouped by topic, not by version — an item sits here until it's ready to become real work, then moves to Completed above with the date it landed.

Each item carries a tag for *intent*, not difficulty or sequencing — Spatium's goal is to become the de-facto standard for computational work on arbitrary mathematical spaces, and that happens by earning the interest of people who know this territory well, not by chasing a release checklist:

- **[course]** — genuinely useful, actively where we're steering the project right now.
- **[want]** — we'd just like this to exist, or it'd be fun/interesting to build; no strong pull yet, no promise of order.
- **[dare]** — big, uncertain, or hard enough that it won't get built by us any time soon — a standing, genuine invitation for whoever wants to take a swing at it.

A PR against any tag here is welcome. So is a PR against nothing here — if you build something we didn't list, that's a real step toward Spatium becoming a standard in its own right, not a detour from this list.

### Manifold applications

- **[dare]** Fiber bundles (tangent/cotangent)
- **[dare]** Geodesic FEM (Laplace-Beltrami, heat equation)
- **[want]** Stochastic processes on Riemannian manifolds — Brownian motion / SDEs via Euler-Maruyama through `exp_map` retraction (drift + diffusion tangent vector, the same shape the existing Lie-group integrators already use). Reuses the `RiemannianManifold` concept directly, no new abstraction needed. Genuinely underserved outside bespoke research code — general manifold-ML libraries don't cover this.
- **[want]** Fisher-Rao information geometry as a `RiemannianManifold` — the manifold of probability distributions, metric = Fisher information, distance relates to KL divergence. Real ML utility on its own (natural gradient descent, distributional RL) for exactly the audience this project already names. Brainstormed 2026-09-09 as "a raytracer inside an ML model's mind" (light traveling through KL-divergence space) — the *space* is the genuine value; the raytracing-visualization framing is optional decoration on top and risks looking hollow if built before the underlying math means something concrete. Land the manifold first, decide separately whether a visualization earns its keep.
- Considered and explicitly parked, 2026-09-09: acoustic raytracing through the Ellis/Morris-Thorne wormhole (send a sound pulse down the throat, solve the wave equation against the curved metric, get a real .wav of "what curved spacetime does to a voice"). Spectacle, not utility — doesn't serve any of the audiences this project actually names, unlike the two items above. The underlying general capability (wave-PDE solving against a curved/non-Euclidean metric, not just the flat `wave_string`/`wave_membrane` domains) might be worth having on its own merits someday; the wormhole framing specifically isn't a roadmap item, just a joke that came up.

### GIS

- **[want]** Ellipsoid (WGS84) as Space — geodesic distance on Earth
