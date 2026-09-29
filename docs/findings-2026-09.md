# Findings, September 2026: audit, RSC, the binary, the root

One session's findings, written down so the next session starts from them
instead of from a transcript. Everything here was measured or read in the
code; each item says which. Two baselines were audited: `c064ad5` (the
library-wide audit, sections 1-3) and `fa12857` (RSC, the binary, sections
4-5). Items from the first baseline that later commits fixed are marked
**fixed upstream** where known; everything else should be re-checked
against current `main` before it is worked on, since the tree moves fast.

How to use this file: pick an item, re-verify it (most carry the probe that
found it), fix it with a test that fails on the old code, and delete the
item here in the same commit. When the file is empty it is deleted too.

Status tags: **fixed (patch)** -- fixed by the commit that added this file
or the one before it; **fixed upstream** -- fixed on main since `c064ad5`;
**open**.

---

## 0. The one pattern behind most of the bugs

Almost every serious bug below passed its tests, because **the test checked
a property the wrong answer also has**: SVD checked by reconstruction (a
non-orthonormal U reconstructs too), the DEC heat solver by the constant
mode and dt = 1 (where the wrong sign still decays), the symplecticity
verifier on diagonal planes only, polygon union with `>=` bounds, orbitals
at 1s only, where the wrong factor is 1. Rule adopted from this: **every
verifier and every numeric test gets a case it must fail.**

The second pattern: **absolute epsilons where the quantity scales**
(section 2). The third: **a comparison on `Dual` that reads the value
only**, so an entry that is zero with a non-zero derivative is treated as
absent -- found in `invert()` by this session, suspected elsewhere
(section 2).

And the pattern that finds them, stated as the project's rule: **consistency
comes only from automatic checks, and the cheapest check is a symmetry** --
two paths to one answer that must agree. None of this month's deepest bugs
was found by a test someone wrote for it; each was two paths disagreeing
(closed-form distance against |log|_g, `Dual²` curvature exposing
`invert()`), and a symmetry that held is equally a result (erased and typed
CCD charts: the same stops, reverts and piercings, section 11). It follows
that the more of the library is reused through its interfaces, the more
paths meet, and the more such checks exist for free: connectivity is not
only economy, it is coverage. The architecture in `architecture.md` is much
larger than the parts worked on; the drift away from it happened where it
was not re-read before new work.

---

## 1. Correctness (library audit at c064ad5)

| # | What | Where | Status |
|---|---|---|---|
| 1 | BVH traversal stack a fixed `std::array` with no bound; SAH depth unbounded (205 on 3000 geometric boxes), ASan overflow | `spatial/bvh.hpp` | **fixed upstream** (#63) |
| 2 | DEC heat solver runs heat backwards: `A = mass - L*dt` assumes a negative semi-definite L; `differential.hpp` says positive. 3.4e36 after 50 steps at dt=1e-3 | `mesh/dec.hpp:288` | open |
| 3 | LGVI divides by `J_j + J_k` instead of `J_i`: a symmetric body turns at half speed; the documented "~3% energy drift of Cayley" is this bug (2.3e-14 after the fix, probed) | `physics/mechanics/lgvi.hpp:94-96, 144-146`; docs/concept-driven-physics.md §6.2; ROADMAP | open |
| 4 | Hydrogenic radial normalisation uses `((n+l)!)^3` (old convention) with the modern Laguerre: ∫R²r² = 0.25 for 2s, 6.9e-5 for 3d | `physics/atomic/orbital.hpp:131-135` | open |
| 5 | `verify_symplecticity_drift` checks only the (q_i, p_i) planes: `p += S q` with S antisymmetric (not symplectic) reads 1e-12, like Verlet. Needs ‖MᵀJM − J‖ | `physics/mechanics/symplectic.hpp:144-163` | open |
| 6 | Segment-segment distance for N != 3 compares endpoints only: crossing segments give 1, overlapping plus-shaped polygons 4, star-of-David triangles 1.28 | `geometry/distance.hpp:81-111` | open |
| 7 | `Triangle<3>::contains` ignores the distance from the plane (a triangle 5 above another is at distance 0); degenerate triangle "contains" everything | `geometry/triangle.hpp:105-114` | open |
| 8 | `ParametricSurface`/`ImplicitSurface` `exp_map` and `log_map` are first-order -- a step in the ambient space projected back, the difference q - p projected to the tangent plane -- not the exponential and logarithm: exp(log) is not the identity, |log|_g is not the distance (the matrix: these cells L1 on every scalar that compiles). A missing pillar, not glue: geodesics on a chart need the geodesic ODE through its Christoffel symbols (plan stage 2.2). They did not compile at all until 2026-09-29 (`auto diff = q - p` captured an expression template) | `spaces/parametric.hpp`, `implicit.hpp` | open |
| 9 | Sphere chart normal points inward: `offset(sphere, +0.5)` has radius 0.5 | `spaces/chart.hpp:105` | open |
| 10 | Generic CCD tunnelled (distance from a local projection is an upper bound) | `narrow_phase.hpp`, `rigid_contact.hpp` | **fixed upstream** (#60, certified floor) |
| 11 | Second `.moving()` wraps both motions in an opaque lambda; exact form kept when it should be cleared; `MotionEnv::origin` lost | `io/build.hpp:900-911` | open |
| 12 | `content_hash` merged different shapes (fixed [0,1]² sampling, 4x4 grid, no equality check) | `io/build.hpp` | **fixed upstream** (#62) |
| 13 | Polygon booleans: [0,2]²−[1,3]² area 2 (should be 3); union of disjoint squares is the hull; CW input gives "no overlap" | `geometry/boolean.hpp:99, 231-292, 318` | open |
| 14 | `Polygon<3>` area/normal by fan from centroid (L-shape 6.33 vs 5); `project` of an interior point returns an edge point | `geometry/polygon.hpp:45, 74, 129` | open |
| 15 | `FiniteSet<FiniteSet>`: `<=` means subset, which is not a strict weak order; {{2},{1},{2},{3},{1}} has size 5; UB in sort | `discrete/finite_set.hpp:78, 84` | open |
| 16 | Schild's ladder transports the full-length vector: 50% error at |v|=1 over a quarter great circle | `mesh/transport.hpp:41` | open |
| 17 | `distance(Circle2, Circle2)` uses disk semantics (concentric r=1, r=3 give 0) | `geometry/distance.hpp:146` | open |
| 18 | `mesh_quality` max angle = π − 2·min (30-60-90 reports 120°) | `mesh/quality.hpp:71` | open |
| 19 | `mesh::transform(AffineTransform)` calls a nonexistent `.apply` with N off by one | `mesh/operations.hpp:39-45` | open |
| 20 | `Trace` not copyable while `is_copy_constructible_v<TraceNode>` says it is | `io/build.hpp:52-67`, `field.hpp:688` | open |
| 21 | JSON: short `\u` eats the closing quote; `1e`, `01`, `-.5` accepted; nan/inf written but not read; 200k nested `[` segfaults | `io/json.hpp` | open |
| 22 | Scene save/load keeps only `base_color` | `io/scene.hpp:410-418` | open |
| 23 | Truncated binary STL read as ok; OBJ face index out of range accepted; OBJ written with 6 fixed decimals | `io/stl.hpp`, `io/obj.hpp` | open |
| 24 | `Morphism::invert()` dereferences an empty optional | `morphism.hpp:40-42` | open |
| 25 | `Hyperbolic` subdivision projects vertically, not to the geodesic midpoint (area 9.42 vs 2.23 after 5 levels); fix: midpoint = `spaces::midpoint` | `spaces/hyperbolic.hpp:77`, `mesh/subdivision.hpp` | open |
| 26 | `Hyperbolic::contains` absolute tolerance rejects 147/256 of its own exp_map outputs at r ≥ 3 | `spaces/hyperbolic.hpp:35` | open |
| 27 | `GeodesicMethod` overload fails `static_assert` without Eigen even for Dijkstra | `mesh/geodesic.hpp:131` | open |
| 28 | `MeshTopology` silently overwrites a third face on an edge | `mesh/topology.hpp:86` | open |
| 29 | Tessellated sphere chart keeps a full vertex row per pole: χ = 0, two boundaries, zero-area triangles hitting the ±1e6 cotangent clamp | `spaces/parametric.hpp:366`, `mesh/primitives.hpp:178` | open |

## 2. Precision and scale (one family)

| What | Where | Status |
|---|---|---|
| `relative_epsilon(s)` falls back to absolute eps for s < 1: eigen_decompose(1e-15 A) stops at once, svd(1e-8 A) wrong and U not orthonormal | `core/epsilon.hpp:35`; `eigen_decomp.hpp:95`, `svd.hpp:113` | open |
| `solve_cubic/quartic` branch on disc vs absolute eps (disc ~ root^6): roots 1e-3..3e-3 become a triple 0.002; SPD<3> wrong on small eigenvalues through it | `algebra/polynomial.hpp:130, 143, 199`; `spaces/spd.hpp:83` | open |
| Absolute singularity thresholds: `(1e-5 I).inverse()` "singular" | `matrix.hpp:161, 175`, `linear_solve.hpp:46, 84` | open |
| `SO3::log` near π: 1.4e-4 error at π−1e-6; with Real50 a π rotation logs to 0 (`T{std::numbers::pi}` is double π) | `algebra/groups/so3.hpp:100` | open |
| acos/acosh distances lose precision at short range: d(p,p) = 1.5e-8..2.1e-8 on Sphere/Hyperbolic (measured by `test_space_access.cpp`); use atan2(|a×b|, a·b) and 2 asinh(|a−b|_L / 2) | `spaces/sphere.hpp:38`, `hyperbolic.hpp:41` | open |
| `epsilon<Dual<T>>()` is 1e-10, `numeric_limits<Dual>::max()` is 0: Dual runs take other branches than the primal | `core/epsilon.hpp:18-26` | open |
| **`invert()` skipped `factor == 0` rows; for Dual that compares the value only, so derivatives of the inverse were lost on every coordinate plane** -- exact Kerr read a Ricci residual of 2.5e-5 there | `algebra/linear_solve.hpp:97` | **fixed (patch)** |
| Same `== T{0}` pattern, same risk for Dual: `solve_linear`, `solve_quadratic`, `solve_cubic`, `solve_quartic` degenerate branches | `algebra/polynomial.hpp:61, 71, 105, 180` | open |
| `Complex::is_real` absolute 2.2e-16 on the imaginary part | `algebra/complex.hpp:59` | open (#60 reports the CCD symptom did not reproduce) |
| ImplicitSurface FD step h = 2.8e-11 absolute: at |x| = 1e6 the normal is (0,0,0) | `spaces/implicit.hpp:99` | open |

## 3. Connectivity: where "any space gets everything" stops

Measured by `docs/connectivity.md` (first run 2026-09-29), beyond items
listed below:
- `Real50` fails in `solve_quadratic` (SPD affine-invariant's eigenvalues)
  and `fmod` on a Boost expression (`ParametricSurface`'s periodic wrap).
- SPD affine-invariant over Dual: the derivative cell runs to a non-finite
  value, alone and as a product's factor -- its eigen-decomposition's
  derivative, not yet looked at.

- Mesh calculus is Euclidean-ambient: cotangents by `cross`, `face_area_3d`, `heat_geodesic` ignores the space. On `Hyperbolic<2>` heat errs 342%, `Mesh::area` 22.0 vs 2.23; Dijkstra on the same mesh is right, so the two methods answer different questions. Fix: intrinsic edge lengths (law of cosines, Heron). Doc `concept-driven-physics.md:147` promises a `metric_at` Laplacian that does not exist.
- SO3/SE3 are LieGroups but not Riemannian spaces: add `LieGroupManifold<G>` (exp_map = p·exp(tv), log_map = log(p⁻¹q)); then `frechet_mean` averages rotations. With `core/access.hpp` this can be ADL functions instead of an adapter type.
- `Dual` does not flow through geometry: 58 unqualified-free `std::abs/std::sqrt` calls; `ray_triangle`, `Triangle::distance`, seg-seg, `ray_torus`, `BVH<Triangle<3,Dual>>` fail to compile with Dual. Dual lacked `T/Dual`, `log`, `floor`, `atan2`, `sinh`, `cosh` (added 2026-09-29; `pow(Dual, int)` still goes through the real-exponent overload).
- SPD limited to N ≤ 3 by a stale reason (general `eigen_decomp.hpp` exists); `distance` costs ~4 eigendecompositions where 2 suffice.
- `riemannian_minimize` requires `HasNormal`, excluding SPD and products.
- `Scalar` requires `totally_ordered`: no `Complex`, no `Interval`.
- BVH only N = 2, 3; `RayHittable` checked against `Ray<3>` for any N.
- Islands in include/: `units.hpp` unused by any integrator; `Point<Space>`, `Morphism`, `lazy()`/`TransformExpr` unused inside the library (only tests).
- Four composition mechanisms without bridges: `Morphism`, `lazy()`, `VecField`/`.moving()`, `Instanced`.
- Duplicates to fold into one owner: Dijkstra (geodesic.hpp = voronoi.hpp; scatter re-runs full Voronoi per site), finite differences where Dual exists (5 places), RK4 (ode.hpp vs integrator.hpp), basis-from-normal (6 copies), Fréchet mean (spd.hpp vs calculus.hpp), parametric tessellation (2 copies), bounding boxes (3 types), factorial (2), thread pools (2).
- `geodesic.hpp`'s `killing_energy`/`killing_angular_momentum` assume BL index 3 = φ; wrong for Cartesian Kerr-Schild.

## 4. RSC at fa12857

**What it is in code.** Three things sharing a registry: an oracle-labelling
factory (candidates run against a reference -- the valuable part), a search
lab (chains over numbers, meshes, CCD primitives), and an MLP dispatcher
(the least important part). **Nothing in include/ consumes a learned choice
at runtime**: `Chooser` has one consumer, the chooser overload of
`surface_ccd.hpp`'s `first_contact`, called only by a test;
`load_embedded_base` is called only by a test.

Accuracy vs the feature-only Bayes ceiling (probe): rootfind 0.966 / 0.965
(0.989 with `a` as a feature); ODE 0.97 / 0.985; integrator 0.84 (CE 0.94) /
0.977; linear 0.8685 / 0.889; mesh 0.78 / **0.983** (exact lookup);
geodesic 0.536 / **0.553** (mesh level is the only feature). The models
sit at the ceilings their features allow: features and labels limit, not
learning.

| # | Breakage | Where |
|---|---|---|
| B1 | ODE/integrator labels are not cheapest-correct: the rule compares to the *last candidate*, not the reference. ODE: in 72% neither Euler nor RK4 is within tolerance; 18.9% labelled Euler with RK4 >10% off | `comparison_task.hpp:93-101` |
| B2 | Geodesic capped by its single feature (0.536 = "Dijkstra at vc ≤ 162, Heat at ≥ 642"); README reads it backwards | `geodesic_task.hpp` |
| B3 | Rootfind plateau: Newton's outcome depends on `a` and sign(x0), which the feature drops | `rootfind_task.hpp:376` |
| B4 | Mesh plateau is the MLP, not the data (lookup 0.983, MLP at the best single threshold 0.77) -> fit trees on labels | `mesh_task.hpp` |
| B5 | Linear test `post > pre` passes on 0.8655 → 0.8685; its premise is inverted (Jacobi dearer at every size tried) | `test_rsc_linear.cpp:113` |
| B6 | One-sided gate exercised only on toy candidates; `choose_gated` falls back to the last candidate silently; `with_cost` untested | `comparison_task.hpp:109, 119` |
| B7 | Generated CCD header: tree reproduces bit for bit, header numbers stale after #90, output-path argument missing from its provenance line, vertex-face worse than default in its own table; regeneration ~53 CPU-minutes | `rsc/include/generated/ccd_chooser_tree.hpp` |
| B8 | "A setting changes cost, never correctness" is false under a finite budget: seed 1 query 64, ball-only policy spends 2^22 pairs and answers hit at 0.9075, default answers miss; 4/720. `ccd_policy.cpp:100` checks the wrong policy's budget | `core/concepts.hpp:202-205` |
| B9 | Reproducibility checks protect nothing: `embedded_base_is_current()` always false, `load_embedded_base` validates nothing, snapshot does not pin feature semantics (`feature_mean` depends on the seed: −1.344/−1.336/−1.369) | `checkpoint.hpp`, `rootfind_task.hpp:364` |
| B10 | Cross-entropy measured better (2000 updates vs REINFORCE's 16000) but every trainer and test still uses REINFORCE | `train_base.cpp:33` |
| B11 | README stale: synthesis "not attempted" (three tools exist), "deployment: a tree per domain" (none exists), no mention of CCD or `Chooser` | `rsc/README.md` |
| B12 | (`surface_ccd.hpp`'s dangling reference to the chooser header: fixed (patch)); geodesic test `post > 0.5` passes for a constant; `ccd_synthesize` picks winners in-sample; its hybrid uses 2^16 where the library uses 2^20; `dispatch_cost` numbers (260x/3.7x) do not reproduce (1102x/8.5x) and it times random nets, not the real ones | various |

**Not attempted:** calibration (no RSC path calls `minimize`), surrogates
(no Padé/Chebyshev/fits anywhere), IR rewriting (the bit-exact oracle
exists), structure building (splits are fixed: `NUM_BINS=12`, `kLeaf=4`).

**Choices in include/ with no hook:** geodesic step (fixed RK4, caller dλ),
black-hole ray step `clamp(0.05·hd, 0.01, 1.5)` ×0.35, capture at hd < 1.02,
escape at 3·|cam|; `render_level()` = "exact if available" (measured wrong
past ~100 objects); tessellation 48×24; `GeodesicMethod` by caller;
rigid_contact budgets 2^20/2^24/64; xpbd 4, lgvi 8, ray_parametric 24,
calculus 200/100; BVH bins 12, ball tree leaf 4.

**Re-framing adopted.** RSC is the library's *error-budget economy*:
subsystems report their error as a certificate; RSC allocates the budget of
the final result (a pixel, a time of impact) among them and chooses where to
refine; learned models only predict which choice will pass, never hold
correctness. Loop: **Refine, Select, Certify** (one possible reading of the
name). Consequences: a `Certified<T>{value, bound, certified}` next to
`Result<T>`; `Chooser` extended to a budget (`choose(features, eps)`) with
a certifier that escalates; labels from oracles, one-sided, measured cost,
a "refuse" class; cross-entropy and trees, deployed as generated headers
held by a test (the #89 pattern). Three searches on one machine: computation
settings (CCD, done), type compositions (the connectivity matrix, below),
metric forms (the binary).

## 5. The binary black hole at fa12857

**What exists.** Metric as ten IR fields (`metric_field.hpp`), Kerr BL and
Kerr-Schild, both vacuum-tested; GLSL with exact first derivatives; one
clock (a motion's time is coordinate 0, so positions are retarded along each
ray); pair = η + Σ f_i l_i⊗l_i, unboosted, spins about z; Newtonian orbit
shrinking by Peters to `stop`, then circling there forever (merger not
modelled, said in the header); dust as test particles on the metric's
geodesics.

**Fixed (patch):**
- Disk and dust redshift numerator was p·e0 evaluated at the emitter; for a
  pair p_t drifts 40-60% along a ray, so g was off 12% at 14 M and 25% at
  8 M, brightness (∝ g⁴) by 1.6-2.4x. Now the camera value, 1 by
  construction. Identical to 6e-6 for one hole. Run on the Iris Xe: `--check`
  unchanged (device against double, worst 2.2e-6); the default binary frame
  at 960x540 barely moves -- mean luminance 33.36 to 33.34, 0.5% of pixels
  change by more than 2 levels, median ratio over lit pixels 1.000 -- since
  the Doppler damping (0.6) and the tone curve compress g⁴. The error is in
  the physics the picture reports, not in the picture as it is shipped.
- `vacuum_residual` is exact (nested Dual; the old one kept as
  `vacuum_residual_fd`): Kerr 1e-17..1e-20 against 1e-11..1e-14.
- `lorentzian_at()`; `vacuum_residual` returns NaN and `residual_at` an
  error where the superposition is not a spacetime (det g ≥ 0 between two
  equal holes at d ≤ 4 M), instead of the 0 a silently zeroed inverse gave.

**Measured, open:**
- Residual (√|R_μν R^μν| / √|Kretschmann|) at the midpoint: 37% at 6 M, 23%
  at 10 M, 17% at 14 M, 7.8% at 40 M. Near a hole it is ≈ 0.78·v -- the
  missing boost dominates; **boosted KS terms cut it 18x at 14 M** (probe,
  ~40 lines of IR); the midpoint keeps a nonlinear floor (0.169 → 0.126).
  The header's "a few percent" for the missing boost is ~10% measured.
- Ray numerics are 10-1000x over-provisioned (worst direction change 6e-5
  rad vs a 1080p pixel of 8.6e-4); the error budget is model error. 2/1296
  rays grazing a horizon at 8 M break the null condition (0.91).
- Exact oracles available: null constraint g(k,k) (almost free); helical
  invariant p_t + Ω L_z for a circular pair (conserved to 2e-5); Ricci
  focusing κ_R along a ray (magnification error ≈ 2κ: median 1.6%, p90 15%
  at 14 M); exact residual per event (~29 µs).
- `christoffel()` still substitutes a zero inverse on failure (now guarded
  in the residual path only; the geodesic path should return a Result).
- Escape at 3·D_cam leaves ~2 px of weak-field deflection; use the analytic
  tail M·b/s².
- Doppler damping 0.6 is an artistic knob presented beside physics: split
  settings into physical and artistic, show which are on.
- Majumdar-Papapetrou (U = 1 + Σ m_i/r_i) is an exact static binary
  expressible in today's IR: a test metric for multi-centre machinery.

**Gap to "inspiral, merger, ringdown, everything".** Fair as ~1/3 for a
convincing picture with a measured, bounded error; ~15% for physical
fidelity. Next, in order: boosted terms (S); arbitrary spin axes (S);
PN + spin-orbit + precession on an adaptive ODE, with an **IR table/spline
op with a derivative** so the metric can read an integrated trajectory (M);
remnant mass/spin/kick fits (S; equal mass: M_f ≈ 0.952 M, χ_f ≈ 0.686) and
QNM ringdown visual ((2,2,0): M_f ω ≈ 0.533 − 0.081 i, τ ≈ 12 M_f) (S-M);
merger as a single-KS morph g = η + H k⊗k with η-null k -- Lorentzian by
construction, residual shown (M-L, research); matched-asymptotics metrics
(Mundim/Ireland) later (L). The merger metric itself exists only in full
numerical relativity; surrogates (NRSur7dq4, EOB, IMRPhenom) give
waveforms and trajectories, not a near-zone metric. Rendering a
certificate overlay (colour = residual or κ_R along the ray) makes the
picture show where it holds.

As synthesis: the metric is a search over compositions of primitives
(terms, boosts, spin rotations, time maps, blends, PN paths, remnant fits)
plus calibration of their parameters (differentiable through one more Dual
level), graded by the exact residual **and anchors** (total mass and
angular momentum at infinity, PN far limit, remnant Kerr late) -- zero
residual alone only says "some vacuum". Events sampled by an adaptive tree
over spacetime, the same refine-earliest loop CCD uses.

Concept gap the scene exposes: spacetime is pseudo-Riemannian (no distance,
log not general); the hierarchy stops at `RiemannianManifold`, so
`MetricField` lives outside `spaces/`. Add `PseudoRiemannianManifold` with
exp via Christoffel.

## 6. The root: from implementations back to interfaces

`architecture.md` (unchanged since the first public commit) describes an
STL-shaped library: concepts as the access vocabulary, free algorithms,
lazy adaptors (Vec expression templates, `lazy()` transforms, `Morphism`,
pipes), materialisation at the edge. What grew instead: laziness replaced by
type erasure (`std::function` in 17 headers: `ParamFn`, `ImplicitFn`,
offset thickness, `project_fn`/`normal_fn`, narrow_phase, surface_ccd,
Morphism) and eager materialisation (36 public signatures take a concrete
`ParametricSurface<T>`/`ImplicitSurface<T>`/`Mesh`), then laziness
reinvented at runtime in `io/field.hpp`. Counter-example that proves the
point: `relativity/geodesic.hpp` stayed templated ("never std::function")
and nested Dual worked there with zero changes.

Five pillars, in migration order:
1. **Customization points with derivations** -- started by `core/access.hpp`
   (member → ADL → derived). Next: move algorithms onto them one at a time;
   derive `exp_map` from any metric callable via Christoffel + ODE for any
   dimension (generalise `geodesic.hpp` off `Vec<T,4>`).
2. **One laziness model, two backends**: `ParametricSurface<F>` /
   `ImplicitSurface<F>` templated on the callable with an explicit `Erased`
   alias at boundaries; the Field IR moved into core as the expression
   language (static templates ↔ IR pool via `lower`), `TransformExpr` and
   `Morphism` folded into it.
3. **Views and materialisers**: offset, product, restrict, pullback,
   quotient (flat torus = R²/Z²), moving, compactify; `to_mesh`, `lower`,
   `sample`, `cook` at the edge; `tessellate` moves to mesh/ (breaks the
   spaces ↔ mesh cycle); mesh algorithms over a `DiscreteSurface` concept.
4. **Policies and certificates** as parameters (section 4).
5. **Hierarchy**: exp as the base, log optional; `PseudoRiemannian`;
   numeric tower Ring → Field → OrderedField → RealLike; `dynamic_extent`;
   semantic requirements written down; archetypes.

CI metrics that keep it honest: share of public templates constrained by
concepts vs taking concrete library types; no `std::function` in
core/algebra/spaces/geometry/spatial; archetype instantiation of every
algorithm; the **connectivity matrix** -- one translation unit per
(scalar × space × algorithm) cell, bodies instantiated, levels L0 compiles
/ L1 runs finite / L2 axioms / L3 matches a reference; green share as the
number behind "any space gets everything". A generator for 180 cells
(5 scalars × 6 spaces × 6 algorithms) was written this session but not run;
it is a small Python script emitting one `main()` per cell -- rewrite into
`scripts/` when picked up. `Dual<Dual>` curvature was found by accident;
this finds such connections on purpose.

## 7. Continuous collision: the quality bar

Never late, never a miss -- certified, not statistical; external benchmark
(Wang et al. TOG 2021 sample queries, TDIB) with false negatives 0 reported
next to false positives and time; an adversarial set (spikes, grazes, speed
≫ size, scales 1e-6..1e6, starts on the surface, concave shells); fuzzing
against an interval oracle; no stalls; long runs with an independent
penetration check (`feat/napkins` does this: 3 545 piercings without the
continuous pass, 0 with -- that branch is unmerged); the same on Sphere and
Hyperbolic via exp_map. Budget exhaustion must surface (B8).

## 8. Tracing on manifolds (design notes; main has implemented much)

Implemented on main since: `Bound` concept, `BVH<Shape, Bound>`,
`GeodesicBall` with closed-form ray intervals on E/S/H, SAH by the boundary
measure (Crofton: 2π sin r / 2π sinh r in 2D, 4π sin² r / 4π sinh² r in
3D -- not the enclosed volume), `geometry_bounds`, `GeodesicBallTree`.
Still open: chunked mode for integrated geodesics -- a tree lower bound on
distance truncated at depth k is still a valid lower bound, so sphere
tracing and BVH fuse into one marcher; `ManifoldCamera` as a frame from
`metric_at` or an isometry; light via `log_map`, shadow rays with `t_max`;
recentering in H (the Minkowski form loses all digits near d ≈ 20). The
injectivity/convexity radius matters for building and merging balls and
for CCD swept volumes, not for pruning, which the triangle inequality
makes correct at any radius.

## 9. Using the project

- Docs split by reader: README (what, install, three 10-line examples,
  stability table), Guide ("your space in 20 lines" -- now three ADL
  functions), Reference (concepts, customization points, derivations), and
  a Journal (today's ROADMAP and rsc/README as they are). ROADMAP shrinks
  to an index.
- Stability table: stable core/algebra/spaces/geometry/mesh; experimental
  io/build, physics/relativity, render/gpu*, rsc.
- `v1.0.0` is far behind main: semver, CHANGELOG from PR titles, monthly
  release.
- **The compiler requirement is overstated**: README says GCC 15+/Clang
  19+, while all of main with tests builds and passes on GCC 14.2 (this
  session). Test it in CI and lower it.
- CI additions: build matrix (GCC 14/15, Clang 19; Debug/Release;
  noboost/noeigen), ASan+UBSan on tests, connectivity matrix, generated
  header freshness, a `std::function` lint.
- Stale: `Point`/`Morphism`/`lazy()` unused; modules build behind; units
  island; dispatch_cost numbers; README claims ("distance for all pairs",
  GCC 15+).
- Missed: Compiler Explorer links for every README snippet (header-only is
  ideal for it); vcpkg/Conan port; a docs site from docs/; issue templates
  ("my space does not fit"); good-first-issues from section 1.
- 28 of 29 remote branches are merged; `feat/napkins` is the only one ahead
  of main.

**Web solver (rsc/README "Site").** The gate is open: WASM is compiled C++
(commit 42c55af), Boost.Multiprecision is header-only. GitHub Pages cannot
set COOP/COEP, so no SharedArrayBuffer and no WASM threads -- single thread
first, or the coi-serviceworker workaround. The honest core is a parser
(the description door), certified computation and a transparent method
choice -- fix B1 before the site explains choices in public. MVP domains
with certificates already: polynomial roots with f64 → Real50 escalation,
distances/geodesics on S²/H² (and user spaces via `core/access.hpp`), "will
these collide" with a certified time of impact. For the gallery, a WGSL
emitter beside `field_glsl`/`metric_glsl` is one more backend of the same IR.

**GPU.** Two paths exist: legacy CUDA in `gpu/` (sympy-transcribed
Christoffels, one BL Kerr) and Vulkan compute from the IR (exact derivatives
without sympy). Keep one; CUDA at most as a cross-check. Device has no
fp64 (Iris Xe); fp32 passes the gate (3.5e-4 at the critical curve). The
ROADMAP's next step, an IR interpreter in GLSL, removes shader recompiles on
scene change.

## 10. Media, as decided

Not an SDK sale: sell the new physics, let engineers find the library under
it. Show the impossible, never list features; every piece: a claim that
sounds impossible, the description beside the result, a measured counter
on screen, a control showing where the usual approach fails, a live knob.
Audit findings stay internal. Candidates: "break the physics" sandbox
(napkins, a public counter of attempts vs piercings = 0), two million
objects from 34 nodes, frame 10^6 without stepping, a scene in a QR code
reproduced bit for bit in a browser, one metric function spinning a black
hole, a chain found on a sphere replayed on a torus, a picture that shows
where it is right (certificate overlay).

## 11. Carried from the conversation, not in the first write-up

Decisions and measurements that were in the pasted reports or this month's
CCD work but not in the sections above. Same rules: re-verify, fix with a
failing test, delete here.

**RSC, cheap experiments on the existing machinery** (from the RSC report):
1. Cross-entropy everywhere and trees fitted on labels -- mesh 0.78 → ~0.98.
2. `a` as a rootfind feature (+2.4%).
3. Per-vertex features for geodesic.
4. An absolute gate and a "refuse" class in ODE and integrator (B1).
5. `with_cost` in linear -- it will likely collapse to "always direct", an
   honest answer.
6. A `Chooser` on `sweep_sphere_surface` from `ccd_synthesize`'s winners.
7. `explicit_contact_window`'s explicit-or-Newton window as a dispatcher.
8. The render level by the measured ~100-object boundary.

**RSC contract, as sketched.** Three layers: *certificates* in the library
(`Certified<T>{value, bound, certified}` or a `Bracketed<T>`, next to
`Result<T>`; CCD's `certified` flag and #90's bracket are the first),
*decisions* by RSC under a budget, *prediction* by learned models only for
the order of refinement and the first choice.

```cpp
concept BudgetedChooser = requires(const C& c, const F& f, double eps) {
    { c.choose(f, eps) };   // a choice under a tolerance
};
// a Certifier checks bound <= eps; if not, Refine
```

One loop for every domain -- the one `first_contact` already runs: refine
the earliest or the most uncertain, discard what a certificate settles, stop
when the bound is within budget. Labels from oracles, one-sided, with
measured cost and a refuse class; deployed as generated headers held by a
test. B8 is why the contract has to change: under a finite budget a setting
does change the answer.

**Certificates per quantity, for the binary** (all cheap next to a frame):

| Quantity | Certificate | Cost |
|---|---|---|
| metric at an event | signature, det g < 0 (hard); exact residual through `Dual²`, over local curvature (soft) | 1 evaluation / ~29 us |
| orbit | energy balance dE/dt = -F; size of the last PN term | cheap |
| ray | null condition g(k,k) (free from the first RK stage); step doubling; helical invariant p_t + Omega L_z for a circular pair (to 2e-5) | <= 2x |
| pixel | Ricci focusing kappa_R along the ray (magnification error ~ 2 kappa); redshift at the emitter against the camera | 2-5x on a preview |
| capture | sensitivity to the threshold; growth of the null condition; later a horizon finder | cheap / M |
| escape | the analytic weak-field tail M b / s^2 | free |

Residual without the boost, near a hole: 12% at 10 M, 10% at 14 M, 4.4% at
80 M; at the midpoint, 5% at 80 M (the rest in section 5).

**The binary uses all six kinds of RSC:** method choice (the metric's
regime per event: superposition, boosted, matched asymptotics, remnant Kerr;
ingoing or outgoing KS); calibration (the ray step constants, blend widths,
capture thresholds, against the host double oracle); surrogates (remnant
mass, spin and kick fits; QNM frequencies and damping; PN coefficients --
the first real use of the kind); synthesis (the chain of regimes in time,
PN → morph → ringdown, graded by the residual); IR rewriting (a smaller
metric pool for GLSL, bit-exact, the oracle exists); structure building (the
pixel/ray refinement tree, LBVH and balls for dust, splits by measured
time).

**Why this scene is the peak.** It is the one scene where space itself
moves, not objects in it; it has no closed form, so it forces search and
certificates; it touches every subsystem at once, the longest path through
the library's graph; infinity is physically in it (asymptotic flatness, the
horizon's infinite redshift, the sky at infinity); and it has an outside
truth to be checked against (NR catalogues, LIGO waveforms, remnant fits).

**Merger as synthesis, its limits.** Zero residual means *some* vacuum, not
ours -- anchors are required (mass and angular momentum at infinity, the PN
limit far away, the remnant Kerr late), or the search finds empty space. A
family of formulas has a residual floor at merger: the search finds the
best in the family and shows where the floor is; the exact solution lives
in another search space, discretisations (numerical relativity). The
single-KS form g = eta + H k⊗k with eta-null k is Lorentzian by
construction and turns the search into one for a function H, represented
finitely as an IR expression. For a picture less is needed than for a
solution: a metric and a certificate on every ray -- an image with a proven
map of its own error.

**Infinities, which the library already has without naming them.** The
infinitesimal is `Dual` (e^2 = 0), `Dual<Dual>` second order, a jet --
which is why `Dual²` simply worked. The infinitely large by
compactification: infinity moved to a finite boundary (conformal
compactification in GR; the render's exit at 3 D_cam is a missing infinity,
to be replaced by the analytic tail; the hyperbolic ideal boundary is where
tiling rays go; projective coordinates). The infinite-dimensional as
function spaces: a `Field` expression is an element of one represented
finitely, a `MetricField` a section of a bundle, the merger search a search
over metrics. "Infinities on a PC" means a lazy or symbolic representation
plus a certified bound.

**CCD's configuration space as a space.** The two-surface search runs in
(u_a, v_a, u_b, v_b, t) -- chart × chart × time, a product of spaces --
implemented as boxes and indices. As a `ProductSpace` with its own metric,
cells would be geodesic balls (the ball tree exists) and the Lipschitz
constant a local norm. The first probe of "CCD on the library's spaces".

**What survived of the original style**, to be made the core's law rather
than chosen per subsystem: `point_to` as ADL overloads; `chart_of` as the
DSL's ADL entry; the `Bound` concept and `BVH<Shape, Bound>`, geodesic
balls; `Chooser`, a policy as a parameter; the Field IR with `lower` and
GLSL, a lazy description with several backends; `Placed<T>`, a mesh only on
request; the templated metrics in relativity.

**Derivations the customization points should supply** (`core/access.hpp`
has the first four):

| Operation | Derived from |
|---|---|
| `distance` | the norm of `log` in the metric |
| `midpoint`, geodesic interpolation | `exp(log / 2)` |
| parallel transport | Schild's ladder through `exp`/`log` |
| `exp_map` for any metric function | Christoffel symbols through `Dual`, then an ODE (works in relativity; to become the default for all) |
| curvature | `Dual<Dual>` on the metric |
| `normal` | a chart's `d_u × d_v` |
| `bounds` | a support function `support(dir)` |
| `project_tangent` | the normal |

Numbers behind the diagnosis: `Vec<T,3>` appears 612 times in include/
(`rigid_contact.hpp` 68, `spherical_polygon.hpp` 51, `narrow_phase.hpp` 45,
`ray_surface.hpp` 39, `field.hpp` 37). A fifth CI metric beside section 6's:
a count of early materialisation -- a `Mesh` built inside a function that is
not a materialiser.

**CCD, open from this month's work** (numbers in ROADMAP):
- *To drop IPC from the dependencies* the CCD is not the missing part;
  missing are friction in cloth-cloth contact, a response other than
  stopping short (it kills the velocity), a resting stack that is still (the
  napkins' kinetic energy dithers at 1-30 without settling), and the
  head-to-head on the same scene -- `cloth_sphere_probe` already runs
  ipc-toolkit behind `SPATIUM_IPC_TOOLKIT`.
- *Napkins, the limit to reach*: no piercing and no revert over 10 s
  including settling; <= 4 ms a substep for 16 napkins of 25x25 on the
  laptop (60 fps, 4 substeps); the ipc-toolkit comparison; energy decay and
  stack height; a count of napkin pairs left entangled. Today: 8 of 19x19
  thrown at 8 m/s, 12 ms a substep on 8 threads, 0 piercings.
- *A ball against a curved surface* is a flat minimum that both CCD paths
  pay for (155 ms exact, or 6.7 ms and 74 answers early); the fix is the
  surface offset by the radius, where the minimum is sharp -- a torus stays a
  torus, a chart needs curvature bounds.
- *The bracket's witness* (two cells' centres within tolerance) is weak on a
  contact along a curve (a torus flat on a plane: 4.5M pairs at 1e-6); a
  centre within tolerance of the other cell's slab would witness at once.
- The search still allocates a cell vector and a queue per query: what is
  left after typed charts (7-15%).
- *The napkins repeat IPC's limitation*: the cloth is a triangle mesh with
  vertex-face and edge-edge queries, as in IPC, whose barrier and CCD are
  defined on a distance but computed only between triangles. Ours are
  defined and computed on charts. Next: the napkin as a smooth patch
  (bicubic, control points as degrees of freedom -- cf. "Simulating
  Parametric Thin Shells by Bicubic Hermite Elements", arXiv 2312.14839),
  contact through `surface_ccd` between patches, where it measured 7-20x
  TDIB-CCD; the triangle scene stays as the comparison.
- *Energy.* IPC's implicit Euler dissipates energy; it is robust, not
  energy-conserving. The library's variational and symplectic integrators
  (`variational.hpp`, `lgvi.hpp`) with a barrier on the true distance are
  the route to both: no intersection and energy held.
