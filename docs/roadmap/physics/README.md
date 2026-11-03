# Physics

## Where it stands

Mechanics on spaces and charts (`physics/mechanics/`): units, integrators, symplectic and variational schemes, certified continuous collision computed on the analytic description, gravity of many sources on closed spaces. ipc-toolkit is no longer a dependency; it is the comparison.

## Now

- nothing is taken up here now

## Later

- the library's own IPC, on smooth patches rather than triangles

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area physics`.

## Reasoning and earlier history

### Physics v2: Mechanics, Symplectic/DEC/Lie-group, Variational, Contact (2026-04-24)

Top-level abstraction: **discrete Lagrangian variational integrator on Riemannian manifold with DEC forms** — not the Bullet-style "RigidBody + Force + Constraint" tree.

- **Mechanics foundations.** Compile-time SI via `Quantity<M,L,T,I,K,N,J>` + `std::ratio`, `PointMass<N,T>`, `RigidBody<N,T>`, composable forces (`UniformGravity`, `PointGravity`, `Spring`, `Damper`), integrators (`euler_step`, `semi_implicit_euler_step`, `verlet_step`, `rk4_step`). Kepler orbit over 10 periods under Verlet: energy + L conserved to 1e-3.
- **Symplectic + DEC + Lie-group.** `SymplecticManifold<S>` concept + `CotangentBundle<M>` wrapper + `verify_symplecticity_drift`, `Yoshida4` composition integrator (KAM 1000-period proof), `lie_rkmk4_cf_step` (commutator-free 4-th order RKMK, Celledoni-Marthinsen-Owren 2003), `PointOnManifold<M>` + analytical `geodesic_step<Sphere<N>>`. DEC primitives in `mesh/dec.hpp`: `Form0/Form1/Form2`, `exterior_derivative_0/1`, `hodge_star_0/1/2`, `laplace_beltrami_dec`.
- **Variational + continuum.** `DiscreteLagrangian<LD,N,T>` concept, `SeparableMidpointLagrangian`, `variational_step_separable` (kick-drift-kick Verlet as closed-form DEL). LGVI on SO(3) via Cayley 1-cut (orientation + Noether + Casimir exact; energy drift ~3 % is a Cayley-1-cut limitation documented as follow-up). `DecHeatSolver` (backward-Euler `(M + dt·L) φ = M·φ_n` pre-factored via Eigen `SimplicialLDLT`). `DeformationMap<MFrom,MTo>`, `deformation_gradient` (FD), `right_cauchy_green`, `green_strain`, `StrainEnergy<W,T,N>` concept, `SaintVenantKirchhoff`.
- **Contact (partial).** `ipc_barrier(d, d̂)` = -(d-d̂)²·log(d/d̂) with gradient + Hessian + default stiffness, `XpbdParticle`/`XpbdDistanceConstraint` + Gauss-Seidel projection + `xpbd_step`, `build_distance_constraints`, `build_bending_distance_constraints`. **`ContactSurface<S,T>` concept + free `point_to(p, surf)` overload set** for `Sphere<2>` and `Torus` (closed-form) + generic `Surface` fallback that picks up `ParametricSurface`, `ImplicitSurface`, and any user type satisfying the Surface concept. `ipc_contact_energy`/`ipc_contact_force_on` concept-constrained helpers. Full implicit `cloth-on-obstacle` pipeline (and Klein-bottle USP demo) deferred pending ipc-toolkit integration — explicit position-based contact ceiling was hit with overhanging corners yanking cloth through the obstacle faster than the contact band reacts. (Continued: "Contact physics build dependency", "Implicit-contact Newton solver" below.)

#### Tests + benchmarks

- **673 test cases, 6 549 Catch2 assertions**, all passing on modules ON build (GCC 15 + Eigen 3.4), at this point in the project's history.
- `benchmarks/bench_narrow_phase.cpp`: `point_to_sphere` 3.7 ns/op, `point_to_torus` 17.7 ns/op, `point_to` on a parametric surface (Newton-on-UV via `Surface::project`) 5 245 ns/op, full IPC pipeline (query + energy + force) 12.6–25.5 ns/op for analytical surfaces. Order-of-magnitude faster than GJK (~100-500 ns/op) on closed-form paths.

#### Reference specification

- `docs/concept-driven-physics.md` — a rigorous specification of the concept hierarchy (Set → Manifold → Surface + `SymplecticManifold`, `DiscreteLagrangian`, `ContactSurface`) and the dispatch pattern. Fixes the current status matrix and limitations.

### Contact physics build dependency (2026-08-28)

`ipc-toolkit` (implicit IPC: Newton + log-barrier + line-search filter, Li-Kaufman 2020 §5) is now an optional CMake dependency: `SPATIUM_IPC_TOOLKIT` (default OFF), `FetchContent`-declared against ipc-toolkit v1.6.0, matching `SPATIUM_EIGEN`'s opt-in pattern. Linked directly onto consuming targets (`spatium_tests` today) rather than the exported `spatium_sdk` INTERFACE target — unlike `Eigen3::Eigen` (an IMPORTED `find_package` target), `ipc::toolkit` is built locally via FetchContent and isn't part of any install export set, so linking it into `spatium_sdk` breaks `install(EXPORT SpatiumTargets ...)`. Verified end-to-end, not just at configure time: builds ipc-toolkit's own dependency tree (libigl, TBB, spdlog, TinyAD, tight_inclusion, scalable_ccd, xsimd, robin-map, abseil, filib, all via its CPM recipes), links, and a `CollisionMesh` smoke test (`tests/test_ipc_toolkit.cpp`) runs and passes. Default build is unaffected — option OFF by default, `SPATIUM_HAS_IPC_TOOLKIT` compiles to 0, all 811 test cases / 80 800 assertions still pass with the option off.

### Implicit-contact Newton solver, mechanism proven (2026-08-28)

`examples/cloth_sphere_probe.cpp` — a hand-rolled implicit-Euler Newton solve (`E(x) = 1/2(x-x̂)ᵀM(x-x̂) + h²Ψ_elastic(x) + h²B(x)`) combining a quadratic-spring cloth energy with ipc-toolkit's `BarrierPotential`, replacing explicit XPBD's Gauss-Seidel projection for a cloth-on-sphere scene (the sphere triangulated and merged into the same `CollisionMesh` as a static obstacle — ipc-toolkit's own collision API has no entry point for Spatium's analytical `ContactSurface<Sphere<2>>` directly). Three real bugs found and fixed via diagnostic instrumentation, not guessed: Newton starting from the unchecked predictor `x̂` instead of the last known-safe state (could tunnel before any CCD check ran); an Armijo line search that silently accepted a non-decreasing-energy step whenever every backtracking halving failed instead of erroring out; a CCD min-distance floor set to a physically-meaningful gap (`1e-4`, comparable to `dhat`) instead of a numerical safety margin, permanently deadlocking Newton once contact settled near that value. After all three fixes: confirmed converging via the diagnostic trace (`gnorm` 3.2e-2 → 2.6e-6 over ~14 iterations, no explosion, no tunnelling) — the mechanism genuinely works, beating XPBD's 0/20 baseline on the one config actually traced. Not done: a full multi-config sweep matching the XPBD investigation's own 20-point discipline (see Backlog → Contact physics / RSC) — even one reduced-scale config (81 cloth vertices, coarse sphere, 90 frames) took multiple CPU-minutes under ipc-toolkit's TBB-based collision detection (thread-pool overhead dominates at this problem size). Deliberately not chased further at the time — mechanism correctness was the goal, not throughput.

### Materials & structures

- **[want]** Porous/microstructure analysis — geodesics on pore surfaces, effective conductivity/strength, and grain-shape/orientation analysis in metals, all as manifold queries on a reconstructed surface rather than a bespoke materials-science pipeline.

### Contact physics / RSC

- **[course]** Full multi-config sweep for the implicit-contact Newton solver, matching the XPBD investigation's own 20-point discipline, plus a performance pass (current single traced config takes multiple CPU-minutes under ipc-toolkit's TBB-based collision detection).
- **[course]** RSC's calibration-search against the implicit-contact pipeline — not yet built.
- **[dare]** RSC "real-time control of complex dynamics" domain (illustrative: an underwater drone) — depends on the base+custom deployment split actually working.
