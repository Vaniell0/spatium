# Concept-Driven Physics: A C++23 Reference Specification

**Status:** draft, v0.1
**Reference implementation:** [Spatium](https://github.com/Vaniell0) (this repository)
**Target audience:** C++ library authors, physics-engine designers

---

## Abstract

This document specifies a concept hierarchy for computational physics on arbitrary smooth manifolds and a dispatch pattern that turns specialisation into a compile-time decision rather than a virtual call. The hierarchy is rooted at `Set → TopologicalSpace → MetricSpace → Manifold → RiemannianManifold → Surface`, extended with physics-specific refinements `SymplecticManifold`, `DiscreteLagrangian`, and `ContactSurface`. We claim that under C++23 concepts, a single generic entry point can dispatch between an analytical closed-form routine (3–25 ns/op on the reference hardware), a parametric Newton solver (µs/op), and a user-supplied specialisation, *all without a single virtual function call, runtime type tag, or branch in the hot path*. The reference implementation — Spatium, ~23 kLOC C++23 — realises this for classical mechanics on SO(3)/SE(3), for discrete exterior calculus on triangulated manifolds, for narrow-phase contact on analytical surfaces, and for variational and Lie-group integrators.

This is a reference specification of *what concept-driven physics looks like in practice* — the boundaries where it succeeds, and the limitations it still has to answer.

## 1. Motivation

### 1.1 Why concepts

A modern physics engine's narrow-phase inner loop fires millions of times per simulated second. Every dispatch in that loop — between a sphere-sphere primitive test, a GJK convex-convex test, and a mesh-triangle test — is a classification problem. The traditional OO answer is a virtual function through a `btCollisionShape*` or `PxShape*` base pointer; the cost is (a) an indirection on every call, (b) an instruction-cache miss on first touch of the vtable entry, and (c) the inability of the optimiser to inline across the boundary. The ECS/DOD answer replaces the virtual call with a runtime type tag and a switch, paying a predictable branch-misprediction penalty but still defeating inlining.

C++20 concepts give a third option: the dispatch *is* the type system. A single generic free function

```cpp
template<ContactSurface<T> S>
ContactQuery<T> point_to(const Vec<T, 3>& p, const S& surf);
```

resolves at compile time to one of its overloads — closed-form sphere, closed-form torus, Newton-on-UV parametric, or a user-supplied one-liner. Each call site gets the implementation its type deserves. No vtable, no branch on a type tag, and the optimiser inlines everything. Bullet/PhysX would route all four shapes through the same GJK path at ~100–500 ns/op to preserve the uniform API; we pay the boilerplate of several overloads and in return get an order of magnitude on the analytical shapes.

### 1.2 Why non-Euclidean

Every mainstream real-time physics engine hard-codes its state space to ℝ³. Bullet, PhysX, Havok, Chaos (Unreal 5), Jolt (Godot 4, Horizon Zero Dawn), Box2D, ODE, Newton Dynamics — all of them. When robotics papers ask for dynamics on SE(3), on the sphere, on the hyperbolic plane, on a parametric surface, or on a non-orientable immersion of a Klein bottle, they reach for Drake (C++17, Toyota Research) or Dojo.jl (Julia) or write Python glue by hand.

The pedagogy is older than the engines. Arnold, Marsden, Ratiu, Holm, Bloch — the canonical references for classical mechanics — all derive Hamilton's equations on a symplectic manifold $T^\ast Q$, not on the state vector $\mathbb R^{6}$. Getting from their notation to Bullet's `btRigidBody::applyTorque` is a one-way translation: you lose the manifold structure and never recover it at runtime.

The niche we want is the intersection: *fast, concept-driven, manifold-native physics in C++*, for scientific computing and serious robotics research, not for games. We are not trying to beat Jolt at ragdolls; we are trying to give Drake something to cite as a reference C++ skeleton.

### 1.3 Why now

C++23 closes the last ergonomic gaps. `requires` expressions in lambdas, explicit `this` parameter, the deducing-this shorthand, `std::ranges` fold utilities, module files (where supported) — these are what make the hierarchy below feel natural rather than template-metaprogramming. GCC 15, Clang 18, and MSVC 19.38 all implement enough of C++23 to compile the reference implementation today.

## 2. Concept hierarchy

The spec is defined top-down, from the most general to the most refined. Every concept names the header file that defines it in the reference implementation (`include/spatium/...`). Syntax is expressed in `requires` clauses; axioms that cannot be expressed at compile time (metric symmetry, Jacobi identity, discrete Noether) are enforced through test harnesses.

### 2.1 Space concepts (`core/concepts.hpp`)

```
Scalar              — std::regular + std::totally_ordered + arithmetic + 0, 1
Set                 — PointType, ScalarType, dimension
TopologicalSpace    — Set + contains(p)
MetricSpace         — TopologicalSpace + distance(p, q)
VectorSpace         — Set + VectorType + linear ops
NormedSpace         — VectorSpace + MetricSpace + norm(v)
InnerProductSpace   — NormedSpace + inner(u, v)
Complete            — opt-in tag (S::is_complete)
BanachSpace         — NormedSpace + Complete
HilbertSpace        — InnerProductSpace + Complete
EuclideanSpace      — HilbertSpace + static dimension
Manifold            — TopologicalSpace + TangentVector + exp_map + log_map
RiemannianManifold  — Manifold + MetricSpace + metric_at
Surface             — Manifold + project + normal
```

The hierarchy matches the mathematical one: `EuclideanSpace<N>` is the canonical example of every concept up through `HilbertSpace`, `Sphere<N>` specialises `RiemannianManifold` and `Surface`, `Hyperbolic<N>` specialises `RiemannianManifold` only, `ParametricSurface<T>` specialises `Surface`. No concept demands `ScalarType == double`; `Real50` and `Real100` from `boost::multiprecision` satisfy `Scalar` and participate in all higher concepts, enabling multiprecision orbit propagators without a single `if constexpr`.

### 2.2 Algebraic concepts (`algebra/concepts.hpp`)

```
Magma, Semigroup, Monoid, Group, AbelianGroup
Ring, Field
LieGroup            — Group with exp, log, adjoint
LieAlgebra          — vector space with bracket, Jacobi identity
```

These refine the algebraic structure side rather than the geometric one. `SO3<T>` is a `LieGroup`; `so3<T>` is a `LieAlgebra`; the exp/log maps bridge the two. `Quaternion<T>` is a `Ring` (and a non-abelian group under multiplication), directly usable anywhere a ring is required.

### 2.3 Shape and contact concepts (`geometry/concepts.hpp`, `physics/mechanics/narrow_phase.hpp`)

```
Shape               — dimension + ambient type
ClosedShape         — Shape + closed()
Measurable          — Shape + measure() (length/area/volume by dim)
Bounded             — Shape + bounding_box
DistanceQueryable   — Shape + distance(p)
BoundedRegion       — Shape + contains(p)
ContactSurface<T>   — free point_to(p, s) returns ContactQuery<T>
```

`ContactSurface` is the newest and is the primary vehicle of the dispatch experiment. A type satisfies it when an overload of `point_to(const Vec<T,3>&, const S&)` exists and returns `ContactQuery<T> = {distance, closest_point, normal, inside}`. No inheritance. No type tag. Users add their own shapes by supplying their own `point_to` overload; the concept catches up automatically by ADL.

### 2.4 Physics concepts (`physics/mechanics/*.hpp`)

```
SymplecticManifold    — Configuration, Momentum, State, dimension
DiscreteLagrangian    — callable (q_k, q_{k+1}, h) → scalar with derivatives
StrainEnergy          — W(F) scalar with PK1 and tangent
```

`SymplecticManifold<S>` is satisfied by the canonical cotangent bundle `CotangentBundle<M>` for any `Manifold` $M$; the symplectic form $\omega = dq \wedge dp$ is synthesised there.

`DiscreteLagrangian<LD, N, T>` captures the Marsden–West discrete Lagrangian: a callable taking two consecutive configurations and the step size, returning the discrete action increment, plus the partial derivatives needed to build the discrete Euler–Lagrange (DEL) equation. `SeparableMidpointLagrangian` is the canonical example; the resulting `variational_step_separable` is Störmer–Verlet in closed form.

`StrainEnergy<W, T, N>` is the hyperelastic interface — the template parameter takes a deformation gradient $F \in \mathbb R^{N \times N}$, returns the energy density $W(F)$, and exposes $\partial W / \partial F$ for the first Piola-Kirchhoff stress. `SaintVenantKirchhoff` ships as the canonical example.

## 3. Dispatch patterns

### 3.1 Static overload resolution via ADL

```cpp
// Three concrete overloads, each the optimal routine for its type:
template<Scalar T>
ContactQuery<T> point_to(const Vec<T,3>& p, const Sphere<2, T>& s);      // ~3.7 ns/op

template<Scalar T>
ContactQuery<T> point_to(const Vec<T,3>& p, const geometry::Torus<T>& t); // ~17.7 ns/op

template<Scalar T>
ContactQuery<T> point_to(const Vec<T,3>& p, const ParametricSurface<T>& s); // ~5.2 µs/op

// Generic fallback for anything that satisfies `Surface` but has no
// specialised overload — goes through the ambient-space Newton projection:
template<Surface S>
ContactQuery<typename S::ScalarType> point_to(
    const Vec<typename S::ScalarType, 3>& p, const S& surf);

// The concept itself:
template<typename S, typename T = double>
concept ContactSurface = requires(const S& s, const Vec<T,3>& p) {
    { point_to(p, s) } -> std::same_as<ContactQuery<T>>;
};
```

Client code uses the bare name:

```cpp
template<ContactSurface<double> S>
Vec<double,3> force_on_particle(const Vec<double,3>& p, const S& surf, double d_hat) {
    return ipc_contact_force(point_to(p, surf), d_hat);
}
```

Compile-time resolution selects the closest overload per concrete `S`. No runtime dispatch. User-defined shapes participate on equal terms — the extension mechanism is "define your own `point_to(p, my_shape)` in your own namespace, and ADL finds it".

### 3.2 Algorithm by concept refinement

`make_dec_heat_solver(topo, dt)` takes a `MeshTopology<S>` for any `Surface S` — it builds the Hodge operators from the surface's `metric_at` if available and falls back to the Euclidean Hodge star otherwise. The fall-back path is not a branch; it is a concept refinement: `spatium::mesh::build_laplacian<S>` overloads on `RiemannianManifold` versus plain `Surface`, and the compiler dispatches.

### 3.3 Verify suite

Concepts cannot encode the axioms of a metric space (symmetry, triangle inequality), of a Lie algebra (Jacobi identity), or of a symplectic manifold (closedness, non-degeneracy of $\omega$). The specification pairs every concept with a runtime verifier:

```cpp
template<MetricSpace S> bool verify_metric(const S& space, std::size_t n_samples);
template<LieGroup G>    bool verify_lie_group_axioms(const G& group, std::size_t n_samples);
template<SymplecticManifold S> double verify_symplecticity_drift(
    const typename S::State& s, auto step, double eps, double dt);
```

These are test helpers, not part of the concept. They compile only when the concept is satisfied.

## 4. Reference implementation status

Numbers and status as of commit `7cd472a` on `main`. Tests run with `ctest`; benchmarks with Google Benchmark at `-O3` on an i5-1235U.

| Area | Status | Evidence |
|---|---|---|
| **Space concepts** (`Set` → `Surface`) | shipping | `core/concepts.hpp`; 12 concrete types in `spaces/` satisfy them |
| **Algebraic concepts** | shipping | `algebra/concepts.hpp`; `SO3`, `SE3`, `Quaternion`, `Complex`, `Real50`/`Real100` |
| **Shape + contact** | shipping | `geometry/concepts.hpp`, `narrow_phase.hpp` |
| **Symplectic geometry** | shipping | `SymplecticManifold<S>`, `CotangentBundle<M>`, Yoshida-4 integrator |
| **Lie-group integrators** | shipping | RKMK commutator-free 4, Lie–Euler, Lie–midpoint |
| **DEC** | shipping | `Form0/1/2`, `d_0/d_1`, `*_0/*_1/*_2`, `laplace_beltrami_dec`, backward-Euler heat |
| **Variational integrator** | shipping for separable flat case | `variational_step_separable` = kick-drift-kick Verlet from DEL |
| **LGVI on SO(3)** | Cayley 1-cut; ~3 % energy drift | `lgvi.hpp` — full LLM Newton is documented follow-up |
| **Continuum** | scaffolding | `DeformationMap`, `SaintVenantKirchhoff`, `StrainEnergy` |
| **IPC kernel** | shipping | `contact.hpp` — barrier, gradient, Hessian |
| **XPBD** | shipping for distance + distance-bending | `xpbd.hpp` |
| **Narrow-phase** | shipping for Sphere, Torus, ParametricSurface | `narrow_phase.hpp` |
| **Cloth visuals** | removed (explicit contact ceiling) | headless Newton solve proven correct (`cloth_sphere_probe.cpp`); no viewer-facing demo, full sweep, or performance pass yet |

### 4.1 Narrow-phase benchmark

| Benchmark | ns/op | Comment |
|---|---|---|
| `BM_PointToSphere_OutsideBand` | 3.70 | closed form |
| `BM_PointToSphere_InBand` | 3.70 | same branch — no band-entry penalty |
| `BM_PointToTorus_OutsideBand` | 17.7 | closed form, one frame change |
| `BM_PointToParametric_TorusSurface` | 5245 | grid seed + Newton UV |
| `BM_IpcContactPipeline_Sphere` | 12.6 | point_to + energy + force |
| `BM_IpcContactPipeline_Torus` | 25.5 | point_to + energy + force |

Typical GJK on convex shapes: 100–500 ns/op. The analytical path beats that by an order of magnitude where the shape is known; the parametric path pays for full generality but stays under 10 µs, which is tolerable for offline scientific simulation and cloth with ~1 k particles.

### 4.2 Test suite coverage

**673 test cases, 6 549 Catch2 assertions, all passing**, running under both the legacy header-only build and the C++23 modules build (GCC 15). Tests that bind the concepts to the implementation:

- `test_concepts.cpp` — static_asserts that the expected types satisfy the expected concepts
- `test_narrow_phase.cpp` — `ContactSurface` static_asserts, FD-verified `ipc_contact_force = −∂E/∂x`
- `test_mechanics.cpp`, `test_lie_integrator.cpp`, `test_variational.cpp`, `test_block_b_finish.cpp` — integrator correctness
- `test_block_c_close.cpp` — LGVI, DEC heat, continuum
- `test_block_d_start.cpp` — IPC barrier, XPBD

## 5. Comparison with existing engines

|  | Spatium | Bullet/PhysX/Havok/Chaos/Jolt | Drake | MuJoCo MJX |
|---|---|---|---|---|
| State space | any `Manifold` | hard-coded ℝ³ | articulated MBS on Lie groups | ℝ³ with joints |
| Narrow-phase on Sphere | 3.7 ns/op | ~100-300 ns/op (GJK) | ~100 ns/op | N/A |
| Contact model | barrier (IPC kernel) + XPBD projection | depenetration impulse | hydroelastic contact | smooth contact |
| No-penetration guarantee | IPC kernel provable; pipeline not yet | best-effort | hydroelastic-proven | smooth, no guarantee |
| Multiprecision | `Real50`, `Real100` via Boost | float32 | double | float32 |
| Symplectic integrator | Verlet, Yoshida-4, LGVI (1st-cut), RKMK CF4 | optional Verlet | Runge–Kutta | semi-implicit |
| Non-orientable surfaces | Klein via Bonan–Jennings in parametric | impossible | impossible | impossible |
| Dispatch | C++23 concepts, ADL overloads | virtual functions | virtual functions | function pointers |
| LOC (core) | ~23 k | 500 k+ | 1 M+ | 200 k |

We are not in their league on throughput or GPU. We are in a different league on generality (`Body<Manifold>`), correctness guarantees (IPC-style), and reference-implementation clarity (concept-driven, header-only path available).

## 6. Limitations and open problems

### 6.1 Non-orientable surfaces

The generic `point_to(p, surf)` fallback on a Klein-bottle immersion (Lawson or Bonan–Jennings) returns a closest UV through `Surface::project` + `Surface::normal`, but `q.inside` and `q.normal` are only *locally* defined: as a particle crosses a self-intersection branch, both flip sign. The demo pipeline cannot prevent this from launching cloth particles sideways. A principled fix needs winding-number or volumetric-SDF contact models — see Jacobson–Kavan–Sorkine-Hornung 2013 *Robust Inside-Outside Segmentation*. Marked as future work.

### 6.2 LGVI energy drift

The Cayley-form step in `lgvi.hpp` achieves exact orientation ($R^TR = I$), exact discrete-Noether momentum ($|L|$ conserved to round-off), and exact Casimir ($|\Pi|^2$ conserved) — but the energy drifts at ~3 % over 5 000 steps at $h = 10^{-3}$, monotonically, which is inconsistent with a symplectic-flow modified Hamiltonian bound. The residual lives in the joint $(F_k, \Pi_{k+1})$ update; a full Newton iteration on the coupled system is the documented next refinement (Lee–Leok–McClamroch 2007, Alg. 1, fixed-point resolution only on $y$ is insufficient).

### 6.3 Cloth-on-obstacle with explicit contact

The removed explicit-contact cloth demos (sphere, torus, Klein) showed the limit of explicit position-based contact + Gauss–Seidel projection on a stiff cloth: structural tension from heavy overhanging corners can yank the cloth through the obstacle in a single substep faster than the IPC activation band can respond, even with contact interleaved into the constraint loop. A working pipeline needs implicit IPC (Newton + log-barrier + line-search filter per Li–Kaufman 2020 §5). `ipc-toolkit` is now wired in as an optional CMake dependency (`-DSPATIUM_IPC_TOOLKIT=ON`), and a headless implicit-Euler Newton solve combining a quadratic-spring cloth energy with ipc-toolkit's `BarrierPotential` (`examples/cloth_sphere_probe.cpp`) is confirmed converging on a cloth-on-sphere scene — no explosion, no tunnelling. A full multi-config sweep and a performance pass (this small a problem is slow under ipc-toolkit's TBB-based collision detection) are deliberately deferred follow-up work, not attempted yet.

### 6.4 Remaining gaps toward a v2.0 reference implementation

- A freestanding subset (no dynamic allocation in narrow-phase).
- ABI stability across the C++23/C++26 boundary; currently no stable ABI on any concept.
- A verified (proof-carrying) version of `verify_*` — at the moment these are test helpers, not theorems.
- **Comparative paper** against Drake (and MuJoCo MJX where comparable) on accuracy under equivalent scenes.
- **Contribute to an existing engine as an optional backend** — Drake or geometry-central — to exercise the ABI and find real adoption friction.

## 7. References

### Concepts, dispatch, and standards

- ISO/IEC 14882:2023 (C++23), clause 13.7 (concepts), 13.4 (constrained templates)
- Stroustrup, *Thriving in a Crowded and Changing World: C++ 2006–2020*, HOPL IV (2020)

### Physics on manifolds

- Marsden, J.E.; Ratiu, T.S. *Introduction to Mechanics and Symmetry*. Springer, 2nd ed. (1999).
- Holm, D.D. *Geometric Mechanics*, vols. I–II. Imperial College Press, 2nd ed. (2011).
- Hairer, E.; Lubich, C.; Wanner, G. *Geometric Numerical Integration*. Springer, 2nd ed. (2006).

### Specific techniques

- Marsden, J.E.; West, M. "Discrete mechanics and variational integrators", *Acta Numerica* 10 (2001), 357–514.
- Lee, T.; Leok, M.; McClamroch, N.H. "Lie group variational integrators for the full body problem in orbital mechanics", *Celestial Mech. Dyn. Astron.* 98 (2007), 121–144.
- Celledoni, E.; Marthinsen, A.; Owren, B. "An introduction to Lie group integrators", arXiv:1207.0069 (2014).
- Müller, A. "Evaluation and implementation of Lie group integration methods for rigid multibody systems", *Multibody Syst. Dyn.* (2024), DOI 10.1007/s11044-024-09970-8.
- Hirani, A.N. *Discrete Exterior Calculus*, Caltech PhD (2003).
- Crane, K. *Discrete Differential Geometry: An Applied Introduction*, CMU lecture notes (2019–).
- Li, M.; Ferguson, Z.; Schneider, T.; Langlois, T.R.; Zorin, D.; Panozzo, D.; Jiang, C.; Kaufman, D.M. "Incremental Potential Contact: Intersection- and Inversion-free, Large-deformation Dynamics", *ACM ToG* 39 (2020), SIGGRAPH 2020.
- Macklin, M.; Müller, M.; Chentanez, N. "XPBD: Position-Based Simulation of Compliant Constrained Dynamics", *ACM MIG* (2016).

### Engines compared

- Drake — <https://drake.mit.edu/>
- MuJoCo / MJX — <https://mujoco.org/>
- Bullet Physics — <https://bulletphysics.org/>
- Jolt Physics — <https://github.com/jrouwe/JoltPhysics>
- geometry-central — <https://geometry-central.net/>
- ipc-toolkit — <https://github.com/ipc-sim/ipc-toolkit>

---

*This document is checked into the reference implementation and is updated as the status matrix in §4 changes. C++23 modules (formerly tracked on `exp/modules`) are folded into `main`'s `CMakeLists.txt` behind the hybrid `SPATIUM_USE_MODULES` flag; individual concepts' source-of-truth is the header file named beside them in §2.*
