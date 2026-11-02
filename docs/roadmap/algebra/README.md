# Algebra

## Where it stands

Scalars and the operations on them (`include/spatium/algebra/`): vectors, matrices, quaternions, octonions, complex, dual numbers, `Series` as arithmetic of limits, calculus and Riemannian descent, quadrature as one door, ODE solvers, polynomial solvers with their steps, roots on a bracket; models with their gradients (softmax, a tanh net, score trees) and step rules; statistics (fitting, information, CART).

## Now

- nothing is taken up here now

## Later

- infinity as an input class across the interfaces
- `Series` as arithmetic of limits, joined to quadrature

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area algebra`.

## Reasoning and earlier history

### Riemannian optimization (2026-08-28)

Gradient descent on manifolds (Riemannian SGD): `riemannian_minimize()` in `algebra/calculus.hpp`, Armijo backtracking line search along the manifold. Retracts via each space's own `exp_map()` — no separate retraction abstraction needed. `raise_gradient()` + `project_tangent()` convert the ambient `gradient()` covector into the actual Riemannian gradient first — index-raising matters once the ambient metric isn't Euclidean, e.g. Hyperbolic's Minkowski form (a real bug caught by a failing test, not review: on Hyperbolic the optimizer silently never moved before this fix, because the uncorrected "gradient" pointed the wrong way and line search rejected every step).

### SO(3)/SE(3) templated on Scalar (2026-08-31)

`SO3`/`SE3` (`algebra/groups/so3.hpp`, `se3.hpp`) moved from hardcoded `double` to `template<Scalar T = double>`, matching every other Scalar-templated type in the tree — `SO3<double>`/`SE3<double>` behave exactly as before (all call sites updated: `physics/mechanics/lgvi.hpp`, `examples/tumbling_body_demo.cpp`, test helpers), and `SO3<Dual<double>>`/`SE3<Dual<double>>` now satisfy `Group`/`LieGroup` (checked via `static_assert`, plus a real differentiation test — rotating a point about Z with a `Dual`-seeded angle recovers the exact closed-form derivative, no hand-derived Jacobian). This is what Sophus/manif already give C++ robotics/SLAM/pose-graph optimization; Spatium didn't until now.

Two real, non-hypothetical things found while templating, not assumed:
- `Dual<T>` gained `acos()`/`tan()` (`algebra/dual.hpp`) — `SO3::log()` needs the former, `SE3::log()` the latter; neither existed before since nothing had exercised `Dual<T>` through a Lie-group log map yet.
- A real latent bug, unmasked by the change itself: `SO3::log()`'s `sin(angle)` call had no `using std::sin;` in scope (unlike `exp()`, which did) and silently worked anyway because it resolved to the global `::sin` leaked in by `<cmath>` — invisible while `so3.hpp` didn't yet include `dual.hpp`. Once it did (needed for the `Dual<double>` `static_assert`), `spatium::algebra::sin<T>` (from `dual.hpp`) became visible via ordinary unqualified lookup inside the same namespace and shadowed the outer `::sin`, breaking the `T=double` case outright until `using std::sin;` was added explicitly. Caught by the compiler on the very first build, not by review.

C++23-modules partitions (`modules/algebra/groups_so3.cppm`, `groups_se3.cppm`) needed `import :dual;` added — the `Dual<double>` static_asserts reference a type those partitions hadn't previously needed to see. 825→827 tests, one pre-existing unrelated failure (`embedded_base_is_current`, a stale CMake-configure artifact from the 2026-08-31 release's git-history squash — `kSpatiumCommitSha` reads "unknown" in the cached `build/`, needs a `cmake` reconfigure to pick up the real HEAD SHA, untouched by this change).

**Follow-up same day: `SO3::exp()`/`log()` rewritten for real differentiability at the origin.** Building an actual rotation-averaging demo (`examples/primitives_demo.cpp --scene rotavg`, below) exposed that `gradient(f, Vec3{0,0,0})` came back a hard `[0,0,0]` on a genuine rotation-averaging objective — confirmed via a standalone probe against the same objective evaluated 0.01 away (a real, nonzero gradient there). Root cause: `SO3::exp()`'s old `if (angle < eps) return identity();` and `SO3::log()`'s old `if (angle < eps) return AlgebraType{};` were both exactly correct in *value* at the origin/identity but silently *v-independent* there — under `Dual<T>` this zeroed every derivative exactly at the single most common optimization starting point, even though `exp()`/`log()` are analytically smooth there (the singularity is entirely in the intermediate `axis = v/angle` division and `acos'(1) = -1/√0`, not in the functions themselves). Fixed by re-deriving both in terms of `θ² = v·v` (a smooth polynomial in `v`, unlike `angle = √θ²`) with Taylor-series coefficients near the origin instead of a constant-value shortcut — `SO3::exp()`'s `R = I + a(θ)K + b(θ)K²` now uses `K = skew(v)` directly (not `skew(v/angle)`) with `a,b` the (now removable-singularity-safe) `sin(θ)/θ`, `(1-cos θ)/θ²`; `SO3::log()`'s general formula was refactored around `vee(R-Rᵀ)` (always linear/smooth in `R`) scaled by a coefficient that's Taylor-expanded near `cos_angle=1` *before* ever calling `acos` there, since `acos`'s own derivative diverges at exactly that point. The `angle≈π` branch is untouched — a genuine, non-removable coordinate singularity of the axis-angle chart itself (two antipodal axes represent the same rotation there), not a Dual artifact, and out of scope for this fix. Two new regression tests (`tests/test_algebra.cpp`): `SO3<Dual<double>>` differentiates correctly *at* `v=0` (checked against the known first-order Taylor expansion `exp(v) ≈ I + skew(v)`) and through a `log(exp(v))` roundtrip at `v=0`. 827→829 tests, same one pre-existing unrelated failure as above.

Re-deriving `SO3::exp()` this carefully surfaced the same *shape* of problem in `SE3::exp()`'s translation part — and cross-checking it against an independent ground truth (a 40-term brute-force Taylor sum of the 4×4 se(3) generator's own matrix exponential, computed independently of `so3.hpp`/`se3.hpp`) turned up a **real, pre-existing, silent correctness bug**, not just a differentiability gap: the old `V` matrix used `skew(ω/angle)` (the unit axis) with coefficients `(1-cos θ)/θ²`, `(θ-sin θ)/θ³` that are only valid for `skew(ω)` (unnormalized) — wrong for any nonzero rotation combined with a translation. Confirmed numerically: the old formula gave `t=[0.897, 0.646, 0.218]` against the ground truth's `[0.975, 0.603, -0.062]` for a representative `(ω,v)` — not close, and the *z* component even has the wrong sign. No prior test caught it: `"SE3 exp/log roundtrip"` is self-consistent under either convention (both sides used the same wrong one), and `"SE3 exp known-good: 90deg Z rotation + translation"` builds `T` via `from_Rt()` directly, bypassing `exp()`'s `V` matrix entirely — a real, honest test-coverage gap, not carelessness caught late. Fixed by extracting a shared `translation_jacobian(ω)` (the corrected, ground-truth-verified `V(ω) = I + b(θ)K + c(θ)K²` with `K = skew(ω)`, Taylor-safe near `ω=0` the same way `SO3::exp()` is) used by both `SE3::exp()` (as `V` itself) and `SE3::log()` — which now inverts it via the library's own general `invert()` (`algebra/linear_solve.hpp`) rather than a second hand-derived `V⁻¹` closed form, deliberately trading a small runtime cost for not risking a *third* manually-derived formula after the first one turned out wrong. One more regression test (`tests/test_algebra.cpp`): `SE3::exp()` checked directly against the brute-force ground truth, closing the exact coverage gap that let the bug through. 829→830 tests, same one pre-existing unrelated failure as above.

**`examples/primitives_demo.cpp --scene rotavg`** — the demo that surfaced all of this: SO(3) rotation averaging (N noisy measurements of one rotation → the estimate that best explains them all, a real SLAM/sensor-fusion building block), solved via `calculus.hpp`'s existing `gradient()`/`minimize()` with zero hand-derived Jacobian, animated one Armijo-backtracking step per frame (the same per-iteration logic `minimize()` uses internally, exposed step-by-step) so the convergence itself is watchable — bright RGB gizmo (live estimate) converging toward a gray reference (ground truth) alongside dim RGB gizmos (the noisy inputs).

### Calculus

- **[course]** Region-aware integration over Spatium's own geometry (Box/Sphere/Polygon volumes, intersection overlap) — `integrate()` in `algebra/calculus.hpp` is deliberately scoped to 1-D; this is the next step the function's own docstring already flags.

### Polynomial solvers: the degenerate leading coefficient

Design settled 2026-09-15. **Shipped 2026-09-17**, and the three
`[!shouldfail]` pins did exactly the job they were written for: the day
the solvers stopped dividing by a vanishing leading coefficient all three
started passing, Catch2 reported *that* as a failure, and the tags had to
come off. A bug recorded as a comment gets skimmed; one recorded as a
test that must keep failing cannot be quietly forgotten, and cannot be
quietly left marked broken after it is fixed either.

The return type is `UpTo<Complex<T>, N>` in a new `core/up_to.hpp` — at
most N values with a count, no allocation — and the name states the thing
that changed. `solve_linear` exists as a named function because a
degenerate quadratic is exactly that, and each solver now collapses to
the one below it. The donut frame is byte-identical afterwards, which is
the expected result rather than a disappointing one: a ray exactly
parallel to a quadric's null direction is measure-zero among camera rays,
and the defect was always about the cases a renderer does not sample.

**Two decisions the design did not contain, recorded because both had a
plausible alternative.**

*The zero test is exact equality, not a tolerance.* A relative test — "a
is negligible beside b and c" — catches more cases and is the wrong trade:
it turns a genuine but ill-conditioned quadratic into a linear one and
drops a real root that is large but honest, replacing one silent wrong
answer with another. Exact zero is the only threshold that says *this is
not a quadratic* rather than *this quadratic is awkward*, and it is also
exactly what the geometry produces — for a ray along a cone's generator
the t² coefficient is one quantity minus itself, zero to the bit. What
this leaves open, stated so it is not mistaken for covered: catastrophic
cancellation when the leading coefficient is tiny but nonzero. The fix
for that is the numerically stable form (`q = -(b + sign(b)·√disc)/2`,
then `x₁ = q/a`, `x₂ = c/q`), it is a separate defect, and bundling it
would have made one change two.

*`ray_quadric_proximity` checks the answer, not the coefficient.* The
`[!shouldfail]` test's own comment predicted the fix would be a
degeneracy check on the t² coefficient *before* `solve_quadratic`, on the
grounds that checking after would catch the symptom. That prediction was
written when the return type could not say "fewer roots", and it is
superseded rather than merely unmet: asking whether the answer has two
roots is what the function's model actually requires, and it leaves one
place deciding what degenerate means. A coefficient test at the call site
would be a second copy of the solver's rule, free to drift from it.

**The defect.** `solve_quadratic`, `solve_cubic` and `solve_quartic` all
divide by the leading coefficient with no guard. `solve_quadratic(0, 2, -1)`
is a perfectly good equation with the root `0.5` and returns `NaN` and
`-inf`. Because `NaN` reports `is_real() == true`, the garbage then walks
through every caller's filter. Concretely: a ray parallel to one of a
cone's own generators meets the surface exactly once, and the hit is
silently lost.

**Where the risk actually lives, stated as a boundary rather than a file
list.** The leading coefficient is dangerous exactly where it is *data*
and safe exactly where it is *structure*.

- **Data.** A ray's direction is user input and can degenerate — down a
  cylinder's axis, along a cone's generator. `geometry/ray_surface.hpp`,
  `geometry/ray_hit.hpp`, `io/scene.hpp`, `physics/mechanics/rigid_contact.hpp`.
  A third live instance sits at `ray_surface.hpp:394` in
  `ray_quadric_proximity`, and it presents worse than the other two. A ray
  down a cylinder's axis returns `Result` **success** with
  `closest_t = NaN`, a `NaN` point, and `miss = 0` — because the NaN root's
  imaginary part is exactly zero, so `abs()` of it is a perfectly plausible
  number in the one field a caller is most likely to branch on. The other
  two hand back NaN, which is at least conspicuous. This one reports a
  clean grazing hit for a ray that runs down the middle of the tube and
  never approaches the wall.
- **Structure.** A characteristic polynomial is `det(A - λI)` by
  definition; its leading coefficient is 1 as a matter of mathematics and
  no input can change that. `spaces/spd.hpp`, `algebra/eigen_decomp.hpp`,
  `core/precision.hpp` are safe not because they were checked but because
  it could not be otherwise. Those three pass `T{1}` literally.

These are different grades of certainty and the second is the stronger
one. An `assert(leading != T{0})` at the structural sites is still worth
its line, but for a different reason than at the data sites: there it
guards a future refactor changing what gets passed, not an input.

**The return type.** A fixed-capacity container with a count, no
allocation — it must be able to say "one root" where the degree collapsed.

- `size()` is the count; `capacity()` is 2/3/4. Range-`for` takes its
  bounds from the count, so every existing `for (auto& r : roots)` call
  site keeps working unchanged.
- `operator[]` is **unchecked in Release** — no branch in a hot path — so
  indexing at or past `size()` is undefined there. This has to be said
  loudly in the docstring, because a reader arriving from `std::array`
  expects a fixed size where every index is always valid, and that is
  precisely what stops being true. In Debug it carries
  `assert(i < size())`.

  That assert was nearly dropped on the theory that it would fire on
  correct code, `spaces/spd.hpp` indexing `roots[2]` being the worry.
  Checked instead of assumed — `git grep "roots\["` finds exactly four
  library sites. The two in `spd.hpp` are monic, so the count is always 3
  and the assert can never fire there. The other two are both inside
  `ray_quadric_proximity`, where the count really can fall short — and
  they are the ones already producing `miss = 0` today. So the assert does
  not break correct code; the only code it breaks is the code that is
  already wrong.
- `at(i)` returns `Result<Complex<T>>` and **never throws**. Checked
  access signals; unchecked access does not; different contracts, so
  different return types. Exceptions appear nowhere else in this library
  and must not appear here.
- In Debug, slots between `size()` and `capacity()` are poisoned with
  `NaN`. Zero instructions under `NDEBUG`, and under the Debug CI job it
  turns reading past the count from undefined into visibly wrong.

**Tests have to cover the case that changed, not the one that already
worked.** `roots[2]` used to be valid always, by virtue of a fixed-size
array; now it is valid only when the count reaches 3. So both need
covering: indexing within the count, and `at()` refusing beyond it.
Testing only in-range indexing tests the old behaviour.

**The one non-library caller, and the convention question it raises.**
`rsc/include/precision_ops.hpp`'s `flatten_cubic_roots` takes
`const std::array<Complex<double>, 3>&` explicitly and copies three roots
unconditionally, so it needs a signature change and a check.

The check is an `assert`, not a `Result`, and the reasoning is worth
keeping because the first answer here was the wrong one. The registry op
receives `std::span<const double> in` with `in[0]` as the leading
coefficient, which *looks* like a boundary with untrusted input. It is
not. A boundary is where data arrives from outside the process and may be
anything — `load_obj`, `json::parse`. Here the data comes from a harness
living in the same binary, generating inputs from a distribution it
defines itself. That is an internal pipeline invariant, and
`conventions.md` already answers those: `assert`, the same as
`physics/mechanics/`, `spaces/` and `mesh/`.

So: `assert(in[0] != T{0} && "harness contract: monic cubic")` in the op,
`assert(roots.size() == 3)` before the copy, and `flatten_cubic_roots`
keeps returning `void`. A monic cubic always has three roots with
multiplicity; if it has two, either the input or the solver is broken,
and both are bugs rather than failure modes. `Result` here would be
signalling a refusal that cannot legitimately happen.

**`Op::Fn` returns `void`, deliberately, until something needs otherwise.**
An op with no failure channel is the right shape for every op the registry
has, because none of them can legitimately fail. A channel is warranted
when an op appears that genuinely can — an external solver, a network
call, a parse, anything with unbounded input — and not before. Recorded
explicitly so the next reading of "but the convention says `Result`"
does not have to re-derive that the convention is about boundaries rather
than about everything.

**NaN poisoning is a debugging aid, not a channel.** Slots between
`size()` and `capacity()` get NaN under Debug so that reading past the
count is visible rather than undefined. It is the same category as
`_GLIBCXX_ASSERTIONS`: a signal while debugging, never part of the
contract, and zero instructions under `NDEBUG`.

The training checkpoint is **not** affected. `rsc/include/precision_task.hpp`
builds coefficients from three sampled roots by Vieta and returns
`{1.0, b, c, d}` — the leading coefficient is literally 1 in every
sample, the degenerate case never occurs in the task's distribution, and
`solve_cubic_f64 4 4 6`'s pinned in/out sizes stay exactly right. The
distribution is not a guarantee, which is why the op still gets the
check; it is the reason the checkpoint does not have to be retrained.

### Native math (working defaults, not dependency removal)

Reframed 2026-09-13. The old heading was "dependency reduction" and the item below claimed native SVD/eigendecomposition was still open — both stale. Dense linear algebra landed natively on 2026-09-06 (`0372794` general-N symmetric eigendecomposition via cyclic Jacobi, `4c77d38` general MxN SVD via eigendecomposition of AᵀA, benchmarked in `cb7978a`), works for any `Scalar` including `Real<Digits>`, and is cross-validated against `Eigen::JacobiSVD` in `tests/test_svd.cpp`.

Corrected 2026-09-15: the sentence here used to add "which passes, contrary to the 'one pre-existing unrelated failure' note repeated in several Completed entries above". Both halves were true and the conclusion drawn from them was wrong. The test did pass, in Release; the earlier note was written from a Debug run where it aborted. The cause was ours, not Eigen's — the test asked `JacobiSVD` for thin U/V on a fixed-size matrix, which Eigen allows only when the column count is dynamic and otherwise asserts. `NDEBUG` removes that assert, and every CI job builds Release, so "it passes now" read as "the note went stale" rather than "the check stopped running". Fixed in #29 by asking for full factors, which the test's singular-value comparison does not care about either way. The general lesson is recorded as an open item below: a Release-only CI matrix cannot see this class of defect at all.

The goal was never to remove dependencies. It is to have a working default for every operation so no capability is gated behind an optional package, while any specialized library can still be absorbed behind the same interface — the pattern `SPATIUM_HAS_EIGEN` / `SPATIUM_HAS_BOOST_MULTIPRECISION` / `SPATIUM_HAS_IPC_TOOLKIT` already implements. Spatium does not need to beat specialized libraries at their own work; it needs to compose with them and still stand alone without them.

- **[course]** Sparse linear algebra — the real remaining gap, and the only place with no native default at all. The heat method, DEC, `differential.hpp` and `mesh/geodesic.hpp` all require `SPATIUM_EIGEN=ON` because there is no sparse matrix type or sparse factorization in the tree. Until that exists, geodesics via the heat method are a capability a default build simply does not have.
- **[want]** `Eigen::Ref<…>` public-API adapters — was blocked on dense SVD/QR, which has since landed, so this is now unblocked rather than waiting.

- **[want]** **Spatial acceleration on a manifold — the one part of "any space" that stops at Euclidean.** Noticed 2026-09-16 while measuring instancing, and it is a gap rather than a defect: nothing is currently wrong, because everything that would hit it routes around it. `examples/hyperbolic_tessellation_demo.cpp:8` says so in its own header — *"not a mesh + BVH::ray_cast() (that stack is Euclidean ray-primitive intersection)"* — and sphere-traces with `space.distance()` instead, linear in marker count, which is why that demo has dozens of markers and not thousands. `blackhole_gr_demo` has no BVH at all.

  An AABB cannot simply be carried over, and not only because it would be loose: a box needs coordinates, a chart does not preserve distances, and a sphere has no global chart without singularities — so a node that misses in the chart can be hit on the manifold. The intrinsic bound is the geodesic ball `B(p,r) = { q ∈ M : d(p,q) ≤ r }`: chart-independent, metric-consistent, and exact.

  **Half of this is thirty-five years old and should be borrowed rather than invented.** Ball trees (Omohundro 1989), vantage-point trees (Uhlmann 1991; Yianilos 1993) and M-trees (Ciaccia 1997) index metric spaces by balls for exactly this reason, with established build heuristics. What they answer is nearest-neighbour and range queries.

  **The open piece is the ray.** Casting a *geodesic* against a metric ball through such a tree is where the literature runs out; the closest relatives are sphere tracing in curved space (Hart / Hawksley / Matsumoto / Segerman, hyperbolic VR, ~2017 — which is what our own hyperbolic demo already does) and GR ray tracing, which brute-forces or adaptively steps rather than using a tree. The ball test has closed forms where the geodesics do: spherical and hyperbolic law of cosines. General manifolds need numerical geodesic integration, which `physics/relativity/geodesic.hpp` already has. SAH needs a volume measure, and on a manifold that is curvature-dependent — `2π(1 − cos r)` on the sphere, `2π(cosh r − 1)` in hyperbolic 2-space, an integral of the metric otherwise.

  Not an optimisation: without it there is no acceleration structure off flat space at all. And explicitly **not** an RSC problem — ROADMAP's own description above has RSC choosing *which existing* call to make, with Spatium as the exact substrate. It selects among methods that exist and cannot invent the ray-ball test; it gets work here once there are two ways to do this, not before.

  **Placed in the order, 2026-09-17, so that "not an RSC problem" does not turn into "not a problem yet".** `[want]` was the right grade for how settled the design is and the wrong signal about when to look at it: an item whose first line is *"the one part of any space that stops at Euclidean"* does not belong in the far drawer. The sequence is the HN post, then RSC, then this — and it comes directly after RSC rather than at some later unnamed point, because RSC is precisely the thing that will have two methods to choose between only once this exists. So the dependency runs the opposite way from what the paragraph above might suggest: RSC does not gate the ball tree; the ball tree is what eventually gives RSC something to select here.
