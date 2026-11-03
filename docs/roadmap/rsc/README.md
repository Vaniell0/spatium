# RSC: synthesis and search

## Where it stands

`rsc/core` is generic: a typed graph of the library and a tree search guided by a learned prior and value, which reads on past a value that is not a number. `rsc/spatium` binds it to the library: particular analytic solutions synthesised as IR fields on an interval, certified on a dense grid, piecewise where one expression does not reach the tolerance, and written as a person writes them where that still holds (`sqrt(-p)`, not the interval's own variable). The definition is `rsc/README.md`.

## Now

- Pow among synthesis's operations, so the root of x^p + p = 0 is (-p)^(1/p) and not a Chebyshev polynomial

## Later

- structural candidates for where to split an interval, before bisection and before a learned policy
- infinity and chaos as classes in the search, after infinity is an input class in algebra
- a model over the IR that names what the library lacks, admitted only through samples, certificates and green PRs

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area rsc`.

## Reasoning and earlier history

### RSC: a trained dispatcher on top of Spatium

`rsc/` — a trained dispatcher deciding *which* Spatium call to make and with *what* parameters, not a learned recurrence pretending to be a computer: Spatium stays the exact substrate. Minibatched policy-gradient (REINFORCE) training, hand-derived backward. Domain pipeline, ordered by how proven the underlying mechanism is (named by content, never by position):

- **Precision-critical dispatch** (`precision_ops.hpp`) — `solve_cubic` f64 vs. Real50 on casus-irreducibilis cubics; ~0.99 held-out accuracy.
- **Geodesic/mesh dispatch** (`geodesic_task.hpp`) — Dijkstra vs. Heat method; proven mechanism, partial accuracy (~0.536 on class-balanced sampling — real signal, not yet near-ceiling).
- **Root-finding dispatch** (`rootfind_ops.hpp`) — Newton vs. bisection around `f(x)=x^3-a`'s `f'(x)=0` inflection.
- **Cauchy/IVP dispatch** (`ode_ops.hpp`) — Euler vs. RK4 across Decay/Oscillator/CircularOrbit families.
- **Mesh-strategy dispatch** (`mesh_ops.hpp`) — uniform vs. anisotropy-adapted UV tessellation; plateaus ~0.77-0.79 over a ~0.26 baseline.
- **General linear-solve dispatch** (`linear_ops.hpp`) — Jacobi iteration vs. direct Gaussian elimination for `Ax=b`.
- **Rigid-body integrator dispatch** (`integrator_ops.hpp`) — five real `PointMass` steppers (`euler`/`semi_implicit_euler`/`verlet`/`rk4`/`yoshida4`), already present in `physics/mechanics/integrator.hpp` with no new algorithm work needed, across three closed-form test families (UniformGravity, Spring, PointGravity). Measured dispatch mix (not uniform): Verlet ~63%, RK4 ~24%, Euler ~8%, semi-implicit Euler ~4%, Yoshida4 ~1% — every candidate wins real problems. 4.9% → 84.5% held-out accuracy over 2500 REINFORCE updates.
- **Contact physics / soft bodies** — investigated and found structurally blocked: a 20-configuration compliance/substep sweep of explicit XPBD contact against a heavy cloth overhang failed on all 20 (stretch ≥23%, often exploding). See "Contact physics build dependency" below for the resolution in progress.
- **Real-time control of complex dynamics** (illustrative: an underwater drone) — see Backlog → Contact physics / RSC; sequenced last, depends on the base+custom deployment split actually working.

New Spatium primitives built specifically to support this pipeline: `algebra/dual.hpp` (forward-mode autodiff), `algebra/calculus.hpp` (`gradient`/`integrate`/`minimize`), `algebra/ode.hpp` (generic `euler_step`/`rk4_step`/`integrate_fixed` over `Vec<T,N>`), `algebra/linear_solve.hpp` (`solve_direct`/`solve_jacobi`/`diagonal_dominance_ratio`), `spaces/parametric.hpp`'s `parametrization_anisotropy`/`normal_at`. Full per-domain writeups, measured accuracy numbers, and postmortems (each domain surfaced at least one real bug during training, not just tuning) live in `rsc/README.md`.

### RSC as search, not classification

Recorded 2026-09-15. Untested, but the pieces it needs already exist and were checked before writing this down.

**The reframing.** RSC today classifies: given a problem, pick the operation. The harder and more valuable question is *synthesis* — find the shortest correct **sequence** of operations. That is not open-ended code generation, it is combinatorial search over a closed, finite action space (a few dozen registered operations, every one already implemented and tested), which is the classical setting for program synthesis over a fixed DSL. No LLM belongs in that loop.

**Why the pieces fit.** `rsc/include/dispatcher.hpp` already emits a softmax distribution over the registry — that is a policy network in the AlphaZero sense, already trained. What is missing is a value head (a second small output estimating how close a partial chain is to the goal, so a tree search prunes instead of enumerating) and the search itself, wrapping the existing `Registry`/`Chain`. Reward shape is the one `rsc/README.md` already cites AlphaDev for: correctness is a hard filter via the existing `grade()`, length is minimised only among candidates that pass it.

**The cheap first experiment, using only what exists.** `rsc/include/chain.hpp` is teacher-forced on add/multiply chains, with a dispatcher measured at 0.312 → 0.972 on 3-step chains (`rsc/README.md:241`), and `tests/test_rsc_chain.cpp` already in place. Replace teacher forcing with free-running search that uses the trained dispatcher as its heuristic, and see whether it finds chains as short as, or shorter than, the teacher-forced path. `chain.hpp:21` records that multi-step credit assignment was never attempted precisely because teacher forcing sidesteps it — which is the gap this closes. Three existing pieces connected, nothing built from scratch; if it fails, the failure points at a specific next thing to diagnose, which is how the rest of this project is already run.

**The part nobody else can do, and the reason it is worth trying.** A chain found this way is expressed over *concepts*, not concrete types. An algorithm discovered by searching on `Sphere` transfers to `Hyperbolic` and to an arbitrary `ParametricSurface` by construction, because it was never written against a particular space. AlphaTensor found faster matrix multiplication, but that discovery is bound to a specific ring and does not generalise itself to another algebraic structure. Space-agnostic discoveries are not something anyone has searched for, because nobody else has a single interface through which such a generalisation could even be *checked*. That is the claim worth testing, and it costs one experiment to start.

**Measured 2026-09-18, and it changes where the search has to happen.** A first attempt built the search over a toy numeric action space -- five integer operations, reach a target -- and the result was a clean negative: a trained policy used to order a best-first search was worse than uninformed breadth-first on every axis at once, solving less often (90.8% against 100% at depth 6), expanding 1.6-1.7x *more* nodes, and returning a shortest chain half as often, while costing 600 ns per expansion on top. See `rsc/tools/search_heuristic.cpp`.

The reason is structural rather than a training failure, and it generalises: breadth-first is the provably optimal order for shortest paths under uniform cost, so a heuristic there competes with an optimal ordering rather than improving a poor one. It has to carry real information merely to break even.

But the deciding fact is about the **state type**, not the method. Numbers collide: `10+3+3` and `2*8` reach the same state, so the search space collapses. Measured on that toy, 9 765 625 paths at depth 10 reduce to 13 836 distinct states -- a 706x collapse, and an effective branching factor of 2.13 rather than 5. A space that small is searched exhaustively for nothing, and no value head changes that.

`Registry::chainable_ops()` selects on arity -- `in_size==2 && out_size==1`, so the output feeds back as the next accumulator -- which for the current registry means `add`/`multiply` over a *scalar* accumulator. That was always labelled a mechanism check, and the collapse is a property of scalars, not of chains.

**Geometries do not collide.** Two different sequences of mesh operations essentially never produce the identical mesh, so generalising the accumulator from a scalar to an object of the substrate leaves the branching factor intact and puts the search in the regime this entry is about. The selection rule does not change at all; only the type of the state does. (**Measured 2026-09-18**, `rsc/tools/mesh_collisions.cpp`, and the claim as first stated was too strong. Geometries *do* collide -- translations and scales compose, so two sequences can land on one mesh -- but roughly four times less. Six mesh operations from an icosahedron keep an effective branching factor of 4.66 out of 6 by depth 6, against 2.13 out of 5 for the integer chain, both converged. On the same depth that is 13 590 distinct meshes against 664 distinct numbers. Counting meshes equal when they agree to 1e-9 rather than bit-exactly lowers it to 4.24, also converged -- worth knowing, since that gap widens with depth as accumulated arithmetic makes more near-coincidences.)

**Measured 2026-09-22, both sides.** An exhaustive pass over the mesh space takes 0.20 s at depth 5, 1.55 s at 6, 10.50 s at 7 — roughly seven times per level, above the 4.66 branching factor because the meshes themselves grow as well. Extrapolated: about 73 s at depth 8, eight minutes at 9, an hour at 10. So **exhaustive search stops being affordable around depth 8–9**, and that is the number that decides when a heuristic earns its place: below it breadth-first wins and a policy is overhead, above it there is nothing to be exhaustive with.

At depth 7 that is 81 540 expansions in 10.50 s — **129 µs per node**. A network evaluation is 608 ns, so the policy costs 0.47% of one expansion, and the distilled per-domain tree 0.0065%. Against the numeric toy, where a policy cost more than ten nodes, that is a two-hundredfold inversion, now measured from both ends rather than argued.

It also inverts the economics, which the toy got backwards. Expanding a node over numbers is five integer operations and a map probe -- tens of nanoseconds -- so a 600 ns policy evaluation costs more than ten nodes and can never pay. Expanding a node over geometry means *running* a mesh operation: milliseconds. Against that a learned heuristic is free, and the distilled 8 ns tree (see `rsc/tools/distill_tree.cpp`) is free several times over.

**What the search is actually for, stated so the ambition is calibrated.** It finds *compositions* of registered operations; it cannot invent an operation. "A new algorithm" here means a sequence nobody wrote down -- the same sense in which AlphaDev found shorter sorting sequences out of existing instructions. What makes that more than a curiosity is the concept-generic part: a chain of concept-constrained ops is a program over the concept, so running it on `Sphere` and on `Hyperbolic` is two instantiations of one algorithm rather than two algorithms.

**And the existing design already predicts which half transfers.** This file's own dispatch/calibration split -- dispatch is the trained model, calibration is classical optimization and not the model -- maps exactly onto it: the *dispatch* transfers, since which operations in which order is space-agnostic by construction, while the *calibration* does not, because a step size tuned on a unit sphere will not hold on a hyperbolic space. Structure travels, constants are re-fitted per space.

**Measured 2026-09-18** (`rsc/tools/transfer.cpp`), and it came out as the design predicted rather than as hoped. A chain was searched for on `Sphere<2>` and then replayed unchanged on `Euclidean<3>` and on a `ParametricSurface` torus -- three spaces that behave differently under the same operation, since subdivision projects new vertices onto the sphere, onto nothing in Euclidean space, and onto the torus.

| specification | transfers |
|---|---|
| relative (10x the starting faces, edge at most 0.45x the starting edge) | **3 of 3** |
| absolute (at least 200 faces, edge in [0.40, 0.50]) | 1 of 3 |

The relative chain is `subdivide -> subdivide` and holds everywhere. The absolute one is `subdivide -> subdivide -> scale_up`, and fails on Euclidean by 6% (0.5303 against a 0.50 ceiling) and on the torus by half (0.2840 against a 0.40 floor).

**Structure transfers, constants do not** -- measured, and landing exactly on the dispatch/calibration split this design already drew for unrelated reasons. Two independent routes to the same seam is the sort of agreement worth trusting.

The check is also the compiler's: the operations are templates over `Surface S` and none of them names a space, so the file compiling for three of them *is* the single interface through which such a generalisation can be checked. That phrase stops being a metaphor here.

**And the follow-up question is answered too.** "The transferred chain still satisfies the spec" is weaker than "it is as good as a local search would have found", so each target space was also searched directly:

- **Relative spec: transfer costs nothing.** Searching Euclidean and the torus directly returns a chain of the *same length* -- in fact the same chain. One caveat: the search returns the first solution at minimum depth and `subdivide` is listed first, so the identity of the chains is partly an artefact of op ordering. The equality of *length* is not.
- **Absolute spec: two different failures, and they must not be conflated.** On Euclidean the transferred chain fails and a local search finds nothing at all -- the specification is unreachable there, which is not transfer's fault. On the torus a local search does find one, of the same length but a *different* composition (`scale_down -> subdivide -> scale_up`). There the structure genuinely did not carry.

Related and already in this file: manifold-native RL (Object model section), Fisher-Rao as a `RiemannianManifold` (Manifold applications) — the closest to buildable, since it needs no new abstraction — stochastic processes on manifolds, and hyperbolic embeddings for hierarchical data (Discrete / graph geometry).

### The labels check the outcome, not the preference

Asked 2026-09-22: how do we know a dispatcher's choice is the *right* one,
and was measured speed used when comparing the distilled tree against the
network? Two different holes, and the answer to the second is no.

**The label scheme is sound where it was examined.** Ground truth is
generated rather than hand-written: candidates are run cheapest-first and
checked against a reference-quality last candidate -- `Real50` for
precision, 50-iteration bisection for root-finding, exact Gaussian
elimination for the linear domain. The cheap one wins when it is accurate
enough. That is the right notion, "cheapest of the correct", and the
accuracy half is genuinely measured, with each two-regime split probed
before being written down.

**Which candidate is the cheap one is not measured.** It is fixed by
registry index. Measured (`rsc/tools/solver_cost.cpp`), dense, 200k
repetitions:

| | direct | jacobi(20) | |
|---|---|---|---|
| N=4, diag 6 | 64.9 ns | 430.9 ns | jacobi **6.6× dearer** |
| N=4, diag 20 | 50.2 ns | 489.5 ns | 9.8× |
| N=8 | 401 ns | 1121 ns | 2.8× |
| N=16 | 1599 ns | 3844 ns | 2.4× |

Jacobi is dearer at every size tried, so a dispatcher trained to prefer it
when it converges is trained to lose whenever it says so. In fairness the
domain is not wrong in principle -- Jacobi wins asymptotically and on
sparse systems -- but `solve_jacobi` here is dense and the training sizes
are 4 to 16, so it is wrong at the sizes it trains on.

**And the distillation numbers were fidelity to that label, not to speed.**
The 0.005 accuracy lost by a per-domain tree measures how faithfully the
tree reproduces the network's agreement with the label. Where a label's
cost premise is inverted, both network and tree reproduce it faithfully.
Those numbers are unaffected in fact -- the base model is Tier-1,
precision, root-finding and Cauchy/IVP, and the linear domain is not in it
-- but the method was not checking what the question assumed.

### One RSC domain is trained on a choice nothing makes

Checked 2026-09-22 by exact name, after a crude grep over guessed names
produced a confident wrong answer for the other domains -- `lgvi` is
called from `tumbling_body_demo` and three tests, and bisection lives in
`ode.hpp` and `spherical_polygon.hpp`, neither of which the first pass
found. Only the row below was verified properly, and only it is claimed.

**`solve_jacobi` has no callers anywhere outside RSC's own training, and
`solve_direct` survives only in comments** -- the code that used to call it
now calls `invert()`, after profiling showed four independent solves at
~1 µs/step against ~0.3 µs for one inverse. So `rsc/include/linear_ops.hpp`
trains a dispatcher to choose between two implementations the library does
not choose between anywhere.

That is not an argument against the domain: it is a reasonable thing to
have trained, and the registry is meant to hold operations whether or not
the rest of the tree calls them. It is an argument about what "connect RSC
to the pipeline" can mean. A dispatcher needs a call site with a real
alternative, and this one has none, so connecting it would mean inventing
the consumer first.

The two places that *do* have the shape are recorded above: the
render-level choice, where the rule reads one feature of three and is
measurably wrong past a hundred objects, and `RenderLevel::Newton`, a
third option the automatic rule cannot produce at all.

### Where Newton is not needed, measured

The argument, recorded before it could be tested: while a Newton step is
in the pipeline the IR does not compress and there is no interactivity,
because a Newton step is `solve_linear(H, g)` then a line search then an
update — a long chain whose control flow depends on its own intermediate
values, which lowers to no kernel and serialises into nothing. Newton is
the largest opaque leaf in the physics path.

**The alternative now exists and is a composition of pieces the library
already had separately:** `point_to`, then `ipc_contact_force`, then
`verlet_step`. No Hessian, no linear system, no backtracking. Every step
is arithmetic, so the chain is expressible as operations — which is the
same property the closed-form section above is about.

So the question is not whether Newton can be replaced but where the
boundary lies, and `rsc/tools/explicit_contact_window.cpp` sweeps it.
Dropping a body onto a unit sphere, reporting penetration or energy
growth over six seconds:

| approach | kappa 1 | 100 | 10⁴ | 10⁶ |
|---|---|---|---|---|
| v = 0 | through at every dt | ok at every dt | ok at every dt | 396× at dt 1e-2, ok below |
| v = 5 | through | through | ok below dt 3e-3 | 240× at 1e-2, ok below 1e-3 |
| v = 50 | through | through | ok only at dt 1e-3 | ok below dt 1e-3 |

**Two failures close in from opposite sides** — a barrier too weak to stop
the body, and a barrier too stiff for the step — and between them is a
window that narrows as the approach speed rises.

**Which reframes what a dispatcher is for here.** It does not replace
Newton. It decides *where Newton is not needed*, and inside that region
the description becomes processable: no linear solve, no data-dependent
control flow, an IR that compresses. The features are cheap and all three
are already to hand — approach speed, stiffness, step size.

One cell is not explained and is left visible rather than smoothed: at
v = 50 and kappa 10⁴, dt 1e-3 holds while the *smaller* 3e-4 penetrates.
Non-monotonic in the step size, which nothing above predicts.

### A closed form is not only faster — it is the only thing that fits in an IR

Named 2026-09-22, after the same decision arrived twice in one session
from unrelated directions.

`ParametricSurface::project` searched for parameters where a sphere's
projection is `p·r/|p|`, and was wrong on 12.9% of points because a search
cannot converge at a degenerate chart point. `sweep_point_surface` advanced
conservatively toward a first contact where a quadric's answer is a root of
a quadratic, `ray_quadric` having solved it exactly for years. Both were
fixed the same way: use the closed form where the shape has one, iterate
only where it does not.

**The speed is the lesser half.** A quadratic solve is arithmetic and can
be written down as a sequence of operations. A convergence loop whose exit
condition depends on its own intermediate values cannot: it lowers to no
kernel, and it serialises into nothing. So an iterative implementation is
**opaque in exactly the sense an opaque leaf is** — the thing this file
already measures when it counts what a scene can export.

Which joins two halves of the project that had been discussed separately.
Dispatch picks an implementation; the implementations that are closed forms
are IR-expressible and the iterative ones are not; so **dispatching toward
a closed form is simultaneously a performance decision and an
exportability decision**, and a search that finds a shorter chain has also
found a more exportable one.

It also gives the exportability measurement a second lever. Until now the
only way to raise the structural count was to rewrite opaque leaves in the
vocabulary. This is another: wherever an operation has a closed form and
takes the iterative path anyway, that is a leaf that did not have to be
opaque.

### Two searchable domains, and why owning the library is what makes them searchable

Recorded 2026-09-18, pulling together threads that were already here separately.

**The reason everything is hand-rolled turns out to be a training reason, not only an independence one.** A search needs to grade its candidates, and grading at the resolution a search needs means running instrumented code you own. On someone else's library an operation is a black box: it can be called, not scored. Here every registered operation is simultaneously a legal action and a measurable outcome, which is the whole substrate. `Vec`, the solvers, the eigendecomposition, `json.hpp`, `wav.hpp` — the argument for those was self-containment, and this is the second one.

There are exactly two places where an oracle is free, exact, and ours.

**Rewriting the IR.** "The trace is the IR" (see the Declarative scene DSL section). A rewritten op graph must compute the identical result, and the oracle for that is already queued rather than new work: the POD interpreter asserting bit-exact agreement with `eval_into`. The test that would demonstrate exportability *is* the grader for a rewrite — one artifact, two uses. This is also, precisely, "arrive at one algorithm's results from a different set of primitives": a rewritten graph computes the same thing by a different composition, which is the shape AlphaTensor's result had.

**Building trees.** The oracle is measured query cost — build the hierarchy, then run queries and time them.

The difference between the two is worth keeping straight, because it decides what each search can be asked for. **IR search can only equal**: bit-exact equality admits no improvement on correctness, so the objective is cost subject to exact identity. **Tree search can beat**: there is no output to match, only a cost to lower, so improvement over the reference is on the table.

**Domain-specific models are not a hypothesis here.** Measured 2026-09-18: one distilled tree per domain loses 0.005 accuracy against the network at 8.43 ns, where a single tree covering all four domains loses 0.008 at 11.18 ns — better and roughly half the cost, because the domain is known at the call site and should not be costing the tree two levels to re-derive. The base-plus-custom split already in the RSC design is the same idea one storey up.

**And the manifold tree is where the transfer claim gets tested cheaply.** A BVH and a ball tree differ only in the bound; construction and traversal are the same. Given a bound concept — bound a set, merge two bounds, test a query against a bound — one hierarchy serves axis-aligned boxes, geodesic balls, OBBs and k-DOPs alike. A split policy learned over that concept on Euclidean boxes can then be run against geodesic balls, and whether it still helps is a measurement, not an argument. That is the space-agnostic discovery claim in its cheapest concrete form, and the compiler does the checking.

**Dependencies, so the order is not confused.** IR search needs the op vocabulary and the POD interpreter, both already queued. Tree search needs the bound concept, which does not exist yet. Neither needs the other, so either can be taken first.
