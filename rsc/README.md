# RSC

A trained dispatcher sitting on top of Spatium, not a WKV-style learned
recurrence pretending to be a computer. Spatium stays the exact substrate;
RSC only decides *which* Spatium call to make and with *what* parameters.
Framed explicitly as a real, serious project ("не игрушка... промышленный
калькулятор"), parallel to noesis, not a replacement for it.

No formal expansion of the name — it started as "Recurrent Solver Core,"
but that stopped being accurate once "no learned hidden state by default"
was settled (below). RSC is just the name now.

Status: implementation started. `F`, the autodiff foundation both Tier-2 ops
and calibration need, is done (`spatium/algebra/dual.hpp` +
`algebra/calculus.hpp`). Build plan, done items marked, named by content —
not by position, nothing should ever refer to one of these by a number:
registry skeleton (done), Tier-1 task generator (done), minimal feedforward
dispatcher (done), hand-derived backward for the dispatcher (done),
policy-gradient training loop (done), validated accuracy actually rises on
held-out tasks (done — the Tier-1 mechanism-sanity-check plan is complete).

## Registry & dispatch

A fixed jump table, two tiers:
- Tier 1: library-independent math primitives (arithmetic, `poly_eval`,
  `dot`/`cross`/`lerp` — already free functions in `algebra/functions.hpp`).
- Tier 2: Spatium's actual versioned calls (`solve_quartic`,
  `symplectic_step`, `geodesic_distance`, the approximation-family ops
  "fast" needs — see below).

Same table drives the model's dispatch head size *and* the training task
generator — one source of truth, not three designs to keep in sync by hand.

**Built and tested:** `rsc/include/registry.hpp` — `rsc::Registry`,
`rsc::Op`, `rsc::OpSignature`. `add()` returns an op's index (the exact value
a trained classification head would need to produce); dispatch is
`registry[index](in_span, out_span)`; arity mismatches throw immediately
instead of reading past a span. This is meant as a real API, not internal
plumbing — normal application code registers/dispatches through it the same
way training or calibration would. `rsc/include/tier1_ops.hpp` wires in
the first four Tier-1 ops (`add`, `multiply`, `dot3`, `lerp3` — the latter
two wrap Spatium's own `dot`/`lerp`) to prove the API end to end; Tier-1
coverage is not exhaustive yet.

`rsc/include/task.hpp` — the task generator, direct analog of noesis's
`train_chain()` (`experiments/A0_state_probe/micro_wkv.py`) but against this
real registry instead of a two-op toy: `TaskGenerator::sample()` picks a
random registered op, generates random operands matching its declared
arity, and computes the correct output via the registry itself (exact by
construction — no model needed to get ground truth). `grade()` is the hard
pass/fail gate from Training, above: wrong op index or output outside
tolerance both fail, nothing in between. Single-op tasks only for now —
chains (composing several ops) come later, once single-op dispatch works;
starting with chains would mean debugging composition and dispatch at once.

`rsc/include/dispatcher.hpp` — the minimal feedforward dispatcher: one
hidden layer, ReLU, linear output to `registry.size()` logits, no
recurrence. **Important scope note, caught while building this, not before:**
its input is `[operands, observed output]`, not operands alone — for
same-arity Tier-1 ops (`add` vs. `multiply`) the operands alone carry no
signal, since `TaskGenerator` assigned the correct op independently at
random. Given the observed output too, "which op reproduces this
input→output pair" is genuinely learnable (`target == a+b` vs.
`target == a*b`), not a coin flip. This is a mechanism sanity check — does
the forward/backward/training loop work at all — not a preview of Tier-2's
real dispatch signal, which reads a property of the input alone (e.g. a
discriminant's magnitude, for polynomial/precision-critical dispatch) with
no "observed output" involved.
Both shapes are legitimate; they're just different problems, and conflating
them would have quietly made the wrong one the template for Tier-2 later.

`rsc/include/dispatcher.hpp` also carries the hand-derived backward
pass now: `forward_cached()` keeps the intermediates, `backward(cache,
dlogits)` gives closed-form gradients for `w1/b1/w2/b2` (matmul/ReLU/matmul
in reverse — no general reverse-mode AD engine). `cross_entropy_loss()` /
`cross_entropy_dlogits()` are the categorical-policy loss shape both a
numerical gradient check and the coming REINFORCE loop use (reward=1
recovers ordinary supervised cross-entropy, which is what gradient-checking
needs). Validated by finite-difference gradient checking across every
weight, 3 random cases — not asserted, actually checked
(`tests/test_rsc_dispatcher.cpp`).

`rsc/include/train.hpp` — the policy-gradient (REINFORCE) loop and its
validation. Tried plain REINFORCE first, exactly as planned, and *measured*
it rather than assuming it would work: eval accuracy rose to ~0.61 then
regressed back to the ~0.20-0.25 starting point over further training, and
learning rates above 0.01 never learned at all. Cause: reward=0 (a wrong
guess) gives *zero* gradient with no baseline, so wrong guesses carry no
corrective signal at all — a known REINFORCE weakness, confirmed real here,
not assumed. Fixed with a running-average reward baseline
(`advantage = reward - baseline`) — wrong guesses now get a negative
advantage instead of none. Result, fully deterministic (every seed fixed):
0.203 pre-training → 0.560 post-training accuracy on 5000 held-out tasks
(different seed from training) after 50k steps at `lr=0.01`
(`tests/test_rsc_train.cpp`). This closes the Tier-1 mechanism-sanity-check
plan end to end — forward, backward, and training all validated, not just
compiling.

`rsc/include/describe.hpp` — friendlier input/output, resolved as a
separate descriptive layer rather than changing the numeric path at all:
`OpSignature` carries optional per-slot names (`"add"` → `a`,`b` → `sum`),
and `describe()` formats any op call as `"add(a=2, b=3) -> sum=5"` instead
of raw indices — the original doc's "показывает шаги решения" goal, and the
layer a future DSL front-end naturally produces/consumes. `features()` (the
padded vector actually fed to the dispatcher) is untouched; this is purely
additive.

No parser: the model's output head *is* the jump-table index (classification
over K) plus per-op continuous parameters (regression head), not text that
gets parsed into a call. Removes a whole class of tool-calling failure modes
(bad arity, hallucinated function names) by construction — the head's range
is exactly K.

## DSL — task input/output formalization

A standing idea: RSC needs a formal language for *describing a problem* on
the way in and *reporting the answer plus the solution steps* on the way
out — "нам нужна инфраструктура не меньше чем сама модель." Not tied to the
original v0.3 WKV draft's own DSL sketch (`Вход (DSL/ASCII Math) →
[Парсер] → ... → Learned Solver Head`) — that architecture, parser
included, is exactly what's since been superseded by the Spatium-registry
design in this file. Follows the *current* principles, not v0.3's.

This doesn't conflict with "no parser" above — that's about the model's
*output dispatch decision* specifically (no free text there, ever). A DSL's
own parser is ordinary deterministic software (a grammar/AST reader, like
any compiler front-end), not an ML model guessing at ambiguous text — no
hallucination risk, different end of the pipe entirely.

Resolved once chains gave it something real to describe (Chains, below): no
text grammar at all — the DSL *is* the structured `Task`/`ChainTask` record
(named fields: which op, what operands, what result), and `describe()`/
`describe_chain()` render it as readable text on the way out
("3 -> add(a=3, b=2) -> sum=5, then ..."). A single Tier-1 op had nothing
to formally describe, which is exactly why this waited for chains instead
of being designed in the abstract beforehand.

## Training

Two genuinely different mechanisms, only one of which is "the model."

### Dispatch — the model, trained

Dispatch (which op/family) is a hard discrete choice through
non-differentiable C++ — no gradient flows through execution, so it trains
via RL/search (AlphaZero/AlphaTensor/AlphaDev-style self-play), not
end-to-end backprop. Spatium's own exactness is the reward signal — no
hand-labeled data needed. Speed/compute is optimized only as a secondary
objective, gated hard by accuracy first (AlphaDev's own reward shape:
correctness is a pass/fail gate, latency is what gets minimized among
passing candidates) — a "catastrophic" accuracy loss disqualifies a
candidate outright, it isn't a point on a soft trade-off curve. The
tolerance itself is task-relative, not one global constant — reuses
Spatium's existing `epsilon<T>()`/`relative_epsilon<T>()`
(`core/epsilon.hpp`) for the type-precision part, but the calling context or
the specific Tier-2 op declares what error is acceptable for that case.

### Calibration — not the model, classical optimization

Once dispatch has picked a family/op, fitting *its* continuous parameters
is an ordinary nonlinear optimization problem (Newton/gradient descent),
not something that needs learning or generalizing across cases — it gets
re-solved fresh, exactly, every time it's invoked.

Concretely (the actual XPBD contact bug): dispatch picks "resolve this
contact with XPBD, parameters (compliance α, friction μ, damping β)."
Define `L(α,μ,β)` = simulate N steps, measure penetration depth. Because
the whole simulation is already `template<Scalar T, ...>`, seeding one
parameter at a time as `Dual<T>::variable(...)` gives `∂L/∂α` (etc.)
straight from `calculus.hpp`'s `gradient()` — one simulation run per
parameter, no backprop through a network, no training data. The result gets
cached per scene-class (see Model shape, below). No neural network involved
anywhere in this path — "the model" only ever covers the discrete choice
above.

Built and tested: `calculus.hpp`'s `minimize(f, theta0)` — gradient descent
with Armijo backtracking line search, calling `gradient()` internally. Not
yet wired to a real Spatium scene (XPBD above is still the target, not done
— tested so far against toy losses with known closed-form minima).

Curriculum is one training run, not phases with a hard switch: general math
(Tier 1) early — required regardless of target domain, no task is solved
correctly without it — Tier 2 mixed in with growing composition depth, Tier
1 kept in the sampler forever at a floor quota (replay) so base competence
doesn't get crowded out. Plain catastrophic-forgetting mitigation, not a
novel mechanism.

Pure C++, no Python pipeline anywhere — hand-rolled (reusing `Dual<T>`'s
machinery for the small policy net), not a libtorch/PyTorch dependency.
Keeps the header-only, build-it-from-concepts story intact end to end.

No learned hidden state by default. Alpha*-style: a stateless evaluator
navigating an explicit search/composition process, not a recurrent net
carrying an invisible hidden vector — Spatium's own numeric state already
carries what a program needs. Revisit only if a concrete task demonstrates
the visible state is insufficient. Expected scale is small either way —
hundreds of thousands to low millions of params, nowhere near noesis's
billion-scale WKV.

Toy-scale precedent for the hard-dispatch mechanism already exists in
noesis: `ModularChainController` in `experiments/A0_state_probe/micro_wkv.py`
(`torch.where` masking by `op_code` between op branches) — built, not run,
two-op toy only. `train_chain()` in the same file is the Phase-1 curriculum
prototype.

## Chains

Multi-step composition, generalizing Tier-1's single-op sanity check
(dispatch trained, calibration classical, all above) to a sequence sharing
one evolving accumulator instead of independent random tasks. Built in
`rsc/include/chain.hpp`.

`ChainTask` is a sequence of ordinary `Task`s where step *i*'s inputs are
`{accumulator_i, operand_i}` and its expected output is `accumulator_{i+1}`
— reuses `Task`/`grade()`/`features()` completely unchanged; a chain step
and an i.i.d. Tier-1 task are graded and trained identically, only how
they're generated differs. `Registry::chainable_ops()` selects which ops
can participate — derived from arity alone (`in_size==2 && out_size==1`,
so the output can feed back as the next step's accumulator), not tagged by
hand; for the current registry that's `add`/`multiply` (`dot3`/`lerp3`
don't have this "state in, state out" shape and are excluded
automatically). Teacher-forced: the accumulator always follows the chain's
own ground-truth trajectory, never the model's actual prediction, so one
wrong guess at a step doesn't corrupt the next step's input — free-running
(the model's own choice determines the next state) is a real next step,
not attempted here, since it needs credit assignment across steps.
Chain length is fixed per call, not adaptive; the model deciding when to
stop is a real next step too, not built yet.

**Real finding, not assumed:** single-sample REINFORCE (`train_chain()`,
same mechanics validated on Tier-1) was tried first and *measured* on
chains — restricted to only add/multiply, with no dot3/lerp3 to make the
overall number look better, it's exactly the hard 50/50 case Tier-1 could
partly hide inside its aggregate accuracy. Accuracy stayed stuck at the
~0.50 chance level over 20k+ steps, never learning `target==a+b` vs.
`target==a*b` — reward=1 samples carry a real but very high-variance
single-sample gradient, too noisy alone to make progress on that
genuinely nonlinear decision. Fixed with minibatching (`train_chains_batch()`
in `chain.hpp`, `train_batch()` in `train.hpp` for the i.i.d. case):
average the gradient over a batch of ~32 samples before one weight update.
Measured result, fully deterministic: 0.312 pre-training → 0.972
post-training per-step accuracy on 500 held-out 3-step chains, 1200
batched updates at `lr=0.1` (`tests/test_rsc_chain.cpp`). Both the failure
and the fix are asserted in tests, not just the happy path.

`describe_chain()` is the DSL's output side (above) — renders a solved
chain as `"3 -> add(a=3, b=2) -> sum=5, then multiply(a=5, b=4) ->
product=20"`, reusing `describe()` per step.

## Comparison-based dispatch: one generic template, not one class per domain

Both domains below independently arrived at the same shape: compare N
candidates against a reference, dispatch to the cheapest one accurate
enough. `rsc/include/comparison_task.hpp`'s `ComparisonTaskGenerator<Problem,
Output>` writes that `sample()` once; `precision_task.hpp`/`geodesic_task.hpp`
are now thin configuration (candidates, a problem sampler, a feature
extractor, a reference, a distance metric) instead of two full bespoke
classes. Expensive setup (geodesic's Heat-method Cholesky factorization)
stays cheap per-sample by precomputing once into closures the template
just calls — it doesn't know or need to know whether those closures are
backed by fresh computation or a cache.

Real regression caught refactoring this, not assumed away: the selection
criterion started as a fixed-tolerance gate (`distance(candidate, reference)
<= tolerance`), which happened to be correct for precision dispatch (the
most-capable candidate *is* the reference there, distance 0) but silently
changed geodesic dispatch's actual ground truth (originally "whichever
candidate is closer to exact," not "is this candidate within an absolute
epsilon of exact") — caught immediately by the existing regression tests
(`test_rsc_geodesic.cpp`'s `coarsest_heat_frac < 0.1` failed). Fixed by
comparing each candidate's error against the *last* candidate's own error
instead of a fixed constant — reduces to the gate for precision dispatch
and to "closest wins" for geodesic dispatch, both correctly, not by luck.

The reusable win this unlocks: a third comparison-based domain now costs
configuration, not a new class.

## Domains, in pipeline order

Not arbitrary — ordered by how proven the mechanism each one needs is, base
math/physics competence (above) being the one prerequisite common to all.
Named by content below, not by position — nothing outside this list should
ever refer to one of these by a number.

- **Polynomial/precision-critical dispatch** — done and tested, the first
  Tier-2 domain built. `rsc/include/precision_ops.hpp` registers
  `solve_cubic_f64`/`solve_cubic_real50` (same math, different precision —
  a genuinely different registry shape than Tier-1's "which op produced
  this" ops). `rsc/include/precision_task.hpp` generates cubics with
  three real roots (casus irreducibilis, always) in two deliberately mixed
  regimes — widely-separated roots (double is fine) and closely-spaced
  ones (double and Real50 measurably disagree) — ground truth comes from
  comparing the two candidates against each other, not from a hand label.
  Getting there required first fixing a real bug in Spatium itself:
  `solve_cubic<Real50>`/`solve_quartic<Real50>` didn't compile at all on
  any branch — qualified `std::` math calls bypassing ADL, `std::numbers::
  pi_v<T>` restricted to standard float types, and Boost.Multiprecision's
  `number<>` arithmetic returning lazy expression types that broke template
  deduction across `solve_quadratic`/`solve_cubic`'s cross-calls (fixed in
  `polynomial.hpp`/`complex.hpp`, covered in `tests/test_polynomial.cpp`'s
  new Real50 cases). **Second real finding, not assumed:** the first
  task-generator version sampled the
  root cluster's center from a wide `[-10,10]` range and the dispatcher
  never learned past chance (~0.50) over 1000+ batched updates — absolute
  root position dominated the raw coefficients' scale and swamped the
  actual separation signal the label depends on. Narrowing that range to
  `[-1,1]` removed the confound; the same raw coefficients as features,
  same minibatched REINFORCE already validated on chains, reach ~0.99
  held-out accuracy (`tests/test_rsc_precision.cpp`). No engineered
  discriminant feature needed once the confound was gone. Smallest,
  cleanest registry — the starting point, as planned.
- **Geodesic/mesh pipelines** (`geodesic.hpp` Dijkstra vs. Heat method) —
  done, tested, and only partially working — the domain the "decision
  features + typed execution" architecture question was actually about,
  proven out but not yet at full accuracy. `rsc/include/geodesic_task.hpp`:
  a mesh doesn't fit `rsc::Registry`'s `span<const double>` shape, so this
  reads a fixed-size summary (vertex count) extracted from the real typed
  `Mesh`/`Surface` objects and calls `geodesic_distances()`/
  `heat_geodesic_distances()` directly — no generic "typed registry"
  abstraction built, deferred until a second domain needs the same thing.
  Ground truth: compare both methods against the sphere's exact closed-form
  (great-circle) distance, same methodology `tests/test_heat_geodesic.cpp`
  already uses for verification, used here to generate labels instead.
  **Two real, separately-found bugs, not one:** (1) same scale confound as
  precision dispatch — raw vertex count spans 12..2562, blew up training;
  fixed with `log(vertex_count)`. (2) a deeper one log-scaling alone didn't
  fix: `log(vertex_count)` is *always positive*, and for a single-feature
  input, `ReLU(w·x)` for `x>0` is decided purely by `sign(w)` at
  initialization — roughly half the hidden units are permanently dead
  before training even starts, and a bad initial sign can't self-correct
  since dead units pass no gradient back. Watched training make this
  *worse* over 500 updates before diagnosing it. Fixed by centering the
  feature (`log(vertex_count) - mean`), which restores real positive/
  negative variation and normal ReLU on/off behavior.

  Third real bug, found after the above two: training then plateaued
  around 0.6-0.67 (the ~64% majority-class/"always predict Heat" baseline),
  never finding the one level (unsubdivided icosahedron, 12 vertices)
  where Dijkstra is correct ~100% of the time — a class-imbalance/local-
  optimum issue minibatch REINFORCE fell into with uniform-across-levels
  sampling. Fixed with class-balanced sampling: sample 50/50 between
  "Dijkstra correct"/"Heat correct" buckets instead of uniformly across
  levels (both partitions precomputed once, at construction, same "closest
  wins" comparison the generic template's own `sample()` uses at grading
  time). **Partial, honestly measured result, not a full fix:** by 2000
  updates the dispatcher correctly favors Dijkstra at the unsubdivided
  icosahedron and correctly favors Heat at the two finest levels; only the
  middle level (vc=162, ~88% Heat) still leans the wrong way, pulling the
  bucket-weighted aggregate accuracy down to ~0.536 (up from ~0.47
  pre-training, a real but partial improvement over the old ~0.64
  majority-only ceiling once the metric itself is measured consistently —
  see `tests/test_rsc_geodesic.cpp` for why the *old* per-level percentages
  aren't directly comparable anymore: balanced sampling deliberately
  distorts the observed label frequencies it trains on).
- **Root-finding dispatch** (`rootfind_ops.hpp`) — Newton vs. bisection
  on `f(x)=x^3-a` (root = `cbrt(a)`, exact closed form). Chosen because
  `f` has a single real inflection at `x=0` where `f'(x)=3x^2` vanishes —
  Newton's actual textbook failure mode, not contrived. Newton reuses
  `Dual<double>` exactly as `calculus.hpp` does (one Dual-seeded call
  gives value + derivative together); bisection over a generous fixed
  bracket plays the same "reference-quality candidate" role Real50 played
  for precision dispatch. `x0` sampled from two mixed regimes (away from
  the inflection vs. near it). Applied the geodesic-domain lesson
  proactively this time instead of rediscovering it: the feature
  (`|f'(x0)|`, always >=0) is log-scaled and centered by a *measured*
  mean (sampled from the real `x0` distribution at construction), not
  guessed. One real surprise: the "risky" regime was assumed to fail
  >20% of the time; measured ~13% (most of `[-0.15,0.15]` still recovers
  inside the 20-iteration budget) — training still succeeds on this
  minority class, so the test thresholds were widened to match the real
  number rather than the generator being changed to force a rounder
  split.
- **Cauchy/IVP dispatch** (`ode_ops.hpp`, subsumes what was originally a
  separate "orbital mechanics" entry) — Euler vs. RK4 across three test
  families sharing one mechanism: Decay (`y'=-ky`), Oscillator
  (`[x,v]'=[v,-w^2x]`), and CircularOrbit (2D gravity, circular initial
  velocity). Needed a real Spatium primitive first, the same two-layer
  pattern precision dispatch needed: `include/spatium/algebra/ode.hpp` —
  generic `euler_step`/`rk4_step`/`integrate_fixed` over `Vec<T,N>`,
  decoupled from `physics/mechanics/integrator.hpp`'s `PointMass`-tied
  integrators entirely. **Correction made while building this:**
  `physics/orbital.hpp` — referenced by the original "orbital mechanics"
  pipeline entry — turned out to be *atomic* orbital theory (Legendre
  polynomials, Slater's rules), not celestial mechanics; there was no
  existing closed-form two-body solution to reuse. Rather than pull in
  eccentric Kepler-equation solving (itself a root-finding problem),
  CircularOrbit's closed form (`w=sqrt(GM/r0^3)`) is the real, if
  restricted, orbital-mechanics family used here. Feature = `log(char_freq
  * dt) - mean` (freq*dt predicts whether a step falls inside Euler's
  bounded stability region), same proactive log+center discipline as
  root-finding, plus two one-hot bits distinguishing which family a given
  freq*dt belongs to (the three families' error curves aren't identical).
  All test thresholds passed on the first measurement this time.
- **Mesh-strategy dispatch** (`mesh_ops.hpp`) — uniform UV tessellation
  vs. an anisotropy-adapted UV step, for one quad cell of a parametrized
  surface. Generalizes the Klein-bottle seam mesh-twist bug (see
  `examples/primitives_demo.cpp`'s `KleinBottle` comment) into a
  measurable, dispatchable quantity: a UV-uniform step produces long,
  thin ("sliver") triangles wherever the parametrization's local fu/fv
  magnitude ratio is high. New Spatium primitive first, same two-layer
  pattern every prior domain needed:
  `ParametricSurface::parametrization_anisotropy(u,v)` (ratio of the two
  eigenvalues of the first fundamental form — 1.0 means locally
  isometric) plus a public `normal_at(u,v)`, both in
  `include/spatium/spaces/parametric.hpp`. Two test families for two
  different real causes of anisotropy: Torus by aspect ratio `R/r` (thin
  ring = severe), Cone by distance from the apex (a genuine degenerate
  point, not just a global ratio).

  **Real, measured correction made while building this** — the first
  design used face-*normal* angular deviation as ground truth (matching
  how `subdivide_adaptive()` already judges curvature). A standalone
  probe showed that error was essentially *identical* across torus
  aspect ratios from 1.2 to 30: sweeping `u` around a torus rotates the
  tangent plane by exactly `du` radians per unit `du` regardless of
  `major_r` — a plain circle in the xy-plane — so normal-angle error
  tracks *curvature*, not the fu/fv *magnitude ratio* anisotropy
  actually measures. Two genuinely different quantities, confused once.
  Switched the ground-truth metric to ambient **edge-length aspect
  ratio** of the UV-uniform quad (`spatium::mesh::face_aspect_ratio()`'s
  own hi/lo convention, computed directly from two surface evaluations)
  — confirmed by the same kind of probe before committing: `uniform_
  aspect` tracks `parametrization_anisotropy`'s own range at every level
  tested (e.g. torus R/r=30 gave anisotropy up to ~31x and measured
  aspect ratio up to ~21x in the same run — same phenomenon, two
  different but correlated exact quantities).

  The "adapted" candidate picks per-axis UV step sizes that equalize
  ambient edge length while keeping the same UV-area budget as the
  uniform quad — a genuinely different triangulation choice, not just
  "more compute" on the same one, the same distinction Newton-vs-
  bisection and Euler-vs-RK4 already draw. `tolerance=0.2`, not the
  domain-wide default `1e-3` — measured: at `1e-3` the "uniform already
  fine" bucket came out empty (every quad needed the adapted strategy),
  which segfaults the stratified sampler's `bucket[idx(rng)]` on an empty
  bucket; `0.2` is the smallest value found, by direct probing, that
  keeps both buckets genuinely populated at every level. Same class-
  imbalance mitigation as geodesic dispatch (50/50 stratified sampling
  between "uniform wins"/"adapted wins" buckets, precomputed once).
  **Measured, not aspirational:** training plateaus around 0.77-0.79
  (checked every 500 updates out to 4000, never crosses 0.8) over a
  ~0.26 pre-training baseline — real, substantial signal, not a
  near-ceiling domain like precision's ~0.99 (`tests/test_rsc_mesh.cpp`).
- **General linear-solve dispatch** (`linear_ops.hpp`) — Jacobi iteration
  vs. direct Gaussian elimination for `Ax=b`, N=4 fixed. Filled a real gap
  surveyed early in the session: Spatium had closed-form `inverse()` only
  for 2x2/3x3 and `determinant()` via LU for any N, but no general solve.
  New Spatium primitive first, same two-layer pattern every domain has
  needed: `include/spatium/algebra/linear_solve.hpp` — `solve_direct()`
  (Gaussian elimination with partial pivoting, exact regardless of
  conditioning, fixed O(N^3)), `solve_jacobi()` (fixed 20-iteration
  budget, no convergence check — same "small enough that a bad case
  measurably fails" discipline as root-finding's Newton budget), and
  `diagonal_dominance_ratio()` (`min_i |A_ii| / sum_{j!=i}|A_ij|`, the
  textbook sufficient condition for Jacobi's guaranteed convergence — the
  real dispatch signal). Direct is always registered last/reference
  (distance to itself is 0, same shape precision dispatch's Real50
  reference has); Jacobi wins whenever it actually converges inside its
  budget. Feature = `log(diagonal_dominance_ratio) - mean`, the same
  proactive log+center discipline every prior always-nonnegative-feature
  domain needed (geodesic's `vertex_count`, root-finding's `|f'(x0)|`,
  Cauchy's `char_freq*dt`) — applied from the start this time, not
  rediscovered as a dead-ReLU bug. Two regimes sampled per matrix (fixed
  diagonal, off-diagonal entries drawn from a "safe" or "risky" scale) —
  measured with a standalone probe before writing the task generator, not
  assumed: at N=4/diag=6.0, off-diagonal scale in `[0.5,2.0]` ("safe") vs.
  `[4.0,7.0]` ("risky") gives a real ~36% needs-direct fraction under the
  exact comparison rule this domain trains on. `tests/test_rsc_linear.cpp`.
- **Rigid-body integrator dispatch** (`integrator_ops.hpp`) — the domain a
  code-level audit flagged as RSC-ready with no new algorithm work: five
  real, independently-implemented steppers already existed in
  `include/spatium/physics/mechanics/integrator.hpp` (`euler_step`,
  `semi_implicit_euler_step`, `verlet_step`, `rk4_step`, `yoshida4_step`)
  for the same `PointMass` state-stepping problem, the identical shape
  Cauchy/IVP dispatch already proved out, just never wired. Three test
  families (`integrator_ops.hpp`), all mass=1/reference-length=1 so a
  single `sqrt(rate)` characteristic frequency applies uniformly across
  all three, unlike Cauchy/IVP's CircularOrbit-only `sqrt()` case:
  UniformGravity (exact closed-form quadratic motion — only Euler and
  semi-implicit Euler have nonzero error here, Verlet/RK4/Yoshida4
  integrate it exactly regardless of step size), Spring (isotropic 2D
  Hookean, reduces to two independent SHM axes), PointGravity (2D
  circular orbit, same shape as Cauchy/IVP's CircularOrbit). Candidates
  ordered cheapest-to-priciest by force evaluations per step (1/1/2/4/6);
  Yoshida4 registered last, playing the same near-reference role RK4
  plays in Cauchy/IVP dispatch. Measured (seed=1, n=20000, not guessed):
  dispatch is genuinely uneven across the five, not uniform — Verlet
  ~63% (exact for the UniformGravity third of samples and the right
  cost/accuracy point across much of the rest), RK4 ~24%, Euler ~8%,
  semi-implicit Euler ~4%, Yoshida4 ~1% (only needed for the hardest
  tail of the sampled range) — every candidate genuinely wins real
  problems, none is dead weight. Same minibatched REINFORCE already
  validated on every other domain: untrained accuracy 4.9% (below the
  20% chance level for 5-way dispatch — a small random-init model isn't
  guaranteed to start at chance), 84.5% after 2500 updates.
  `tests/test_rsc_integrator.cpp`.
- **Contact physics / soft bodies** (`contact.hpp`, `xpbd.hpp`) —
  **investigated and found structurally blocked, not just next in line.**
  `git log` on the deleted `examples/cloth_*_demo.cpp` files
  (`3df4681`->`1b3ce2c`) documents why they were removed: explicit XPBD
  position-based contact races against stiff structural constraints and
  either tunnels or stretches. Built `examples/cloth_sphere_probe.cpp`
  (headless, measures worst-penetration *and* max-stretch-ratio, which
  the removed demos never tracked) and swept 10 compliance/substep
  combinations x 2 contact-interleaving schemes — all 20 failed (stretch
  >=23%, kinetic energy not settling, often exploding, even at
  compliance=0). This confirms the project's own prior postmortem with a
  real sweep instead of one previously-tried point: calibration-search
  cannot fix this until the contact backend itself is swapped for
  `ipc-toolkit` (implicit IPC + Newton) — not a mechanism this pipeline
  entry can attempt yet. `ipc-toolkit` is now a build dependency
  (`-DSPATIUM_IPC_TOOLKIT=ON`) and the implicit-contact Newton solve
  (`E(x) = 1/2(x-x̂)ᵀM(x-x̂) + h²Ψ_elastic(x) + h²B(x)`, rebuilt into
  `cloth_sphere_probe.cpp`) is wired up and confirmed converging on a
  cloth-on-sphere scene — three real bugs found and fixed getting there
  (unchecked-predictor Newton start, a line search that silently
  accepted energy-increasing steps, an oversized CCD min-distance floor
  that deadlocked Newton). Not done: a full multi-config sweep matching
  this entry's own 20-point XPBD discipline — even one reduced config
  took multiple CPU-minutes under ipc-toolkit's TBB-based collision
  detection (thread-pool overhead dominates a problem this small).
  Deliberately deferred, not chased — calibration-search against this
  pipeline (the rest of this domain) waits on that follow-up.
- **Real-time control of complex dynamics** (illustrative example: an
  underwater drone) — hand-tuning breaks across the operating envelope.
  Mechanism: calibration-search *plus* "fast" (below) to replace
  expensive online integration with a cheap surrogate for onboard use.
  Depends on the base+custom split (below) actually working — sequenced
  last for that reason, not because its own mechanism is unproven.

"Fast" (below) is not its own tier in this ordering — once its
approximation-family ops exist in the registry it's dispatch+calibration
like the domains above, so it can be built alongside them rather than
deferred.

## "Fast"

The user's term (from the original RSC doc's embedding note): a closed-form
surrogate for an expensive iterative process, evaluated cheaply many times
without going back through the model — CfC (Hasani 2022, replaced a liquid
time-constant network's per-step ODE solver with an analytical closed form,
~20x speedup) is the concrete precedent.

Not a research problem: classical approximation theory already has the
needed toolbox — Padé approximants, Chebyshev/spectral fits, exponential-sum
(Prony) fits, perturbation/averaging methods, harmonic balance, model-order
reduction. None are implemented in Spatium yet, but implementing them is
ordinary Tier-2 registry work, the same category as `solve_quartic`, not new
architecture. Once registered, "fast" collapses into the exact same
dispatch+calibration+gate mechanism as everything else: dispatch picks which
family, calibration fits its coefficients, the gate compares against exact
Spatium ground truth. If no family in the menu passes the gate for a given
case, RSC doesn't claim "fast" there and falls back to the exact/iterative
call — graceful degradation built into the mechanism, not a coverage gap.

A surrogate's natural representation is a `Function` (already
built in `calculus.hpp`) — cheap to evaluate/differentiate afterward.

## Model shape: base + custom

Same frozen-backbone-plus-adapter pattern noesis and fleeb83's own
state-interface work already use — not a new technique. A base ("болванка")
trained on Tier 1 + broad Tier 2 dispatch competence, plus a custom layer
calibrated per deployment on top of the frozen base:

**First real base checkpoint, not just the design (2026-08-25).**
`rsc/include/base_task.hpp` merges Tier-1, precision, root-finding, and
Cauchy/IVP into one combined classification head (10 ops) and feature
space (22 slots: 4-way domain one-hot + each domain's own native
features in non-overlapping slots) — before this, every domain had its
own separately-initialized `Dispatcher` that never saw another domain's
tasks, which isn't what "base" was meant to be. Geodesic isn't included
yet (depends on `SPATIUM_HAS_EIGEN`; a `Dispatcher`'s dimensions are
fixed at construction, so conditionally changing them per build flag is
a real complication, left as a follow-up). Curriculum is a fixed 0.25
floor per domain (the design doc's "growing depth" schedule isn't
implemented — plain fixed mixing was the simplest thing to try first).
**Honestly partial, not equalized across domains:** precision (~0.98),
root-finding (~0.86), and Cauchy/IVP (~0.70) all train well; Tier-1 lags
at ~0.60 — its add-vs-multiply pair is the one genuinely nonlinear
decision here (same operand pattern, only the input->output relationship
tells them apart) and gets diluted to roughly 1/8 of total training
samples under the fixed curriculum, far less than `chain.hpp`'s own
dedicated experiment spent isolating that exact pair.

**Two real, checked fixes, not guesses, resolving this to a genuinely
good result — user's call ("нет модель первичней") to fix the model
before touching infrastructure.**

(1) A single shared REINFORCE `baseline` scalar (the running-average
reward `train.hpp` subtracts to form the advantage) sits near the
*pooled average* reward across all four domains — this should
miscalibrate the advantage for whichever domain's own difficulty sits
far from that average. Confirmed across a 4-seed sweep at both hidden=32
and hidden=64: per-domain baselines (`base_task.hpp`'s `domain_of()`
routes each sample to its own baseline) are a strict improvement or a
wash everywhere tried, never worse — but only part of the story, since
hidden=64 still trailed hidden=32 afterward (mean tier1 accuracy ~0.37
vs ~0.53 across the same 4 seeds), so capacity was still hurting
something the baseline fix alone didn't explain.

(2) Ruled out "just single-seed noise" with that same multi-seed sweep,
then checked a specific, principled hypothesis: bigger hidden layers let
the easy domains' shared representation converge to an overconfident
policy faster under REINFORCE, starving exploration on Tier-1's harder
decision before it's seen enough signal (premature policy-entropy
collapse, a known RL failure mode) — not a representational-capacity
problem, so not an argument for changing the `Dispatcher`'s architecture.
Added `entropy_beta` to `reinforce_gradient()` (`train.hpp`, default 0.0,
every existing caller unaffected) — a standard entropy bonus, gradient
derived and sign-checked directly (pushes toward a less-peaked softmax
when beta>0). A beta=0/0.05/0.1/0.15/0.2 x 3-seed sweep at hidden=64
showed a clean, monotonic improvement, not a one-off: mean tier1 accuracy
0.39 -> 0.53 -> 0.65 -> 0.73 -> 0.87 as beta rose, every seed better at
each step, never worse.

(3) Given that the real blocker was exploration collapse rather than
capacity, hidden=32->64 was retried (the user's own read: the model is
still small enough that trying more capacity again cost nothing once
the actual cause was fixed) with entropy_beta=0.2. Result: **every**
domain improved, not just Tier-1 — precision 0.98->0.99, rootfind
0.86->0.97, ode 0.70->0.97, Tier-1 0.60->0.81 (measured directly through
the shipped test, not just the diagnostic sweep). `Dispatcher`'s
architecture (one shared MLP, flat softmax over all ops) never changed
through any of this — every gain came from fixing *how the model
trains*, confirming the user's "модель первичней" call was right, and
that the instinct to suspect architecture next was reasonable to check
but wasn't where the actual fix lived. `tests/test_rsc_base.cpp`'s
per-domain thresholds reflect this final, measured result with real
margin.

**Separately investigated while diagnosing why a diagnostic sweep felt
slow, corrected after checking properly:** a standalone `-O0` compile of
the training loop measured ~130s against a standalone `-O2` compile's
~23s for the same workload, which looked like a clean Debug-vs-Release
story. Set up `build-release/` (`CMAKE_BUILD_TYPE=Release`,
`-DSPATIUM_EIGEN=ON`, modules left off) to capture that win for real —
but timing the *actual* project build both ways (Debug via `build/`,
Release via `build-release/`) on the identical current test showed
**no real difference** (21.7s vs 22.7s). The standalone `-O0` result
doesn't reproduce inside this project's real build (which uses C++20
modules; the standalone timing script used plain headers) — cause not
tracked down, flagged honestly rather than papered over. `build-release/`
stays (a Release config is reasonable to have either way), but the
~22s this specific test costs is apparently just the real cost of
256,000 sequential REINFORCE steps through this codebase's type-erased
`std::function` dispatch (~86us/step) — not something Debug vs. Release
explains. The actual lever, if this ever needs to be faster, is the
sweep-level parallelism discussed but not built (independent
seed/hyperparameter configs have zero shared state and could run as
concurrent processes).

- Custom starts as cached coefficients only (the calibration-search-with-
  cache mechanism above) — no change to the dispatcher's weights at all,
  just a lookup table keyed by scene/deployment class. Cheap enough to run
  "on the spot," in the field, not just once at development time.
- A small adapter on top of the frozen base (LoRA-style, or given the base
  is already small, possibly just retraining a head) is a later step, added
  only if plain coefficient-caching proves insufficient for some deployment
  — not built up front.

**Deployment split built and tested (2026-08-25), the question raised
earlier ("вплавить веса в бинарь" vs. a live file) now confirmed by the
user: both, split base vs. custom exactly along the reproducibility
boundary.** Base weights are embedded into the binary at CMake configure
time — `rsc/tools/embed_checkpoint.py` reads a checkpoint (same
plain-text format `checkpoint.hpp` writes) and emits a `constexpr`-array
header (`generated/embedded_base_data.hpp`); the top-level
`CMakeLists.txt` runs it against `rsc/checkpoints/base_v1.checkpoint`
every configure, same pattern as the git-SHA capture, with a
`SPATIUM_HAS_EMBEDDED_BASE 0` stub if no checkpoint exists yet (fresh
checkout before `train_base` has run). `rsc/include/embedded_base.hpp`:
`load_embedded_base()` (zero runtime file I/O), `embedded_base_is_current()`
(compares the embedded checkpoint's commit SHA against the current
build's — real check, not assumed, catches the checked-out code having
moved on since the embedded base was trained). The custom layer stays
genuinely separate, real file I/O on purpose — it has to update in the
field without recompiling, the base doesn't. `rsc/include/
custom_layer.hpp`: `CustomLayer` (scene-class key -> calibrated
coefficients, plus which base commit it was calibrated against),
`save_custom_layer`/`load_custom_layer`, `custom_layer_matches_base`
(the "base change invalidates every custom layer built on it" check
from this section, above). **Honestly, not yet exercised by a real
scenario:** tested against a synthetic entry only — no concrete
calibration case is wired to this yet, since contact physics (the
original motivating case) is blocked pending the ipc-toolkit backend
swap. 759->764 tests, zero regressions.

## Reproducibility

Ground truth *is* Spatium's own code here (the reward compares against
exact Spatium computation) — so unlike typical ML, changing Spatium changes
both the environment and the definition of correct at once. Every base
build records the exact Spatium git commit SHA *and* a full registry
snapshot (op name → index → parameter arity — index reshuffles must be
caught too, not just behavioral changes). Loading a checkpoint whose
recorded commit/registry doesn't match the current build doesn't silently
proceed — it's flagged for revalidation. Custom layers are versioned
independently of the base: a base change invalidates every custom layer
built on it; a single deployment's real-world drift invalidates only that
one custom layer.

**Built and tested (2026-08-25), not just designed.** `rsc/include/
checkpoint.hpp`: `RegistrySnapshot`/`RegistryOpEntry` (name, global
index, arity — built fresh from `build_tier1_registry()`/`build_
precision_registry()`/`build_rootfind_registry()`/`build_ode_registry()`
at their `kBaseOpOffset`s, not a literal merged `Registry` object, since
`train.hpp` never needs one — see the base-model section above),
`save_checkpoint()`/`load_checkpoint()` (plain-text, `max_digits10`
precision for exact double round-trip — checked directly, not assumed),
`validate_checkpoint()` (`Match`/`CommitMismatch`/`RegistryMismatch`,
the index-reshuffle case checked explicitly by swapping two ops' indices
in a test and confirming it's caught). The commit SHA is captured at
CMake configure time via `execute_process(git rev-parse HEAD)` into a
generated `build_info.hpp` (top-level `CMakeLists.txt`) — baked into the
binary, not computed at runtime (a deployed binary shouldn't need the
git repo present). `rsc/tools/train_base.cpp` runs the actual winning
config (hidden=64, entropy_beta=0.2, 8000 batched updates) and writes a
real, pinned checkpoint to disk, self-validating immediately after
saving. 755→759 tests, zero regressions.

## Site

The reason RSC exists beyond "a trained dispatcher is interesting on its
own": a Wolfram-Alpha-like static site (GitHub Pages), built *after* RSC
works in genuinely sufficient volume — not before, and not something
being built yet. Two zones:

- **Решатель (Solver)** — Wolfram-like in usage: input a problem, see the
  dispatch decision *transparently shown* (which method/precision RSC
  picked and why — a real differentiator vs. Wolfram's black box), get
  the exact answer plus step-by-step (`describe()`/`describe_chain()`
  already render this, at zero extra cost), plus optional visualization
  for geometry/mesh domains.
- **Галерея (Gallery)** — showcase demos (Klein bottle, black hole
  ray-tracing, geometric deformations) demonstrating Spatium's actual
  differentiator: arbitrary M->N-dimensional immersions, concept-based
  generalization across Euclidean/Hyperbolic/Sphere, near-miss detection
  via complex roots (`ray_quadric_proximity`). Not solving user tasks —
  drawing contributors, who then deepen the geometry side.

Real Wolfram Alpha limitations this design turns into non-tradeoffs by
construction, not generosity: free-tier complexity/time limits (server
cost, since Wolfram runs centrally) and step-by-step gated behind Pro.
RSC runs client-side (once WASM-ported, see below) — no shared server to
protect, so both are free by construction.

Task-menu categories, spanning what Spatium actually supports (not "all
of math"): algebra/equations (polynomials, `Complex`, `Dual`, small
linear algebra), geometry (primitives, intersection/distance, boolean
ops, convex hull, analytical ray-quadric), mesh/computational geometry
(triangulation, subdivision, geodesics, DEC, Voronoi), physics
(rigid-body integrators, contact/XPBD once unblocked, orbital-style ODE
dispatch), precision (Real50/Real100 — a cross-cutting quality axis, not
its own category), spatial acceleration (BVH).

MVP criterion: 2-3 reliably-working domains with transparent dispatch and
shown steps — not coverage of everything above.

**Not started, flagged plainly:** GitHub Pages has no backend, so
Spatium+RSC need a WASM build to run client-side at all, and the Vulkan
viewer needs a separate WebGL/WebGPU port — not just a recompile. This is
the concrete gate on "RSC works in sufficient volume" turning into real
site work, not a date.
