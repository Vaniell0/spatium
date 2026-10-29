# RSC: Route Synthesis & Certification

RSC finds, chooses and certifies alternative routes through Spatium, against
Spatium's own definitions. A route is a chain of the library's operations or
an assembly of its mechanics; the definition of the same result that the
library already has is the oracle; a certificate says where the route agrees
with it. Nothing outside the library is needed to grade a candidate, which is
why the search can be trusted and the models in it are allowed to be wrong.

This page states what a test or a measurement holds and marks everything
else as design. The prototype's notes (domain by domain, with the debugging
history) are in git: `git show 924b98f:rsc/README.md`.

## The idea

Treat the library as a graph. Nodes are types, spaces and representations;
edges are operations, derivations and adapters (`lower`, `Real50`, `Dual`,
GPU emission). Mathematics says that many pairs of paths between the same two
nodes must agree, and every such pair is a free test
([`docs/open-problem-rsc.md`](../docs/open-problem-rsc.md)). Two of the
library's deepest bugs were found that way, by paths that disagreed, not by a
test anyone wrote.

Every RSC job has three roles:

- **Reference**: a function or a relation the library already defines. The
  library is the oracle.
- **Search**: a way to find another route to the same result. Enumeration, a
  trained policy, minimisation of continuous knobs, a decision tree fitted
  to labels: any of them, and the job declares which.
- **Certificate**: where, and how strongly, the route agrees with the
  reference.

A job without a reference or without a certificate is not an RSC job. The
search is for a cheaper or independent path to an answer that is already
known, not for the answer.

Four verbs, one currency:

| verb | meaning |
|---|---|
| choose | among routes that exist, by features of the query |
| find | build a route that is missing (a red cell of the connectivity matrix) or cheaper than the reference |
| witness | an independent second route becomes a test, even when it is not cheaper |
| lower | the same route in a cheaper form (straightened, IR, GLSL, specialised to a scene), checked by a pair against the unlowered form |

The currency is cost (counters, not wall time) against tolerance (the
certificate). A search that finds nothing within its budget leaves a record
too: the class of the gap, "no glue" or "no footing", is a measured border.

Online, the same machine runs as the loop recorded in
[`docs/findings-2026-09.md`](../docs/findings-2026-09.md) section 4: refine,
select, certify. Subsystems report their error as a certificate, a result's
budget is divided between them, and a choice is made where to refine.

## The contract

**Job.** A declaration with typed fields: the reference (function or
relation, the input region, the tolerance); the components a route may use;
the cost measure; the search organisation; a budget; the certificate level
required; where the artifact goes; the consumer.

**Artifact.** A tuple: reference, route, region of inputs, certificate, cost,
fallback to the reference. The region is a partition of the input manifold
into cells, each with a route and a certificate; outside it the fallback runs.
It ships as a generated header under a test, the pattern of
[`rsc/include/generated/ccd_chooser_tree.hpp`](include/generated/ccd_chooser_tree.hpp).

**Certificate.** Two axes, not one scale. The kind of reference: a function
(agreement with it) or a relation (an axiom such as `exp` after `log`). The
strength of the check: on a sample, or proved over a box of inputs (`Interval`,
not built yet). A relation is not stronger or weaker than a function; it is a
different thing, and an artifact says which it has.

**Consumer.** A call site in `include/` on a path that a test in CI or a
measured binary executes. No consumer, no job: a solution nothing calls is a
note. `tests/symmetry/` counts as a consumer only for the outcome "witness".

**Independence of a witness.** The transitive dependency graphs of the
reference and of the route share no node outside a declared trusted base: a
list of primitives (arithmetic, `exp`, `sqrt`, `acos`, `cos`), each held by
its own pair against `Real50`. The graph comes from the `uses` that components
declare; an instrumented build counts calls by name and fails when the
declaration and the measurement differ. A witness that fails the check stays a
test, marked dependent, and is not counted.

**Freshness.** The pair test is the freshness check: a stale route fails its
comparison with the reference. What was proved once over a box carries a hash
of its dependency closure and drops to "sampled" when the closure changes,
until the proof runs again.

**One discipline.** Production and research are the same standard here.
Nothing enters without a test, a certificate, a consumer and a roadmap entry.
The only asymmetry is the direction of dependency: `rsc` depends on the
library; the library does not depend on `rsc`; a generated artifact is an
ordinary header of the library under the library's tests.

## What exists, and how it is held

| part | state | evidence |
|---|---|---|
| Labels by running candidates against a reference | built | `include/comparison_task.hpp`, `tests/test_rsc_*.cpp` |
| Connectivity matrix: every (scalar, space, probe) cell built and graded L0-L3 | built | `docs/connectivity.md` (1115 of 1168 cells at L2 or above), `tests/test_connectivity.cpp`, regenerated in CI |
| Registry of pairs of paths with certificate tolerances | designed | stage 1 of the plan; `tests/test_space_access.cpp` is the shape |
| A learned choice in a pipeline: the CCD chooser | built; consumer is a test only | `core/concepts.hpp` (`Chooser`), `surface_ccd.hpp` (chooser overload), `include/generated/ccd_chooser_tree.hpp`, `tests/test_rsc_ccd.cpp`; held out, time against the default 0.80 (torus-patch), 0.93, 0.94, 0.96, 0.97, 1.00; a choice per query is worth 5-20 %, the large factors came from the floors |
| Seven dispatch domains (precision, geodesic, root finding, ODE, mesh strategy, linear solve, rigid integrator) | built as models | `include/*_ops.hpp`, `*_task.hpp`, `tests/test_rsc_*.cpp`; the models sit at the ceilings their features allow, except mesh and integrator ([findings](../docs/findings-2026-09.md) section 4); `solve_jacobi` has no caller outside RSC |
| Search over chains: numbers, meshes, CCD primitives | built as tools | `tools/toy_search.hpp`, `tools/mesh_collisions.cpp`, `tools/ccd_synthesize.cpp` |
| Structure transfers between spaces, constants do not | measured | `tools/transfer.cpp`: a relative spec 3 of 3, an absolute one 1 of 3 |
| Calibration of a constant by an exact residual | tool built, not wired to a scene | `tools/metric_calibrate.cpp` (binary-merger boost; result in ROADMAP, 2026-09-29) |
| Cross-entropy against REINFORCE on the base dispatcher | measured | `tools/supervised_vs_reinforce.cpp`: 0.955 against 0.551 at 2000 updates; the trainers still use REINFORCE (finding B10) |
| Distillation of a dispatcher into a decision tree | built | `tools/distill_tree.cpp`; the cost figures did not reproduce exactly (finding B12), the order of magnitude did |
| Where an explicit contact chain replaces Newton | measured, window mapped | `tools/explicit_contact_window.cpp`; one non-monotone cell unexplained |
| `Certified<T>`, `Interval`, proof over a box | designed | plan stages 1b and 4 |
| Jobs, the generator, the grown registry | designed | this page |
| Mechanics assemblies, straightening, the error budget, IR passes | designed | this page |
| Surrogates, IR rewriting, structure building | not attempted | the catalogue in the prototype's notes |
| Browser solver | not started | gated on a WASM build |

## Measured borders

These are recorded in the ROADMAP and the findings, with the tool that
measured each.

- Over numbers a chain space collapses: 9 765 625 paths at depth 10 reach
  13 836 distinct states (706 times). A learned heuristic there is worse than
  breadth-first search: 1.6 to 1.7 times the nodes, and it solves 90.8 % of
  problems against 100 % (`tools/search_heuristic.cpp`).
- Over meshes the branching stays at 4.66 of 6 and exhaustive search ends at
  depth 8 or 9: 0.20 s at depth 5, 10.5 s at depth 7, a node costs 129 µs
  against a policy evaluation of 608 ns.
- Over continuous-collision primitives, 819 chains are 22 to 109 distinct
  algorithms per kind of query once compared by their answers.
- A label must check cost, not only correctness: the first labelling of the
  linear-solve domain preferred an algorithm that is dearer at every size
  tried (finding B5).
- A dispatcher needs a call site with a real alternative.
- A setting changes cost and never correctness only without a budget: under a
  finite budget a policy can change the answer (finding B8).

## The design, in the order it will be built

Everything below is design. Each item names the experiment that will keep or
drop it.

**Where space enters.** The input is a point of a manifold (a contact
configuration, a metric, a chart), so anchors, geodesic nearest neighbours
(`GeodesicBallTree`) and Fréchet means work without a model. The continuous
part of a move is a point of a manifold (a split direction on a sphere, a
rotation), so a policy is a distribution there, through `exp` and `log`. The
parameters of a model are points of a space, so training can follow the Fisher
information instead of the raw gradient, with exact Hessians through nested
`Dual`. The architecture can be built from the problem's symmetry: invariant
features, and a piecewise model whose gate is a partition. The search space
itself collapses by equal results, a quotient of the space of paths.

**The collapse number.** Paths divided by distinct states, measured on short
depths and written into the job, chooses the organisation of the search:
where it is large, exhaustive search over states is cheap and neither a tree
nor a model is needed; where it is small, a value estimate pays.

**A ladder of hypotheses.** A threshold on a derived invariant (a condition
number, a scale ratio) → a decision tree → a boundary (a level set of the
cost difference in one to three features, with a margin band and a fallback)
→ nearest anchors → a small network → search with a depth estimator. The
simplest rung that passes a held-out measure is what ships. Reinforcement
learning, supervised learning with expert iteration, and distillation are
ways to build any rung; reinforcement learning is for rewards that are a cost
with no label.

**Training the proposers.** The trainers are plain SGD on flat vectors and one
network shape for every domain (`include/dispatcher.hpp`). The first
experiment is diagnostic: exact full-batch training, then natural gradient,
against SGD and REINFORCE at equal oracle calls on all seven domains. It
separates domains whose ceiling is their features from those limited by the
optimiser or the architecture. It is dropped if no domain moves. The further
experiments, each with its threshold, are natural gradient where the reward is
a cost, invariant against raw features for transfer between spaces, a
piecewise model against one network and one tree, a policy on a manifold of
actions for building structures, and a hyperbolic embedding of chain prefixes
against a Euclidean one.

**Mechanics.** The ladder of components runs from an expression, through a
library function, to a mechanic (energy, constraint, integrator, nonlinear
and linear solver, ordering, contact) and a scene. A job over mechanics is a
sketch: the skeleton of a step with holes, and the holes are choices. The
references are a slow accurate trajectory over a window and relations (no
penetration, bounded energy drift, momentum). Certified contact is the spine of
every assembly and is never bypassed by a proposal, so the periphery is
searched freely and a wrong proposal is rejected, not turned into wrong
physics. Proof over a box is for primitives and the spine; a whole scene gets
a sample over its class of scenes plus the spine.

**Lowering.** An iteration with a data-dependent exit cannot be an IR
kernel. When a certificate gives a number of iterations that is enough on a
region, the loop becomes a straight chain of fixed length with a guard that
falls back outside the region. Then IR passes (common subexpressions, constant
folding, evaluation schemes; bit-exact first, within a tolerance against
`Real50` once intervals exist), fusion and layout, specialisation to a scene
instance, and placement and precision (a lower precision on the device with a
margin, the certificate on the host). Each step is a pair against the step
before it (`interpret` against `eval_into`, host against device), so the
registry of pairs checks the lowering.

**The error budget of a stack.** Dividing a result's tolerance between layers
is minimising the sum of the layers' costs under a bound on the sum of their
errors weighted by the sensitivity of the result to each; the sensitivities are
exact derivatives through `Dual`. `integrate_volume` already weighs the
statuses of inner integrals by how much they matter. The search runs from the
top down, family of mechanics first, then straightening, then IR passes, then
placement, re-measuring after each, with the spine fixed.

**Growth and shrinking.** A certified route that has shipped becomes a
component with its declared cost, region, certificate and `uses`, so the next
search is shorter. It is admitted when it reduces the size of the library
plus everything expressed through it; two implementations that always agree
leave one of them redundant. The weights of a model are a memory of searches
that makes the next one cheaper to propose; the artifacts are what is
certified. The two grow together, and a pinned feature contract (finding B9)
keeps weights from going silently stale.

**Using it.** Four entries, none of them a public interface yet. A library
user gets `Auto` at a choice point and, on request, the route as steps
(`describe_chain`). A job author writes a declaration and gets an artifact, a
pair in the registry and a roadmap entry. A query through the scene language
becomes a task whose route is composed step by step with explanations and a
certificate. A route query over the graph, how to get B from A, uses the same
machine over the library's own capability index.

## Order

1. The foundation: the registry of pairs, the levels of the matrix,
   `Interval` and `Certified<T>`, cost counters and `uses` on components.
2. First jobs, each with a consumer named: `Method::Auto` of the integration
   door (reference: tanh-sinh at a strict tolerance, with a witness); a derived
   `exp` against the closed forms; the CCD policy in the job format; the render
   level; the explicit contact window; the choices in `include/` that have no
   hook yet (the geodesic step, the ray constants of the black hole, the
   budgets of `rigid_contact`, the bins of the BVH, the leaf of the ball tree).
3. The diagnostic experiment on the weights, which blocks nothing.
4. Filling red cells and witnesses in place of hand-written pairs.
5. Surrogates through `L2Interval`, structure building, IR rewriting.
6. A query as a route, baking, the browser solver after a WASM build.
7. Physics as the consumer: the replacement of IPC is where lowering and the
   budget pay, and what the line shows to be redundant or missing there
   decides the physics concept.

**The removal rule.** What has a consumer is a solution; what has none is a
note. When the contact path of the physics stage is written, the CCD
chooser's tree moves into it as its consumer, and the rest of this directory
leaves the tree: the numbers stay in the findings and the ROADMAP, the code
stays in history (a tag marks the last commit that has it), and a job
re-creates anything that earns a consumer.

## Layout today

The removal rule above was applied on 2026-10-20: the prototype listed below
left the tree, and the tag `rsc-prototype-final` marks the last commit that
has it. The tree now holds `rsc/core` and `rsc/spatium`, described in the
root `CLAUDE.md`.

- `include/`: the registry, tasks, dispatcher, trainer, chains, checkpoint
  and the per-domain ops and tasks; `include/generated/` holds the one shipped
  artifact.
- `tools/`: the searches, trainers and measurements named in the table above.
  Build them with `-DSPATIUM_BUILD_RSC_TOOLS=ON -DSPATIUM_BOOST=ON`.
- `checkpoints/`: the base checkpoint, pinned to a commit and a registry
  snapshot.
- Tests are `tests/test_rsc_*.cpp`.

