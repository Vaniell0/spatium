# RSC: Route Synthesis & Certification

RSC is the glue of the library: within what Spatium defines, it can find a route that is missing, choose
among the routes that exist, turn a second route into a test, and lower a route into a cheaper form, and it
learns from every search it runs. Its line goes back to AlphaDev and AlphaTensor on one side and Wolfram on
the other: a search over a library's own operations, answering with the steps it took and how far they can be
trusted. Put to a project, it is a set of tools that keeps the project in order and helps it grow: it finds
where the project is wrong or says the same thing twice, fills the gaps between implementations, and takes no
change without a check.

The first consumer it is built toward is the black-hole scene: light bent by a real metric, through analytic
pieces RSC synthesised and certified, fast enough to draw because what was learned is baked into the scene.

## Why it can be trusted

The library is a graph. Its nodes are types, spaces and representations; its edges are operations,
derivations and adapters (`Real50`, `Dual`, lowering to GLSL). Many pairs of paths between the same two nodes
must agree, and every such pair is a test nobody had to write. Two of the library's deepest bugs were found
that way, by paths that disagreed. So the library is its own oracle: a search is for a second path to an
answer already defined, never for the answer.

Every job has three roles:

- **reference**: what the library already defines, a function or a relation;
- **search**: a way to another route to the same result;
- **certificate**: where, and how strongly, the route agrees with the reference.

A job without a reference or a certificate is not an RSC job. Because nothing is taken without a check,
whoever proposes may be wrong, and a model may propose before it is any good.

Research here is production: what a search finds enters the library only as a generated, reviewable,
regenerable artifact held by tests, under the same rules as code a person writes. The only asymmetry is the
direction of dependency: `rsc/` depends on the library, the library never on `rsc/`.

Who owns what. RSC holds the solutions and the tools: search, synthesis, certification, lowering, the loop and
its judge. The library holds the mathematics they run on, models included: a softmax over a linear score, a
tanh net and a score tree (`include/spatium/algebra/models/`), the step rules AdaGrad, Adam and Muon
(`include/spatium/algebra/analysis/`), fitting, entropy, cross-entropy and KL (`include/spatium/statistics/`),
each held by two paths like any other domain of the library. RSC works because the library can train models
by its own means.

## One loop, four quarters

Each quarter is a verb, not a method. Synthesis can enumerate, search a tree, fit, or follow a model;
certification can sample, bound over a box, or set two paths against each other; lowering can straighten,
rewrite the IR or specialise; learning can fit a table, a tree, a net or a loop. Which method answers a job is
chosen by measurement on that job, the cheapest one the check accepts, and not fixed in advance. The tree
search is both the teacher of what is learned and its second check.

| quarter | what it does | built | where |
|---|---|---|---|
| synthesise | finds a route that is missing, or cheaper than the reference | routes over the typed graph, run and certified; particular analytic solutions on an interval, in pieces where one expression does not reach the tolerance; the laws of the spaces, found by enumeration | `rsc/core/graph.hpp`, `rsc/spatium/rsc_find.py`, `rsc/spatium/rsc_laws.py`; the synthesis of particular solutions in `include/spatium/answer/synthesis.hpp`, in the library because the site's page calls it |
| certify | says where and how strongly a route agrees with the reference | pairs of paths with the bound each path promises; certificates on a dense sample; numbers whose proved digits are counted by `Interval`; a judge that takes a change only on checks | `tests/symmetry/`, `include/spatium/core/interval.hpp`, `rsc/spatium/rsc_judge.py` |
| lower | the same route in a cheaper form, checked by a pair against the form before | IR fields lowered to POD and GLSL, the interpreter against the lowered form | `include/spatium/io/field_glsl.hpp` |
| learn | keeps a memory of searches that makes the next one cheaper | the tree search makes the targets (the distribution of its visits, the result of a game); a student is fitted to them and distilled down to the cheapest rung that still agrees with the search, its cost counted in bits; the derivations of the laws written out as training items | `rsc/core/tree_search.hpp`, `rsc/spatium/rsc_laws.py`; the models it fits are the library's |

A search ends in one of four outcomes, the verbs of RSC; lowering is both a quarter of the loop and the outcome
that puts a route into a cheaper form:

| verb | meaning |
|---|---|
| find | build a route that is missing, or cheaper than the reference |
| witness | an independent second route becomes a test, even when it is not cheaper |
| choose | among routes that exist, by features of the query, learned and then distilled into a tree |
| lower | the same route in a cheaper form (straightened, IR, GLSL, specialised to a scene), checked by a pair |

The currency is cost, counted in calls rather than time, against tolerance, given by the certificate. A
search that finds nothing within its budget leaves a record too: the gap it hit is a measured border. Every
decision a verb makes leaves a row to learn from: what was there, what was chosen, what the check said.

## The game

Two players and a judge.

- The **breaker** looks for an input, or a sequence of calls, on which a law or a contract fails: degenerate,
  extreme, coincident, aligned, not finite. Every case it finds is shrunk to the smallest one and becomes a
  regression test, so the tests grow from real violations.
- The **fixer** changes the code until the case holds and nothing else got worse.
- The **judge** decides by running checks, and cannot be talked round.

A player is scored against the other player, never on an absolute scale.

## Who runs the loop

An agent owns it: it reads the state of the project, picks the next move (break here, fix there, a test for
this clause, an abstraction in that cluster), calls the tools and the judge, and records the game. The tools
are its hands and the judge its conscience; its goal is the gain the judge measured, and its memory is the
map, the laws, past verdicts and the loopholes found. Today the agent is a general model run headless; the
first fix it made alone, the hyperbolic distance of a point from itself, was taken by the judge.

On a project, the part of the loop that runs today is the loop of tasks over the laws:

1. **Survey.** `rsc_laws.py discover` enumerates the laws and leaves them as data; it is dear and is run when
   a branch is grown.
2. **Regression.** `rsc_laws.py check` evaluates that data again on every change, a column per object.
3. **Task.** A finder writes a task: an object that breaks what its neighbours keep, a contract no test
   reaches, a meaning written twice.
4. **Proposal.** `rsc_propose.py --task N` hands it to the fixer in a work tree of its own, with nothing it can
   push.
5. **Verdict.** `rsc_judge.py BRANCH --task N` decides; a branch it accepts becomes a pull request.

The judge checks, in order: the size of the change, every test of the base kept by name, the build, the
tests, the task done, no law lost where the change reached. Each way a branch slipped through becomes a rule
and a case the rule must refuse. A finder earns its place by the past: replayed on the trees before the
project's fixes, it must catch some of them.

## The place of a model

This is settled and mostly not built: today the tree search and the models it trains exist, the store of
weights and the baking do not.

A model is the memory of searches. Its weights hold what the searches over the library found, are fetched
when a project is built rather than retrained, and are cached; they accumulate, they are not a disposable
cache. One body of weights serves many policies, each a way of doing one verb on one kind of task. It is an
expectation, what should hold before it is checked: it says where to look and what to propose, and it is never
the judge. It takes a slot only by beating the rule that holds it, on the ladder threshold, decision tree,
boundary, nearest anchors, small net, search with a value estimate; the simplest rung that passes a held-out
measure is what ships. The teacher is the certificate, and tree search turns an exact, dear referee into
targets a cheap model can learn: the full distribution of the search's visits, not one chosen move.

Cross-entropy and KL are the one measure across the loop. Cross-entropy in bits is the length of a route's
description under a model; KL is how far a student is from its teacher. Synthesis keeps the shorter of two
expressions that agree, distillation the cheaper rung at the same KL, growth the component that shortens the
library with everything written through it: the same rule each time.

What is learned is baked: a model assembles the choice or the chain, and the scene or the program then runs
without the model in the loop.

## How it is used

- **From the library.** A choice point takes a chooser (`Chooser` in `include/spatium/core/concepts.hpp`); what
  RSC learns there ships as a generated header under the library's tests, as the prototype's policy for
  continuous collision once did.
- **As a question.** A formula or a question goes in; the answer comes back with its steps, its certificate
  and pieces where no one expression holds (the answer door in `include/spatium/answer/`, the site's calculator).
- **On a project.** The loop of tasks above, on Spatium and on GoodNet, a plugin network kernel.
- **As a route query.** How to get B from A over the typed graph: `rsc_route`.

## What is measured, good and bad

- Over numbers a chain space collapses: 9 765 625 paths at depth 10 reach 13 836 states, and a learned
  heuristic there needs 1.6 to 1.7 times the nodes of breadth-first search and solves 90.8 % against 100 %.
  The ratio of paths to distinct states picks how to search.
- A chain found on one space, replayed on others, meets a relative specification three times in three and an
  absolute one once in three: structure transfers, constants do not.
- Laws found on samples are mostly chance until a number's degree in length is part of its type and the
  breaker moves the inputs: 970 classes of spaces without them, 25 with.
- A fix that removes a defect shortens the description: the hyperbolic distance made exact took the laws from
  834 to 637. That fix was proposed by a general model run headless and taken by the judge.
- Learning from the search's answer beats learning from reward alone: cross-entropy 0.955 against REINFORCE
  0.551 at 2000 updates on the same dispatcher.
- On chess a linear policy is cheaper than any tree that imitates it: 416 bits against 13 403 at nearly equal
  agreement, so a tree is for distilling what is not linear.
- An attempt by the fixer is bounded by building, not by thinking: two builds from cold are about 25 of its 50
  minutes.
- The law check would not have caught the library's past fixes: 0 of the 13 it could build, and its tolerance
  of 1e-6 cannot see a fix at 1e-10. One finder covers a narrow slice.
- The shape of the call map says little about where bugs were: AUC 0.65-0.72 here and 0.57-0.58 on GoodNet
  against size, as the literature finds; history (churn, files changed together) is the stronger signal.
- Why bugs happened, in 113 past fixes of Spatium and GoodNet: a contract no code checked at a boundary, and one
  meaning written twice (one shutdown fault fixed three times, in three sibling plugins).
- A low-rank model of the law table is better calibrated than the rules and no better at anything else: the
  table holds about 2 kbit, too shallow for a model.
- On GoodNet, with sanitizers: 1523 tests, a leak at shutdown (a task posted before shutdown that never runs
  keeps its state), and the project's own sanitizer build no longer compiling with its compiler; the shared
  teardown suite runs on 4 of its 9 link plugins.

## Where things are

`rsc/core` is generic: a typed graph with the query for routes, and the tree search; it includes nothing of the
library. `rsc/spatium` binds it to this library: the graph (`gen_graph.py`), the routes and their witness pairs
(`rsc_find.py`), the laws, the judge and the fixer, each described in the root `CLAUDE.md`. The answer door the
site calls, with the synthesis it runs, lives in the library, `include/spatium/answer/`: the boundary is who waits,
not whether there is a search. What runs while a person waits, bounded and reproducible, is the library's; open-ended
search, learning and the loop over a project are `rsc/`, and what they find reaches the door as artifacts. The design not built yet, the
contract of a job, how RSC is brought into another project, and the prototype of 2026-09 are in
[`docs/roadmap/rsc/README.md`](../docs/roadmap/rsc/README.md); where the idea came from, in
[`docs/open-problem-rsc.md`](../docs/open-problem-rsc.md).
