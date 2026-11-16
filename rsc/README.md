# RSC: a project that checks itself

RSC (route synthesis and certification) is the tools that make a project check
itself: they find where it is wrong or says the same thing twice, have a fix
proposed, and take the fix only on a check, never on anyone's word. The checks
come from the project itself -- two paths to one answer must agree, a law that
holds on a class of objects must hold on each of them -- so whoever proposes may
be wrong. Spatium is where it is built and measured first; GoodNet, a plugin
network kernel, is the second project it runs on. Where the idea came from is
[`docs/open-problem-rsc.md`](../docs/open-problem-rsc.md).

One rule above the rest: **structure says where to look; only a check gives a
verdict.**

## The parts

| part | what it does | here | state |
|---|---|---|---|
| map | who calls whom, who owns what | `rsc/spatium/graph.json`, the call map | built |
| laws | contracts read off behaviour: what holds, on which objects | `rsc/spatium/rsc_laws.py`, `laws.json` (637 laws, 25 spaces); each law it derives with its chain (`rsc_laws <size> <file>`) | built; consequence pruning to take from QuickSpec |
| trace invariants | what a stateful component always does across its calls: once per session, never after shutdown | a shim on the plugin table's calls | to build (Daikon's idea) |
| breaker | moves inputs, or sequences of calls, to where a contract breaks, and shrinks what it finds to the smallest case | the tuples `rsc_laws.cpp` evaluates on | inputs built; sequences and shrinking to take (RapidCheck) |
| siblings | implementations of one interface run on the same inputs or calls; the one that differs is the task | spaces by the laws; plugins of one table by a shared conformance suite | built for spaces; plugins to build |
| scheduler | the same interleaving of asynchronous work run twice, so an ordering bug can be checked | an event loop on a fake clock, interleaving by seed | to build; the one heavy part |
| structural queries | who may call what, which declared surface is empty, which size field goes unchecked | Joern or CodeQL | to take |
| sanitizers, fuzzing, mutation | the cheapest finders, and the strength of a test | ASan, TSan, libFuzzer, Mull | taken |
| history | which past fixes a finder would have caught; files that change together | `scripts/replay_fixes.py`; co-change from git | replay built; co-change to build |
| fixer | whoever writes the fix: a person or a general model run headless | `rsc/spatium/rsc_propose.py` | built |
| judge | takes a fix only on checks, and cannot be talked round | `rsc/spatium/rsc_judge.py` | built; rules grow from loopholes |
| task format | one row for every kind of task: the check, the witness, the neighbours that hold it, the verdict's command | `rsc/spatium/tasks.json` | laws' tasks built; one schema for all kinds to build |
| axioms | what is meant, kept apart from what is observed: each norm with its source, yes or no from a person | `laws-triaged.json`; a norms file | triage built; norms to build |

The breaker and the fixer are two players; the judge decides between them.
Every case the breaker finds, shrunk, becomes a regression test.

## What a project needs

A build in one command, tests addressable by name, a sanitizer configuration,
and fix commits that can be found in its history. GoodNet adds one more: its
bugs are about order (shutdown, completion after disconnect), and a check of
order needs a scheduler that runs the same interleaving twice.

## The loop

1. **Survey.** `rsc_laws.py discover` enumerates the laws and leaves them as
   data; it is dear and run when a branch is grown.
2. **Regression.** `rsc_laws.py check` evaluates that data again on every
   change in CI, a column per object.
3. **Task.** A finder writes a task: an object that breaks what its neighbours
   keep, a contract clause no test reaches, a meaning written twice.
4. **Proposal.** `rsc_propose.py --task N` hands it to the fixer in a work tree
   of its own, on its own branch, with nothing it can push.
5. **Verdict.** `rsc_judge.py BRANCH --task N` decides; a branch it accepts
   becomes a pull request.

## Kinds of task and their check

| kind | found by | the judge accepts when |
|---|---|---|
| a bug | siblings disagree (implementations of one interface on the same inputs), a breaker's case, a sanitizer | the law or clause now holds and nothing else got worse |
| a clause no test checks | a clause or law no test reaches, an attribute with no single owner | the new test passes on the code and fails on mutants the judge picks |
| a missing abstraction | a type used in many places and owned by few, one technique repeated, files that always change together | tests and laws hold and the score below is positive |

## What a person does

Four things, each a file with a format and a check, never a procedure in prose:
how a plugin is set up for the shared conformance suite; the table from a
contract's clauses to the tests that check them; the triage of what a finder
raised and was no defect (`laws-triaged.json`, each with its reason); and yes
or no on each norm harvested from past fixes, tests and contracts, each with
its source.

## Task format

One row of `tasks.json` per task: its kind, the object, the check it names
(a law, a clause, an invariant), the witness (an input or a sequence of
calls), the neighbours that hold it, and the command that gives the verdict.

## The judge's rules

In order: the size of the change, the tests no weaker (every test case of the
base kept by name, no fewer assertions), it builds, the tests pass, the task is
done against the base's snapshot, no law stops holding where the change
reached, the generated checks pass. Each way a branch slipped through becomes a
rule and a case it must refuse. Holes known and not yet closed: an assertion
weakened to a tolerance, an asynchronous test made to pass by a wait.

A change is scored as a stake: its price grows faster than its size, and it is
paid only for gain the judge measured -- a shorter description of the project,
new laws or clauses checked, a lower counted cost -- growing faster than the
gain, so a bold change that holds earns more than the small ones adding up to
it. A refused change loses its stake. A reformat gains nothing and pays its
size.

## How a finder earns its place

By the past: replayed on the trees before each of the project's fixes, it must
catch some of them. A finder that catches none is a note, however good it looks.

## What is measured

Good and bad alike.

- Over numbers a chain space collapses (9 765 625 paths at depth 10 reach 13 836
  states); a learned heuristic there needs 1.6 to 1.7 times the nodes of
  breadth-first search and solves 90.8 % against 100 %.
- Laws found on samples are mostly chance until a number's degree in length is
  part of its type and a breaker moves the inputs: 970 classes of spaces
  without them, 25 with.
- A fix that removes a defect shortens the description: the hyperbolic
  distance made exact took the laws from 834 to 637. That fix was proposed by a
  general model run headless and taken by the judge.
- An attempt is bounded by building, not by thinking: two builds from cold are
  about 25 of its 50 minutes.
- The law check would not have caught the library's past fixes: 0 of the 13 it
  could build, and its tolerance of 1e-6 cannot see a fix at 1e-10.
- The shape of the call map says little about where bugs were: AUC 0.65-0.72
  here and 0.57-0.58 on GoodNet against size, as the literature finds;
  history (churn, files changed together) is the stronger signal.
- Why bugs happened, in 113 past fixes of the two projects: a contract no code
  checked at a boundary, and one meaning written twice (one shutdown fault
  fixed three times, in three sibling plugins).
- A model of the law table (low rank) is better calibrated than the rules
  (log-loss 0.19 against 0.23) and no better at anything else: not at
  predicting a hidden space's laws, not at choosing which to check first (AUC
  0.966 against 0.969). The table holds about 2 kbit.
- On GoodNet, with sanitizers: 1523 tests, a leak at shutdown (a task posted
  before shutdown that never runs keeps its state), and the project's own
  sanitizer build no longer compiling with its compiler; the shared teardown
  suite runs on 4 of its 9 link plugins.

## Where things are

`rsc/core` is generic: a typed graph and the query for routes, and a tree search
guided by a learned prior and value (`tree_search.hpp`); it includes nothing of
the library. `rsc/spatium` binds it to this library: the graph (`gen_graph.py`),
the routes and their witness pairs (`rsc_find.py`), the laws, the judge and the
proposer, each described in the root `CLAUDE.md`. The prototype of 2026-09 left
the tree on 2026-10-20 (tag `rsc-prototype-final`); its parts, the contract of a
job and the design not built yet are in `docs/roadmap/rsc/README.md`.
