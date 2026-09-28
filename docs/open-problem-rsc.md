# Open problem: can a library find its own connections?

This is an unsolved problem, stated as plainly as we can, because we think
it is more interesting than anything in this repository that is finished.

## The setting

Spatium has two halves. **Concepts** say what a space is -- a metric, an
exponential map, a surface with a normal. **Pillars** are hand-written
implementations of the best-known mathematics for particular systems:
geodesics on a sphere, Kerr's metric, a certified time of first contact,
Christoffel symbols by forward-mode differentiation. Everything is templated
on the space and the scalar, so a pillar written for one system can in
principle run on another.

"In principle" is the problem. Between the pillars there is a space of
compositions nobody wrote: a pillar applied to a type it was not written
for, one operation derived from others, one scalar nested inside another.
Some of those compositions work. Some work and disagree with a pillar that
computes the same thing another way. Nobody knows which, because nobody
enumerates them.

## What made us think this is a real question

Two bugs were found in one afternoon, neither of them by looking for it.

1. A generic `distance` derived from a space's own `log` and metric,
   |log_p q|_g, was compared against each space's closed form. They agreed
   everywhere except at p = q, where the closed forms (acos, acosh) return
   1.5e-8 to 2.1e-8 and the derived form returns exactly 0.
2. Curvature was taken by running the Christoffel symbols -- already
   computed with `Dual<T>` -- on `Dual<Dual<T>>`, a composition nobody had
   written and which compiled unchanged. Exact Kerr, which is a vacuum
   solution, came out with a Ricci residual of 2.5e-5 on every coordinate
   plane. The cause was three layers down: the matrix inverse skipped rows
   whose entry compared equal to zero, and for a dual number that
   comparison reads the value alone, so an entry that was zero but had a
   derivative was dropped.

In both cases **two different paths to the same answer disagreed**, and the
disagreement was the bug report. Neither path was a test anybody wrote.

## The question

Treat the library as a graph: types and concepts as nodes, operations,
derivations, adapters and representation changes (`lower`, GPU emission,
`Real50`, `Dual`) as edges. Mathematics says that many pairs of paths
between the same two nodes must agree. Every such pair is a free test.

- **Can the library enumerate its own compositions and find the ones that
  work, the ones that disagree, and the ones that should exist by symmetry
  and do not?** Output: a connectivity matrix, a list of disagreements
  (bugs), a list of gaps.
- **Can it choose, per region of inputs, which path to use** -- a closed
  form or a derivation, `double` or fifty digits, analytic or tessellated --
  from measured cost and certified error rather than a rule someone wrote?
- **Can it shrink itself?** If two implementations always agree and one is
  slower, one of them is dead weight. The objective in its strongest form
  is description length: the size of the library plus the size of
  everything expressed through it.
- **Can it grow?** A chain of pillars found once to meet a specification
  could become an operation in its own right, so the next search is shorter.

## What exists

The module that does this is called RSC. It began as "a dispatcher that
picks the implementation", and turned out to be stranger than that.
Measured so far:

- Ground truth is generated, not labelled: candidates are run against a
  reference. This is the part that works best.
- Searches over chains of operations: over numbers, 9.7 million paths at
  depth 10 collapse to 13 836 distinct states; over meshes the branching
  stays at 4.66 of 6; over continuous-collision primitives, 819 chains are
  22-109 distinct algorithms per kind of query once compared by their
  answers bit for bit.
- A chain found on a sphere, replayed unchanged on a plane and a torus,
  still meets a relative specification, and a direct search there finds a
  chain of the same length. Structure transferred; constants did not.
- A learned heuristic on small spaces is *worse* than blind breadth-first
  search (1.6-1.7x the nodes).
- One learned choice ships in the library, as a generated header held by a
  test: a policy for certified continuous collision, 5-20% cheaper.

## What is hard, and open

- **The oracle.** A disagreement needs a tolerance. For exact paths it is
  zero; for a closed form against a derivation it depends on the input
  (the acos case above). Where does the tolerance come from -- a certified
  error bound carried by each path?
- **The measure.** "Agreed on these inputs" means little unless the inputs
  cover the space. Sampling by the space's own measure (Haar on a sphere or
  a group, area on a surface) gives "agreed on all but a set of estimated
  measure". Is that the right notion of checked?
- **The size.** Paths grow combinatorially with length. Which pairs are
  worth checking first? That is a value function over a graph, and it is
  the one place a learned model looks necessary.
- **Equivalence.** Two programs that compute the same thing should count
  once. Bit-exact fingerprints work for small searches; for rewriting an
  expression, equality saturation over the scene IR looks like the right
  tool. Nobody has tried it here.
- **Growth without rot.** Whatever the library learns has to enter it as
  generated, reviewable, regenerable artifacts held by tests, never as
  edits nobody can reproduce. How much can grow that way?
- **The limit.** It composes; it cannot invent a pillar. Where exactly does
  "a composition nobody wrote" end and "new mathematics" begin?

## Where to start

Small, real entry points, each useful on its own:

- Write one more pair of paths that must agree, and see whether it does.
  `tests/test_space_access.cpp` shows the shape.
- Generate the connectivity matrix: one translation unit per (scalar,
  space, algorithm), bodies instantiated, graded compiles / runs / passes
  axioms / matches a reference.
- Bring a space that does not fit the concept hierarchy. That is the most
  valuable thing anyone can send.

The detailed state -- every breakage in RSC as it stands, measured -- is in
[`findings-2026-09.md`](findings-2026-09.md), section 4.
