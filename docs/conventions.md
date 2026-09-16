# Conventions

Spatium grew domain by domain, each one solving whatever problem was in
front of it at the time. That's fine for getting things working, but left
several conventions that were never written down — different domains
independently arrived at different answers to the same question (which
namespace, which error-handling mechanism, when to subdivide a directory).
This document names the actual rule where one already exists in practice,
and states it as a requirement going forward. It is not retroactive
amnesty for every existing file — known violations are called out
explicitly below, with the fix tracked separately.

## Namespace

Every subdirectory of `include/spatium/` is a namespace: `geometry/` is
`spatium::geometry::`, `mesh/` is `spatium::mesh::`, and so on, matching
`io/`, `spatial/`, `render/`, `physics/`, `viewer/`. Bare `spatium::` is
reserved for the small set of umbrella types that live directly at
`include/spatium/` and don't belong to any one domain (`Point`,
`Morphism`).

`algebra/` is the one deliberate exception to "qualify to see it": it
uses `inline namespace algebra` rather than a plain one, because `Vec`,
`Matrix`, `Complex`, `Quaternion`, and `Dual` are load-bearing building
blocks for every other domain, referenced unqualified from inside
`geometry::`, `mesh::`, `physics::`, `spaces::`, and `io::` throughout the
tree. `inline namespace` makes `spatium::algebra::Vec` and `spatium::Vec`
the same entity — `algebra/` gets the same discoverable, named home every
other domain has, without forcing a qualification change on every
existing consumer. This is not extended to other domains: `Triangle`,
`Mesh`, and similar domain-specific types still require full
`spatium::geometry::`/`spatium::mesh::` qualification, because unlike
`Vec` they aren't used as raw building blocks everywhere.

**Known violation (being fixed by this change):** most of `algebra/`
(`vector.hpp`, `matrix.hpp`, `quaternion.hpp`, `complex.hpp`,
`calculus.hpp`, `dual.hpp`, `functions.hpp`, `ode.hpp`,
`linear_solve.hpp`, `polynomial.hpp`, `vec_expr.hpp`,
`eigen_interop.hpp`, `format.hpp`) currently sits directly in bare
`spatium::` with no `algebra::` qualification available at all, while
`algebra/groups/`, `algebra/concepts.hpp`, and `algebra/verify.hpp`
already use plain (non-inline) `spatium::algebra::`. Target: the first
group moves into `inline namespace algebra`, so it's addressable both
ways. The second group keeps its plain `namespace algebra {}` block and
is documented as `spatium::algebra::`-qualified — `Group`, `SO3`, `SE3`
are concepts/types you opt into, not ambient building blocks, so that's
the form to use even though C++'s inline-namespace rule (inline-ness is a
property of the namespace name as a whole, once established anywhere)
means the bare form incidentally compiles too in any translation unit
that also includes the inline group.

## Directory subdivision

A domain gets subdivided into subdirectories once either condition holds:

- it contains two or more clearly separable thematic clusters (different
  audiences, different dependency footprints, or different maturity), or
- it exceeds roughly 15 files.

Subdivision is not automatic once a domain is merely large — `geometry/`
(21 files) and `mesh/` (15 files) are cohesive single-purpose domains and
stay flat. It fires when there's a real seam to name, the way
`algebra/groups/` (SO(3)/SE(3), Lie-group specific) and
`physics/mechanics/` + `physics/relativity/` already separate from each
other.

**Known violation:** `physics/` mixes `mechanics/` and `relativity/`
(core physics) with five flat top-level files that are visualization
support, not physics — `atom_model.hpp`, `atom_palette.hpp`,
`atom_svg.hpp`, `bohr_model.hpp`, `orbital.hpp`. Target: these move to
`physics/atomic/`. `physics/elements.hpp` stays at the top level — it's
the data backing `atomic/`'s models, and is separately called out below
as the project's one compiled-translation-unit exception.

## Error handling: `Result<T>` vs `assert`

`Result<T>` (`std::expected<T, Error>`) is required at any function that
is a domain's boundary with arbitrary, untrusted input — geometric
construction and queries (`geometry/`), file parsing (`io/`), and any
`mesh/` entry point that builds a mesh from raw data. These are places
where a real, expected, recoverable failure mode exists (degenerate
triangle, malformed file, parallel lines) and the caller needs to be able
to observe and handle it.

`assert()` is for internal invariants over state the domain has already
validated — a physics integrator stepping a state it constructed itself,
an internal `mesh/` operation working on a mesh that's already known
well-formed, `spaces/` manifold math operating on coordinates already
inside the valid domain. These aren't recoverable-by-the-caller failure
modes; they're bugs if they trigger, and `assert` is the right tool.

This isn't a new idea — `docs/architecture.md`'s "three principled
exceptions" section already documents `physics/mechanics/` not using
`Result<T>` as a known, deliberate gap rather than a silent
inconsistency. What was missing is that `spaces/` and `mesh/` follow the
exact same assert-only convention (54 and 6 `assert` sites respectively,
zero `Result<T>`) while `architecture.md` counts them as part of the
uniform "spine" — which is accurate for their concept/API design, not for
their error-handling. That list is corrected to name them explicitly (see
`architecture.md`'s exceptions section).

This does **not** mean retrofitting every `assert` in `spaces/`/`mesh/`
into `Result<T>` — most of them are genuinely internal invariants and
converting them would be exactly the kind of validation-for-scenarios-
that-can't-happen this project should avoid. It means the convention is
now named, so a new fallible *boundary* function in either domain has an
unambiguous answer: `Result<T>`.

## One name serving two contracts

Named 2026-09-16, at the second instance. Once is a case; twice is a
pattern worth catching before the third.

**The class: one operation's name covers two contracts that differ in what
the caller is promising, and code written against the wrong one works —
quietly doing something else.** It is not a bug that shows up as a wrong
answer; it shows up as a capability that mysteriously cannot be built.

| the name | contract A | contract B | how it surfaced |
|---|---|---|---|
| `offset(base, thickness)` | a closed surface offset into a closed surface | a shell over a surface *with a rim*, where the rim needs a rule | the rule "the base must be closed" was true for A and impossible for B, and the demo did B |
| `.moving(f)` | a **placement**: this whole object goes there | a **deformation**: every vertex maps somewhere | instancing needs A and can never be built on B; `cook()` measured 19 800 of 19 801 objects blocked, all of them meaning A |

Both were resolved the same way and it is the resolution worth
remembering: **not by forbidding the second meaning, but by naming the
two operations separately.** `offset` / `offset_shell` was the first. The
motion slot is the second, and the distinction is load-bearing: a
placement can be read, shared and turned into one geometry plus N
transforms; a point map cannot, because instances that deform differently
have nothing to share.

The tell, in both cases: a capability that "should obviously work" keeps
failing to be implementable, and the reason is never the missing feature —
it is that the name already promised something weaker. When that happens,
ask what two things the name is covering before adding the feature.

## A flat array addressed by index, not a tree of pointers

Named 2026-09-15, at the fourth independent application. The first three
were each decided on their own merits; naming it here is so the fifth
decision is quicker rather than re-derived.

| where | what it replaced | what the flatness bought |
|---|---|---|
| `Trace` | a tree of `std::function` closures | walkable, dumpable, re-interpretable by a second pass written later |
| `Field` | a linked expression tree | evaluation in one linear pass, children already computed — the shape a GPU runs |
| `Cooked` | operations resolved on demand | one expansion of operations into objects, done once |
| `Runtime` | state owned by each node | the mutable part is small and contiguous, which is what a frame loop wants |

The recurring reason is not cache friendliness, though that comes along.
It is that **an index is a name that survives**: it can be stored, compared,
serialized, handed to another pass, and used by code that does not own the
structure — none of which a pointer into a tree can do. `type_index` does
not survive serialization and neither does a closure; an index does.

The two obligations that come with it, both learned the hard way:

- **Append-only.** If the array is in dependency order, inserting in the
  middle breaks that order silently — nothing renumbers, and the first
  read gets an uncomputed slot. Check the ordering at the append, per
  element, where it can break; not afterwards over the whole array, which
  is both more expensive and a weaker statement.
- **Say who owns the array.** "Append, never mutate" means one thing when
  each value owns its array and a completely different thing when many
  values index a shared one. `Field` owns its pool, so its `+=` is
  ordinary value mutation and there is no aliasing question at all —
  writing the shared-pool guarantee into its docs would have promised
  something it does not do. See `io/field.hpp`.

And one thing flatness does not decide: **which array.** An index only
means something inside one space of indices, so operations and objects
need separate spaces when one operation can produce many objects — a
`Scatter` node is one operation and 600 sprinkles, and there is no node
index that names the 37th.

## Failure to signal, wearing the costume of an answer

Named 2026-09-15, after finding the third instance and realising the
first two had been treated as unrelated bugs.

The class: **a value that looks like a valid answer, in the place where
the algorithm should have said "I cannot answer".** Not a wrong number —
a wrong number gets noticed. A *plausible* one.

Instances found so far, all in one week:

| what came back | what it meant | why it passed |
|---|---|---|
| a `NaN` root through a filter phrased as "skip if `t < 0`" | no root exists | every comparison against `NaN` is false, so a filter that lists what to *skip* skips nothing |
| a test passing under `NDEBUG` | the assertion never ran | a check that does not execute reads as a check that succeeded |
| `miss = 0.0` from `ray_quadric_proximity` | the model does not apply here | the imaginary part is a literal zero, so `abs()` of it is a clean number |

The third has a mechanism worth stating exactly, because it is the same
mechanism as the first. `Complex<T>`'s one-argument constructor sets
`im` to a literal `T{0}`; every real-root path in all three polynomial
solvers uses it; and `is_real()` tests only `abs(im) <= eps`. So when a
degenerate leading coefficient turns `re` into `NaN`, `im` stays exactly
zero and `is_real()` answers **true** — correctly, for the question it
was asked ("is the imaginary part negligible"), and disastrously for the
question the caller meant ("is this a usable real root"). One mechanism,
three symptoms, not three bugs.

What to do with the name: when a caller reads a field as a signal, ask
what that field holds when the algorithm had nothing to say. If the
answer is "whatever it was initialised to", the signal is not a signal.
Phrase filters as what to **keep**, not what to skip; refuse at the point
where the model stops applying rather than patching the value it
produced; and prefer a type that can say "no answer" over a sentinel that
has to be recognised.

## A check needs the thing that runs it, first

Added 2026-09-15, after getting the order right by accident.

An `assert` guarding an invariant is worth nothing until something
actually executes it. Every job in this repository's CI built Release
until `Debug build (asserts actually run)` was added, so for the whole
life of the project an `assert` in library code was decoration in CI —
and a misuse of Eigen that Eigen itself asserts on rode through green for
eleven days because of it. Worse than being missed: when someone later
saw the test pass in Release, `ROADMAP.md` gained a confident wrong
explanation ("the old failure note went stale"). A check that never runs
does not read as absent. It reads as passing.

The same shape recurs at every level:

| the mechanism | is decoration without |
|---|---|
| `assert` | a build that keeps `NDEBUG` off |
| poisoning unused slots in Debug | a Debug test run |
| a `[!shouldfail]` pin on a known bug | a suite that reports it as expected-failed rather than skipping it |
| an `is_structural()` flag | something that reports the *consequence*, not just the fact |

So the rule: **when adding a mechanism that makes a class of error
visible, land what executes it first, or in the same change.** Not after.
The infrastructure is the part that can be silently missing, because a
mechanism with nothing running it looks exactly like a mechanism that
passes.

### `reserve(size() + k)` inside an append is an anti-pattern

Added 2026-09-15, after it cost a factor of 1522.

`reserve` carries two different contracts in most people's heads: "make
room for exactly this much" and "make room for at least this much, and
keep growing geometrically". `std::vector::reserve` is the first. Asking
for exactly `size() + k` therefore sets capacity to exactly that, and the
*next* append reallocates — so a hand-written append that reserves what
it needs on every call reallocates on every call, and a linear
accumulation becomes quadratic.

This is not a subtle case. `push_back` gets it right by doubling
internally; the moment you write your own append and reach for `reserve`
to be tidy, you have opted out of that and have to restore it by hand:

```cpp
const std::size_t need = size() + incoming + 1;
if (need > v.capacity()) v.reserve(std::max(need, v.capacity() * 2));
```

Measured in `io/field.hpp`: building a 19 800-term field took **12 514 ms**
with `reserve(size() + k)` and **7.95 ms** with the guarded geometric
form. Same algorithm, same asserts, one line apart.

The general rule: **reserving exactly what you need is only right when
you will not append again.** Once before filling a fresh container, yes.
Inside the append itself, never.

### A baseline is a number *and* the configuration that produced it

Added 2026-09-15, after a false alarm that cost a worktree to settle.

`donut_demo --photo` is compared by hash to prove a change is
geometrically neutral, and that works — it is how the `offset` /
`offset_shell` split was shown not to touch a pixel. But a hash carried
in someone's head, or written down beside a change, is not yet a
baseline. Reconfiguring a build directory — modules on, Eigen on, a
different optimization level — moves the hash with no change to a single
line of geometry, and the mismatch then looks exactly like a regression.

A recorded hash whose configuration is not recorded with it is not
"stale". It is **not a value at all**, and keeping it leaves a false
explanation waiting for the next reader. Either record the configuration
beside it, or do what settles it in minutes: build the baseline commit in
the *same* configuration — a `git worktree` at `main` costs one command —
and compare the two renders you just produced, rather than one you
produced against one you remember.

Same family as the row above it: the comparison is the mechanism, and the
configuration is the thing that has to run it.

## `measure()` / `area()` / `length()` / `volume()`

Already stated in `CLAUDE.md`: `measure()` is the dimension-generic name
backing the `Measurable` concept; `area()`/`length()`/`volume()` are
convenience aliases that must forward to `measure()`, never reimplement
the formula. This is a correctness rule, not a style preference — two
implementations of the same formula drift.

**Known violations:** `geometry/triangle.hpp` and `geometry/circle.hpp`
(the `Disk` case) both invert the rule — `area()` holds the real
per-dimension formula and `measure()` is the one-line alias.
`geometry/box.hpp` is the compliant reference implementation. Target:
`triangle.hpp` and `circle.hpp` swap which function holds the body, to
match `box.hpp`.

## Documentation comments

Every public function and class gets at least a one-line comment stating
its purpose. This is a floor, not a target — dense derivations like
`physics/relativity/kerr.hpp`'s are a bonus where the math genuinely
warrants it, not something every file needs to match. What isn't
acceptable is the current gradient, where comment density tracks how
novel or difficult a file felt to write rather than any judgment about
what a reader needs.
