# Contributing to Spatium

Spatium is a C++23 library for arbitrary mathematical spaces, geometric
primitives and mesh operations.  The library is concept-driven: most types
are templated on a `Scalar` parameter and on space concepts, so new
contributions should follow the same pattern instead of hard-coding `double`
or Euclidean assumptions.

## Where to start

`docs/ROADMAP.md`'s Backlog is the honest list, including about what is
deliberately *not* going to happen and why. Entries marked `[want]` are
things that would be good to have and that nobody is working on; several
are described in enough detail to implement directly.

Two kinds of contribution are worth more than their size suggests:

- **A number that disagrees with one of ours.** Most of this repository's
  interesting history is plausible claims that did not survive being run.
  If a benchmark on your machine says something different, that is a
  finding, not a complaint.
- **Arguing with a ROADMAP entry.** Several changed after someone pushed
  back on the reasoning. "This is listed as impossible and here is why it
  is not" is a welcome issue.

## Building

```bash
nix develop                         # GCC 15, Vulkan, Eigen, Catch2, gbenchmark
cmake --preset release              # or `cmake --list-presets` for the full set
cmake --build --preset release
ctest --preset release --output-on-failure
```

If you do not use Nix you need GCC ≥ 15 (or Clang ≥ 19) and CMake ≥ 3.28
(modules build requires GCC ≥ 14 and CMake ≥ 3.28; the legacy header tree
also requires 3.28 because some downstream FetchContent dependencies pin it).

## Workflow

1. Branch off `main` with a short topic name: `feat/<topic>`, `fix/<topic>`,
   `docs/<topic>`.
2. Keep each commit atomic — one logical change per commit.  Tests must be
   green after every commit, not just at the tip.  CI rebases will rerun
   `ctest` per commit on Linux/GCC 15.
3. Commit messages follow a flat-prose style with a 50-character subject line
   and a 72-character wrapped body (no bullet points, no emoji trailer).
4. Merge with `--no-ff` so each topic stays as a recognisable bubble in the
   history.

## Coding conventions

- Namespace `spatium::` is the top level; sub-domains live in
  `spatium::geometry`, `spatium::mesh`, `spatium::physics::mechanics`, etc.
  `algebra/` is `inline namespace algebra` (Vec/Matrix/etc. stay reachable
  as bare `spatium::X` since they're used unqualified everywhere else in
  the tree) — see `docs/conventions.md` before adding a new domain
  namespace of your own.
- Public types are `PascalCase`, free functions are `snake_case`.
- Errors travel as `Result<T> = std::expected<T, Error>`; reserve exceptions
  for exceptional cases such as OOM or unsupported configurations.
- New space types should `static_assert(spatium::Surface<MySpace>);` (or the
  applicable concept) directly in their header so a regression collapses the
  build, not a runtime test.
- Operator overload convention (do not break it without a docs update):
  `|` = intersect / pipe, `&` = boolean intersect, `+` = union,
  `-` = difference.
- Add `std::formatter` for any new public type that has a natural textual
  representation; `std::println("{}", thing)` should "just work".
- Internal helpers go into `namespace detail { … }` and are documented as
  unstable.  Do not expect them to be picked up by downstream users.

## Tests

- Unit tests live under `tests/`, one Catch2 file per module
  (`test_<module>.cpp`).  New code without an accompanying test will not be
  merged.
- Regression tests for a bug must fail on the parent commit and pass after
  the fix in the same series.
- Numerical code that depends on `Real50`/`Real100` should include a
  precision sweep so we notice when a routine quietly silently relies on
  `double` rounding.

## Examples

Each example under `examples/` is a stand-alone reproducer for a feature
(geometry, fractal, geodesic, viewer scene).  When adding a new one:

- Provide `--help` and `--force`, and write outputs only when `--force` or
  the file does not yet exist.
- Wire it into `examples/CMakeLists.txt` and, if it is a CPU-only renderer
  worth shipping, into `flake.nix` (`installPhase` plus an `apps.<name>`
  entry).
- Avoid pulling in Vulkan unless the demo actually needs it.

## Reporting issues

The issue templates ask for what a maintainer would otherwise have to
come back and ask you for: compiler and version, the CMake invocation,
and the smallest reproducer. For a numerical or geometric bug, the
expected value *and where it comes from* — a closed form, a reference
implementation, a paper — is far more useful than a screenshot.

One thing that catches people out and is worth stating here: several
capabilities are absent rather than broken without their optional
dependency. The heat method, DEC and geodesics need `SPATIUM_EIGEN=ON`;
`Real50`/`Real100` need `SPATIUM_BOOST=ON`. A default build genuinely
does not have them, and that is a documented state rather than a bug.

## What review will ask

Not a gate, just the questions that come up every time, so you can answer
them in the PR and save a round trip:

- **Does the test fail without the fix?** A regression test that passes
  either way is not testing anything. The same standard applies to a new
  feature's test against the behaviour before it.
- **Did anything get measured, and what did the measurement say?**
  Including when it said the opposite of what was expected.
- **Is anything in the docs now untrue?** Three CI scripts check that
  documentation names files that exist, that `CLAUDE.md` covers every
  directory, and that the dependency graph matches the real includes —
  but none of them can tell that a sentence became wrong.
