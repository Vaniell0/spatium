---
date: 2026-10-08
title: The consumer rule, counted
stage: 1
files: [scripts/gen_consumers.py, docs/consumers-baseline.json, .github/workflows/ci.yml, scripts/ci_changes.py]
open:
  - the match is by identifier, so a name shared with another thing reads as used and a type used only through `auto` or a helper reads as unused; 170 names are on the baseline and some of them are used
  - 154 names are used only by tests; they are reported, not held, because holding them would block an interface whose caller lands in the next change, and making them a ratchet is for after the existing ones are cleared
  - the baseline lists the names that have no use today; it does not say which of them should be deleted and which given a caller, and nobody has gone through the list
  - the match does not look at the header's own other functions: a name used only by a sibling in the same header reads as unused
---

The rule that a module with no real caller counts as no code and no tests is in CLAUDE.md now, and until today nothing could say how far the library was from it. `scripts/gen_consumers.py` reads the names that `gen_capabilities.py` lists as declared by each header (not the ones in headers that call themselves internal) and asks whether anything else in include/, src/, tests/, examples/, benchmarks/ or rsc/ mentions them. Of 920 names, 170 are mentioned nowhere else, 154 only by tests, the rest by something that is not a test; physics has the most names with no use (45), then io (27), render (24), algebra (23) and geometry (22).

It is a ratchet. `docs/consumers-baseline.json` holds the 170 by name, one per line so that it merges, and `--check` (run by the capability job in CI) fails with the name of any that is added or that loses its last use; removing one from the baseline, which `--update-baseline` does for the ones that were fixed, is a visible change. That makes the rule something a contributor, and an agent writing a pair or a fix, is told on the pull request rather than something a reviewer has to remember.
