---
date: 2026-11-13
title: A judge of whether a branch does a task of the law check, and nothing gets worse
files:
  - rsc/spatium/rsc_judge.py
  - rsc/spatium/test_rsc_judge.py
  - rsc/spatium/rsc_laws.cpp
  - rsc/spatium/rsc_laws.py
  - .github/workflows/ci.yml
open:
  - one profile, this library's: the build, the tests, the law check and the generated checks are spatium's; which parts are a core any project shares and which a project's own is not drawn yet
  - a task of tasks.json is one of two kinds, a law a space breaks and an operation a tangent lacks; the other finders (the call map, the matrix's red cells) give none yet
  - nothing proposes yet: a branch is handed in by whoever wrote it
---

`rsc_judge.py BRANCH --task N` takes a task of `rsc/spatium/tasks.json` and a branch that claims to do it, and decides in order, stopping at the first step that fails: the size of the change against the base (files touched, lines added and removed, S = touched / (net growth + 1), reported, never a reason); that the tests are no weaker -- every TEST_CASE of the base still there by name, and no fewer test cases or assertions under `tests/`; that the test binary builds; that every test case but the connectivity matrix's passes; that the task is done -- for a law, `rsc_laws.py check` on the task's space reports it as now holding and no law there stops holding, for an operation, `rsc_laws.py missing` (new: the arithmetic a space's tangent lacks, read by `requires` with no enumeration, 40 s) no longer lists it; that no law stops holding on a space whose headers the change touched; and that the generated checks pass. It works in a work tree of its own with a build directory kept between verdicts, so the checkout it runs from is untouched, and prints every step with what it saw and how long it took.

Two branches held it. One gave `Matrix` a unary minus, the task the law check had raised for SPD2-affine's tangent: accepted, every step passed -- 2 files, +21 lines, S = 0.09; 1556 test cases, 5698 assertions against 1555 and 5696; built in 879 s the first time; 1501 test cases passed in 144 s; the minus present; 834 laws on the 16 spaces whose headers include the matrix held as before. The other made the same change and removed the test case "Matrix identity": refused at the second step. The counts alone had not refused it, since the new test case kept their total, which is why the names of the base's test cases are compared and not only their number. `test_rsc_judge.py` holds the rules on text: S, the counts, a test case swapped for another, an assertion dropped, what the law check reports, a task done or undone, a regression.
