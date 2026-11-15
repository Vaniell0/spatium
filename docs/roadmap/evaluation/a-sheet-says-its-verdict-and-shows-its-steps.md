---
date: 2026-11-05
title: a sheet says its verdict in words and shows its steps
files: [include/spatium/answer/respond.hpp, tests/test_answer_respond.cpp, web/src/view/answer.ts, site/style.css, site/gallery/gallery.json]
open:
  - only the taken route has steps; Vieta's formulas, the proof by intervals and bisection do not yet write their own chains
  - the sheet does not yet say what else the question can become (an inequality of it, its integral, a parameter in it)
  - a step's label is still its whole caption; a bold lead and a note beside it are not split
---

A sheet said how strongly its answer holds in one line mixing the strength with how it was found ("Closed form, checked by Vieta's formulas; no other real root", "4 roots proved, and no other on the box"), and a number's computation came folded under "18 steps", so a page showed the value with no steps in sight. The sheet now carries a `verdict`: its `strength` ("exact", "proved to 16 digits", "relative error below 10⁻¹⁵") and `notes` on how ("closed form", "checked by Vieta's formulas", "no other real root"); `label` is the two in one line, so there is one source. A number's steps are a section titled Steps, open, like every other chain, and a section says in a few words what is inside (`gist`: "18 steps", "closed form, 4 steps", "4 roots", "3 routes"). Input is the question as the library read it and what it asks for, as facts, instead of the rows again. The page draws the strength in the accent where the answer holds, each note beside it, and the gist on the right of a section's title. The check: eleven questions, each sheet's label equal to its verdict in one line, every Steps section open with a gist (`test_answer_respond`).
