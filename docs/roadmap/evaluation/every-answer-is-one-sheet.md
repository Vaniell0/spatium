---
date: 2026-11-04
title: every answer is one sheet, and the gallery is read by how strongly its answers hold
files: [include/spatium/io/steps_text.hpp, web/src/view/cells.ts, include/spatium/answer/respond.hpp, include/spatium/answer/evaluate.hpp, web/src/view/answer.ts, web/src/view/dom.ts, web/src/page.ts, site/style.css, examples/gallery.cpp, site/gallery/gallery.json, site/data/reference.json, site/wasm/check.mjs, tests/test_answer_respond.cpp, tests/test_answer_evaluate.cpp]
open:
  - a route's own fields stay in the answer beside the sheet (the graph and the checks read them); they go when routes are a table
  - the closed form's steps are still written by hand (`to_steps`); the proof of roots shows its roots, not its cuts and Newton steps
  - a route's cost differs between libms by a few evaluations in a hundred thousand (a synthesis: 117953 on glibc, 117945 on musl); a cost a model learns from must not
  - the card asked in code (prism) and the scene card carry no sheet
---

Every answer now carries one sheet (`sheet`, `include/spatium/answer/respond.hpp`), whatever route gave it: what it is (`outcome`: a number, roots, a function, a parameter, refused), how strongly it holds (`guarantee`: exact, proved digits, every root proved and no other, sampled, estimated, computed, none) and the line that says so, the result in one line of LaTeX, and its sections in one order -- Input, Steps, Roots, Certificate, Check, Pieces, Properties, How it was found -- each a list of steps, of facts or a table, with whether it starts open. It is the sheet of the mock the owner chose, filled by the core.

The page draws a sheet and nothing else (`web/src/view/answer.ts`, 197 lines to 91): fourteen branches on a route's name became one, the scene card's, which has no sheet; the words under "How it was found" come from the core with the route (`route_caption`). A closed form's check now shows its second route: each root of the formula inside an interval proved to hold one (2x^2 − 3x − 5 = 0: −1 to 15 digits, 2.5 to 16), and where the line's far ends overflow a double, that it is undecided there.

The gallery is grouped by guarantee, strongest first, under titles the gallery program writes: exact (the arithmetic card, the exercise from a page of problems), every root proved and no other (2^x = 8, sin x = 0.5 on [0, 10]), proved digits (2x^2 − 3x − 5 = 0), computed (the complex pair, Vieta's formulas checked), certified on a grid (the six syntheses), not answered. A card keeps its routes and what they gave, not their counts.

The pairs: ten questions, one of each kind, give the outcome, guarantee and label measured, with their sections in the one order and each section holding something (`test_answer_respond`); WebAssembly and host agree on all 51 questions, sheets included, the cost column of a table compared as a candidate's cost is, not at all, since a route's count of evaluations differs with the libm's last digits.

The page was drawn by headless Chromium and driven over its DevTools protocol, which found what the type checks could not. A long question's result pushed its answer off a phone's screen, since it repeated the expression: the result is now the value alone, or a name and its value (the exercise's sheet reads 12). The on-screen keyboard's own keys -- delete, the arrows, enter, shift -- were blank and looked disabled: one style rule for the page's icons, `svg { fill: none }`, erased every glyph MathLive draws by fill. It is the icons' rule now, and the keyboard is a panel at the bottom edge with room made above it, where before a window below the row covered the answer. Typed keys reach the formula (7, +, 8, an arrow, a delete, read back from the field). On a phone the bar is two rows, the name and its tools, then the page's switch and its panels; a row written with no name no longer shows `answer =` beside it, and a power made of a chain of products is labelled power.
