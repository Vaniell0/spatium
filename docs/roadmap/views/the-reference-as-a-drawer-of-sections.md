---
date: 2026-11-05
title: the reference as a drawer of sections, written by the library, and docs/language.md beside it
files: [examples/site_reference.cpp, examples/gallery.cpp, examples/CMakeLists.txt, include/spatium/answer/respond.hpp, site/wasm/check.mjs, site/data/reference.json, docs/language.md, web/src/page.ts, web/test/registry.test.ts, site/index.html, site/style.css]
open:
  - the canvas still shares the calculator's page; the drawer puts a canvas example at the end of the cards, not at the centre of a view
  - the logarithm in base 10 is proved to 15 digits, not exact: the core takes log(100) / log(10) in intervals
  - shapes are described in the reference and drawn by no page yet
---

The reference was a flat list beside the page. It is now a drawer over the page, which does not
move: 480 px on the right from 1024 px up, a sheet from the bar down below that. It has a filter,
an index of seven sections, and remembers its scroll for the session. Esc, the book button or the
scrim closes it.

Every section is written by `examples/site_reference.cpp` from the library's own tables, and the
same run writes `docs/language.md`, so the page and the docs cannot differ (`--check` is a ctest):

- **Writing**: what is typed in the formula editor, the formula it becomes (the library's
  `question_latex`), and the line of the language it reads as.
- **Questions**: seven kinds, now with an inequality and a system, each with the answer's label.
- **Functions**: the editor's notation first. `\ln`, `\lg`, `\log_b`, `\log_2`, `\tan` and
  `\sqrt[n]` have rows of their own. Each row gives the line of the language the notation reads as.
  Every MathJSON head the reader takes has a row here or is a form the Writing and Operators
  sections show. Otherwise the program fails.
- **Operators**, **Constants**: as before.
- **Language**: the canvas's text, a line per kind of statement. The program fails if a kind has no
  line. Each example is read by the library and shown with its rows. Shape operations follow, and
  then what the library says of a document it does not take, in its own words.
- **Answers**: the seven guarantees and what each promises. The table moved into
  `include/spatium/answer/respond.hpp`, so the gallery's groups and the reference read one list (the gallery's
  file is unchanged).

On the calculator an entry shows LaTeX first; on the canvas the language first. An example puts its
rows in the line, the old line going into the history.

Checks: `site/wasm/check.mjs` reads every Writing row and every hand-written notation row back
through the editor's parser and asks that it answers as the language line beside it does: 20 of 20.
Swapping `\lg(100)` for `\lg(1000)` fails it ("value 2.9999999999999996"). The host and the
WebAssembly core agree on 56 of 56 questions; 23 of 23 gallery and reference questions read back.
In headless Chromium at 1280 and 400 px, light and dark: the filter "base" leaves 3 rows, all in
Functions; the `\lg(100)` example lands in the line and is answered (≈ 2, proved to 15 digits); Esc
closes the drawer; there is no horizontal scroll.
