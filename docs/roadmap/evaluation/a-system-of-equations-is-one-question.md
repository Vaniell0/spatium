---
date: 2026-11-04
title: a system of equations is one question
files: [include/spatium/io/document.hpp, include/spatium/answer/closed_form.hpp, include/spatium/answer/evaluate.hpp, include/spatium/answer/query.hpp, include/spatium/answer/respond.hpp, include/spatium/algebra/linear_solve.hpp, include/spatium/core/epsilon.hpp, include/spatium/core/rational.hpp, tests/test_document.cpp, tests/test_answer_respond.cpp, examples/gallery.cpp, site/wasm/questions.txt]
open:
  - a system that is not linear is refused by name; its roots on a box, each proved by Newton's interval step in several unknowns, are the next route
  - at most four unknowns, the sizes solve_direct is taken at
  - the steps are the matrix and the solution, not the elimination itself
  - a system with a parameter is not synthesised: it is refused as reading more than its unknowns
  - complex coefficients: the language has no imaginary unit yet
---

A card of two rows, x + 2y = 5 above 3x - y = 1, was refused: each row above the last was read alone, and an equation in two unknowns has no root of its own. Rows that are all equations, reading as many unknowns as there are rows, are now one system, and a document writes it on one line, its unknowns and then an equation each: `x, y: x + 2y = 5, 3x - y = 1`. It is a root with the others beside it (`Statement::also`), so `find` gives each unknown, the line prints and reads back, through JSON too, and the page's rows come back an equation a row. `evaluate` solves it by a new route, a linear system: each equation's coefficients are read from its own values at zero and at each unit vector, and that it is linear is checked at twice each, at each pair together and at one point more, so a square or a product of two unknowns is refused by name; in fractions where every coefficient is one, by the library's own elimination (`solve_direct`, now with a threshold of zero for an exact arithmetic, which says so), else in doubles. x/2 + y/3 = 1, x - y = 1/4 is x = 13/10, y = 21/20 exactly. The check is a second route: each equation, as a field, at the solution. The pairs: the rows from the editor, the rows typed and the document's line give one sheet; sqrt(2) as a coefficient goes to doubles and agrees with 1/(√2 + 1) to 1e-15; a nonlinear and a dependent system refused (`test_answer_respond`); the line read, printed, through JSON, and three refusals (`test_document`); WebAssembly against the host 53 of 53, the gallery's new card read back from its rows.
