---
date: 2026-11-04
title: a card is one question
files: [web/src/doc.ts, web/src/answers.ts, web/src/view/cells.ts, web/src/view/answer.ts, web/src/page.ts, site/style.css, include/spatium/answer/respond.hpp, tests/test_answer_respond.cpp, web/test/doc.test.ts, web/test/notation.test.ts]
open:
  - a card's digits and its route are not settings yet: the sheet writes its result in one form, and the engine's routes are not chosen from the page
  - a system of equations is written as one card but not solved: the core reads the rows above the last as a document, not as equations to hold together
  - the underline marks the row the core could not read, not the place in it
  - in fractions a quadratic's roots are exact but the check by Vieta's formulas is still computed in double (x1 + x2 = -2.2e-16)
---

The calculator was a column of cells, each a row with a Formula/Text switch, the name it defined in a gutter and the settings in a bar. Each cell is now a card: the rows of a formula -- an equation, a system, names and what is asked of them, a new row with Enter, a new card with Enter on an empty last row -- or a text in the language on the canvas, with the answer under it and the graph it names inside it. The calculator asks each card alone, so a card that does not read refuses only itself; the canvas still reads every row above. The gear opens the card's settings: numbers decimal or exact, steps shown or hidden. Exact is the engine in fractions: `respond_rows` takes an `engine` as `evaluate_request` does (`engine_of`, one reader for both), and x^2 - 2 = 0 is answered x = ±√2 where the double gives ±1.414213562373095. A row the core could not read is underlined, the word before the caret is told from the reference under the formula (sin(x): sine, ...), and a pasted aligned or cases block is rows of one card. A number's steps come folded, titled by their count. The Formula/Text switch, the completion list, the gutter and fifteen unused style rules are gone. The pairs: the rows in fractions against `evaluate_request` on the same document, and an arithmetic the engine lacks refused (`test_answer_respond`); a link of the old one-row cells and of the windows page read back as cards (`doc.test.ts`); the word before the caret (`notation.test.ts`); WebAssembly against the host, 52 of 52.
