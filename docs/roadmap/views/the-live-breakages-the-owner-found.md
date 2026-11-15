---
date: 2026-11-05
title: the breakages the owner found live: the keyboard, Ctrl+Z, a question offered next, infinity, where systems are
files: [web/src/page.ts, web/src/view/answer.ts, web/src/view/cells.ts, web/src/view/history.ts, site/style.css, site/wasm/check.mjs]
open:
  - the integral offered next for |x^2 - 4| - |x - 1| = 1 does not converge in the core: an integrand with kinks
  - the hint under the caret covers the line's error message
  - a graph loses its smoothing while zoomed, and joins across a pole; both wait for the graph work
---

Five things the owner found on the live page, each reproduced in headless Chromium first.

- **The keyboard.** Its window was 236 px tall, and the fifth row of keys (0, the point, +, the
  arrows, Enter) was cut off. A press there landed on nothing, the formula lost the focus, and no
  key typed again. The window now holds every row; no keycap is covered (each one is under the
  pointer at its centre). A press anywhere on the window gives the focus back to the formula. 1, +,
  2, ×, 3 pressed with answers arriving between them type 1+2×3 = 7.
- **Ctrl+Z.** It undoes what the Undo strip offers: a removed or cleared history, and now the line
  cleared by Esc twice. Inside a formula with something in it, the editor's own undo goes first.
- **A question offered next.** Next put its text in the language into the line, which takes
  formulas, and the editor's text conversion broke it ("dxon", "'At' is not in the language"). It
  now puts the row the library wrote. `check.mjs` reads every row a sheet offers back through the
  editor's parser: 24 of 24 ask what their text asks.
- **Infinity.** `inf`, as the reference names it, typed in the editor became three letters. It is
  now ∞, which the reader takes as infinity; `in` still gives ∈.
- **Where systems are.** A + button in the line adds a row below, as Shift+Enter does. Its label
  says that rows of equations are a system. x + y = 3 and x − y = 1 give x = 2, y = 1, exact.
