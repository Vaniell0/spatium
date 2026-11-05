---
date: 2026-11-05
title: the calculator is one line and the history of what was asked
files: [web/src/calc.ts, web/src/view/history.ts, web/src/view/cells.ts, web/src/doc.ts, web/src/page.ts, web/src/registry.ts, web/test/calc.test.ts, site/index.html, site/style.css]
open:
  - first visit shows an empty line; three gallery questions as examples under it are not there yet
  - Up arrow in the first row, Ctrl+K, Esc Esc and Ctrl+Z for the undo strip are not wired
  - the busy rule under the line and keeping the previous answer dimmed while a new one is asked are not there yet
  - a calculator link is the line as a document (#w=); the shorter #q= is not written yet
  - exact numbers for roots in closed form (x = ±√2) came from a card's Numbers setting, which the calculator no longer has; the core should show both forms on its own
---

The calculator was a stack of cards with an adder, a gear per card and settings no one needed. It is now one line, plus the history of what was asked.
- Enter puts the line's question on top of the history, shown open so the result does not vanish, and empties the line. Shift+Enter writes a row below, for a system. The x on the line does the same as Enter.
- A history entry is the question with its result and verdict in one line, kept in this browser (`spatium.calc`, at most 50, newest first). Asked again, a question moves up instead of being kept twice, compared without spacing. Opening an entry asks the core for its whole sheet. Up puts the entry back in the line, committing the line first. x removes it, and Clear history empties the list, each with an Undo strip under the bar for five seconds.
- A question from the gallery, the reference or a sheet's Next goes into the line, and what was there goes into the history.

The canvas keeps its document of cards, apart from the calculator's state; the bar switches between the two.

Checked in headless Chromium over the DevTools protocol at 400 px:
- Enter empties the line and opens the entry, "proved, closed form, checked by Vieta's formulas";
- after a reload the history is kept;
- x removes the entry, and Undo brings it back;
- Up puts the question back in the line, answered again.

The history's rules have a node test (`calc.test.ts`).
