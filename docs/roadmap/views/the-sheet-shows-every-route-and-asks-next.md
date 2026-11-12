---
date: 2026-11-05
title: the page shows every route's chain, asks what comes next, and plots on demand
files: [web/src/view/answer.ts, web/src/view/cells.ts, web/src/page.ts, site/style.css, include/spatium/answer/respond.hpp, site/gallery/gallery.json]
open:
  - a Next question is asked as a new card; the one-line calculator with its history, where it would replace the line, is not there yet
  - the tabs are the calculator's; the canvas solve node does not reuse them yet
  - the steps setting on a card still hides the steps
---

The core sent every route's chain and what else a question can become, but the page drew only the taken route's steps, and a graph hung under every equation whether or not anyone wanted it. Now:
- Steps carry a tab per route that answered, with its cost. The route a person read last is shown first, else the one taken, and the section's gist follows the tab.
- Next is a list of buttons. Each asks its question as a new card.
- The graph is a section of its own, Plot, closed by default and drawn the first time it is opened.
- A section stays as the person left it. Only a click on its title is remembered (`spatium.open`, `spatium.route`), not the core's default opening it.
- The proof's cut is a step written as a sum, `40 = 37 + 2 + 1`, with what each term counts in its label, so it fits a phone.
- A card's text area no longer styles KaTeX's `\text{}`, which the same class name had turned into a monospace box.

Checked in headless Chromium over the DevTools protocol at 400 px, light and dark, on `2x^2 - 3x - 5 = 0`:
- three tabs, closed form first;
- the third tab shows four steps;
- Plot is not built until opened;
- a Next click adds a card answered by the interval method;
- after a reload, the route and the open Plot are kept, and nothing else is stored.
