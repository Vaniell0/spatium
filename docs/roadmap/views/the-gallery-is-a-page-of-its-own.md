---
date: 2026-11-05
title: the gallery is a page of its own
files: [site/gallery/index.html, web/src/gallery.ts, web/src/page.ts, site/index.html, site/style.css]
open:
  - the hero tile, Riemann zeta in four views with a still of Z(t), waits for the zeta work
  - a question that is a document of several statements opens in the calculator's line, not yet on the canvas
---

The gallery was a panel squeezed beside the calculator. It is now its own page, `site/gallery/`, with the same bar; the Gallery tab is a link to it, and no core is loaded there. The tiles are grouped by how strongly each answer holds, in the order `gallery.json` gives its groups, two columns from 960 px and one below. A tile shows the question as the library writes it, the result and the verdict, and is a link that opens the question in the calculator's line (`../#w=`). The bar's Canvas and Reference links open the main page on the canvas (`#canvas`) or with the reference open (`#ref`). Checked in headless Chromium at 1280 and 400 px: six groups and sixteen tiles, no horizontal scroll at 400 px, and a tile opens `2^x = 8` in the calculator answered "proved, 1 root, no other on the box".
