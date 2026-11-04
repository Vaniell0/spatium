---
date: 2026-11-04
title: the page is a document of cells on its own package, under two policies
files: [web/src/doc.ts, web/src/core.ts, web/src/answers.ts, web/src/page.ts, web/src/view/cells.ts, web/src/view/graph.ts, web/src/view/answer.ts, web/test/doc.test.ts, rsc/spatium/respond.hpp, site/index.html, site/bundle.sh]
open:
  - cells sit in a column; on the canvas they are not placed freely or wired yet, and groups do not exist
  - a parameter's slider moves the graph that holds it, not the answers of the cells that read it
  - a cell offers digits and steps; exact or decimal and the method come when the core gives them
  - scenes have no cell yet: the TypeScript builder, the scene in WebAssembly and its WebGL2 view
  - the gallery is still the formulas it was, not the renders and the scenes as LaTeX
---

The page stands on a package of its own, `web/`, in TypeScript: the document a page holds is cells, each a line written in LaTeX or in the language, or the graph of a name another cell defines, each with its settings. The calculator and the canvas are one document under two policies: what a person may do with the cells, and what the page does to what is typed on its own -- split a pasted text of several lines into a cell per line, follow an answer with its graph, slide a free name. Cells the page adds are derived from the answers, not kept: a link carries only what was typed. The link is `[2, preset, cells...]`, deflate-raw, base64url; a link of the windows page is read as its cells.

A cell is read in the document the cells above it make: its question is every row up to it, a LaTeX row as MathJSON and a row in the language as `{"text": ...}`, which `respond_rows` now takes beside MathJSON (the mixed document answers as the text one). Writing a cell as a formula or as a line carries what is written over: to text as the core read it, to LaTeX by MathLive's ascii-math reader.

A graph is a cell of its own. The core's shader computes the value at every pixel from the view and the held parameters it is given, so panning (drag), zooming (wheel, two fingers) and sliding a parameter are new uniforms, never a new question. The page is marked beta.

Checked: six tests of the document under node's runner, the link also inflated by node's zlib as a second path; 74 assertions of `test_rsc_respond`; the page laid out by `site/bundle.sh` and shot in headless Chromium, the calculator and the canvas answered by the core in the browser, no sideways scroll at 1280 px.
