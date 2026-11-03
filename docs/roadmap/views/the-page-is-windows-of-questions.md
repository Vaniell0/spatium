---
date: 2026-11-03
title: the page is windows of questions, asked as rows of a LaTeX editor, with their graphs
files: [site/index.html, site/style.css, site/app.js, site/core-worker.js, site/wasm/build.sh, rsc/spatium/respond.hpp, tests/test_rsc_respond.cpp]
open:
  - the workspace is not in the link yet (#w=); #solve?q= still opens one question as text
  - the page is not published on GitHub Pages yet
  - windows sit in an auto grid: they are not dragged or resized, and a phone has no sheets with a switch
  - a graph is drawn with the other parameters at the middles of their intervals; there are no sliders for them
  - a synthesised answer of several pieces, or of a Chebyshev series, has no graph, since the shader draws one expression
  - each graph has its own WebGL2 context; a browser keeps about sixteen alive at once
  - the old page-side text preview (spatium_preview) is still in the core, though the page no longer calls it
  - the answer is shown in double precision; the exact roots the rational arithmetic gives are not on the page
---

The page is a workspace of windows in an auto grid, opened from the bar: a question, the reference, the gallery. Tabs and the Focus button are gone. It has one look, Instrument: warm graphite with one amber accent, Geist and Geist Mono, no borders on panels, depth from a lighter surface, a rail through the chain of steps. Both themes are designed, and the dark one is not blue. The scale is the site's: type 14 / 17 / 21 / 27, step 8, controls 40, radius 8 for panels and 6 for controls.

A question is written as rows of a LaTeX editor (MathLive 0.111.1, vendored in `site/vendor/mathlive/`; its fonts are the KaTeX fonts the page already serves, byte for byte, so they are not copied). Each row goes through the editor's parser, `latex-syntax` from the Compute Engine 0.151.0 (`site/vendor/latex-syntax/`, 470 KB against 4 MB for the whole engine), to MathJSON, and the rows go to `spatium_rows`. Enter opens a row below and Backspace in an empty row closes it. A window can switch to text in the language instead, which `spatium_respond` answers; a question from the gallery or an example from the reference opens as text, answered live. Under the rows, the text the core read them as is shown: the second view of the same question.

Answers are live. The core runs in a Web Worker (`core-worker.js`; the WebAssembly build now names `worker` among its environments), each window keeps only its newest question, and a question made stale while the worker is on it is given 400 ms, then the worker is ended and a new one started: that is how a synthesis typed past is cancelled. A refusal names the row it came from (`row`), and the row is marked.

The keyboard is MathLive's own, in a window that opens beside the row being written (below it when there is room, above when not, along the bottom on a phone).

An answer carries the graph a page draws of it (`graph`, a request for `plot_request`, written by the core): an equation solved in closed form as its left side less its right along the unknown, over its real roots and a margin (around the vertex when there are none); a synthesised answer of one expression along its parameter; a value of one parameter along it. `plot_request` also gives the heights the graph reaches (`range`, 256 samples in the shader's fp32), where the y axis goes. The page compiles the core's fragment shader in WebGL2 and draws axes over it; a window's graph expands it to the width of the workspace. A new name on the left that the right side does not read is a definition: `f = a^2 + 1` under `a in [0.5, 2]` is a value of a, drawn along a, not an equation in f.

Checked in a headless browser: the page at 1280 and 400 px in both themes, with no horizontal scroll at 400 (formulas scroll in their own box); the rows of `2x^2-3x-5=0` give the closed form with its six steps, and `p \in [0.5,4]` over `x^2-p=0,\ x\in[0,10]` gives √p certified at 10⁻¹⁵ on 2001 points, computed in the browser.
