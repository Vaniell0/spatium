---
date: 2026-11-04
title: zoomed past fp32, the curve comes in doubles
files: [include/spatium/answer/respond.hpp, web/src/view/graph.ts, web/src/answers.ts, web/src/view/cells.ts, web/src/page.ts, tests/test_answer_evaluate.cpp, tests/test_answer_respond.cpp]
open:
  - a synthesis of a Chebyshev piece, or of several pieces, has no graph: the plot request draws a value of the language, which has no T_k; it would draw pieces as it draws a net
  - the curve is a value a column, so a function that swings within one column is drawn by its centre alone
---

A graph zoomed in near 0.99999 was a staircase, and its ticks all read 0.99999. The shader computes the value at every pixel in fp32, whose step near 1 is 6e-8: on a view 1e-4 wide around f = 1 - x^2/100000, 64 columns hold fewer than 8 of its values. Past what fp32 tells apart -- a view with fewer than four of its steps a pixel, in x or in y -- the graph asks the core for the curve in doubles, the value at the centre of each column (`plot_request`'s `curve`, the sliders' values put in the program's slots by name), and draws it on the 2D canvas, marked "computed in double"; zoomed out, the shader draws again. A tick's number has as many digits as the step between ticks needs (0.99999748, not 0.99999). A synthesis of one piece shows that piece's formula -- for the error function, its Chebyshev sum -- where it showed "1 pieces"; more pieces are "piecewise, k pieces". The pairs: the fp32 grid and the double curve on the same view, fewer than 8 distinct values against 64, each of the 64 the expression computed directly, and a slider's value reaching the program (`test_answer_evaluate`); the erf synthesis's result its sum (`test_answer_respond`); the page shot zoomed in.
