---
date: 2026-11-02
title: abs in the field IR, and |x| on the page
files: [include/spatium/io/field.hpp, include/spatium/io/field_pod.hpp, include/spatium/io/field_glsl.hpp, include/spatium/io/field_text.hpp, include/spatium/io/field_parse.hpp, site/data/reference.json, site/wasm/questions.txt, site/wasm/plot_check.html, tests/test_field_parse.cpp, tests/test_field_pod.cpp]
open:
  - the reader takes abs(x), not |x|; bars belong to the interactive input that writes LaTeX
  - a metric field refuses Abs (its derivative jumps at zero); the synthesis does not search over it
---

A calculator asked for the size of a number and the language had none: `max(x, -x)` computed it and printed as `\max(x, -x)`. `Op::Abs` is appended after `Symbol`, since an op's number is its code in lowered data. It is total and one instruction on every target. The field's interpreter, the POD interpreter and the GLSL emitter carry it; it prints as `abs(x)` in text and `\left|x\right|` in LaTeX; the reader's table has it, so the page's reference lists it with an example the library answers.

Two paths hold it. The test that interprets every op bit for bit includes it, and a WebGL2 graph of `abs(x - q) - 1` drawn in the browser matches the core's fp32 twin with zero error at every pixel. The WebAssembly core and the host agree on 33 of 33 questions, two new ones among them: `abs(2 - 7) + abs(-1/4)` is 5.25, and `solve abs(x - 1) - p = 0 for x in [1, 10]` is synthesised as `p + 1`.
