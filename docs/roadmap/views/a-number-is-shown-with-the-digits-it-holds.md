---
date: 2026-11-05
title: a number is shown with the digits it holds
files: [include/spatium/core/interval.hpp, include/spatium/io/field_text.hpp, include/spatium/io/steps_text.hpp, tests/test_rsc_respond.cpp, tests/test_rsc_evaluate.cpp, tests/test_scene_symbols.cpp]
open:
  - long digits are not yet grouped in threes, and a fraction's decimal is not yet written beside it everywhere
  - the page's own digits setting, which would cap what the core shows, is not there yet
---

A page showed sin(pi/6)*2 as 0.500000000000000 \cdot 2 \approx 1.0000000000000, a root as -1.4142135623730951, a check of Vieta's formulas as -2.220446049250313 \cdot 10^{-16} = -0, and 2^x = 8 as x = 3.000000000000000. Proved digits now drop the zeros after their last that tells (`without_zeros`): 3, 0.5, 1 -- how many are proved stays in the answer's label; a step writes its operands the same way, so the step reads 0.5 \cdot 2 \approx 1. A double in LaTeX shows at most 15 significant digits and -0 as 0; its text, which reads back to the same bits, is unchanged. A check of Vieta's formulas that holds writes a value off by the last bits as the one it is near, with \approx: x_1 + x_2 \approx 0 = -b/a. The check: nine questions, their results and every step, held against a trailing run of zeros, sixteen digits or more, and -0 (`test_rsc_respond`); the step's operand is the value's own digits; the scene's LaTeX against its h written to 15 digits (`test_scene_symbols`); WebAssembly against the host 53 of 53.
