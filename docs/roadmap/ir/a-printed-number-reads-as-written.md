---
date: 2026-10-30
title: a printed number reads as written, and a factor of exactly 1 is left out
files: [include/spatium/io/field_text.hpp, tests/test_field_symbols.cpp]
open:
  - nothing else is simplified: x + 0, 2 * 2 and p + p stay as they are, and a fitted coefficient near 1 stays a number
  - a piece of a piecewise answer is written in its interval's variable, (p - c) * scale, which reads poorly on a page
---

Two things the gallery data showed, both changes of spelling that change no value. A power of ten in LaTeX lost its padding and its sign: 1.8e-05 is now `1.8 \cdot 10^{-5}`, 1e-26 is `10^{-26}`, -2.5e+20 is `-2.5 \cdot 10^{20}`. And a factor of exactly 1 is left out of a product in both spellings, since x * 1 is x to the last bit for every x: the root of x^2 - p reads `sqrt(p)` instead of `1 * sqrt(p)`, and the prism's volume `0.16666666666666663 * p + 1.0000000000000002` instead of ending in `* 1`. A coefficient of 1.0000000000000002 is a number and stays one.
