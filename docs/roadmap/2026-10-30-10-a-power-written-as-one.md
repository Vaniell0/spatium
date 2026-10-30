---
date: 2026-10-30
title: a product of an op with itself is written as a power
files: [include/spatium/io/field_text.hpp, tests/test_field_symbols.cpp, tests/test_scene_symbols.cpp, site/gallery/gallery.json]
resolves: [2026-10-30-09#2]
---

`x^2` is read into the product of x with itself, and the pool merges the two factors into one op read twice, so the printer can tell a square from any other product exactly: `to_text` writes `x^2` and `to_latex` `x^{2}`, and a third factor of the same op on the right makes `x^3`. A square of a square stays `(x^2)^2`, since `x^4` would be read back as x * x * x * x, another order of multiplication; a negative base is grouped, `(-2)^2`; in LaTeX a fraction or a base with a superscript of its own is grouped, `\left(x^{1}\right)^{2}`. The round trip of `parse` against `to_text` on 400 random fields holds as before. The gallery's arithmetic card reads `2 \cdot \left(3 + 4\right)^{2} - \frac{1}{4} = 97.75` and the asymptotic answer `2 * (n / log(n)^2)`.
