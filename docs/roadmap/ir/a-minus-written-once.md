---
date: 2026-11-03
title: a minus is written once, in text, on the page and in the steps
files: [include/spatium/io/field_text.hpp, include/spatium/io/steps_text.hpp, tests/test_field_parse.cpp, tests/test_solution_steps.cpp]
open:
  - a negative number after a minus keeps (-1) * 3 in text, since -3 would be read back as one number
  - u - (-2) and u - (-v) are written as typed, not as u + 2 and u + v, which are other ops
---

A minus is now written once, wherever it appears in the output.

**Text.** `-v` was written `(-1) * v`, and `2 * -x` was written `2 * ((-1) * x)`. Now they are `-v` and `2 * (-x)`, as typed, and `parse` reads them back to the same ops. A power or an atom follows the minus bare (`-u^2`, `-sqrt(u)`), and anything looser goes in parentheses (`-(u - v)`, `-(-u)`).

**On a page.** A minus that leads the top or the bottom of a fraction now stands in front of it. `-1/4` is `-\frac{1}{4}`, not `\frac{-1}{4}`; `u / -2` is `-\frac{u}{2}`; and `v + u / -2` is `v - \frac{u}{2}`.

**In the steps.** -c/b has its sign worked out, so x - 3 = 0 is `\frac{3}{1}` and not `-\frac{-3}{1}`. The roots' fraction writes no `-0` when b = 0. A negative 2a turns the fraction over, since (b ± √D)/|2a| is the same pair: -x^2 + 4 = 0 is `\frac{\pm \sqrt{16}}{2}` and not `\frac{-0 \pm \sqrt{16}}{-2}`. A test had pinned the old `-\frac{-1}{2}`, and it now holds the new form.

Every new case failed before the change. The round trip of 50,000 random fields through text still holds. The gallery and the reference did not change, and the WebAssembly core and the host agree on 36 of 36 questions.
