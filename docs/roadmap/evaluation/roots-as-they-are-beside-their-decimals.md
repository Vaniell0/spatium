---
date: 2026-11-05
title: a quadratic's roots as they are beside their decimals, and a step's lead and note
files: [include/spatium/answer/closed_form.hpp, include/spatium/answer/evaluate.hpp, include/spatium/answer/respond.hpp, include/spatium/io/steps_text.hpp, tests/test_answer_respond.cpp, tests/test_answer_evaluate.cpp, site/gallery/gallery.json]
open:
  - a number's own computation (sum, power, sine) has no lead or note beyond its label
  - a fraction from an evaluation is written as it is beside its decimal already; one from a linear system is not yet
---

A quadratic's roots as they are, fractions or (p ± q√m)/d, were found only when the engine was asked to compute in fractions, which a card's Numbers setting did and the calculator no longer has. So the calculator showed 2.5 where the answer is 5/2. They are now sought whatever the engine, and the roots' step writes both forms:
- a fraction is followed by its decimal, `x_1 = \frac{5}{2} = 2.5`, or by twelve digits when the decimal does not end;
- a pair with a root in it is followed by its decimals to twelve digits, `\frac{3 \pm \sqrt{17}}{2} \approx 3.56155281281,\ -0.561552812809`;
- a pair about zero is written as one, `\pm\sqrt{2} \approx \pm 1.41421356237`.

The values stay the double's, and the guarantee is exact where the roots are.

A step now carries what a page sets in bold and what follows it (`lead`, `note`), beside the label code keys on. The quadratic's steps read "Read the coefficients from the equation's values at 0, 1 and −1", "Discriminant b² − 4ac", "Case D > 0: two real roots", "Roots", and "Check by Vieta's formulas". The sheet writes them, while the answer's own steps are unchanged.

The check (`test_answer_respond`, `test_answer_evaluate`):
- `2x^2 - 3x - 5 = 0` is exact, its roots `\frac{5}{2} = 2.5` and `-1`;
- `x^2 - 2 = 0` reads the same in double and in fractions;
- with no engine, `x^2 - 3x - 2 = 0` still finds `\frac{3 \pm \sqrt{17}}{2}` while its values stay doubles.
