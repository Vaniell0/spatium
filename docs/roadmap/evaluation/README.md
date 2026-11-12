# Evaluation

## Where it stands

One door for a statement's value, `evaluate(doc, name, chosen, engine)` (`include/spatium/answer/evaluate.hpp`): a number when no parameter is left free, a function of the one it reads otherwise. Every route -- evaluated, closed form with its steps, every root proved on a box, a root on a bracket, an integral, synthesised -- is kept as a candidate with its cost or refusal. The `Engine` says how it computes: double, float, a 64-bit fraction, decimal64 where the compiler has it, Real50; in fractions a quadratic's roots are exact (5/2, (3 ± √17)/2). A number's digits are proved over the reals: it is computed again in intervals (`core/interval.hpp`), each constant held as written, every end rounded outward, and the digits both ends share are proved (sqrt(2) + 1 to 16 digits); where every operation is one on fractions the value is the fraction itself (0.1 + 0.2 is 3/10). A number shows how it was computed, an operation a step as the program ran, and is exact where every operation is one on fractions (an exercise's 12, not 12.000000000000002). An equation's roots are proved on its bracket, or on the doubles' whole range without one: the box is cut where undecided (`prove_on`), a piece is excluded where the interval leaves out zero, and holds exactly one root where Newton's interval step lands inside it (`roots_on`); 2^x = 8 is 3 to 16 digits, with no other root on the line. An inequality's answer is the set where it holds, by the interval method: the axis cut where it is undecided until each piece is proved to hold it or not, its ends roots proved alone (x^2 - 4 > 0 is (-inf, -2) U (2, inf)), drawn on the page as the method is on paper. A system of equations, as many as its unknowns and linear in them, is solved together by elimination, in fractions where its coefficients are (x + 2y = 5, 3x - y = 1 is x = 1, y = 2), and checked by putting the solution back; one that is not linear has every root on every double proved (Krawczyk's operator in `roots_on`). `respond.hpp` writes it as the JSON a page reads, and every answer as one sheet: what it is, how strongly it holds, the result in one line and its sections in one order, which is all the page draws.

## Now

- refusals that say why: a bracket ending where the equation is infinite (0^-1) is refused as having no value
- the ends of the doubles, where a polynomial overflows, by a chart at infinity (x = 1/t, `Series`)

## Later

- integrals and synthesis in the engine's arithmetic, and proved: a synthesis proved on its pieces; roots at fifty digits
- exact roots of a cubic and a quartic
- steps for every route, written from the step data by templates

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area evaluation`.
