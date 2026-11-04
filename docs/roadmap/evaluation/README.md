# Evaluation

## Where it stands

One door for a statement's value, `evaluate(doc, name, chosen, engine)` (`rsc/spatium/evaluate.hpp`): a number when no parameter is left free, a function of the one it reads otherwise. Every route -- evaluated, closed form with its steps, every root proved on a box, a root on a bracket, an integral, synthesised -- is kept as a candidate with its cost or refusal. The `Engine` says how it computes: double, float, a 64-bit fraction, decimal64 where the compiler has it, Real50; in fractions a quadratic's roots are exact (5/2, (3 ± √17)/2). A number's digits are proved over the reals: it is computed again in intervals (`core/interval.hpp`), each constant held as written, every end rounded outward, and the digits both ends share are proved (0.1 + 0.2 is 0.300000000000000 to 15 digits). A number shows how it was computed, an operation a step as the program ran, and is exact where every operation is one on fractions (an exercise's 12, not 12.000000000000002). An equation's roots are proved on its bracket, or on the doubles' whole range without one: the box is cut where undecided (`prove_on`), a piece is excluded where the interval leaves out zero, and holds exactly one root where Newton's interval step lands inside it (`roots_on`); 2^x = 8 is 3 to 16 digits, with no other root on the line. `respond.hpp` writes it as the JSON a page reads.

## Now

- refusals that say why: a bracket ending where the equation is infinite (0^-1) is refused as having no value
- inequalities on a box, the other consumer of `prove_on`: where a statement holds, fails, or is undecided
- the ends of the doubles, where a polynomial overflows, by a chart at infinity (x = 1/t, `Series`)

## Later

- integrals and synthesis in the engine's arithmetic, and proved: a synthesis proved on its pieces; roots at fifty digits
- exact roots of a cubic and a quartic
- steps for every route, written from the step data by templates

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area evaluation`.
