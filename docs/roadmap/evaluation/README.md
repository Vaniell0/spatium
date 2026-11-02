# Evaluation

## Where it stands

One door for a statement's value, `evaluate(doc, name, chosen, engine)` (`rsc/spatium/evaluate.hpp`): a number when no parameter is left free, a function of the one it reads otherwise. Every route -- evaluated, closed form with its steps, a root on a bracket, an integral, synthesised -- is kept as a candidate with its cost or refusal. The `Engine` says how it computes: double, float, a 64-bit fraction, decimal64 where the compiler has it, Real50; in fractions a quadratic's roots are exact (5/2, (3 ± √17)/2). `respond.hpp` writes it as the JSON a page reads.

## Now

- refusals that say why: a bracket ending where the equation is infinite (0^-1) is refused as having no value

## Later

- roots, integrals and synthesis in the engine's arithmetic
- exact roots of a cubic and a quartic
- steps for every route, written from the step data by templates

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area evaluation`.
