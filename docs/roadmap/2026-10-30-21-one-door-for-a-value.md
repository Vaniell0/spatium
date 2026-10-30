---
date: 2026-10-30
title: one door for a definition's value, with every route it was chosen among
files: [rsc/spatium/evaluate.hpp, rsc/spatium/respond.hpp, rsc/spatium/query.hpp, rsc/spatium/closed_form.hpp, include/spatium/io/document.hpp, include/spatium/io/document_trace.hpp, site/wasm/core.cpp, tests/test_rsc_evaluate.cpp]
open:
  - a root solved in closed form is a quadratic in numbers; one whose coefficients read the parameter is synthesised, not solved as sqrt(p) by the formula
  - a definition reading another root or integral is still refused; one definition is evaluated at a time
  - the route taken is the first in a fixed order that gives a value; nothing learned chooses yet
---

`evaluate(doc, name, chosen)` (`rsc/spatium/evaluate.hpp`) gives a definition of a document its value. The parameters in `chosen` are put in as numbers, and what the definition still reads decides what its value is: no parameter, a number; one, a function of it on its interval; more, a refusal that names them.

    p in [0.5, 4]
    x: x^2 - p = 0 on [0, 10]

With nothing chosen, x is synthesised as `sqrt(p)` with its certificate, as the calculator shows it. With `p = 2` it is the number, solved in closed form and kept where it lies on the bracket, as a scene draws it. With two parameters and one chosen it is a function of the other. A value that reads a free parameter is its own expression, exact.

A root may now be written without a bracket, `x: 2x^2 - 3x - 5 = 0`. It is solved in closed form, every root; one that is not a quadratic is refused with the reason and that it needs a bracket. A question, and a scene, refuse it by name.

Every route that can give the value is tried and kept as a candidate, with its cost in evaluations of the definition's field or the reason it gives none; the value is the first, in a fixed order, that gives one (`take`). So a closed form and a root on the same bracket are both computed, and the record says what was chosen among and at what cost, not only the choice. The closed form counts 8 evaluations (three for the coefficients, five to check the degree); a synthesis of `sqrt(p)` about 118,000.

`respond` is now the calculator's shorthand written out as a document and evaluated: `2 + 2` is the value `answer = 2 + 2`, an equation the root without a bracket, a `solve ... for ...` question its parameter and root. Its answer carries the document it was read as and the candidates. The gallery's earlier fields are unchanged to the character. A card keeps the routes but not their costs: a root found to the last bit takes as many evaluations as the libm's last digits lead it to (117,953 on glibc, 117,945 in WebAssembly on musl).

The WebAssembly core takes a document: `spatium_evaluate` reads one JSON object, the document as text or as its JSON, the name and the chosen parameters, and the check holds it against the host build on six such requests.

The pairs: the synthesised function at five values of p against the number with p chosen, to the tolerance; the closed form against the root on the same bracket; a value written as a function against the same value chosen, to the bit; the shorthand against the document line it stands for, as the same JSON; a document sent as text against the same document as JSON.
