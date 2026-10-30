---
date: 2026-10-30
title: a document of named statements, and a question read through one
files: [include/spatium/io/document.hpp, rsc/spatium/query.hpp, tests/test_document.cpp, include/spatium/io/field_parse.hpp]
open:
  - a document holds values, parameters, roots and integrals; shapes and scene operations (torus, offset) are not statements yet
  - a root or an integral that reads another root or integral is refused: one question is answered at a time
---

`io/document.hpp` is what a person writes before anything is computed: statements a line, each reading only the names above it.

    r = 0.5
    p in [0.5, 4]
    x: x^2 - p = 0 on [0, 10]
    I = integral exp(-p * t^2) dt on [0, 1]

A product is written with `*`, or left out only after a number or a closing parenthesis (`2x`, `2(x + 1)`): `pi * r^2`, never `pi r^2`, since a space between two names would make names and functions ambiguous. A parameter is chosen from an interval (`in`); a root and an integral are found on one (`on`), the integral's variable named by its `dt`. `read_document` reads the text, skipping blank lines and `#` comments; `check` holds the names (each defined once, none a word of the language, every name read defined above, a parameter's ends numbers, a root's equation reading its unknown) and says the line, and a parse error its column in that line. Each field keeps the text it was read from, so `to_text` prints a document back as it was written; `to_json` carries the same texts and `document_from_json` reads them through the same reader, so the two forms cannot mean different things. Both round trips are tested, with a case that must fail.

A typed question is now read through its document. `question_document` reads `solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]` into the parameter and the root it asks for, from the question's own text, so an error still names a column of what was typed; `query_of` reads a root or an integral of one parameter into the Query `answer` synthesises, with the values it reads put in place. `ask` is the two together. The word reader it replaced stays as the witness: on ten questions (those of the gallery and the reference, and the forms around them) both give the same Query, field for field, and both refuse the same broken ones. One thing the document route does that the word reader did not: a question reading a name nothing defines is refused when it is read, not when it is answered. The gallery, the reference and the WebAssembly core give the same answers as before (22 of 22 agree).
