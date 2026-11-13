---
date: 2026-11-04
title: a question is the rows of a formula, and the answer gives them back
files: [include/spatium/io/document.hpp, include/spatium/io/field_parse.hpp, include/spatium/answer/respond.hpp, include/spatium/answer/query.hpp, examples/gallery.cpp, examples/site_reference.cpp, site/gallery/gallery.json, site/data/reference.json, site/wasm/core.cpp, site/wasm/build.sh, site/wasm/check.mjs, site/wasm/questions.txt, site/style.css, web/src/page.ts, tests/test_document.cpp, tests/test_answer_query.cpp, tests/test_answer_respond.cpp, tests/test_answer_evaluate.cpp]
resolves: [document/the-rows-of-a-latex-editor-are-a-document#1]
open:
  - a cell typed as text still takes the document's own words (`p in [0.5, 4]`, `x: ... on [0, 10]`); the formula editor is the input, and ∈ typed as text is not read
  - a parameter's interval written backwards, p in [4, 1], is read, and refused only when a question is asked of it
  - the prism card is asked in code and described in words; it has no rows
  - the steps of an answer are not yet written from what was computed
---

A question has one form: the rows of a formula editor. The English grammar the calculator had beside them, `solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]` and `integrate f for t in [a, b], p in [c, d]`, is gone, with its reader (`ask`, `question_document`, the word reader `ask_by_words`), the `for` branch of `respond`, and the preview the old page used. The same question is the rows `p ∈ [0.5, 4]` and `x² − p = 0, x ∈ [0, 10]`, which the core already read as the document lines `p in [0.5, 4]` and `x: x^2 - p = 0 on [0, 10]`; a row reading a parameter above is answered as a function of it, as before. Three ways to write a question became two views of one document, and the code (library, page, tests) is 93 lines shorter.

`to_latex(Statement)` (`io/document.hpp`) writes a statement as the row an editor reads back as it: `p \in \left[0.5, 4\right]`, `x^{2} - p = 0, x \in \left[0, 10\right]`, `I = \int_{0}^{1} \exp\left(-p \cdot x^{2}\right)\,dx`, a name of more than one letter upright, a row written with no name with none. Every answer carries its question so (`question_latex`), the gallery's cards and the reference's examples show it as formulas, and a click opens it as those rows.

The pair: each of the 12 cards of the gallery and the 5 examples of the reference is written as rows by the library, parsed by the page's own formula parser in node, and asked of the WebAssembly core again; all 17 give back the same rows, route and answer (`site/wasm/check.mjs`). It found a name an editor does not keep: `\mathrm{erf}` is read as the error function, not as a name, so the card's integral is called E. A line in a document's own words that does not read is refused with the document's reason (`x reads 'p', which no line above defines`), not the shorthand's. Messages lost the `ask:` prefix of the reader they came from.
