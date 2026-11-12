---
date: 2026-11-12
title: a statement has one spelling, read and written from one table, and MathJSON is read as a tree
files: [include/spatium/io/document.hpp, include/spatium/io/field_mathjson.hpp, include/spatium/io/field_parse.hpp, include/spatium/io/field_text.hpp, include/spatium/answer/respond.hpp, tests/test_io_roundtrip.cpp, tests/test_field_mathjson.cpp, site/wasm/check.mjs]
open:
  - a row of the editor that a form does not fit is refused as the furthest form saw it; a person may want the nearest form named instead
  - the MathJSON a form is matched against is a list of patterns per row; a pattern for a statement the editor writes another way is added by hand
---

A statement reached the library in four forms -- a line of a document, the document's JSON, the shorthand a calculator takes, a row of a formula editor as MathJSON -- and left it in three: text, JSON, the LaTeX row a page shows. Each pair was written by hand. The grammar of a statement was in nine places (four readers: `read_statement`, the shorthand, the system's reader, `document_from_json`, which rebuilt a line of text; five writers: `to_text`, `to_latex`, `to_json`, `mathjson_row_text`, and the questions the answer sheet offers next, put together with `std::format`), and the order of operations in three (the reader, the printer's `Binding`, MathJSON's `Tight`), with nothing holding them to each other.

First the measurement: `tests/test_io_roundtrip.cpp` reads back what every writer writes, for every kind of statement and every form a root takes -- the text, the JSON, the LaTeX rows through MathJSON, the shorthand, and the document an answer says it read. Its first run was the list of breaks: a system's LaTeX printed one equation and the page wrote the others beside it, and `less`, which LaTeX writes as an Iverson bracket, did not read back.

Then the language once. A function of the language is one row of `functions` with its MathJSON head; a MathJSON tree is read straight into the field the text reader builds for the same formula, op for op (`field_mathjson.hpp`), instead of being printed as text and read again, and where a document needs text the library's printer writes it. `["Boole", ["Less", p, 1]]` reads as `less(p, 1)`; a whole power is a power, `p - p^2`, as typed.

Then the statement once: `statement_forms` (`io/document.hpp`), one row a form -- the line as typed, the editor's LaTeX row, the tail a bracket adds, the MathJSON trees the editor writes for it. Reading a line, printing it, its LaTeX, its JSON and a row of the editor all walk these rows. A form with no name is a shorthand (`read_shorthand`: `2 + 2` is `answer = 2 + 2`, `x^2 = 2` a root in x), and lines that are each an equation are a system (`read_system`) unless each defines a new name by the names above it, `a = 2` over `f = a^2 + 1`. A row that fits no form is refused where it departs from the nearest pattern, by the node's name.

Held by the round trips (every form reads back one of the cases, every writer's output reads back as what it started from), by the reference's examples read back from their rows on the page's own parser (35 of 35), and by the WebAssembly core against the host on every question of the check.
