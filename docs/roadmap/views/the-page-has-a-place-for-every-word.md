---
date: 2026-11-05
title: the page has a place for every word of the library's language
files: [examples/site_reference.cpp, web/src/gen/language.ts, web/src/registry.ts, web/test/registry.test.ts, include/spatium/io/document.hpp, rsc/spatium/evaluate.hpp, tests/test_document.cpp]
open:
  - the sheet's shapes and section titles are string literals in respond.hpp, not a table, so they are not in the contract yet
  - the requests of the WebAssembly core (rows, evaluate, plot) are not in the contract; core.ts names them by hand
  - the registry says where a word is shown, but nothing yet checks that the place really draws it (a keyboard key, a canvas node)
---

The page's code named the library's words by hand: a relation, a route and an arithmetic were strings in TypeScript, and a word the library gained reached the page only if someone remembered it. `site_reference`, which already wrote `site/data/reference.json` from the reader's tables, now also writes the statement kinds, the relations, the routes with their captions, the arithmetics and the shape operations, each from its own C++ table (`statement_kinds` and `relations` in io/document.hpp, `routes` in rsc/spatium/evaluate.hpp, beside `arithmetics` and `shape_operations`), and writes the same words as TypeScript, `web/src/gen/language.ts`: a constant list each and its union type. `--check` holds both files. `web/src/registry.ts` gives every word a place on the page, a record over each union, so a new word is a type error there until it has a place or is named a gap with what would show it. The coverage test reads the words from reference.json, fails on a word the registry lacks or has that the library does not, and prints the gaps: today 11 of 50 words (a shape and its seven operations, which no page draws yet; float, decimal and real50, which no control computes in). The relations table is the reader's: each is read and printed back (`test_document`).
