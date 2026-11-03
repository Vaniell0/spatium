---
date: 2026-11-03
title: the language read from MathJSON, the output of a LaTeX editor
files: [include/spatium/io/field_mathjson.hpp, include/spatium/io/document.hpp, tests/test_field_mathjson.cpp, site/wasm/questions.txt]
open:
  - the page does not send MathJSON yet: MathLive in the question line is the page's step
  - the text reader still reads what a person types on the page; once MathLive writes the input, the reader's page-only paths can go
  - MathJSON that says more than the language -- a derivative, an integral inside an expression, a sum, a complex number -- is refused by name; each is a function of the language to add, not a rule of the translation
  - a repeating decimal ({"num": "1.(3)"}) is refused, though the exact arithmetics could read it as 4/3
---

A page that takes its input from a LaTeX editor (MathLive, and the Compute Engine behind it) gets MathJSON from it: a number, a symbol, or a function as an array whose first item names it, `["Add", ["Square", "x"], ["Negate", 2]]`. A field of a document's JSON may now be MathJSON in place of text. It is read as the text of the language first (`mathjson_text`, `io/field_mathjson.hpp`), by a table of the MathJSON functions the language has -- Add, Subtract, Negate, Multiply, Divide, Rational, Half, Power, Square, Sqrt, Root, Exp, Ln, Log, Sin, Cos, Tan, Abs, Max, Min, Delimiter -- and its named constants (Pi, ExponentialE, PositiveInfinity, NaN). That text then goes through the language's own reader, so a formula means one thing whichever way it came in, is refused for the same reasons, and prints back as the language. No grammar of our own is added: the LaTeX is read by the editor, and what crosses is a table.

A root's equation may stand whole as its left side, `["Equal", a, b]`, which is what the editor writes for `2x^2-3x-5=0`. The ends of an interval may be JSON numbers. A subscripted name, `x_1`, is a name. The text is written with parentheses only where they are needed and a minus where a sum has one: MathLive's canonical form of `2x^2-3x-5` reads as `2 * x^2 - 3 * x - 5`.

What MathJSON can say and the language cannot is refused with the name it was given and where it stood: `document: a, value: at [2]: 'Integrate' is not in the language`. The editor's own error, a complex unit, the wrong number of arguments and a name the language cannot read are refused the same way.

The test holds each function of the table two ways: the document read from MathJSON is the same document as the one written by hand (`same_document`), and the text it is read as is that hand-written text; the table and the cases are checked against each other, so a function added to the table without a case fails. 599 assertions with the document tests. The WebAssembly core and the host agree on 46 of 46 questions, among them an exact root, (3 ± √17)/2, a synthesised √(−p) and a refusal, each sent as MathJSON.
