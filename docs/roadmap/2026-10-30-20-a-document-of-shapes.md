---
date: 2026-10-30
title: a document's shapes, built as the scene its calls would build
files: [include/spatium/io/document.hpp, include/spatium/io/document_trace.hpp, include/spatium/io/build.hpp, include/spatium/io/trace_text.hpp, examples/gallery.cpp, tests/test_document.cpp]
resolves: [2026-10-30-19#1]
open:
  - an integral is refused in a scene: a scene takes values, parameters, roots and shapes
  - a shape's motion, colour and material are not statements; a document builds the rest shape
---

A document now has shapes: a line `S0 = torus(R, r)` names the shape one operation makes, from numbers and the shapes above it.

    R = 2
    r = 0.5
    h: (r + h)^2 - 2 * r^2 = 0 on [0, 1]
    S0 = torus(R, r)
    S1 = offset(S0, h)

The operations are a table, `shape_operations`, with the kind of each argument (a number, a field, a shape by its name) and how many may be left out: `torus`, `cylinder`, `sphere`, `flake`, `offset`, `scatter`, `union`. Their names are words of the language. `check` refuses a call with too few or too many arguments, a shape where a number is taken and a number where a shape is, and says the line; a shape line is printed back as written and carried through JSON as its operation and the texts of its arguments.

`to_trace(doc, chosen)` (`io/document_trace.hpp`) builds the `Trace` a program would build by calls, line by line: a value is `let`, a root `let` with its bracket, a shape the factory or operation it names, its numbers evaluated with every name above bound and a field (an offset's thickness) handed over as written, so the scene's LaTeX says it as the document did. A parameter is drawn at the value `chosen` gives it, through the new `Trace::choose`, which keeps its interval and refuses a value off it; LaTeX writes it `p = 1, p in [0.5, 4]`. A builder's refusal (an offset of an open surface, a count that is not whole) comes back naming the line's definition.

The pair: the gallery's shell with a scatter of spheres on it and the union of both, built from the document and by calls, gives the same nodes (kind, recipe, base, item, target, count, seed, children), the same definitions to the bit and the same LaTeX; a wider tube is a different recipe and a different line. The gallery's scene card is now built from its document; its LaTeX reads `2r^2` where the hand-built field read `2r * r`.
