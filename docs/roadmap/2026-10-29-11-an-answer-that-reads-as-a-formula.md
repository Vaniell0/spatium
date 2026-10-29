---
date: 2026-10-29
title: a field as text, and an answer that reads as a formula in its own variable
files: [include/spatium/io/field_text.hpp, rsc/spatium/synthesis.hpp, tests/test_rsc_synthesis.cpp]
open:
  - the text is the lowered ops written out, with no simplification: a coefficient of 1 and a sum of two copies of a term stay as they are
  - plain text only; LaTeX or MathML for a page is the same walk with another spelling, not written
  - noise, gather and the geometric inputs are printed by name and component, not as a formula
---

`io/field_text.hpp` writes a field as text: its lowered ops in the usual notation, the inputs under names the caller gives, a constant in the shortest digits that read back to the same double; a field that does not lower has no text and says why. It is how a synthesised answer is read. The volume of the prism under the triangle (1, 0, 2), (0, 1, 4), (0, 0, p) as a function of p, computed for each p by the integration door over its shadow, comes back from the synthesis on p in [0, 10] as `(0.16666666666666663 * p + 1.0000000000000002 * 1)`, certified to 1.8e-15 on 2001 points: 1 + p / 6, the shadow's area times the mean height.

It did not read like that at first: `(0.119 * ((p - 5) * 0.2 + ((p - 5) * 0.2 + p)) + 1.238 * 1)`, the same line written through the interval's own variable. Three rules of description length in the synthesis made it read in p. At equal error the set of fewer IR nodes wins, and the interval's variable is the five nodes (x - centre) * scale, not one. Every pair is tried among the first terms (the shortest) when there are too many to try them all, rather than none. And below rounding there is nothing to improve, so a swap that only moved the noise no longer trades a short form for a long one.
