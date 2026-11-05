---
date: 2026-11-04
title: every root of a system on a box, proved
files: [include/spatium/algebra/roots.hpp, modules/algebra/roots.cppm, rsc/spatium/evaluate.hpp, rsc/spatium/respond.hpp, tests/test_rsc_respond.cpp, examples/gallery.cpp, site/wasm/questions.txt]
open:
  - a root where the Jacobian is singular (two curves that touch) stays undecided, as a double root does on a segment
  - at most four unknowns; the budget is 65536 evaluations of the system
  - no steps yet: the proof's boxes and Krawczyk's steps are what a page would show
---

A system that is not linear was refused by name. `roots_on` now takes a box in N unknowns as it takes a segment: `prove_on` cuts the box where it is undecided, a piece where some equation's interval leaves zero out holds no root, and Krawczyk's operator K(X) = m - Y f(m) + (I - Y J(X))(X - m) stands for Newton's step -- J the Jacobian on the piece, N passes of dual numbers of intervals, Y the inverse of its middle: K(X) missing X leaves no root, K(X) inside X holds exactly one, and a root on a cut is retried on K(X) widened, as on a segment. `evaluate` asks it where the linear route says a system is not linear, on every double of each unknown. x^2 + y^2 = 4, x - y = 1 has its two roots proved, x = (1 ± √7)/2 to 1e-14, and no other on the plane, in 1100 evaluations; a circle and a line that do not meet are proved to have no real root. The pairs: each root against the quadratic the substitution y = x - 1 leaves, 2x^2 - 2x - 3 = 0; the linear route asked first and refusing by name (`test_rsc_respond`); WebAssembly against the host 56 of 56, the gallery's circle and line among them.
