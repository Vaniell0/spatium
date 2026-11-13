---
date: 2026-11-13
title: A matrix has a unary minus, which the law check asked of SPD's tangent
files:
  - include/spatium/algebra/matrix.hpp
  - tests/test_matrix.cpp
  - rsc/spatium/laws.json
  - rsc/spatium/tasks.json
  - rsc/spatium/laws-report.md
resolves: [rsc/the-laws-of-the-spaces-found-and-held#2]
---

`Matrix` had sum, difference and multiples, and no `-m`. A space whose points are matrices -- `SPDAffineInvariant`, whose tangent is a symmetric matrix -- therefore had no negation of a tangent, which `rsc_laws.py missing` reported, and a product with it as a factor (`ProductSpace<SPDAffineInvariant<2>, Sphere<2>>`) declared a negation whose body did not compile, which no `requires` could see. `operator-()` is now there, entry by entry; the test holds it against the matrix times -1 and against itself applied twice, and the product's tangent negates and comes back. With it the law check's list of missing arithmetic is empty, and its laws are the same 834. It was the first change a judge of branches accepted against a task the law check raised.
