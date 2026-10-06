# What depends on the dimension

Most of what the library does is the same at every dimension: `Euclidean<N>`,
`Sphere<N>` and `Hyperbolic<N>` answer the same questions through the same
code, and the connectivity matrix (`docs/connectivity.md`) runs every probe
over each family at N = 1, 2, 3, 4, 7 and 8 to hold that. Where it is not
the same, mathematics decides, not the implementation: some structures exist
at certain dimensions only. This page names them, with the theorem behind
each, and says where the library stands.

| Capability | Exists at N | Theorem | Library today |
|---|---|---|---|
| Cross product: bilinear, orthogonal to both factors, `\|a×b\|² = \|a\|²\|b\|² − (a·b)²` | 3, 7 (and trivially 0, 1) | Brown and Gray, 1967 (Eckmann, 1943) | `Vec<T,3>::cross` only; N = 7 would need octonion structure constants |
| Normed division algebra, `\|ab\| = \|a\|\|b\|` | 1, 2, 4, 8 | Hurwitz, 1898 | reals (1), `Complex` (2), `Quaternion` (4); N = 8, the octonions, is not implemented and is not associative |
| Group structure on the sphere S^N | 1, 3 | Hopf; Serre | none: S¹ (unit complex numbers) and S³ (unit quaternions) are not exposed as groups; `SO3`, `SE3` are separate types |
| Global tangent frame on S^N (parallelizable) | 1, 3, 7 | Bott and Milnor, 1958; Kervaire, 1958 | no code |
| A nowhere-vanishing tangent field on S^N | N odd | hairy ball, Brouwer, 1912 (even N) | no code; `unit_tangent` in tests is a pointwise choice, not a field |

## Two directions, one of them checkable

`tests/test_dimension_theorems.cpp` holds the library to the direction a
machine can check: **what the library implements lies inside what the
theorem allows.** A cross product at N = 4, or a group structure on S², would
be a defect however plausible it looked, and the test fails on it. The other
direction -- allowed and not implemented -- is a gap of glue, listed in the
last column and not a failure.

A red cell of the connectivity matrix can therefore be of three kinds:

- **no glue** -- the pieces exist and nothing joins them (an ADL cylinder
  as a product factor, before `ProductSpace` took one);
- **no support** -- a missing basis (geodesics on a chart, before
  `MetricChart`);
- **impossible by theorem** -- no code can turn it green. The first
  capability of that kind the library has a type to ask about is the cross
  product at N ≠ 3, 7: asking `Vec<T,4>` for one does not compile, and
  should not.

## What the dimension axis found

Run over N = 1 … 8, the probes were uniform except in one place: the Fréchet
mean on `Hyperbolic<8>` returned NaN. The Karcher iteration took a step of 1,
exact in a flat space and right near the mean on a sphere, but in negative
curvature it overshoots once `d coth d` passes 2 -- points more than about
1.5 from their mean -- and random tangents are longer in eight dimensions
than in two. `frechet_mean` now halves the step until the variance falls
(`algebra/calculus.hpp`), and a test holds spreads of 1, 2 and 3 on
`Hyperbolic<8>`.
