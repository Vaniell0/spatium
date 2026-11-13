---
date: 2026-11-13
title: The laws of the spaces, found by enumeration with no concept told, and held column by column
files:
  - rsc/spatium/rsc_laws.cpp
  - rsc/spatium/rsc_laws.py
  - rsc/spatium/test_rsc_laws.py
  - rsc/spatium/laws.json
  - rsc/spatium/laws-report.md
  - rsc/spatium/tasks.json
  - rsc/spatium/laws-triaged.json
  - .github/workflows/ci.yml
open:
  - Hyperbolic::distance(x, x) is about 4e-5, not 0, for x eight away from the origin -- acosh of 1 plus e^{2d} of rounding; the open task of H1, H2, H4 and S2xH2
  - a tangent of SPD2-affine, a matrix, has no unary minus, and a product with it as a factor declares one whose body does not compile, which a requires-check cannot see
  - two classes the laws draw have no concept in the library -- log global (E, H, SPD, H2xE3, and S2-derived, a chart that is no sphere) and the symmetric spaces (where -log(exp u) = log(exp(-u)) holds) -- and the member-based MetricSpace, Complete and RiemannianManifold are told apart from Distanced by no law
  - the hand-written probes of the matrix are laws too; which of them the enumeration finds, and which it misses at size 7, is not yet set side by side
---

`rsc_laws.cpp` builds terms bottom up over the operations of `core/access.hpp` -- log, exp, the metric, the distance -- with the arithmetic of tangents and numbers, typed by sort (a number, a point, a tangent at a point term) and by a number's degree in length, so that numbers of different degree are never added or compared and a fact about where a sampler puts its points ("d <= 4") cannot be written. Every term is evaluated on the 25 spaces of the connectivity matrix on the same 732 tuples: 192 drawn by each space's sampler, the rest moved by a breaker -- each subset of q, r, u, w drawn toward p, coincident points, tangents along log(p, q) at six scales, and tangents made 8 long in the space's own metric, past any cut locus here. Without the breaker and the degree, 970 classes of spaces came out, most of them inequalities that held by chance; with them, 25. A term equal to an older one wherever both are defined is merged, and the merge is a law of every space; a law that is an instance of a smaller one, or follows from smaller ones by equality, congruence or a constant factor, is dropped. At size 7 that is 790 terms kept, 357 merged and 834 laws, in 182 s on eight E-cores beside another job.

The enumeration is the survey, run when a branch is grown; what it leaves is data. `laws.json` keeps each law with the spaces it holds and is defined on, and `rsc_laws.py check` evaluates it again, one column per space -- all of them in CI, 130 s of which the compile is most, or those whose headers a change touched (`--changed REV`, by each space's include closure) -- failing on any law that starts or stops holding. Held against `discover` itself, the check finds no change on any of the 25 spaces. `test_rsc_laws.py` holds what the comparison reports and what becomes a task, on a snapshot of ten spaces.

The classes the laws draw are set against the concepts of `graph.json` on the 15 spaces that are nodes of it. `Distanced` comes out exactly: 330 laws, exp(log(q)) = q among them, hold where the library gives a distance, and neither erased surface. `tasks.json` lists a space that breaks laws holding on at least 80% of where they can be asked, and arithmetic a tangent does not have; `laws-triaged.json` keeps the ones looked at and found to be no defect, each with its reason. The first run found that Hyperbolic lost digits far from the origin (fixed in `spaces/hyperbolic-far-from-the-origin`), and after that the distance of a far point to itself, still open.
