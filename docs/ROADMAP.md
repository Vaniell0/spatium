# Roadmap

## Vision

Computational geometry on arbitrary Riemannian manifolds. The only C++ library that unifies abstract space theory and practical geometry through concepts. Target: geodesic Voronoi/Delaunay, shortest paths, remeshing — on any Surface, not just Euclidean.

This document tracks what's done and what's open, by area — not by version number and not by date. Version-gated milestones invite the failure mode they're meant to prevent: negotiating what counts toward which release, uneven development where some gates get stuffed and others starve, and a "two parallel tracks" split that has to be reconciled later — the same fragmentation branches/forks cause when they diverge instead of staying one line. Completed work below is dated where the date is known, named by content where it isn't. Open work is grouped by topic in **Backlog**; an item moves up into **Completed** the day it actually ships, not on a schedule.

---

**Open findings from the September 2026 audit live in [`findings-2026-09.md`](findings-2026-09.md)** -- bugs with their probes, RSC's breakages, the binary's measured errors, the root-level direction. Pick an item there, fix it with a failing test, delete it there.

---

## The layers

The library is read bottom up, and each layer reads only the one below it. This is the map of what each layer holds now; what each takes next is in its area below.

| Layer | What it holds | Where | Area |
|---|---|---|---|
| 1. IR | a quantity as one `Field`: numbers, the chart's `u` and `v`, named symbols, operations; `bind`, `free_symbols`, `simplified` | `io/field.hpp`; as text and LaTeX `io/field_text.hpp`; read back `io/field_parse.hpp`, and from MathJSON `io/field_mathjson.hpp` | [ir](roadmap/ir/) |
| 2. Document | what a person writes before anything is computed, a line each: a value `r = 0.5`, a parameter `p in [0.5, 4]`, a root `x: x^2 - p = 0 on [0, 10]`, an integral, a shape `S1 = offset(S0, h)`; text and JSON read into one; its shapes built as a `Trace` | `io/document.hpp`, `io/document_trace.hpp`, `io/build.hpp` | [document](roadmap/document/) |
| 3. Evaluation | one door for a statement's value, `evaluate(doc, name, chosen)`: a number when it reads no parameter left free, a function of the one it reads otherwise; the route evaluated, solved in closed form with its steps, a root on a bracket, an integral, or synthesised and certified, with every route tried and its cost; `respond` the calculator's shorthand written out as a document, and the WebAssembly core taking a document | `rsc/spatium/evaluate.hpp`, `respond.hpp`, `query.hpp`, `closed_form.hpp`, `synthesis.hpp`; `algebra/roots.hpp`, `algebra/quadrature.hpp` | [evaluation](roadmap/evaluation/) |
| 4. Lowering | what a computed document becomes: text, LaTeX, GLSL for Vulkan compute (std430) and for WebGL2 (GLSL ES 3.00, a symbol read from a uniform slot, so a moved number is new data and not a new shader), the graph of a function as a WebGL2 shader, a picture from the CPU ray tracer | `io/field_text.hpp`, `io/trace_text.hpp`, `io/field_pod.hpp`, `io/field_glsl.hpp`, `io/plot_glsl.hpp`, `render/`, `render/gpu/` ([`gpu.md`](gpu.md)) | [lowering](roadmap/lowering/) |
| 5. Views | the page: windows over one document, showing what the core answers and composing nothing | `site/` | [views](roadmap/views/) |

A question typed into the calculator is a document of one statement: `solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]` is read as the parameter and the root it asks for, and answered by layer 3. There is no second language for questions.

What crosses out of the library is the document, as text or JSON; there is no C ABI. The boundary and the rules of lowering (structure is code, numbers are data) are [`gpu-abi-design.md`](gpu-abi-design.md).

---

# Completed

---

## The areas

The roadmap is kept by area, in [`roadmap/`](roadmap/): each area is a folder whose README says where it stands, what is taken up now, and what would be good later, above its reasoning and earlier history. Each change adds one file to its area, named for what it did -- `python3 scripts/gen_roadmap.py --area <area>` prints an area's entries by date and what they left open, `--open` what is still not done in every area.

| Area | Where it stands | Now | Later |
|---|---|---|---|
| [IR: the field](roadmap/ir/) | A quantity is one `Field` (`io/field.hpp`): numbers, the chart's `u` and `v`, named symbols, and operations -- arithmetic, min/max, less, sin/cos, sqrt, log/exp, abs, any power (`Pow`), noise, gather, instance inputs, and a net of one hidden tanh layer whose weights are data (`Net`). | the page's input through MathLive, sending MathJSON; then the text reader's page-only paths removed; refusals as data: what, where, why, and what to write instead -- for a person editing and for automation retrying | vectors and matrices as values, once the matrices-and-systems fork is decided; interval arithmetic, for a certificate that bounds rather than samples |
| [Document](roadmap/document/) | What a person writes before anything is computed, a statement a line (`io/document.hpp`): a value, a parameter on an interval, a root with or without a bracket, an integral, a shape. | — | a TypeScript builder that writes the document's JSON; matrices and systems of equations as statements; a workspace for documents and scenes, edited as a canvas |
| [Evaluation](roadmap/evaluation/) | One door for a statement's value, `evaluate(doc, name, chosen, engine)` (`rsc/spatium/evaluate.hpp`): a number when no parameter is left free, a function of the one it reads otherwise. | refusals that say why: a bracket ending where the equation is infinite (0^-1) is refused as having no value | roots, integrals and synthesis in the engine's arithmetic; exact roots of a cubic and a quartic; steps for every route, written from the step data by templates |
| [Lowering](roadmap/lowering/) | What a computed field becomes on a device: GLSL for Vulkan compute (std430) and for WebGL2 (GLSL ES 3.00), a symbol read from a uniform slot and a net's weights from a uniform array, so a moved number or a trained net is new data and not a new shader (`io/field_glsl.hpp`). | — | slots and a net's weights in Vulkan; CUDA C from the IR, for those who port the library; a net wider than three inputs, which needs vectors in the IR |
| [Views: the page](roadmap/views/) | The page (`site/`) shows what the core answers and composes nothing: the solver, the gallery written by `examples/gallery.cpp`, the reference written from the reader's tables. | the page, step one: windows over one document in a grid, input through MathLive, the graph in its window, the workspace in its link, published on GitHub Pages; testing the page by many questions, from a generator shared with the synthesis data | free windows; on a phone, sheets with a switch; scenes in a window, from a TypeScript builder, drawn in WebGL2 |
| [RSC: synthesis and search](roadmap/rsc/) | `rsc/core` is generic: a typed graph of the library and a tree search guided by a learned prior and value, which reads on past a value that is not a number. | Pow among synthesis's operations, so the root of x^p + p = 0 is (-p)^(1/p) and not a Chebyshev polynomial | structural candidates for where to split an interval, before bisection and before a learned policy; infinity and chaos as classes in the search, after infinity is an input class in algebra; a model over the IR that names what the library lacks, admitted only through samples, certificates and green PRs |
| [Algebra](roadmap/algebra/) | Scalars and the operations on them (`include/spatium/algebra/`): vectors, matrices, quaternions, octonions, complex, dual numbers, `Series` as arithmetic of limits, calculus and Riemannian descent, quadrature as one door, ODE solvers, polynomial solvers with their steps, roots on a bracket; models with their gradients (softmax, a tanh net, score trees) and step rules; statistics (fitting, information, CART). | — | infinity as an input class across the interfaces; `Series` as arithmetic of limits, joined to quadrature |
| [Spaces](roadmap/spaces/) | The concept hierarchy (`core/`) and the spaces on it (`spaces/`): Euclidean, sphere, hyperbolic as one constant-curvature family, Minkowski as a Lorentzian space with no distance, L^2 of an interval, Lie groups as spaces, products, parametric and implicit surfaces with exact partials, a space from its metric alone. | — | — |
| [Geometry, mesh, spatial](roadmap/geometry/) | Primitives and their intersections and distances (`geometry/`), meshes as the materialiser at the edge (`mesh/`: subdivision at the space's own midpoint, geodesics, transport, Voronoi, tessellation, marching cubes), and spatial indices (`spatial/`: a BVH over any bound, geodesic ball trees). | — | — |
| [Physics](roadmap/physics/) | Mechanics on spaces and charts (`physics/mechanics/`): units, integrators, symplectic and variational schemes, certified continuous collision computed on the analytic description, gravity of many sources on closed spaces. | — | the library's own IPC, on smooth patches rather than triangles |
| [Relativity](roadmap/relativity/) | Metric-agnostic geodesic integration (`physics/relativity/`): Schwarzschild, Kerr, a metric as ten IR fields, its GLSL, the accretion disk, and the black-hole renderer (`examples/blackhole_live`). | — | two black holes merging, reporting on itself |
| [Discrete](roadmap/discrete/) | Finite sets and combinatorics, and the rules of chess as an oracle that reads no text (`discrete/`): the board, the moves, the position and its outcomes, a board of any number of axes with pieces as orbits of the cube's symmetries; the agents that learn from it are `examples/chess_agent` and `examples/chess_zero`. | — | what remains of the chess line, as the synthesis work needs it |
| [Sound](roadmap/sound/) | The synth engines and scores; the reasoning is below. | — | the work on its local branch, when it is taken up |
| [Tooling](roadmap/tooling/) | How the repository keeps itself true: CI with path-filtered jobs, the generators and their `--check` modes (connectivity, capabilities, consumers, dependency graph, roadmap), the C++23 modules build, and the conventions in `docs/conventions.md`. | the roadmap by area (this change) | a map of the architecture as one page: every layer and how they connect |
| [History](roadmap/history/) | The steering index of September 2026 and the reasoning of that time, kept as they were written. | — | — |
