---
date: 2026-10-30
title: a scene names its numbers, an equation among them, and reads as LaTeX
files: [include/spatium/algebra/roots.hpp, include/spatium/io/build.hpp, include/spatium/io/field.hpp, include/spatium/io/trace_text.hpp, rsc/spatium/query.hpp, tests/test_scene_symbols.cpp]
resolves: [2026-10-29-12#2]
open:
  - a factory keeps its parameters as values, so the description prints torus(2, 0.5) where the scene said torus(R, r)
  - a definition is solved once, when it is made; changing a name means building the scene again, and a slot a device reads (a slider) is still the next step
  - an equation is one unknown on a bracket; a system, or a name defined as a function of time, is not taken
  - the description has no motion, colour or emission beyond saying a node moves
---

`Trace::let(name, value)` names a number of a scene, and `let(name, equation, lo, hi)` names the root of an equation in the names defined before it, solved when it is made by `find_root` (`algebra/roots.hpp`: bisection to two adjacent doubles, refused where the sign does not change or a value is not finite). Each returns a `Named`, a number where a factory wants one and `field()` where a field reads it by name. Every builder that takes a field -- the thickness of an offset, a motion, a colour, an emission -- binds the defined names on the spot, and a field reading a name nothing defined is refused there, by name, so no symbol reaches `cook()` or `materialize()` to evaluate to NaN on the way to a pixel. A motion's factors are bound and the closures that evaluate them rebuilt from the bound fields (`VecField::bind`, `free_symbols`). The query's root (`rsc/spatium/query.hpp`) is now `find_root`, the same door.

The consumer is a torus whose shell holds as much as its tube: R = 2, r = 0.5, and h the root of (r + h)^2 - 2 r^2 on [0, 1]. The root is r (sqrt 2 - 1) within 4 ulp, every vertex of the shell is at r + h from the tube's circle within 1e-12, and the node keeps its thickness written as `h` beside the bound one. A sphere scaled by a defined name cooks to one object of that scale; scaled by an undefined one, or coloured by one, it is refused at `moving()` and `colored()`.

`to_latex(trace, root)` (`io/trace_text.hpp`) writes the scene for a page: the definitions (r = 0.5; the equation as written, its bracket and the root it came to), then the nodes the root reaches, `S_{1} = offset(S_{0}, h)`, in an `aligned` block.
