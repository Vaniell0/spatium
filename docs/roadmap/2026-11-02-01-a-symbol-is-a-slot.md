---
date: 2026-11-02
title: a symbol lowers to a slot, and the graph of a function is a WebGL2 shader
files: [include/spatium/io/field_pod.hpp, include/spatium/io/field_glsl.hpp, include/spatium/io/plot_glsl.hpp, include/spatium/render/gpu/gpu_instances.hpp, rsc/spatium/respond.hpp, site/wasm/core.cpp, site/wasm/plot_check.mjs, site/wasm/plot_check.html, tests/test_field_symbols.cpp, tests/test_rsc_evaluate.cpp]
open:
  - a Vulkan kernel binds no slot yet; there a symbol is refused by name and bound to a number first
  - a WebGL2 field takes no noise, gather, instance, point or coordinate; each is refused by name until a consumer brings it
  - a graph is drawn of a value; a root's or an integral's synthesised answer is not drawn yet
---

A named symbol no longer stops a field at the lowering. `lower()` gives each name a slot (`PodField::symbols`), and the interpreter takes the slots' values: the field lowered once and run with a value in its slot gives the bound field's value to the bit, in double and, through `narrowed`, in fp32. A run with no value for a slot reads NaN, as the unbound field does.

So structure is code and numbers are data. A GLSL module made for WebGL2 (`glsl_es_module`, `emit_scalar_es`, GLSL ES 3.00) reads its slots from `uniform float field_symbols[N]`, and moving a parameter sets a uniform instead of compiling a shader. A Vulkan module binds no slot and refuses a symbol by name, as before.

The first consumer is the graph of a function (`io/plot_glsl.hpp`). y = f(x) is drawn where F = f(x) - y is near zero, at |F| / fwidth(F) pixels from the curve. The line keeps one width however steep the curve is, with no grid of points and no derivative of f. `plot_request`, and `spatium_plot` in the WebAssembly core, draws a value of a document along one parameter. The other parameters become slots, returned with their intervals for sliders.

The check runs in a browser. `site/wasm/plot_check.mjs` opens a page in headless Chromium, on its software WebGL2. The shader writes F into a float target, `readPixels` reads it back, and the core's fp32 twin computes F at every pixel's centre for comparison. On three functions and eight frames the shader is compiled three times, and the worst difference is 0.006 of a pixel (sin; 0.004 for cos, sqrt and log; 8e-6 for exp), against a tolerance of a sixteenth. A frame drawn with its slot a hundredth off the twin's differs by 0.16 of a pixel, and the check catches it.
