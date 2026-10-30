# The boundary: how a description leaves the library

What crosses out of Spatium is a description, and what comes back is an
answer or a shader with its data. The library itself stays C++. This file
was first written as the design of a C ABI and of CUDA kernels; both are
decided against (2026-10-30), and what follows is the boundary as it is
built and the rules its next steps follow. The GPU path itself is
[`gpu.md`](gpu.md).

## Decided: a description, not an ABI

**No C ABI.** A C ABI carries none of the DSL's syntax: what reads well in
C++ (`t.offset(t.torus(R, r), h)`) becomes handles and integer codes, and a
wrapper made convenient on top of it would be used less than the JSON or
than C++ itself. A DSL rewritten in another language over such an ABI would
be a second library, not a binding. A wrapper can still be added when a
caller asks for one; it would sit on the door below, not beside it.

**One door: text or JSON in, JSON out.** The outside sends what a person
writes -- a question, or a document of statements a line (`io/document.hpp`:
`r = 0.5`, `p in [0.5, 4]`, `x: x^2 - p = 0 on [0, 10]`,
`S1 = offset(S0, h)`) -- and one entry point reads it with the same reader
the host uses. A document's JSON is read back through the line reader, so
the two forms cannot mean different things. The page's core is this door
built to WebAssembly (`site/wasm/core.cpp`: `spatium_respond`,
`spatium_preview`), and CI holds it against the host build of the same
calls. A TypeScript builder writes the same JSON and builds nothing itself.

The price is stated rather than hidden: anything complicated has to be
automated inside the program that reads the description -- the core --
because the caller only sends text. That is where the library is anyway,
and it costs the syntax nothing.

**A parser decides what; search decides how.** A description such as
`torus(2, 0.5)` determines its result, so it is parsed: deterministic,
testable, refused at the column where it stops being the language. What a
description leaves open -- evaluate or synthesise, which root finder, how
fine a grid -- is chosen by the evaluation layer (`rsc/spatium/respond.hpp`
picks the route: a value, a closed form, a synthesised answer with its
certificate). A trained model belongs only on the second side.

## The IR is the trace

`Field` is a flat pool of tagged operations, children before parents, and
`Trace` a flat vector of scene nodes over it. Nothing turns a closure into
IR, in C++ or anywhere else: a field built from the combinators (`point()`,
`scaled`, `rotated`, the arithmetic operators) is structural and crosses; a
field holding a lambda is opaque and runs only where C++ runs. The DSL
keeps the lambda as a legal escape hatch that does not export.

Where the IR is needed:

- **To answer in a browser: not at all.** WebAssembly is compiled C++, so
  the core runs the real library, closures included.
- **To draw on a device: always.** A shader cannot run C++, so whatever a
  device evaluates per pixel or per frame is the IR lowered: to GLSL for
  Vulkan compute today (`io/field_glsl.hpp`, `physics/relativity/metric_glsl.hpp`,
  `render/gpu/gpu_instances.hpp`), and to GLSL for WebGL2 for the page.
  Traversing a cooked scene needs no IR (`Cooked<T>` is flat arrays with no
  pointers); evaluating its fields on the device does.
- **To write a scene to a file or a wire:** the document and its JSON.

## Decided: structure is code, numbers are data

A shader is compiled when the description's structure changes -- a line
added, an operation changed -- and not when a number in it moves. A
parameter (`Trace::choose`), time, and the weights of a network are
numbers: each lowers to a read from a uniform or a buffer, std430 storage
for Vulkan compute and uniforms or textures for WebGL2. So compilation is
automatic -- document, IR, GLSL, compiled by the driver -- and a slider or a
retrained network changes data, not code.

Today the lowering refuses a symbol (`lower()` in `io/field_pod.hpp`,
`io/field_glsl.hpp`), so every name is bound to its number before a field
reaches a device, and a moved number is a new shader. Lowering a symbol to
a uniform read is the step that changes this; the interpreter in
`field_pod` is the host twin it is checked against.

## Decided: a network is one operation of the IR

A dense layer with its activation is a single IR operation whose weights
live in a buffer. It is not expanded into scalar units: a net of tanh units
written as fields is H x F nodes with every weight a constant, and
retraining it would mean recompiling the shader. With the layer as an
operation, a neural shader -- an image filter, a learned field in a scene --
is lowered like any other field, and `TanhNet<double>` (`algebra/models/network.hpp`)
on the host is its twin.

## The vocabulary, and the proof that a field crosses

The structural operations were found by reading leaves that could not be
written without them, not by filling out a list. The first pass over the
donut demo's opaque leaves named a small vocabulary -- noise, a hash from
`origin`, smoothstep, trig, lerp, clamp, dot and cross, select, normalize,
an indexed gather from a small constant table -- with nothing doing I/O or
iterating a solver.

One hole decided most of it (closed 2026-09-23): there was no structural
way to write a scale or rotation factor, and no `Min` or `Max`, so a
smoothstep (a clamp, then a cubic) could not be written. The overloads took
the spelling written down before they existed, `scaled(point(),
smoothstep(Field::t()))`, and with `Min` and `Max` the donut went from 32
opaque leaves to 10 and from 77 structural fields of 107 to 98, with a
byte-identical frame. Two of the remaining ten are deliberate: a hard step
needs a comparison, and the dust's tumble axis is an integer hash.

**The proof is a test, not a backend.** `io/field_pod.hpp` lowers a field
to op codes, child indices and numbers, and its interpreter runs beside
`eval_into` over every field of the real scene, compared on the bits
(`tests/test_field_pod.cpp`). No epsilon: a tolerance would pass under a
different order of operations, which is exactly what a SIMT lane does not
guarantee.

## What crosses is small

34 trace nodes describe the donut's 2,021,984 objects; geometry that would
be 64,654,768 vertices with a copy per object is stored as 39,272. Whatever
crosses a boundary -- to a device or to a page -- is the description and the
flat instance arrays, not two million meshes.
