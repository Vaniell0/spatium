# GPU / C-ABI export — design direction

Two problems live under this heading. They are still distinct and still
deserve different answers, but the second one has moved a long way since
this was first written, and most of what follows is new.

1. **`gpu/` hand-transcribes physics it should share.** Narrow, concrete,
   unchanged, and not blocked on anything. Below.
2. **Getting the scene description out of this process** — onto a GPU, into
   a browser, in front of JavaScript. This was written as speculative and
   gated on `rsc/`. Both of those were wrong. It is not speculative, the
   gate named was the wrong gate, and roughly half of what the earlier
   draft sketched as unbuilt is in the tree today.

---

## Problem 1: the gap between the derivation script and the header

**A correction to an earlier draft** which proposed sharing the generic
`Dual<T>` headers via `__host__ __device__` and treating that as the fix.
That contradicts a decision made earlier for real reasons: `Dual<T>`'s
forward-mode machinery is branch- and indirection-heavy in a way that is
specifically bad for SIMT, where every lane in a warp wants the same path.
The existing choice — hand-derived closed-form Christoffel symbols for the
two metrics the GPU path needs — is a deliberate CPU/GPU split, not an
oversight. Withdrawn, and it stays withdrawn.

**The real gap is one step earlier and much narrower.**
`gpu/derive_christoffel.py` already does the hard part: it derives
`Γ^λ_{μν}` symbolically from the exact metric forms in `schwarzschild.hpp`
and `kerr.hpp`, self-checks that Kerr at `a = 0` reduces term by term to
Schwarzschild, and prints every nonzero component in the form that appears
in `christoffel_closed_form.hpp`. What it does not do is *emit* the header.
A human reads the printed output and retypes it, and nothing catches a
transcription slip except a numeric test noticing downstream.

The fix is to make the script print valid C++ (sympy's `ccode()`/
`cxxcode()`, not `str()`) and write the header itself. Adding a metric then
becomes "add the metric, run the script" instead of "add the metric, do the
tensor algebra by hand, hope the retyping was faithful". No CUDA-side
template sharing is involved: the symbolic work and the emission both
happen at derivation time, on the host, in Python, and the emitted header
stays architecturally what it is today.

Not built. It is the next step whenever the formulas change again, rather
than a system to build ahead of that.

---

## Problem 2: getting a scene out of this process

### The blocker named here before was the wrong one

The earlier draft deferred this behind `rsc/README.md`'s statement that its
Site feature waits on RSC "working in genuinely sufficient volume". That
conflates two different kinds of blocker. "No consumer is asking yet" is a
scheduling fact about one feature. "The artifact is not structurally
exportable yet" is a technical fact about the DSL, and it is the one that
decides whether any of this is possible. Only the second is measurable, and
this memo never measured it.

It is measurable, by `field_report()`, which walks a trace and counts what
it finds. On the donut demo, on `main`:

```
fields: 107 total, 77 structural, 30 opaque; leaves 32
        (32 recognized, 0 unknown), 11 distinct types, 3208 B captured
motions: 33 placements (shareable geometry), 0 deformations
34 trace nodes
```

A scene is exportable exactly as far as its fields are **structural** —
built from tagged ops the evaluator knows — because an opaque leaf is a C++
lambda and a lambda cannot cross an ABI boundary. So the number that
matters is the 32.

**Two cautions about that report, both learned the hard way.** It used to
print `73 fields, 66 structural, 7 leaves`, which read as "nearly all of it
is structural" and was wrong in the flattering direction: it never walked
the emission slot, and it never counted the closures a `Scale` or a
`Rotate` carries in `scale_fn`/`rot_fn`. Both are fixed and the numbers
above are the honest ones. And `unknown = 0` is weaker than it sounds — it
says every leaf carries a `type_index`, not that every leaf is a *known
operation*. A lambda that happens to be a smoothstep and a lambda that is
arbitrary are indistinguishable by that number.

### The open question, stated so it can be answered rather than assumed

**Do those 32 leaves fall into a small closed vocabulary?** A first pass
over the demo's fields suggests they do — noise, a hash from `origin`,
smoothstep, trig, lerp, clamp, dot/cross, a select/step, a normalize, and
an indexed gather from a small per-scene constant table (the BOOM
letterform points, the noise permutation tables) — with nothing in the list
doing I/O, iterating a solver, or touching mutable external state.

**That pass has not been verified leaf by leaf and should not be quoted
until it has been.** It was made by reading call sites, and reading call
sites and walking the trace disagree by one, which is itself unresolved:
either a lambda was missed, or "opaque leaf" is wider than "lambda". The
work is to walk the 32 and name each with one word from the vocabulary.
Eight matching and twenty-four not is a completely different conclusion
from the reverse, and only one of them supports the word "exportable".

### The IR, if the answer is yes

Very little needs inventing, which is the encouraging part.

**Value types:** `Scalar<T>` and `Vec3<T>`. Nothing else is needed —
rotations enter as an axis and an angle and are built by ops, rather than
being first-class IR values.

**Scalar ops**, extending today's `Op`: `Const, U, V, T, Add, Sub, Mul,
Div, Neg, Min, Max, Clamp, Lerp, Smoothstep, Select, Sin, Cos, Sqrt,
Noise3(table, x, y, z), Hash(salt, v) → [0,1)`.

**Vector ops**, extending today's `VecOp`: `Point, Origin, ConstVec, Add,
Sub, Scale(scalar), Dot, Cross, Normalize, Gather(table, index),
RotateByAxisAngle(axis, angle)`.

**Representation:** exactly today's `FieldOp`/`VecFieldOp` array — flat,
topologically ordered, children before parents — with the erased callable
and its `type_index` replaced by a closed op tag and operand indices. That
is a deletion, not a redesign.

**`MotionEnv` crosses unchanged.** Its three named fields — `p`, `t`,
`origin` — become three fixed input registers rather than being baked into
op signatures, which is precisely the hedge the type was given a name for.
`origin` is per *instance*, one value per object, which is already the
shape a GPU wants per thread.

**Side tables:** two kinds of small flat constant buffer travel with the op
array — noise permutation tables (512 B each) and point lookup tables. Both
are static per render and neither depends on scene structure.

### What each target needs — and they do not need the same thing

Earlier drafts wrote "CUDA + WASM + JS" as one item. That grouping is
wrong and it hides the most useful fact in this document: **the browser
does not need any of the export work.**

**CUDA needs the IR.** Device code cannot run `std::function`, so a field
holding an erased callable cannot be evaluated on a GPU at all. This is
what the POD mirror of `FieldOp`/`VecFieldOp` and the switch-based
`eval_into` loop are for. CUDA then uploads that op array, `Cooked<T>`'s
object and shape arrays, the mesh buffers and the two side tables as device
buffers, and compiles the interpreter as `__device__` code. Byte buffers
cross, not live objects; `Shape<T>::exact` is a `std::any`, is diagnostic
only, and is dropped at the boundary.

**WASM needs none of it.** WASM is compiled C++. Emscripten handles
`std::function`, RTTI and lambdas, so an opaque leaf in the browser is just
a compiled lambda that runs. Structurality is a constraint imposed by
device code and by serialisation — writing a scene to a file or sending it
over a network — and the browser is neither. Nothing about the vocabulary
gates a WASM build.

That makes the browser the *nearest* target rather than the furthest, which
is the opposite of how this document used to read.

**JS** is where the earlier version of this section was wrong, and the
error is worth keeping visible because it is easy to make again. It said JS
"needs no separate IR — the op array serialises directly to a typed array".
That describes a *consumer* of the buffer. It says nothing about where the
buffer comes from, and if the answer is "from C++", then the JS side is a
viewer and not a DSL.

### Nothing turns a closure into IR. Not in JS, and not in C++ either

This is the part to be clear about, because the C++ side can look like
magic from the outside and isn't.

In C++ you can write a motion two ways. As a closure:

```cpp
.moving([](const Vec<double,3>& p, double t) { return p * (1.0 - smoothstep(t)); })
```

or structurally, out of combinators — `point()`, `scaled`, `rotated`,
`operator+` — which build the op array by construction. The first is an
opaque leaf and does not cross an ABI. The second does. Nothing traces the
first into the second; `is_structural()` exists precisely because that
conversion does not happen.

So the C++ DSL is already the "write structural expressions" model, with a
soft refusal rather than a hard one: the closure is legal, it simply does
not export. **A JS binding should mirror exactly that, and then it is not a
special case at all** — the same combinators, the same legal escape hatch
that runs locally and does not serialise.

The objection to that model is that an expression builder is an IR
constructor with a human face rather than a DSL. It is a fair objection and
the C++ side answers it by naming things well — `scaled(point(), ...)`
reads fine. Where it stops reading fine is exactly where the vocabulary has
a hole, because then there is nothing to reach for but a lambda.

### The hole, found while writing this — closed 2026-09-23

**There was no structural way to express a scale or rotation factor.** Both
`scaled` overloads took a callable, both `rotated` overloads took a
callable, and no overload took a `Field<T>`. So even the structural path
bottomed out in a closure the moment anything changed over time — which is
why the donut's report showed the factor closures it did, and why
`scaled`/`rotated` had to be taught to declare themselves opaque at all.

The fix proposed here was one overload family:
`scaled(point(), smoothstep(Field::t()))`. That is now what it is, spelling
included — which is worth noting, because the spelling was written down
before anything existed to spell.

What the proposal did not see is that the family alone converts almost
nothing. The factors it was aimed at are smoothsteps, a smoothstep is a
clamp followed by a cubic, and the vocabulary had no `Min` and no `Max`:
the cubic was expressible and the clamp was not. With the overloads and
without the two ops, exactly one leaf in the scene converts. With them,
the donut goes from 32 opaque leaves to 10 and from 77 structural fields
of 107 to 98, with a byte-identical frame.

So this section's claim — that the missing factor form is a precondition
for the export question rather than a nicety — held. Its implied estimate
of what the fix costs did not, and the gap was the vocabulary rather than
the signature. The remaining ten leaves are the subject of "name each leaf
against a vocabulary", and two of them are deliberate: a hard step needs a
comparison, and the dust's tumble axis is an integer hash that no
arithmetic vocabulary expresses.

### If a JS binding is built: handles, not an IR

Deprioritised 2026-09-18 — the browser is reachable and stays reachable,
and the description-based door below is the cheaper boundary to build
first. Recorded here so the design does not have to be re-derived when it
comes back round.

Since WASM runs the real library, JS does not have to build an op array at
all, and therefore does not have to solve "how does a JS closure become
IR". **The trace lives in WASM memory and JS holds integer handles into
it.**

```js
const donut = scene.torus(1.0, 0.4);        // -> 7   index into the C++ trace
const iced  = scene.offset(donut, 0.08);    // -> 8
scene.scatter(iced, sprinkle, 11000, 42);   // -> 9
```

Each call reaches an `extern "C"` entry point, appends to the same `Trace`
the C++ API appends to, and returns an integer. Only integers and doubles
cross. There is no second implementation of anything: the trace is built by
the same code that has always built it, so the JS API matches the C++ one
because it *is* the C++ one. It also disposes of a problem that would
otherwise need solving — `VecField` is move-only, and ownership simply
never leaves C++.

**The cost is a wrapper per builder**, which is mechanical but numerous and
should be generated from a list rather than typed. This is where the old
"manifest" idea finally earns its place: it was a solution without a
problem when it was about exporting algorithms, and it is the right shape
here. `T` is pinned to `double` at the boundary.

**The limit is that a JS lambda cannot be passed.** A callback from WASM
into JS is technically possible and would be ruinous per vertex, so in
practice the browser has no escape hatch and the vocabulary *is* the
language there. That makes the missing structural factor above hurt more in
JS than in C++, where a lambda is at least available.

**And performance, without illusions:** WASM runs perhaps 1.5–3x slower
than native, threads need `SharedArrayBuffer` with COOP/COEP headers, and
the renderer here is CPU ray tracing. A two-million-instance frame would
take minutes in a browser. A live browser demo is a small scene, not this
one.

### Two alternatives, named so they are chosen rather than drifted into

**Trace, the way JAX does** — call the user's function with proxy objects
so `p.scale(...)` records into a graph instead of computing. Ergonomic and
real, and *not JS-specific*: the C++ equivalent is a symbolic scalar
substituted for `T`, which is the roadmap's "IR evaluable in `Dual<T>`"
note. If tracing is built it should be built for both rather than bolted
onto one.

**Let JS construct, cook and render alone.** Then WASM is not needed, and
the JS thing is a reimplementation rather than a binding, sharing only a
scene format. That is a legitimate architecture — roughly what glTF is —
and it is a different project.

## The other door: a description, not an ABI

Everything above assumes the boundary is a *binding* — a caller in another
language holding handles or buffers and issuing calls. There is a second
shape for the same boundary, it is cheaper, and half of it is already
built: **the outside sends a description, and one entry point reads it.**

One door instead of a wrapper per builder. No handle protocol, no ABI
surface that grows every time the DSL gains a verb.

### The precedent is `io/scene.hpp`, and it already works

A scene loads from JSON today, and `ShapeRegistry` maps a `shape_kind`
string to a factory. The property worth copying is the one its header
comment leads with: **an object with an unregistered `shape_kind` still
loads and saves correctly**, carrying its string and its raw parameters
through untouched. The format does not break on what it does not know.

That covers *placed objects*. Extending the same idea from placed objects
to the **trace** — `torus`, `offset`, `scatter`, `moving`, as records
rather than as calls — is the whole job, and it needs no IR, no vocabulary
and no structural factor, because the description is read on the C++ side
and turns straight into the builders that already exist.

### But a dispatcher is not a parser

RSC is the obvious thing to reach for here and it is the wrong tool for
half of the job. The two halves must stay separate:

**A precise description → a parser.** `torus(1.0, 0.4)` determines its
result completely. Parsing is deterministic, testable, and fails cleanly on
input it does not understand. A trained model put in that position is
strictly worse: non-deterministic, unverifiable, and unable to say "I don't
know" — it always selects something.

**An under-determined goal → search.** A model earns its place only where
the input genuinely does not determine the output.

### Where RSC does belong: it decides *how*, not *what*

The description says "offset this torus by 0.08". It does not say whether
to stay analytic or tessellate, at what resolution, with which root-finder,
at what precision. **That is exactly what RSC is already trained to
decide** — precision dispatch, root-finding method, mesh strategy are three
of its existing domains.

So the division is clean and both halves have machinery today:

| | decides | mechanism | state |
|---|---|---|---|
| Parser | *what* was asked | `io/scene.hpp`'s string registry, extended from objects to trace nodes | half built |
| RSC | *how* to carry it out | `rsc/include/registry.hpp`, dispatcher trained per domain | built, unapplied |

`rsc/include/registry.hpp`'s `add()` returns the dispatch-head index — the
value the model's classification head emits to select an op. That is the
far side of the same door.

This is also the honest answer to "is the IR overkill". For *this* path it
is not needed at all. The IR earns its place in exactly two places —
evaluating fields on a device, and serialising a scene to a file or a wire
— and neither of those is what a language binding needs.

### The cheap GPU layer, unchanged by any of this

Worth restating because it keeps getting buried under the IR discussion:
**uploading a cooked scene and traversing it on the GPU needs none of this
work.** `Cooked<T>` is already flat POD arrays of objects and shapes with
no pointers. The CPU cooks, the device traverses and shades. The IR becomes
necessary only when fields must be re-evaluated *on the device* every
frame, which is a real-time-animation requirement and not a rendering one —
cooking two million instances takes 2.6 s on the CPU, which is fine offline
and not fine at sixty frames a second.

### The first step is a proof, not a backend

Not CUDA. Not WASM. **A test.**

Write the POD interpreter, run it over the donut's fields alongside
`eval_into`, and assert the two agree **bit for bit** across the scene.

No epsilon. A tolerance would pass under a different order of operations,
and order of operations is exactly what a SIMT lane does not guarantee — so
the tolerance would hide the one class of failure the test exists to catch.
A failure on ordering is a finding, not noise.

If it passes, exportability is demonstrated for a real scene and can be
written down as fact. If it fails, the failure names the operation the
vocabulary is missing. Either outcome is worth more than another round of
discussion, and it needs no GPU toolchain, no Emscripten and no JavaScript
to obtain.

---

## What the earlier draft treated as distant and is in fact built

This section exists because the memo's errors were all of one kind: they
described things as far off that had been built in the meantime, which made
the remaining work look larger than it is.

**"Trace a specific chain once, at concrete types, compile that to a flat
kernel, cache it."** `Trace` *is* that trace — a flat, inspectable record
of operations over real spaces, built once and read many times — and
`cook()` *is* the compile step, deduplicating shapes, folding transforms
and expanding operations into instances in a single pass. JAX/XLA, TVM,
Triton and Numba remain worth naming as prior art for the *model*, but the
question is no longer whether to build a trace/cache design. It is what is
missing to lower `Trace` and `Cooked` to a device kernel, which is a much
smaller question.

**"A manifest naming which concrete instantiations get an entry point,
producing a flat-POD struct."** For the scene half of the problem,
`Cooked<T>` is already close to being exactly that. It is flat `Object[]`
and `Shape[]` arrays, and `Object<T>` is 152 B of `size_t`, `Vec3`, `T`,
`Quaternion`, `bool` and a POD `Material` — no pointers. The per-algorithm
numeric export the manifest idea also covered (closed-form Christoffel
symbols and the like) is a separate and still-open question; only the
scene-manifest half is stale.

**The compression that motivates all of this is measured, not projected.**
34 trace nodes describe 2,021,984 objects; geometry that would be
64,654,768 vertices if every object carried its own copy is stored as
39,272, and the frame renders in under a gigabyte. Whatever crosses the
boundary, it is the 34 nodes and the flat instance arrays — not two million
meshes.
