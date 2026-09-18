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

### What each target needs on top of the IR

**CUDA and WASM need the same thing first:** a POD mirror of
`FieldOp`/`VecFieldOp` with no `std::function` and no `std::type_index`,
plus a switch-based loop reproducing `eval_into`. That loop is already the
shape a SIMT lane wants — one linear pass, no recursion, no branching on
depth — so this is a port, not a design.

**CUDA** then uploads the op array, `Cooked<T>`'s object and shape arrays,
the mesh buffers and the two side tables as device buffers, and compiles
the interpreter as `__device__` code. Byte buffers cross, not live objects.
`Shape<T>::exact` is a `std::any`, is diagnostic only, and is dropped at
the boundary.

**WASM** is the same interpreter under Emscripten, reading the same buffer
out of `WebAssembly.Memory`.

**JS** needs no separate IR at all: the op array is index-addressed and
small — 107 fields for this whole scene — so it serialises directly to a
typed array. A JS caller either runs a short native interpreter for tooling
or hands the identical buffer to the WASM module.

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
