# GPU / C-ABI export — design direction

Two related but distinct problems live under this one heading, and they
deserve different answers, not one hand-wave:

1. **`gpu/` hand-transcribes physics it should share.** This is the real,
   current problem — `schwarzschild_render_kernel.cu` re-derives the same
   geodesic/disk math `examples/blackhole_gr_demo.cpp` already has,
   independently, in a separate file. That's not a standard mechanism,
   it's a one-off port, and it's exactly what creates the desync risk
   `gpu/verify_cuda_render.cpp` exists to catch after the fact instead of
   preventing structurally. This one has a concrete fix, below, and
   doesn't need to wait on anything else.
2. **Getting a stable C ABI out of the generic DSL for WASM/Mojo/JS.**
   This is the more speculative, cross-language problem — still real,
   still worth a target shape, but genuinely gated by `rsc/`'s own stated
   WASM timeline, not something to build ahead of a consumer.

## Problem 1: the real gap is between the derivation script and the header

**Correction to an earlier draft of this memo**, which proposed sharing
the generic `Dual<T>`-based templated headers directly via
`__host__ __device__` and treated that as the fix. That contradicts a
decision already made, for real reasons, before this direction was
written: `Dual<T>`'s general forward-mode autodiff machinery is
branch/indirection-heavy in a way that's specifically bad for SIMT
execution (every lane in a warp needs to take the same path for full
throughput; autodiff's general architecture doesn't guarantee that).
The already-reasoned choice was to hand-derive closed-form Christoffel
symbols for exactly the two metrics the GPU path needs (Schwarzschild,
Kerr) instead of porting the generic engine — a deliberate CPU/GPU
architecture split, not an oversight. Sharing `__host__ __device__`
headers as a blanket fix is wrong for this specific pipeline; withdrawn.

**The real, concrete gap is one step earlier and much narrower than a
generic C-ABI system: `gpu/derive_christoffel.py` already exists and
already does the hard part.** It's a real sympy script that (1) derives
`Gamma^lambda_{mu nu}` symbolically straight from the exact metric forms
in `schwarzschild.hpp`/`kerr.hpp`, (2) self-checks that Kerr at `a=0`
reduces term-by-term to Schwarzschild, (3) prints each nonzero component
in the form transcribed into `christoffel_closed_form.hpp`. What it
doesn't do — and this is the actual "automate it" gap, not an abstract
one — is *emit* that header. A human reads the printed sympy output and
hand-copies it into the `.hpp` file; correctness against the CPU engine
is then checked separately, numerically, in `tests/test_relativity.cpp`.
Two manual steps (transcription, and trusting the transcription was
faithful) sit between a correct derivation and a correct header, and
nothing catches a transcription slip except the numeric test noticing
downstream — after the fact, the same shape of problem the CPU/GPU
per-ray verification has.

The concrete fix: make `derive_christoffel.py` print valid C++
expressions (sympy's own `ccode()`/`cxxcode()` printer, not `str()`) and
write `christoffel_closed_form.hpp` directly instead of printing for a
human to retype — closing the loop from "metric defined in `schwarzschild
.hpp`/`kerr.hpp`" to "correct GPU-side closed-form code" without a
hand-transcription step in between. This is *the* concrete instance of
"the engine automates creation of these calculations, inlining
everything" — not a hypothetical, an existing script one code-emission
step away from doing it. Extending the same script to a new metric (if
one ever gets added) becomes "add the metric, run the script," not
"add the metric, do the tensor algebra by hand, hope the transcription
is right." This doesn't need CUDA-side template sharing at all — the
symbolic math and the code emission both happen at derivation time, on
the host, in Python; the emitted header stays exactly what it is today
architecturally (hand-optimized closed-form CUDA code), just no longer
hand-*transcribed*.

Not built yet — the render currently in flight uses today's hand-
transcribed header, and this doesn't change or need to change for it.
It's the concrete next step for whenever the Christoffel formulas change
again (a new metric, a correction) rather than a general system to build
ahead of that need.

## Problem 2: the C-ABI boundary itself, for WASM/Mojo/JS — still deferred

Once the physics is shared instead of duplicated, a separate question
remains: a `spatium::` algorithm is templated/concept-constrained, and a
C ABI is monomorphic POD-only — crossing that boundary always means a
specific, already-instantiated `T`/`N`, chosen up front. The answer
here, sketched but not built: a manifest — a small, explicit list living
next to the algorithm it exports — naming which concrete instantiations
get a C-ABI entry point, from which a generated (not hand-written)
`extern "C"` shim and flat-POD struct are produced, with failure
signaled as a status code or sentinel (the same shape `render_kernel.h`'s
`captured_flag` already uses informally, since `Result<T>` itself can't
cross a C boundary). This turns "write a new GPU/WASM entry point" into
"add a line to a manifest and regenerate" — the same relationship
`scripts/gen_dependency_graph.py` now has to `docs/dependency-graph.dot`.

Once a C ABI exists this way, WASM itself is the easy part — Emscripten
consumes a C ABI directly, and a thin TS wrapper on top is a solved
problem elsewhere. `rsc/README.md`'s Site section states plainly that
this is gated behind RSC "working in genuinely sufficient volume — not
before, and not something being built yet." That gate is real and this
memo doesn't override it: unlike Problem 1, there's no current consumer
forcing this, so it stays a target shape, not a task.

## Later, bigger idea: fusing a *chain*, not just one function

A harder version of Problem 2, raised alongside it: several algorithm
steps chained together (ray → geodesic integration → disk lookup →
compositing) currently compile to whatever control flow the C++ source
has — real branches (adaptive step control, CCD checks, line-search
backtracking) that don't collapse into the uniform, branch-free control
flow a GPU/SIMD lane wants. The ask is tracing a *specific* chain once,
at known concrete types, compiling *that* to a flat kernel, cached after
the first trace.

Python JIT/tensor compilers do exactly this and are worth naming as
prior art for the *model*, not as dependencies to adopt (Spatium's asset
is the C++ DSL; none of these would have anything to compile):
**JAX/XLA** (`jax.jit` traces once per input-shape signature, XLA fuses
the traced graph into one kernel, caches the executable — the closest
match to "cache after one cycle"), **Apache TVM** (explicit operator
fusion as a compiler pass plus hardware-specific auto-tuning — the
closest match to "chain → one tuned kernel"), **Triton** (a DSL for
hand-writing the fused kernel directly, JIT'd to PTX), **Numba**
(`@cuda.jit`, same trace-cache mechanics as JAX, no graph-level fusion).
What's worth borrowing is the trace-once/fuse/cache-by-concrete-shape
model, most likely as a C++ constexpr/metaprogramming mechanism, not a
dependency on any of the above. This is further out than Problem 1 or 2
and genuinely unresolved — worth having named, not worth planning a
build date for yet.
