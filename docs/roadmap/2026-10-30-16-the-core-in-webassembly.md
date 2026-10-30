---
date: 2026-10-30
title: the core the page calls, in WebAssembly, answering as the host does
files: [site/wasm/core.cpp, site/wasm/build.sh, site/wasm/check.sh, site/wasm/check.mjs, site/wasm/questions.txt, include/spatium/io/field.hpp, flake.nix, .github/workflows/ci.yml]
open:
  - the page itself is not written yet
  - the build is 545 KB of wasm at -O2, not measured against a budget
---

`site/wasm/core.cpp` is what the page calls: `spatium_respond(text)` returns the library's answer as the JSON the page reads, through the same door the gallery is written by, and `spatium_preview(text)` the question as LaTeX while it is typed, or the parse error with its column. `build.sh` builds it with Emscripten 5.0.6 (`nix develop .#wasm`, a new opt-in shell with Emscripten and node) to 545 KB of wasm and 14 KB of JavaScript.

`check.sh` builds the same file for the host too and asks both the 21 lines of `questions.txt` -- arithmetic, equations with their steps and checks, questions answered by synthesis and by Chebyshev polynomials, refusals, previews and parse errors: the answers must be the same JSON, numbers to 1e-9 of their size and a row of coefficients at the scale of its largest (WebAssembly links musl's libm, whose last digits differ from glibc's), and the certificate must say it was computed in wasm. All 21 agree; it is a CI job.

The field IR held its motion closures in `std::move_only_function`, which libc++ does not have, so nothing that included `io/field.hpp` built against libc++ -- WebAssembly, and clang hosts that use it. It is now `std::move_only_function` where the standard library has one, and otherwise an owner of the same shape (a callable held by `unique_ptr`, moved, called const).
