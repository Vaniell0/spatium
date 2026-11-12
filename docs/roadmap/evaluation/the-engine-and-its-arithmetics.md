---
date: 2026-11-03
title: how a value is computed is the engine's, and the arithmetic is one of its settings
files: [include/spatium/core/rational.hpp, include/spatium/core/decimal.hpp, include/spatium/io/field_pod.hpp, include/spatium/answer/evaluate.hpp, include/spatium/answer/respond.hpp, tests/test_answer_evaluate.cpp, site/wasm/questions.txt, site/gallery/gallery.json]
open:
  - only a value with nothing left to read is computed in the engine's arithmetic; a root, an integral and a synthesis still compute in a double and say so
  - the roots of a closed form are not yet exact (5/2, or (3 ± √17)/4) where its coefficients are fractions
  - a decimal exists where the compiler has decimal64 (GCC); the WebAssembly core, built by Clang, says it has none
  - the page does not show the engine yet
---

What is asked and how it is computed are now apart. A document says what is asked. An `Engine` beside it says how it is computed: the arithmetic, whether an arithmetic that cannot compute something refuses (`strict`) or hands it to a double, and a synthesis's tolerance. Its defaults are today's behaviour, so the gallery's numbers did not move; each candidate now says the arithmetic it computed in.

The arithmetics are a double; a float, as a device computes; `Rational`, a fraction of two 64-bit integers, in which 1/3 + 1/6 is 1/2 and 0.1 * 3 is 3/10; `Decimal`, IEEE decimal64 through GCC's `std::decimal` (TR 24733 -- the C++ standard itself has no decimal type, and Clang has none at all), in which ten times 0.1 is 1; and Real50 where Boost is built in. A literal is read through its shortest decimal form, so the 0.1 that was typed is 1/10, not the binary fraction nearest it.

No interpreter was written for them. The one in `field_pod.hpp` was already a template over the scalar. Its ops past + - * / and order moved under `if constexpr`, so it compiles for an arithmetic with no square root and computes only the ops `exact_op` names; `beyond<T>` refuses the first other op by name before it runs ("rational has no Sqrt in it"). `narrowed`, the fp32 twin, became one case of `converted<T>`. A fraction that does not fit in 64 bits, or divides by zero, is spoiled rather than wrapped, and the evaluation says so.

Each arithmetic is held against the double: as a double, a fraction or a decimal is the double's value within 4e-16 relative, and a float within 2e-7. The cases they must fail are that 0.1 * 3 is not 0.3 in a double and is not the double's in a float. The fp32 and double interpreters still agree bit for bit with the field. The WebAssembly core and the host agree on 36 of 36 questions, three of them with an engine.
