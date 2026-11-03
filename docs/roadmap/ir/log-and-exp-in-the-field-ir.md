---
date: 2026-10-29
title: Log and Exp in the field IR, and what the logarithm buys the synthesis next to the critical impact parameter
files: [include/spatium/io/field.hpp, include/spatium/io/field_pod.hpp, include/spatium/io/field_glsl.hpp, tests/test_field_pod.cpp, rsc/spatium/synthesis.hpp, tests/test_rsc_synthesis.cpp]
open:
  - GLSL's log and exp are checked as emitted text here; their numbers on a device are not compared with the host's in this change, and fp32 on a device will differ from fp64 on the host
  - a metric field still refuses Log and Exp; no metric has asked for them
  - the net of tanh units is not yet lowered through Exp
---

`Op::Log` and `Op::Exp` are appended after `Coord`, since an op's number is its code in lowered data and inserting one would renumber what is saved. The field's interpreter, the POD interpreter and the GLSL emitter carry them; the test that every op interprets bit for bit now includes both, and a field of `log` and `exp` lowers to GLSL's own and interprets as the field does, bit for bit.

The consumer is the synthesis. Next to the critical impact parameter the deflection of light is -log(b / b_c - 1) and a constant to leading order; in the job's own variable y = b / b_c - 1 from 1e-6 to 1e-2 it takes, without and with the logarithm among the operations: at 1e-4, eight pieces of 120 nodes against one of 14; at 1e-6, 24 of 537 against four of 75; at 1e-8, 62 of 1555 against six of 161. The exponential has its consumer waiting in a net of tanh units, tanh x = 1 - 2 / (exp(2x) + 1).
