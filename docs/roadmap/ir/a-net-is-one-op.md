---
date: 2026-11-03
title: a net is one op of the field IR, and its weights are data
files: [include/spatium/io/field.hpp, include/spatium/io/field_pod.hpp, include/spatium/io/field_glsl.hpp, include/spatium/io/field_text.hpp, include/spatium/io/plot_glsl.hpp, rsc/spatium/respond.hpp, tests/test_field_net.cpp, tests/test_field_pod.cpp, site/wasm/plot_check.html, site/wasm/check.mjs, site/wasm/questions.txt]
open:
  - a net reads one to three inputs, the arity of an op; a wider input needs a vector the IR does not have
  - a Vulkan kernel binds no buffer for the weights yet; there a net is refused by name
  - one hidden layer, as TanhNet has; a deeper net is a net of nets' outputs, not yet a consumer
  - no document line writes a net; a page draws one through plot_request's `net`
---

A net of one hidden tanh layer, `TanhNet<F, H>`, is now one op of the field IR, `Op::Net`. Its weights are numbers it reads, in the net's own order (`net<F, H>(tanh_net.weights(), inputs)`). Lowered, they are a buffer (`PodField::weights`), and in a WebGL2 module a `uniform vec4 field_weights[]` read by one function, `field_net`. A net trained further is therefore new numbers for the same program, and its shader is not compiled again. Written out as scalar ops, each weight would be a constant in the code. Nothing folds a net of constants, for the same reason.

What the field computes is what the net computes, to the bit: the interpreter follows `TanhNet::forward` op for op, and the plain-data interpreter shares the same function (`Field::net_at`). The fp32 twin is within 1e-6 of it. Two nets of one shape emit the same shader text and different weight buffers; a net used twice is one block of weights; two nets that differ in one weight by 1e-12 are not merged.

The consumer is the graph. `plot_request` takes `{"net": {"hidden": H, "weights": [...]}}` and draws y = net(x) with the weights beside the shader. In the browser check, one program draws two nets, eight units each, with only `field_weights` set between them. One shader is compiled for both, the worst difference from the core's fp32 twin is 1e-5 of a pixel, and a frame drawn with one net's weights against the other's twin is caught.

Graph numbers are fp32, so the host-against-WebAssembly check now compares them to four float steps of their scale. musl's tanhf and glibc's differ in a float's last bit: 0.301456094 against 0.301456034 at scale 1.39. Other questions keep 1e-9.
