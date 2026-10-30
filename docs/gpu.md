# The GPU path

One path: Vulkan compute, with shaders written as GLSL strings and compiled at
run time. A device sees exactly what the host computes, and each shader has a
host twin it is checked against. The CUDA kernels that came first (`gpu/`,
hand-transcribed Christoffel symbols for one Kerr metric) are removed; the
Vulkan path renders the same holes from the metric's IR with exact
derivatives and runs on any device with a Vulkan driver, a rented NVIDIA card
included.

## What is where

| Layer | Files | What it is |
|---|---|---|
| A field as data | `io/field_pod.hpp` | a `Field` lowered to op codes, child indices and numbers, and the interpreter that runs it; it agrees with `eval_into` bit for bit |
| A field as GLSL | `io/field_glsl.hpp` | the same field as straight-line GLSL, one statement per op (`emit_scalar`, `emit_vector`, `emit_placement`) |
| A metric as GLSL | `physics/relativity/metric_glsl.hpp` | a `MetricField` with its partials, and the geodesic step built on it |
| Layouts | `render/gpu/gpu_types.hpp` | the fp32 std430 structs a shader reads a scene through, sizes asserted |
| A scene for a device | `render/gpu/gpu_scene.hpp` | a cooked scene (`render/cooked_scene.hpp`) packed into those arrays, and the traversal again on the host in the shader's precision |
| Trees over moving instances | `render/gpu/lbvh.hpp`, `render/gpu/gpu_lbvh_glsl.hpp` | a linear BVH (Karras 2012) on the host, and its build as compute kernels with the sort left on the host |
| Instances moved on the device | `render/gpu/gpu_instances.hpp` | a Scatter's motion as a kernel generated from its fields |
| Shaders | `render/gpu/gpu_trace_glsl.hpp`, `render/gpu/gpu_splat_glsl.hpp`, `render/gpu/blackhole_glsl.hpp` | the scene's trace and shade, sub-pixel particles drawn by projection, and a spacetime traced ray by ray through its metric |
| Device | `viewer/compute.hpp`, `viewer/gpu_lbvh.hpp` (`src/viewer/`) | headless Vulkan compute: a device, host-visible storage buffers, a kernel compiled from GLSL; the LBVH build on it |

## Generated from the IR, and written by hand

Generated: a field as GLSL (`field_glsl`), a metric with its partials and
the geodesic step (`metric_glsl`), and a Scatter's motion as a kernel
(`gpu_instances`). Written by hand, each a transcription of the host
function it is named after: the scene's trace and shade (`gpu_trace_glsl`),
the splatting of sub-pixel particles (`gpu_splat_glsl`), and the black-hole
frame (`blackhole_glsl`), whose loop is written by hand around the metric
`metric_glsl` generates. A new scene reaches the device through the
generated half; the hand-written half is the frame around it.

Numbers in a generated shader are constants today: the lowering refuses a
symbol, so every name is bound before a field reaches a device and a moved
number is a new shader. The rule the next step follows is in
[`gpu-abi-design.md`](gpu-abi-design.md): structure is code, numbers are
data.

`render/` itself is the CPU ray tracer (camera, supersampling, sky, textures,
PNG); `render/gpu/` is what a device reads. `viewer/app.hpp` is the
rasterising window, kept apart from compute so that the compute path runs on
a machine with no display.

## How it is checked

Every shader is a transcription of a host function, which keeps the name of
the function it transcribes, and the two are compared:

- `field_pod` against `eval_into`: bit for bit over a real scene
  (`tests/test_field_pod.cpp`).
- `donut_live` prints three comparisons, each isolating one layer: fp64 trees
  against the fp32 host trace (what fp32 costs), the fp32 host trace against
  the shader (whether the shader is the spec), and frame time.
- `blackhole_live --shadow-check` measures the Schwarzschild shadow from the
  pixels a horizon ended, against sin(alpha) = 3 sqrt(3) M / D sqrt(1 - 2M/D).
- The LBVH host build is held by `tests/test_lbvh.cpp`.

## Precision and devices

Shaders compute in fp32. The device this is written on (Iris Xe) has no
fp64, and fp32 passes the shadow check at 3.5e-4 at the critical curve.
Buffers are host-visible and coherent: on an integrated GPU device memory is
host memory, so there is no staging copy.

## Programs

- `donut_live`: the donut scene traced on the device; headless comparisons by
  default, `--live` for a window, `--photo` for a frame.
- `blackhole_live`: one hole or a pair described in the scene language;
  `--live`, `--video DIR`, `--frame PATH`, `--shadow-check`.

Both need the viewer (`SPATIUM_BUILD_VIEWER`) and a Vulkan driver.
