---
date: 2026-11-27
title: headless compute is a library of its own and links no window
files: [CMakeLists.txt, cmake/SpatiumConfig.cmake.in, src/viewer/compute.cpp, src/viewer/present.cpp, src/viewer/vk_check.hpp, include/spatium/viewer/compute.hpp, tests/test_compute.cpp, tests/CMakeLists.txt, docs/gpu.md]
open:
  - CI builds no Vulkan target, so the device test runs where a GPU is, not in CI; a software device (lavapipe) in the Nix job would run it there
  - the matrix product in the test is the test's own kernel; a library kernel waits for its consumer
---

`Context`, `Buffer`, `Kernel` and the device BVH were compiled into the viewer library, which links GLFW and, with `IMGUI_DIR`, ImGui. A program that only computes had to bring a windowing library it never opens. They are now `Spatium::compute` (`spatium_compute`), which needs Vulkan and shaderc only and is built wherever those two are found, with or without GLFW. `Spatium::viewer` links it and adds the half that needs a window: the constructor of a `Context` that can present, and the `Presenter` (`src/viewer/present.cpp`). Linking the compute library pulls no GLFW: `ldd libspatium_compute.so` names none, while `libspatium_viewer.so` names it. An installed Spatium offers `find_package(Spatium COMPONENTS compute)`.

The context now takes the instance extensions and a function that makes the surface, instead of a `GLFWwindow`, so the device code does not know which window library made the surface.

`tests/test_compute.cpp` is the first test that runs a kernel on a device. It multiplies two float matrices on the device and on the host, and holds their difference to the bound a sum of k float products can drift by in any order of rounding, gamma_k times the sum of |a||b| (Higham, 3.1), entry by entry. The sizes (33 x 17 x 65, and 1 x 1 x 1000) are not multiples of the 16 x 16 workgroup, so the edge invocations that must write nothing are exercised. A kernel that reads B transposed is compiled beside the right one and must fail the same bound, and does. On the Iris Xe all three checks pass. Where no device is found, the cases are skipped, not passed.

`donut_live` and `blackhole_live`, the existing consumers of the headless half, give the same comparisons as before: shader and host differ by 7 of 691 200 pixels on the donut, and the black-hole device by 1.27e-6 of path length from the host in float.
