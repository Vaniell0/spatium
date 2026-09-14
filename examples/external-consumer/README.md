# External Consumer Example

A minimal, completely standalone CMake project showing how to pull Spatium
into your own project with `FetchContent` -- no Nix, no manual
`cmake --install`, no `find_package(Spatium REQUIRED PATHS ...)` pointed at
a hand-built scratch prefix.

This directory is **not** part of Spatium's own build: it is not listed in
`examples/CMakeLists.txt` and does not use `spatium_add_example()`. It is
its own independent CMake project that happens to live inside this repo as
a worked example. Copy `CMakeLists.txt` and `main.cpp` into your own
project and the `FetchContent_Declare` call pulls in the real Spatium
repository the same way.

## Requirements

- CMake >= 3.28
- GCC >= 15 or Clang >= 19 (C++23)
- A CMake generator (Ninja recommended, Make also works)

Nothing else. No Nix, no Boost, no Eigen, no Catch2, no Vulkan, no Google
Benchmark. `spatium/core.hpp` used to include Boost.Multiprecision
unconditionally, which made Boost headers a hard requirement even here;
since `SPATIUM_BOOST` exists that is opt-in, and this example is built in
CI on a runner with no optional package installed at all.
Pulling Spatium in via `FetchContent` (rather than configuring Spatium's
own repository standalone) builds only the header-only `Spatium::sdk`
interface target -- Spatium's viewer, examples, tests, and RSC tools all
default OFF when Spatium is not the top-level CMake project.

## Build and run

```bash
mkdir build && cmake -B build && cmake --build build && ./build/sphere_distance_demo
```

The first configure clones Spatium (`main` branch) into
`build/_deps/spatium-src`; later configures reuse that checkout instead of
re-cloning.

Expected output:

```
Sphere<2> geodesic distance (north pole -> equator): 1.570796
Expected: pi/2 = 1.570796
```
