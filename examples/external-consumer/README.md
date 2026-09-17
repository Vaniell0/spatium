# External Consumer Example

A minimal, completely standalone CMake project showing the two ways an
unrelated project picks Spatium up: `FetchContent`, which needs no install
step and nothing on the system, and `find_package` against an installed
tree, which is how a system or distro package is consumed.

One project and one `main.cpp` for both. They are the same program asking
the same question of the same library, and a second copy would be a second
thing to keep in step -- the kind that drifts unnoticed, because a build
nobody runs is not a check.

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

The one-command path, pulling the sources at configure time:

```bash
cmake -B build && cmake --build build && ./build/sphere_distance_demo
```

The first configure clones Spatium (`main` branch) into
`build/_deps/spatium-src`; later configures reuse that checkout instead of
re-cloning.

Against an installed Spatium instead:

```bash
# from the Spatium checkout
cmake -B build-lib -DCMAKE_INSTALL_PREFIX=/tmp/spatium-prefix \
      -DSPATIUM_BUILD_TESTS=OFF -DSPATIUM_BUILD_VIEWER=OFF
cmake --build build-lib && cmake --install build-lib

# here
cmake -B build-fp -DSPATIUM_CONSUME=find_package \
      -DCMAKE_PREFIX_PATH=/tmp/spatium-prefix
cmake --build build-fp && ./build-fp/sphere_distance_demo
```

Both are gated in CI. The `find_package` path is the newer of the two, and
it was added because it was the *untested* one: Spatium has shipped
`install(EXPORT SpatiumTargets)`, a generated `SpatiumConfig.cmake` and a
version file since before its first release, and no job ever installed
them. A package config nothing consumes is a mechanism with nothing running
it, which `docs/conventions.md` has a whole section about.

Expected output:

```
Sphere<2> geodesic distance (north pole -> equator): 1.570796
Expected: pi/2 = 1.570796
```
