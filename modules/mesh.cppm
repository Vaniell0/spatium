// Single partition for spatium.mesh: 13 headers covering Mesh<S>, subdivision,
// topology, primitives, operations, lod, quality, geodesic, voronoi, transport
// and (when SPATIUM_EIGEN=ON) differential + heat_geodesic.
//
// One TU keeps every implicit Vec / Matrix / Mesh instantiation under one
// module attachment, avoiding the "conflicting deduced return type" trap
// (pitfall #8 — see lesson_cpp_modules_dual_mode.md).
//
// Phase 4 supersedes the Phase 2 single-partition `:types` shim — this file
// also brings `mesh/mesh.hpp` itself.
//
// Include order respects forward dependencies:
//   mesh → subdivision → topology → primitives → operations → lod → quality
//        → geodesic_types → geodesic → voronoi → transport
//        → differential? → heat_geodesic? (Eigen-only)
module;
#include <cstdint>  // UINT32_MAX macro used by topology.hpp
export module spatium.mesh;
import std.compat;
import spatium.core;
import spatium.algebra;
// NOTE: four mesh headers are DELIBERATELY outside this module; all stay
// header-only and module consumers `#include` them directly:
//   - `mesh/primitives.hpp` uses spatium.spaces; spaces imports mesh → cycle.
//   - `mesh/operations.hpp` uses spatium.geometry; geometry imports mesh → cycle.
//   - `mesh/differential.hpp` / `mesh/heat_geodesic.hpp` / `mesh/dec.hpp`
//     pull <Eigen/Sparse>, whose SSE intrinsics re-declare the x86 ones
//     already carried by spatium.algebra :vec_simd BMI → "conflicting
//     language linkage" (same pitfall that kept `:eigen` out of algebra).
// Clean fix for all four needs C++23 header units for the system headers;
// scheduled for Phase 7+.
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/mesh/topology.hpp>
#include <spatium/mesh/lod.hpp>
#include <spatium/mesh/quality.hpp>
#include <spatium/mesh/geodesic_types.hpp>
#include <spatium/mesh/geodesic.hpp>
#include <spatium/mesh/voronoi.hpp>
#include <spatium/mesh/transport.hpp>
