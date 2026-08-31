// Single partition for spatium.geometry: 20 headers with heavy mutual
// includes plus a forward decl in line.hpp (`template<...> struct Box;`).
// Splitting along header boundaries produced separate module-attached
// `Box`/`Vec` instantiations and "conflicting deduced return type" errors
// — same pattern as Phase 1's vec_expr→:vector and Phase 2's spaces.
//
// Headers must follow include-order (defs before uses):
//   concepts → box → hyperplane → line → triangle → simplex → circle
//          → polygon → convex_hull → boolean → clip → intersection
//          → distance → make → transform → surface_adapter → ray_surface
//          → ray_parametric → format
module;
export module spatium.geometry;
import std.compat;
import spatium.core;
import spatium.algebra;
import spatium.mesh;
import spatium.spaces;
import spatium.point;
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/geometry/concepts.hpp>
#include <spatium/geometry/box.hpp>
#include <spatium/geometry/hyperplane.hpp>
#include <spatium/geometry/line.hpp>
#include <spatium/geometry/triangle.hpp>
#include <spatium/geometry/simplex.hpp>
#include <spatium/geometry/circle.hpp>
#include <spatium/geometry/polygon.hpp>
#include <spatium/geometry/convex_hull.hpp>
#include <spatium/geometry/boolean.hpp>
#include <spatium/geometry/clip.hpp>
#include <spatium/geometry/intersection.hpp>
#include <spatium/geometry/distance.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/geometry/transform.hpp>
#include <spatium/geometry/surface_adapter.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/geometry/ray_hit.hpp>
#include <spatium/geometry/ray_parametric.hpp>
#include <spatium/geometry/format.hpp>
