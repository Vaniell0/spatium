// spatium.spatial — single header (BVH). Depends on geometry for shape concepts.
module;
export module spatium.spatial;
import std.compat;
import spatium.core;
import spatium.algebra;
import spatium.geometry;
import spatium.spaces;
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/spatial/bound.hpp>
#include <spatium/spatial/bvh.hpp>
#include <spatium/spatial/geodesic_ball.hpp>
