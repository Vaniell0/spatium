// spatium.point — Point<S> wrapper + Morphism<From,To> with pipe composition.
// Single module with two partitions (point, morphism) since morphism depends
// on point and they ship together.
module;
export module spatium.point:point;
import std.compat;
import spatium.core;
import spatium.algebra;
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/point.hpp>
