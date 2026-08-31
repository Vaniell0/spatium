// Single partition for spatium.spaces: each header has static_assert blocks
// that implicitly instantiate `Vec<T, N>` operators. Splitting into per-header
// partitions made gcc15 see independent template instantiations across BMIs
// and reject them as "conflicting deduced return type". Folding all 6 spaces
// into one TU keeps every instantiation under a single attachment.
module;
export module spatium.spaces;
import std.compat;
import spatium.core;
import spatium.algebra;
import spatium.mesh;
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/spaces/product.hpp>
#include <spatium/spaces/parametric.hpp>
#include <spatium/spaces/implicit.hpp>
