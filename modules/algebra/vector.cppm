// vector partition embeds vec_expr.hpp directly: a forward declaration of
// `Vec` in a separate :vec_expr partition + the full definition in :vector
// would create two distinct module-attached entities and ambiguity. Sharing
// the same partition merges them.
module;
#include <cassert>  // assert is a macro — must come from a real header, not import std
export module spatium.algebra:vector;
import std.compat;
import spatium.core;
import :vec_simd;
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/algebra/vec_expr.hpp>
#include <spatium/algebra/vector.hpp>
