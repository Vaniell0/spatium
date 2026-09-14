module;
// Guarded the same way core/precision.hpp itself is: with SPATIUM_BOOST
// off there is no Boost on the include path, so the global module
// fragment must not reach for it either. The partition still exists and
// is still `export import`-ed by core.cppm -- it just exports nothing,
// which is legal and keeps the module graph unchanged either way.
#if defined(SPATIUM_HAS_BOOST_MULTIPRECISION) && SPATIUM_HAS_BOOST_MULTIPRECISION
#  include <boost/multiprecision/cpp_dec_float.hpp>
#endif
export module spatium.core:precision;
import std.compat;
import :concepts;
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/core/precision.hpp>
