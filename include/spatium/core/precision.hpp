#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <boost/multiprecision/cpp_dec_float.hpp>
#endif

SPATIUM_EXPORT namespace spatium {

// Arbitrary precision decimal floating-point types.
// These satisfy the Scalar concept and can be used with all Spatium types:
//   Euclidean<3, Real50>, Vec<Real100, 3>, Triangle<3, Real50>, etc.

using Real50  = boost::multiprecision::cpp_dec_float_50;
using Real100 = boost::multiprecision::cpp_dec_float_100;

// Verify Scalar concept satisfaction
static_assert(Scalar<Real50>);
static_assert(Scalar<Real100>);

} // namespace spatium
