#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <boost/multiprecision/cpp_dec_float.hpp>
#endif

SPATIUM_EXPORT namespace spatium {

// Arbitrary precision decimal floating-point, parametrized on digit count.
// `Real50`/`Real100` used to be the only two options (hardcoded
// `cpp_dec_float_50`/`_100` typedefs) -- `Real<Digits>` is the same
// underlying Boost.Multiprecision machinery generalized to any compile-time
// digit count, so "how much precision do you need" is an actual parameter,
// not a choice between two fixed presets. Satisfies the Scalar concept and
// can be used with all Spatium types the same way the old aliases were:
//   Euclidean<3, Real<50>>, Vec<Real<200>, 3>, Triangle<3, Real<500>>, etc.
template<unsigned Digits>
using Real = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<Digits>>;

// Kept as plain aliases for Real<50>/Real<100> -- every existing call site
// (Real50{...}, Vec<Real50,N>, solve_cubic<Real50>, ...) keeps compiling
// unchanged; Real<Digits> is additive, not a replacement.
using Real50  = Real<50>;
using Real100 = Real<100>;

// Verify Scalar concept satisfaction, including a digit count neither
// preset uses, so this isn't just re-checking Real50/Real100 under a new
// name.
static_assert(Scalar<Real50>);
static_assert(Scalar<Real100>);
static_assert(Scalar<Real<24>>);
static_assert(Scalar<Real<250>>);

} // namespace spatium
