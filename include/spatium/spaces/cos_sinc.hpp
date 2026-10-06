#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium {

// cos(sqrt(s)) and sin(sqrt(s)) / sqrt(s) -- cosh and sinh for a hyperbolic
// space -- as functions of s itself.
//
// This is the form in which a geodesic step is smooth in its tangent. The
// exponential map of a sphere is p cos(theta) + v t sinc(theta) with
// theta^2 = t^2 |v|^2 / r^2, and both cos and sinc are even in theta, so
// they are analytic in theta^2 straight through zero. Written with |v| (a
// square root, and a division by it) they are not: at v = 0 the square
// root has no derivative and the quotient is 0/0, so the code kept a
// `speed == 0` shortcut that returned p -- and for a Dual tangent of value
// zero, whose `== 0` compares the value alone, returned p with a zero
// derivative where d/dv exp_p(v) = t. Found by the scout of the "== 0 on a
// Dual" class, after the connectivity matrix found it for t.
//
// Near zero a Taylor series in s takes over -- six terms, which is exact to
// rounding below a threshold that follows the scalar's epsilon, so the
// same function is right on float, double and Real50. Both branches are
// smooth in s and meet at the threshold to rounding, so a Dual's derivative
// is continuous across it; the branch is chosen on the value, which is all
// a choice may read.
template<Scalar T>
struct CosSinc {
    T cos;    // cos(sqrt(s))            (cosh when hyperbolic)
    T sinc;   // sin(sqrt(s))/sqrt(s)    (sinh when hyperbolic)
};

template<Scalar T>
CosSinc<T> cos_sinc_of_square(const T& s, bool hyperbolic = false) {
    using std::cos; using std::sin; using std::cosh; using std::sinh; using std::sqrt;

    // error of the sixth-order series is s^6 / 13!, below eps for
    // s < 0.1 eps^(1/6) with room to spare.
    const double small = 0.1 * std::pow(machine_epsilon<T>(), 1.0 / 6.0);
    if (primal_double(s) < small) {
        // C = sum w^k / (2k)!, S = sum w^k / (2k+1)!, w = -s circular, +s hyperbolic.
        const T w = hyperbolic ? s : T(-s);
        T term_c{1}, term_s{1}, C{1}, S{1};
        for (int i = 1; i <= 5; ++i) {
            term_c = T(term_c * w / T((2 * i - 1) * (2 * i)));
            term_s = T(term_s * w / T((2 * i) * (2 * i + 1)));
            C = T(C + term_c);
            S = T(S + term_s);
        }
        return {C, S};
    }
    const T r = sqrt(s);
    if (hyperbolic) return {T(cosh(r)), T(sinh(r) / r)};
    return {T(cos(r)), T(sin(r) / r)};
}

} // namespace spatium
