#include "../pair.hpp"

#include <spatium/algebra/calculus.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/algebra/vector.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

namespace symmetry {
namespace {

using spatium::Vec;

// ── gradient: the exact derivative through Dual against central differences ──
//
// The bound is the central difference's own: a step of cbrt(eps) leaves a truncation of
// order h^2 f''' and a rounding of order eps / h, both near 1e-11 for this f; 1e-9 of the
// size of the gradient is the promise, and the worst distance measured is 6e-11 of it.
template<class T>
T gradient_test_function(const Vec<T, 3>& x) {
    using std::sin; using std::exp;
    return sin(x[0] * x[1]) + exp(x[2]) * x[0];
}

inline Entry gradient_dual_against_differences(std::size_t n = 2000) {
    return {"gradient: Dual against central differences", n, [n] {
        return check<Vec<double, 3>>(
            "gradient: Dual against central differences", Kind::Function, n,
            [](std::size_t i) { return Vec<double, 3>{-2 + 4 * halton(i, 0), -2 + 4 * halton(i, 1), -2 + 4 * halton(i, 2)}; },
            [](const Vec<double, 3>& x) {
                return spatium::gradient([](const Vec<spatium::Dual<double>, 3>& v) { return gradient_test_function(v); }, x);
            },
            [](const Vec<double, 3>& x) {
                const double h = std::cbrt(std::numeric_limits<double>::epsilon());
                Vec<double, 3> g{};
                for (std::size_t k = 0; k < 3; ++k) {
                    Vec<double, 3> hi = x, lo = x;
                    const double step = h * std::max(1.0, std::abs(x[k]));
                    hi[k] += step; lo[k] -= step;
                    g[k] = (gradient_test_function(hi) - gradient_test_function(lo)) / (2 * step);
                }
                return g;
            },
            [](const Vec<double, 3>& a, const Vec<double, 3>& b) { return (a - b).norm(); },
            [](const Vec<double, 3>& x) { return 1e-9 * (1.0 + spatium::gradient([](const Vec<spatium::Dual<double>, 3>& v) { return gradient_test_function(v); }, x).norm()); });
    }};
}

}  // namespace

const Registrar registered_gradient_dual{"gradient_dual_against_differences", 2000, [](std::size_t n) { return gradient_dual_against_differences(n); }};

}  // namespace symmetry
