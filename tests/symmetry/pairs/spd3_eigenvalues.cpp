#include "../pair.hpp"

#include <spatium/core/precision.hpp>
#include <spatium/spaces/spd.hpp>

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

// ── SPD(3): the eigenvalues through the characteristic cubic, double against Real50 ──
//
// S = Q diag(1, 1 + gap, 3) Q^T with Q from three Euler angles and the gap log-uniform in
// [1e-10, 1], rounded once to double and handed to both scalars, so the two paths see one
// matrix. The bound is the conditioning of a cubic root: the coefficients carry a rounding of
// eps S^3, a root with gap g to its neighbour moves by that over its slope there, about g, and
// the slope cannot be smaller than sqrt(eps) S before the root is a double one, where the error
// saturates at sqrt(eps) S. Thirty times that is the promise (about four times the worst measured, 2.2e-7 at gaps
// near 1e-7, against a bound of 4.3e-6 there).
//
// First run: `eigen_sym<Real50>` did not compile (a unary minus on a Boost number is a lazy
// expression and solve_cubic<T> cannot deduce T from four different types).
struct SpdCase { std::array<double, 9> m; double gap; };

inline SpdCase spd_case(std::size_t i) {
    const double a = 2 * std::numbers::pi * halton(i, 0), b = std::numbers::pi * halton(i, 1),
                 c = 2 * std::numbers::pi * halton(i, 2);
    const double gap = std::pow(10.0, -10.0 * halton(i, 3));
    const double D[3] = {1.0, 1.0 + gap, 3.0};
    const double ca = std::cos(a), sa = std::sin(a), cb = std::cos(b), sb = std::sin(b), cc = std::cos(c), sc = std::sin(c);
    const double Rz1[9] = {ca, -sa, 0, sa, ca, 0, 0, 0, 1}, Rx[9] = {1, 0, 0, 0, cb, -sb, 0, sb, cb},
                 Rz2[9] = {cc, -sc, 0, sc, cc, 0, 0, 0, 1};
    const auto mul = [](const double* A, const double* B, double* C) {
        for (int r = 0; r < 3; ++r)
            for (int k = 0; k < 3; ++k) {
                double sum = 0;
                for (int j = 0; j < 3; ++j) sum += A[r * 3 + j] * B[j * 3 + k];
                C[r * 3 + k] = sum;
            }
    };
    double T1[9], Q[9];
    mul(Rz1, Rx, T1);
    mul(T1, Rz2, Q);
    SpdCase out{};
    out.gap = gap;
    for (int r = 0; r < 3; ++r)
        for (int k = 0; k < 3; ++k) {
            double sum = 0;
            for (int j = 0; j < 3; ++j) sum += Q[r * 3 + j] * D[j] * Q[k * 3 + j];
            out.m[r * 3 + k] = sum;
        }
    for (int r = 0; r < 3; ++r)
        for (int k = r + 1; k < 3; ++k) out.m[r * 3 + k] = out.m[k * 3 + r] = 0.5 * (out.m[r * 3 + k] + out.m[k * 3 + r]);
    return out;
}

template<class T>
std::array<double, 3> spd_eigenvalues(const SpdCase& s) {
    spatium::Matrix<T, 3, 3> M;
    for (int r = 0; r < 3; ++r)
        for (int k = 0; k < 3; ++k) M(r, k) = T(s.m[r * 3 + k]);
    const auto e = spatium::detail::eigen_sym(M);
    std::array<double, 3> v{spatium::primal_double(e.values[0]), spatium::primal_double(e.values[1]),
                            spatium::primal_double(e.values[2])};
    std::sort(v.begin(), v.end());
    return v;
}

inline Entry spd3_eigenvalues_precision(std::size_t n = 4000) {
    return {"SPD(3) eigenvalues: double against Real50", n, [n] {
        return check<SpdCase>(
            "SPD(3) eigenvalues: double against Real50", Kind::Function, n, spd_case,
            spd_eigenvalues<spatium::Real50>, spd_eigenvalues<double>,
            [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
                double worst = 0.0;
                for (int k = 0; k < 3; ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const SpdCase& s) {
                const double eps = std::numeric_limits<double>::epsilon(), S = 4.0;
                return 30.0 * eps * S * S * S / std::max(s.gap, std::sqrt(eps) * S);
            });
    }};
}

}  // namespace

const Registrar registered_spd3_eigenvalues{"spd3_eigenvalues_precision", 4000, [](std::size_t n) { return spd3_eigenvalues_precision(n); }};

}  // namespace symmetry
