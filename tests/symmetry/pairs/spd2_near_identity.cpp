#include "../pair.hpp"

#include <spatium/core/precision.hpp>
#include <spatium/spaces/spd.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

namespace symmetry {
namespace {

using spatium::Vec;

// ── SPD(2): a function of a matrix at a distance h from a scalar matrix, double against Real50 ──
//
// S = c I + h M with M a fixed-shape symmetric matrix of angle phi and h log-uniform in [1e-9, 1e-1],
// c in [0.5, 2]. Both scalars take the function sqrt of S (and its inverse, and the log) through the
// eigendecomposition and are held to each other. The bound is a few ulps: a function of a symmetric
// matrix is well conditioned everywhere, however close the eigenvalues are, and it is the
// decomposition that is not; so a bound of 1e-13 is what the function promises and what a
// decomposition has to keep.
//
// First run: the 2x2 eigen solver took the eigenvalues from the characteristic quadratic and the
// vectors as (b, lambda - a), and V diag(sqrt lambda) V^T missed by 1e-3 at h = 1e-6 and by 0.15 at
// h = 1e-7 (error ~ eps / h^2). It uses the half-angle formulas now, and the error is 1e-16 at every h.
struct NearIdentity { double c, h, phi; };

NearIdentity near_identity(std::size_t i) {
    return {0.5 + 1.5 * halton(i, 0), std::pow(10.0, -1.0 - 8.0 * halton(i, 1)),
            2 * std::numbers::pi * halton(i, 2)};
}

template<class T>
std::array<double, 12> functions_of(const NearIdentity& s) {
    using Aff = spatium::SPDAffineInvariant<2, T>;
    // The matrix is rounded once, to double, and both scalars see the same one.
    const double a = s.c + s.h * std::cos(s.phi), d = s.c - s.h * std::cos(s.phi), b = s.h * std::sin(s.phi);
    spatium::Matrix<T, 2, 2> m;
    m(0, 0) = T(a); m(0, 1) = m(1, 0) = T(b); m(1, 1) = T(d);
    const auto r1 = Aff::sqrt_sym(m);
    const auto r2 = Aff::inv_sqrt_sym(m);
    const auto r3 = spatium::detail::apply_eigen_sym(m, [](T x) { using std::log; return log(x); });
    std::array<double, 12> out{};
    std::size_t k = 0;
    for (const auto* r : {&r1, &r2, &r3})
        for (int i = 0; i < 2; ++i)
            for (int j = i; j < 2; ++j) out[k++] = spatium::primal_double((*r)(i, j));
    // 9 used, pad the rest with the traces so every slot means something
    out[9] = out[0] + out[2]; out[10] = out[3] + out[5]; out[11] = out[6] + out[8];
    return out;
}

Entry spd2_near_identity(std::size_t n = 4000) {
    return {"SPD(2) functions of a matrix near a scalar one: double against Real50", n, [n] {
        return check<NearIdentity>(
            "SPD(2) functions of a matrix near a scalar one: double against Real50", Kind::Function, n, near_identity,
            functions_of<spatium::Real50>, functions_of<double>,
            [](const std::array<double, 12>& a, const std::array<double, 12>& b) {
                double worst = 0.0;
                for (std::size_t k = 0; k < a.size(); ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const NearIdentity&) { return 1e-13; });
    }};
}

}  // namespace

const Registrar registered_spd2_near_identity{"spd2_near_identity", 4000, [](std::size_t n) { return spd2_near_identity(n); }};

}  // namespace symmetry
