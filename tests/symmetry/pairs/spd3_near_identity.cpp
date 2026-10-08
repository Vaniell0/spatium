#include "../pair.hpp"

#include <spatium/core/precision.hpp>
#include <spatium/spaces/spd.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace symmetry {
namespace {

// ── SPD(3): a function of a matrix at a distance h from a scalar matrix, double against Real50 ──
//
// S = c I + h M with M a symmetric matrix of six numbers in [-1, 1] and h log-uniform in [1e-9, 1e-1],
// c in [0.5, 2]; sqrt, inverse sqrt and log of S through the eigendecomposition on both scalars, held
// to each other to 1e-13. A function of a symmetric matrix is well conditioned however close the
// eigenvalues are; the decomposition is what must keep up.
//
// First run: the 3x3 solver took its eigenvalues from the characteristic cubic and its vectors from
// cross products of rows of S - lambda I. On I + hM the reconstruction V diag(l) V^T missed by 1.7e-3 at
// h = 1e-2 and by 0.95 at h = 1e-5, the vectors no longer orthonormal. Cyclic Jacobi rotations replaced
// it; the error is 1e-15 at every h.
struct NearIdentity3 { double c, h; std::array<double, 6> m; };

NearIdentity3 near_identity3(std::size_t i) {
    NearIdentity3 s{0.5 + 1.5 * halton(i, 0), std::pow(10.0, -1.0 - 8.0 * halton(i, 1)), {}};
    for (std::size_t k = 0; k < 6; ++k) s.m[k] = 2 * halton(i, 2 + k) - 1;
    return s;
}

template<class T>
std::array<double, 18> functions_of3(const NearIdentity3& s) {
    using Aff = spatium::SPDAffineInvariant<3, T>;
    // Rounded once, to double, and both scalars see the same matrix.
    const double a00 = s.c + s.h * s.m[0], a01 = s.h * s.m[1], a02 = s.h * s.m[2],
                 a11 = s.c + s.h * s.m[3], a12 = s.h * s.m[4], a22 = s.c + s.h * s.m[5];
    spatium::Matrix<T, 3, 3> m;
    m(0, 0) = T(a00); m(0, 1) = m(1, 0) = T(a01); m(0, 2) = m(2, 0) = T(a02);
    m(1, 1) = T(a11); m(1, 2) = m(2, 1) = T(a12); m(2, 2) = T(a22);
    const auto r1 = Aff::sqrt_sym(m);
    const auto r2 = Aff::inv_sqrt_sym(m);
    const auto r3 = spatium::detail::apply_eigen_sym(m, [](T x) { using std::log; return log(x); });
    std::array<double, 18> out{};
    std::size_t k = 0;
    for (const auto* r : {&r1, &r2, &r3})
        for (int i = 0; i < 3; ++i)
            for (int j = i; j < 3; ++j) out[k++] = spatium::primal_double((*r)(i, j));
    return out;
}

Entry spd3_near_identity(std::size_t n = 3000) {
    return {"SPD(3) functions of a matrix near a scalar one: double against Real50", n, [n] {
        return check<NearIdentity3>(
            "SPD(3) functions of a matrix near a scalar one: double against Real50", Kind::Function, n, near_identity3,
            functions_of3<spatium::Real50>, functions_of3<double>,
            [](const std::array<double, 18>& a, const std::array<double, 18>& b) {
                double worst = 0.0;
                for (std::size_t k = 0; k < a.size(); ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const NearIdentity3&) { return 1e-13; });
    }};
}

}  // namespace

const Registrar registered_spd3_near_identity{"spd3_near_identity", 3000, [](std::size_t n) { return spd3_near_identity(n); }};

}  // namespace symmetry
