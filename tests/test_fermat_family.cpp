// The Fermat family sum x_i^n = 1 -- Lame's superellipsoids -- as typed
// ImplicitSurfaces: one lambda, and the invariants that stay while the shape
// does not.
//
// n = 2 is the unit sphere; as n grows the surface flattens toward the cube
// [-1, 1]^3 and its curvature gathers at the eight corners. Through all of it
// the surface is a sphere topologically, so the total Gaussian curvature is
// 4 pi (Gauss-Bonnet) and the Euler characteristic of a mesh of it is 2.
//
// Two routes to 4 pi that share nothing:
//   smooth  K from the implicit function (exact, through nested Duals), dA from
//           a radial chart of the same surface, by the library's quadrature;
//   mesh    the sum over a closed triangle mesh of the angle defects
//           2 pi - (sum of the angles at a vertex), which equals 2 pi chi for
//           ANY closed mesh -- a combinatorial identity, no curvature in it.
// When they disagree, the quadrature is the one to doubt: at n = 8 the
// integrand is a set of sharp peaks, the rule reports its status Suspicious
// and an estimate of 3e-9 over a true error of 1e-3, and it is the known 4 pi
// -- a second witness -- that says so.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/calculus.hpp>
#include <spatium/mesh/tessellate.hpp>
#include <spatium/core/access.hpp>
#include <spatium/spaces/implicit.hpp>
#include <spatium/spaces/parametric.hpp>
#include <spatium/spaces/sphere.hpp>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <set>
#include <utility>
#include <vector>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

constexpr double pi = std::numbers::pi;

// x^n for an even n, as products: smooth, and callable on any Dual.
auto even_power(int n) {
    return [n](auto t) {
        const auto s = t * t;
        auto r = s;
        for (int k = 2; k < n; k += 2) r = r * s;
        return r;
    };
}

// F_n = x^n + y^n + z^n - 1, with bounds that are not symmetric: a mesh whose
// grid never passes exactly through the surface has no coincident vertices.
auto fermat(int n) {
    return make_implicit<double>(
        [pw = even_power(n)](auto x, auto y, auto z) { return pw(x) + pw(y) + pw(z) - 1.0; },
        {-1.2345, 1.1987, -1.2111, 1.2376, -1.1923, 1.2201});
}

// The same surface as a radial graph over the sphere: p = r(w) w with
// r = (sum w_i^n)^(-1/n). u is the azimuth, v the polar angle.
auto radial(int n) {
    return make_parametric<double>(
        [pw = even_power(n), n](auto u, auto v) {
            using std::sin; using std::cos; using std::exp; using std::log;
            const auto a = sin(v) * cos(u), b = sin(v) * sin(u), c = cos(v);
            const auto r = exp(-log(pw(a) + pw(b) + pw(c)) * (1.0 / n));
            return Vec<decltype(u), 3>{r * a, r * b, r * c};
        },
        {0.0, 2 * pi, 0.0, pi}, true, false);
}

// Total curvature by quadrature: periodic trapezoid in the azimuth (spectrally
// accurate for an analytic periodic integrand), adaptive Simpson in the polar angle.
IntegralResult<double> total_curvature_smooth(int n) {
    const auto level = fermat(n);
    const auto chart = radial(n);
    constexpr int steps = 96;
    const double dphi = 2 * pi / steps;
    return integrate_with_error<double>(
        [&](double theta) {
            double acc = 0;
            for (int k = 0; k < steps; ++k) {
                const double phi = (k + 0.5) * dphi;
                acc += level.gaussian_curvature(chart.evaluate(phi, theta)) * chart.area_element(phi, theta);
            }
            return acc * dphi;
        },
        1e-4, pi - 1e-4, 1e-8);
}

struct MeshInvariants { long euler; double angle_defect; };

// V - E + F, and the sum of the angle defects, over a marching-cubes mesh.
MeshInvariants mesh_invariants(int n, std::size_t resolution) {
    const auto mesh = marching_cubes(fermat(n).erased(), resolution);
    std::set<std::pair<unsigned, unsigned>> edges;
    std::vector<double> angle_sum(mesh.vertices.size(), 0.0);
    for (const auto& f : mesh.faces)
        for (int i = 0; i < 3; ++i) {
            const unsigned a = f[i], b = f[(i + 1) % 3];
            edges.insert({std::min(a, b), std::max(a, b)});
            const auto& p0 = mesh.vertices[f[i]];
            const Vec<double, 3> e1 = mesh.vertices[f[(i + 1) % 3]] - p0, e2 = mesh.vertices[f[(i + 2) % 3]] - p0;
            angle_sum[f[i]] += std::atan2(e1.cross(e2).norm(), e1.dot(e2));
        }
    double defect = 0;
    for (const double a : angle_sum) defect += 2 * pi - a;
    return {static_cast<long>(mesh.vertices.size()) - static_cast<long>(edges.size()) +
                static_cast<long>(mesh.faces.size()),
            defect};
}

}  // namespace

TEST_CASE("n = 2 is the sphere: curvature 1 and the geodesics of Sphere<2>", "[fermat][implicit]") {
    const auto level = fermat(2);
    const Sphere<2> closed;
    const Vec<double, 3> p{0.6, 0.0, 0.8};
    CHECK_THAT(level.gaussian_curvature(p), WithinAbs(1.0, 1e-13));
    const Vec<double, 3> w{0.0, 1.0, 0.0};
    CHECK((level.exp_map(p, w, 1.1) - closed.exp_map(p, w, 1.1)).norm() < 1e-8);
}

TEST_CASE("The curvature gathers at the corners as n grows, and stays positive", "[fermat][implicit]") {
    // On the diagonal the surface is at x = y = z = 3^(-1/n); K there grows
    // with n while K on an axis, at (1, 0, 0), falls to 0 -- the flat faces.
    double corner_before = 0;
    for (const int n : {2, 4, 6, 8}) {
        const auto level = fermat(n);
        const double c = std::pow(3.0, -1.0 / n);
        const double corner = level.gaussian_curvature(Vec<double, 3>{c, c, c});
        const double face = level.gaussian_curvature(Vec<double, 3>{1.0, 0.0, 0.0});
        CHECK(corner > 0.0);
        CHECK(face >= 0.0);
        CHECK(corner > corner_before);
        if (n > 2) CHECK(face < corner);
        corner_before = corner;
    }
}

TEST_CASE("As n grows the surface approaches the cube: the diagonal point goes to 1", "[fermat][implicit]") {
    double previous = 0;
    for (const int n : {2, 4, 8, 16}) {
        const auto level = fermat(n);
        const auto hit = level.project(Vec<double, 3>{1.1, 1.0, 0.9});     // onto the surface, near the corner
        const double expected_diagonal = std::pow(3.0, -1.0 / n);
        CHECK_THAT(std::abs(level(hit)), WithinAbs(0.0, 1e-12));
        const auto on_diagonal = level.project(Vec<double, 3>{1.0, 1.0, 1.0});
        CHECK_THAT(on_diagonal[0], WithinAbs(expected_diagonal, 1e-9));
        CHECK(on_diagonal[0] > previous);
        previous = on_diagonal[0];
    }
    CHECK(previous > 0.8);          // 3^(-1/16) = 0.934
}

TEST_CASE("Gauss-Bonnet by the mesh: the angle defects sum to 4 pi at every n", "[fermat][implicit]") {
    // An identity of any closed triangle mesh of a sphere, whatever its shape.
    for (const int n : {2, 4, 6, 8}) {
        const auto inv = mesh_invariants(n, 48);
        INFO("n = " << n);
        CHECK(inv.euler == 2);
        CHECK_THAT(inv.angle_defect, WithinAbs(4 * pi, 1e-9));
    }
}

TEST_CASE("Gauss-Bonnet by the smooth path: 4 pi, and where the quadrature cannot see it", "[fermat][implicit]") {
    // Where the integrand is smooth the two routes agree to quadrature accuracy.
    for (const auto& [n, tolerance] : {std::pair{2, 1e-6}, std::pair{4, 1e-6}, std::pair{6, 1e-4}}) {
        const auto smooth = total_curvature_smooth(n);
        INFO("n = " << n << ", estimate " << smooth.error_estimate);
        CHECK_THAT(smooth.value, WithinAbs(4 * pi, tolerance));
    }
    // At n = 8 the curvature is eight sharp peaks. The rule's estimate says 1e-8;
    // the true error is 1e-3, which only the other witness -- the mesh, or the
    // known 4 pi -- shows. The status is the rule's own warning, and it is given.
    const auto peaked = total_curvature_smooth(8);
    CHECK(peaked.status == IntegralStatus::Suspicious);
    CHECK(std::abs(peaked.value - 4 * pi) > 100 * peaked.error_estimate);
    CHECK_THAT(peaked.value, WithinAbs(4 * pi, 1e-2));
}
