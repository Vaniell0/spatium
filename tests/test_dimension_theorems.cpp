// What the library offers at each dimension against what mathematics allows
// there.
//
// Some capabilities exist only at certain dimensions, by theorem:
//   a cross product (bilinear, orthogonal to both factors, with
//   |a x b|^2 = |a|^2 |b|^2 - (a.b)^2)       N = 3 and 7        Brown-Gray 1967
//   a normed division algebra                  N = 1, 2, 4, 8     Hurwitz 1898
//   a group structure on the sphere S^N        N = 1, 3           Hopf, Serre
//   a global tangent frame on S^N              N = 1, 3, 7        Bott-Milnor, Kervaire 1958
// The test holds the library to the one direction a machine can check: what
// it implements must lie inside what the theorem allows. A capability the
// library offers at a dimension the theorem forbids is a defect of the
// library or of the theorem's statement as coded -- a 4-dimensional `cross`
// that returned something would be wrong however plausible it looked. The
// other direction, allowed and not implemented (the groups S^1 and S^3), is a
// gap of glue and is listed, not failed. Where the two meet the gap is closed:
// the cross product at 3 and 7, the division algebras at 1, 2, 4 and 8 and the
// tangent frames on S^1, S^3 and S^7 are all implemented, no more and no less.
//
// So a red cell of the connectivity matrix is one of three kinds: no glue (the
// pieces exist and nothing joins them), no support (a missing basis), or
// impossible by theorem -- asking `Vec<T,4>` for a cross product, or
// `Sphere<2>` for a tangent frame, does not compile, and should not.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/complex.hpp>
#include <spatium/algebra/concepts.hpp>
#include <spatium/algebra/octonion.hpp>
#include <spatium/algebra/quaternion.hpp>
#include <spatium/algebra/vector.hpp>
#include <spatium/spaces/sphere.hpp>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <set>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

// Does Vec<double, N> offer a cross product: the member of R^3, or the free
// function of R^7 (algebra/octonion.hpp)?
template<std::size_t N>
constexpr bool has_cross() {
    return requires(const Vec<double, N>& a, const Vec<double, N>& b) { a.cross(b); } ||
           requires(const Vec<double, N>& a, const Vec<double, N>& b) { cross7(a, b); };
}

// Does Sphere<N> offer a global tangent frame?
template<std::size_t N>
constexpr bool has_frame() {
    return requires(const Sphere<N, double>& s, const Vec<double, N + 1>& p) { s.tangent_frame(p); };
}

// Is Sphere<N> a group in the library's own sense?
template<std::size_t N>
constexpr bool sphere_is_group() { return algebra::Group<Sphere<N, double>>; }

template<std::size_t... Ns>
std::set<std::size_t> cross_dimensions(std::index_sequence<Ns...>) {
    std::set<std::size_t> found;
    ((has_cross<Ns + 1>() ? (void)found.insert(Ns + 1) : void()), ...);
    return found;
}
template<std::size_t... Ns>
std::set<std::size_t> frame_dimensions(std::index_sequence<Ns...>) {
    std::set<std::size_t> found;
    ((has_frame<Ns + 1>() ? (void)found.insert(Ns + 1) : void()), ...);
    return found;
}
template<std::size_t... Ns>
std::set<std::size_t> group_dimensions(std::index_sequence<Ns...>) {
    std::set<std::size_t> found;
    ((sphere_is_group<Ns + 1>() ? (void)found.insert(Ns + 1) : void()), ...);
    return found;
}

}  // namespace

TEST_CASE("The cross product exists only where a theorem allows one", "[dimensions]") {
    const std::set<std::size_t> allowed{3, 7};
    const auto implemented = cross_dimensions(std::make_index_sequence<12>{});
    for (const auto n : implemented) CHECK(allowed.count(n) == 1);
    // Both dimensions the theorem allows are written, and no other: the
    // member of R^3 and `cross7` of R^7. The line above is what keeps a
    // third from being added.
    CHECK(implemented == allowed);
}

TEST_CASE("The three-dimensional cross product is the one the theorem describes", "[dimensions]") {
    std::mt19937_64 g(3);
    std::uniform_real_distribution<double> u(-1, 1);
    for (int i = 0; i < 50; ++i) {
        const Vec<double, 3> a{u(g), u(g), u(g)}, b{u(g), u(g), u(g)};
        const auto c = a.cross(b);
        CHECK_THAT(c.dot(a), WithinAbs(0.0, 1e-14));
        CHECK_THAT(c.dot(b), WithinAbs(0.0, 1e-14));
        const double lhs = c.dot(c), rhs = a.dot(a) * b.dot(b) - a.dot(b) * a.dot(b);
        CHECK_THAT(lhs, WithinAbs(rhs, 1e-13));
    }
}

TEST_CASE("The normed division algebras the library has are at dimensions 1, 2, 4 and 8", "[dimensions]") {
    // |ab| = |a||b|: the defining property, at N = 2 (Complex), N = 4
    // (Quaternion) and N = 8 (Octonion). N = 1 is the reals. No other
    // dimension has a type to test, and none may.
    std::mt19937_64 g(5);
    std::uniform_real_distribution<double> u(-2, 2);
    for (int i = 0; i < 100; ++i) {
        const Complex<double> a{u(g), u(g)}, b{u(g), u(g)};
        const auto ab = a * b;
        CHECK_THAT(std::hypot(ab.re, ab.im), WithinAbs(std::hypot(a.re, a.im) * std::hypot(b.re, b.im), 1e-12));

        const Quaternion<double> p{u(g), u(g), u(g), u(g)}, q{u(g), u(g), u(g), u(g)};
        CHECK_THAT((p * q).norm(), WithinAbs(p.norm() * q.norm(), 1e-12));

        const Octonion<double> x{u(g), u(g), u(g), u(g), u(g), u(g), u(g), u(g)},
                               y{u(g), u(g), u(g), u(g), u(g), u(g), u(g), u(g)};
        CHECK_THAT((x * y).norm(), WithinAbs(x.norm() * y.norm(), 1e-11));
    }
}

TEST_CASE("A global tangent frame exists on S^1, S^3 and S^7 and on no other sphere", "[dimensions]") {
    // Parallelizable spheres (Bott-Milnor, Kervaire, 1958): the three of the
    // division algebras. The member `tangent_frame` exists at those N only, so
    // asking S^2 for one does not compile -- which is the theorem, not a gap.
    const std::set<std::size_t> allowed{1, 3, 7};
    const auto implemented = frame_dimensions(std::make_index_sequence<12>{});
    for (const auto n : implemented) CHECK(allowed.count(n) == 1);
    CHECK(implemented == allowed);
}

TEST_CASE("No sphere is a group in the library; two could be", "[dimensions]") {
    // S^1 (unit complex numbers) and S^3 (unit quaternions) are Lie groups
    // and the only spheres that are; the library exposes neither as a Group.
    // The list is what a future adapter fills in -- the guard is that nothing
    // outside {1, 3} ever appears here.
    const std::set<std::size_t> allowed{1, 3};
    const auto implemented = group_dimensions(std::make_index_sequence<10>{});
    for (const auto n : implemented) CHECK(allowed.count(n) == 1);
    CHECK(implemented.empty());
}
