// The octonions (algebra/octonion.hpp): the normed division algebra of
// dimension 8, built by Cayley-Dickson doubling of the library's quaternions.
//
// There is no second table of structure constants to compare the product
// with, and none is needed: the properties below determine the octonions up
// to isomorphism (a unital composition algebra of dimension 8 is the
// octonions, Hurwitz), so an error in a sign of the doubling breaks at least
// one of them. Every identity is checked on random elements, over double,
// over Real50 where Boost is present, and over Dual, where the derivative
// reaches the same answer by two routes.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/complex.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/algebra/octonion.hpp>
#include <spatium/algebra/quaternion.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/spaces/sphere.hpp>
#include <cmath>
#include <random>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

Oct random_oct(std::mt19937_64& g) {
    std::uniform_real_distribution<double> u(-1.5, 1.5);
    return Oct{u(g), u(g), u(g), u(g), u(g), u(g), u(g), u(g)};
}

double distance(const Oct& x, const Oct& y) { return (x - y).norm(); }

}  // namespace

TEST_CASE("The basis units square to -1 and anticommute", "[octonion]") {
    for (std::size_t i = 1; i < 8; ++i) {
        CHECK(distance(Oct::basis(i) * Oct::basis(i), Oct::basis(0) * -1.0) < 1e-15);
        for (std::size_t j = 1; j < 8; ++j) {
            if (i == j) continue;
            const Oct ij = Oct::basis(i) * Oct::basis(j), ji = Oct::basis(j) * Oct::basis(i);
            CHECK(distance(ij, ji * -1.0) < 1e-15);
            // The product of two distinct imaginary units is a signed unit.
            double biggest = 0;
            for (std::size_t k = 0; k < 8; ++k) biggest = std::max(biggest, std::abs(ij[k]));
            CHECK_THAT(biggest, WithinAbs(1.0, 1e-15));
            CHECK_THAT(ij.norm(), WithinAbs(1.0, 1e-15));
        }
    }
}

TEST_CASE("The norm is multiplicative: a normed division algebra", "[octonion]") {
    std::mt19937_64 g(8);
    for (int n = 0; n < 500; ++n) {
        const Oct a = random_oct(g), b = random_oct(g);
        CHECK_THAT((a * b).norm(), WithinAbs(a.norm() * b.norm(), 1e-12));
        CHECK_THAT((a * b).norm_squared(), WithinAbs(a.norm_squared() * b.norm_squared(), 1e-11));
    }
}

TEST_CASE("Every nonzero octonion has a two-sided inverse", "[octonion]") {
    std::mt19937_64 g(9);
    for (int n = 0; n < 200; ++n) {
        const Oct a = random_oct(g);
        const auto inv = a.inverse();
        REQUIRE(inv.has_value());
        CHECK(distance(a * *inv, Oct{}) < 1e-12);
        CHECK(distance(*inv * a, Oct{}) < 1e-12);
    }
    CHECK_FALSE(Oct{0, 0, 0, 0, 0, 0, 0, 0}.inverse().has_value());
}

TEST_CASE("The octonions are alternative, and satisfy the Moufang identities", "[octonion]") {
    std::mt19937_64 g(10);
    for (int n = 0; n < 300; ++n) {
        const Oct a = random_oct(g), b = random_oct(g), c = random_oct(g);
        // a(ab) = (aa)b and (ab)b = a(bb): any two elements generate an
        // associative subalgebra.
        CHECK(distance(a * (a * b), (a * a) * b) < 1e-11);
        CHECK(distance((a * b) * b, a * (b * b)) < 1e-11);
        // Moufang: (ab)(ca) = a((bc)a), a(b(ac)) = ((ab)a)c.
        CHECK(distance((a * b) * (c * a), a * ((b * c) * a)) < 1e-10);
        CHECK(distance(a * (b * (a * c)), ((a * b) * a) * c) < 1e-10);
    }
}

TEST_CASE("The octonions are not associative", "[octonion]") {
    // The defect, and a witness: (e1 e2) e4 = -e1 (e2 e4).
    const Oct lhs = (Oct::basis(1) * Oct::basis(2)) * Oct::basis(4);
    const Oct rhs = Oct::basis(1) * (Oct::basis(2) * Oct::basis(4));
    CHECK(distance(lhs, rhs * -1.0) < 1e-15);
    CHECK(distance(lhs, rhs) > 1.9);
    // And on random elements the associator is nonzero almost everywhere.
    std::mt19937_64 g(11);
    int nonzero = 0;
    for (int n = 0; n < 100; ++n) {
        const Oct a = random_oct(g), b = random_oct(g), c = random_oct(g);
        if (distance((a * b) * c, a * (b * c)) > 1e-3) ++nonzero;
    }
    CHECK(nonzero > 95);
}

TEST_CASE("Quaternions and complex numbers sit inside the octonions", "[octonion]") {
    std::mt19937_64 g(12);
    std::uniform_real_distribution<double> u(-1.5, 1.5);
    for (int n = 0; n < 100; ++n) {
        const Quat p{u(g), u(g), u(g), u(g)}, q{u(g), u(g), u(g), u(g)};
        const Oct product = Oct{p, Quat{0, 0, 0, 0}} * Oct{q, Quat{0, 0, 0, 0}};
        const Quat pq = p * q;
        CHECK_THAT(product[0], WithinAbs(pq.w, 1e-14));
        CHECK_THAT(product[1], WithinAbs(pq.x, 1e-14));
        CHECK_THAT(product[2], WithinAbs(pq.y, 1e-14));
        CHECK_THAT(product[3], WithinAbs(pq.z, 1e-14));
        for (std::size_t k = 4; k < 8; ++k) CHECK(product[k] == 0.0);

        const Complex<double> z{u(g), u(g)}, w{u(g), u(g)};
        const Oct zw = Oct{z.re, z.im, 0, 0, 0, 0, 0, 0} * Oct{w.re, w.im, 0, 0, 0, 0, 0, 0};
        const Complex<double> c = z * w;
        CHECK_THAT(zw[0], WithinAbs(c.re, 1e-14));
        CHECK_THAT(zw[1], WithinAbs(c.im, 1e-14));
    }
}

TEST_CASE("The seven-dimensional cross product has the identities of a cross product", "[octonion][dimensions]") {
    std::mt19937_64 g(13);
    std::uniform_real_distribution<double> u(-1.5, 1.5);
    const auto vec7 = [&] { Vec<double, 7> v; for (std::size_t i = 0; i < 7; ++i) v[i] = u(g); return v; };
    int jacobi_fails = 0;
    for (int n = 0; n < 200; ++n) {
        const auto x = vec7(), y = vec7(), z = vec7();
        const auto c = cross7(x, y);
        CHECK_THAT(c.dot(x), WithinAbs(0.0, 1e-12));
        CHECK_THAT(c.dot(y), WithinAbs(0.0, 1e-12));
        CHECK_THAT(c.dot(c), WithinAbs(x.dot(x) * y.dot(y) - x.dot(y) * x.dot(y), 1e-11));
        // Antisymmetric.
        const auto d = cross7(y, x);
        CHECK((c + d).norm() < 1e-12);
        // Unlike R^3's, it does not satisfy Jacobi: x(yz) + y(zx) + z(xy) != 0.
        const auto jac = cross7(x, cross7(y, z)) + cross7(y, cross7(z, x)) + cross7(z, cross7(x, y));
        if (jac.norm() > 1e-6) ++jacobi_fails;
    }
    CHECK(jacobi_fails > 190);
}

TEST_CASE("S^1, S^3 and S^7 have a global orthonormal tangent frame, polynomial in the point",
          "[octonion][dimensions]") {
    std::mt19937_64 g(14);
    std::normal_distribution<double> n01(0, 1);
    const auto check = [&](auto sphere, auto sample_point) {
        constexpr std::size_t N = decltype(sphere)::dimension;
        double worst = 0;
        for (int k = 0; k < 100; ++k) {
            const auto p = sample_point();
            const auto frame = sphere.tangent_frame(p);
            for (std::size_t i = 0; i < N; ++i) {
                worst = std::max(worst, std::abs(frame[i].norm() - 1.0));           // unit
                worst = std::max(worst, std::abs(frame[i].dot(p)));                 // tangent
                for (std::size_t j = i + 1; j < N; ++j)
                    worst = std::max(worst, std::abs(frame[i].dot(frame[j])));      // orthogonal
            }
            // Continuous: a nearby point has a nearby frame.
            auto q = p;
            q[0] += 1e-6;
            q = sphere.project(q);
            const auto frame2 = sphere.tangent_frame(q);
            for (std::size_t i = 0; i < N; ++i) worst = std::max(worst, (frame[i] - frame2[i]).norm() - 1e-5);
        }
        return worst;
    };
    const auto unit = [&](auto make) { return [&, make] { auto v = make(); return decltype(v){v * (1.0 / v.norm())}; }; };   // a Vec, not the expression
    CHECK(check(Sphere<1>{}, unit([&] { return Vec<double, 2>{n01(g), n01(g)}; })) < 1e-13);
    CHECK(check(Sphere<3>{}, unit([&] { return Vec<double, 4>{n01(g), n01(g), n01(g), n01(g)}; })) < 1e-13);
    CHECK(check(Sphere<7>{}, unit([&] {
        return Vec<double, 8>{n01(g), n01(g), n01(g), n01(g), n01(g), n01(g), n01(g), n01(g)}; })) < 1e-13);
}

TEST_CASE("A derivative through an octonion product reaches the same value by two routes", "[octonion][dual]") {
    // d/dt |a(t) b(t)| computed through the product, and through |a| |b|.
    using D = Dual<double>;
    std::mt19937_64 g(15);
    std::uniform_real_distribution<double> u(-1.5, 1.5);
    for (int n = 0; n < 50; ++n) {
        Octonion<D> a, b;
        for (std::size_t k = 0; k < 8; ++k) {
            a[k] = D(u(g)); a[k].deriv = u(g);
            b[k] = D(u(g)); b[k].deriv = u(g);
        }
        const D through_product = (a * b).norm();
        const D product_of_norms = a.norm() * b.norm();
        CHECK_THAT(through_product.value, WithinAbs(product_of_norms.value, 1e-12));
        CHECK_THAT(through_product.deriv, WithinAbs(product_of_norms.deriv, 1e-11));
    }
}

#if SPATIUM_HAS_BOOST_MULTIPRECISION
TEST_CASE("The octonion identities hold to fifty digits", "[octonion][precision]") {
    using R = Real50;
    std::mt19937_64 g(16);
    std::uniform_real_distribution<double> u(-1.5, 1.5);
    const auto make = [&] { return Octonion<R>{R(u(g)), R(u(g)), R(u(g)), R(u(g)), R(u(g)), R(u(g)), R(u(g)), R(u(g))}; };
    const R tol = R("1e-45");
    for (int n = 0; n < 30; ++n) {
        const auto a = make(), b = make(), c = make();
        CHECK((abs(R((a * b).norm_squared() - a.norm_squared() * b.norm_squared())) < tol));
        CHECK(((a * (a * b) - (a * a) * b).norm() < tol));
        CHECK((((a * b) * (c * a) - a * ((b * c) * a)).norm() < tol));
    }
}
#endif
