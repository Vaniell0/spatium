// The proved zeta: its values where they are known in closed form, two routes to one number that
// must agree, and the zeros on the critical line below 100 against those mpmath gives to 25 digits.

#include <catch2/catch_test_macros.hpp>
#include <spatium/answer/zeta_request.hpp>
#include <spatium/algebra/zeta.hpp>

#include <cmath>
#include <numbers>

using namespace spatium;
using I = Interval<double>;
using C = Complex<I>;

namespace {

// the ordinates of the first 29 zeros, from mpmath's zetazero to 25 digits
constexpr double mpmath_zeros[] = {
    14.13472514173469379045725, 21.02203963877155499262848, 25.01085758014568876321379,
    30.4248761258595132103119, 32.93506158773918969066237, 37.58617815882567125721776,
    40.91871901214749518739813, 43.32707328091499951949612, 48.00515088116715972794247,
    49.77383247767230218191678, 52.9703214777144606441473, 56.44624769706339480436776,
    59.34704400260235307965365, 60.8317785246098098442599, 65.11254404808160666087505,
    67.07981052949417371447883, 69.54640171117397925292686, 72.06715767448190758252211,
    75.70469069908393316832692, 77.14484006887480537268266, 79.33737502024936792276359,
    82.91038085408603018316484, 84.73549298051705010573531, 87.42527461312522940653167,
    88.80911120763446542368235, 92.49189927055848429625973, 94.65134404051988696659793,
    95.87063422824530975874103, 98.83119421819369223332442};

int sign(const I& x) { return x.lo > 0 ? 1 : x.hi < 0 ? -1 : 0; }

} // namespace

TEST_CASE("zeta holds its closed forms, narrowly", "[zeta]") {
    using std::numbers::pi;
    const struct { double s, value; } known[] = {{2, pi * pi / 6}, {4, pi * pi * pi * pi / 90}, {0, -0.5}, {-1, -1.0 / 12}, {0.5, -1.4603545088095868}};
    for (const auto& k : known) {
        const C z = zeta(C{I(k.s)});
        CHECK(z.re.contains(k.value));
        CHECK(z.im.contains(0.0));
        CHECK(z.re.width() < 1e-10);
    }
    // and not a number next to it: pi^2/6 differs from 1.64493406685 in the twelfth digit
    CHECK_FALSE(zeta(C{I(2.0)}).re.contains(1.64493406685));
    // atan, which theta reads, holds pi/4 at 1
    CHECK(atan(I(1.0)).contains(std::numbers::pi / 4));
}

TEST_CASE("zeta's two routes agree: conjugate symmetry, and Z real on the critical line", "[zeta]") {
    const C a = zeta(C{I(0.3), I(20.0)}), b = zeta(C{I(0.3), I(-20.0)});
    CHECK(a.re.lo <= b.re.hi);
    CHECK(b.re.lo <= a.re.hi);
    CHECK(a.im.lo <= -b.im.lo);
    CHECK(-b.im.hi <= a.im.hi);
    // e^{i theta} zeta(1/2 + it) is real: its imaginary part, from theta's Stirling series and
    // zeta's Euler-Maclaurin sum, holds zero, which a theta off by 1e-9 would not
    for (double t : {10.0, 100.0, 1000.0}) {
        const I th = riemann_siegel_theta(I(t));
        const C z = zeta(C{I(0.5), I(t)});
        CHECK((sin(th) * z.re + cos(th) * z.im).contains(0.0));
        CHECK_FALSE((sin(th + I(1e-9)) * z.re + cos(th + I(1e-9)) * z.im).contains(0.0));
    }
}

TEST_CASE("left of the line zeta is the functional equation's: its values at 1 - 2n, and the trivial zeros", "[zeta]") {
    // zeta(1 - 2n) = -B_2n / 2n, from the table of fractions, which the sum to the right reads too;
    // Euler-Maclaurin alone gave 392.98 at -13 and a NaN past -58
    for (int n = 1; n <= 30; ++n) {
        const I value = -zeta_detail::bernoulli_number<double>(n) / I(2.0 * n);
        const C z = zeta(C{I(1.0 - 2 * n)});
        INFO(n);
        CHECK(z.re.lo <= value.hi);
        CHECK(value.lo <= z.re.hi);
        CHECK(z.re.width() < 1e-9 * std::abs(value.mid()));
    }
    // each trivial zero a change of sign along the real axis, proved
    for (int n = 1; n <= 30; ++n) {
        INFO(n);
        CHECK(sign(zeta(C{I(-2.0 * n - 1e-3)}).re) == -sign(zeta(C{I(-2.0 * n + 1e-3)}).re));
        CHECK(sign(zeta(C{I(-2.0 * n + 1e-3)}).re) != 0);
    }
    // both routes where both run, 0 > Re s > -58: they meet, and the reflected one is the narrower
    for (const double sigma : {-0.5, -2.0, -5.0})
        for (const double t : {0.3, 20.0, 1000.0}) {
            const C s{I(sigma), I(t)};
            const C a = zeta_detail::euler_maclaurin(s), b = zeta_detail::reflected(s);
            INFO(sigma << " " << t);
            CHECK(a.re.lo <= b.re.hi);
            CHECK(b.re.lo <= a.re.hi);
            CHECK(a.im.lo <= b.im.hi);
            CHECK(b.im.lo <= a.im.hi);
            CHECK(b.re.width() <= a.re.width());
        }
    // the case it must fail: -1/12 is zeta(-1), not zeta(-1.001)
    CHECK_FALSE(zeta(C{I(-1.001)}).re.contains(-1.0 / 12));
}

TEST_CASE("zeta in doubles, the colour map's, lies in the proved interval", "[zeta]") {
    for (const auto& [x, y] : {std::pair{-0.5, 3.0}, {0.5, 14.0}, {1.5, 40.0}, {2.0, 700.0}}) {
        const C z = zeta(C{I(x), I(y)});
        const Complex<double> d = zeta(Complex<double>{x, y});
        const double slack = 1e-12 * (1 + std::hypot(d.re, d.im));
        CHECK(z.re.lo - slack <= d.re);
        CHECK(d.re <= z.re.hi + slack);
        CHECK(z.im.lo - slack <= d.im);
        CHECK(d.im <= z.im.hi + slack);
    }
}

TEST_CASE("every zero below 100 proved, each holding mpmath's", "[zeta]") {
    const auto z = critical_zeros(0.0, 100.0);
    REQUIRE(z.complete);
    CHECK(z.expected == 29);
    REQUIRE(z.zeros.size() == 29);
    for (std::size_t i = 0; i < 29; ++i) {
        CHECK(z.zeros[i].contains(mpmath_zeros[i]));
        CHECK(z.zeros[i].width() < 1e-10);
    }
}

TEST_CASE("Gram's law first fails at g_126; Lehmer's pair is two changes of sign", "[zeta]") {
    const auto g = gram_points(0.0, 300.0);
    REQUIRE(g.front().n == -1);
    long first = -1;
    for (const auto& p : g)
        if (p.good == Truth::No) { first = p.n; break; }
    CHECK(first == 126);
    // Lehmer's zeros near 7005.0629 and 7005.1006: Z has one sign between them and the other outside
    const int out = sign(hardy_z(I(7005.04))), in = sign(hardy_z(I(7005.08)));
    CHECK(out != 0);
    CHECK(in == -out);
    CHECK(sign(hardy_z(I(7005.13))) == out);
}

TEST_CASE("the page's zeta requests answer, and refuse a window that is none", "[zeta]") {
    const auto zeros = answer::zeta_request(R"({"what": "zeros", "t0": 0, "t1": 30})");
    CHECK(zeros.find("zeros")->as_array().size() == 3);
    CHECK(zeros.find("complete")->as_bool());
    CHECK(answer::zeta_request(R"({"what": "hardy", "t0": 5, "t1": 1})").find("refused") != nullptr);
    // below 0 the mirror: the zeros at -t for those at t, the window across 0 counted as a whole
    const auto both = answer::zeta_request(R"({"what": "zeros", "t0": -30, "t1": 30})");
    const auto& z = both.find("zeros")->as_array();
    REQUIRE(z.size() == 6);
    CHECK(both.find("complete")->as_bool());
    // counted past the window on both sides alike: five zeros each below 33
    const auto& counted = both.find("counted")->as_array();
    CHECK(counted[0].as_number() == -counted[1].as_number());
    CHECK(both.find("expected")->as_number() == 10);
    CHECK(z[0].find("lo")->as_number() == -z[5].find("hi")->as_number());
    CHECK(z[2].find("hi")->as_number() < 0);
    CHECK(z[3].find("lo")->as_number() > 0);
}

TEST_CASE("zeta past t = 1e5 is still proved, and the double still lies in it", "[zeta]") {
    // each Euler-Maclaurin term from the one before: the rising product and the power of N, kept
    // apart, left the doubles here, and every value was NaN
    for (const double t : {1e5, 3e5}) {
        const auto proved = zeta(C{I(0.5), I(t)});
        const auto drawn = zeta(Complex<double>{0.5, t});
        INFO(t);
        REQUIRE(std::isfinite(proved.re.lo));
        CHECK(proved.re.hi - proved.re.lo < 1e-5);
        CHECK(std::abs(drawn.re - proved.re.mid()) < 1e-6);
        CHECK(std::abs(drawn.im - proved.im.mid()) < 1e-6);
    }
}

TEST_CASE("Z by Riemann-Siegel meets the proved Z of Euler-Maclaurin", "[zeta]") {
    // Psi at 1/4, where its quotient is 0 / 0, and at 1/2; symmetric about 1/2, its odd derivatives not
    CHECK(zeta_detail::psi_derivatives(0.25)[0] == 0.5);
    CHECK(std::abs(zeta_detail::psi_derivatives(0.5)[0] - std::cos(3 * std::numbers::pi / 8)) < 1e-15);
    const auto l = zeta_detail::psi_derivatives(0.3), r = zeta_detail::psi_derivatives(0.7);
    for (int n = 0; n <= 12; ++n) CHECK(std::abs(l[n] - (n % 2 ? -r[n] : r[n])) <= 1e-12 * (1 + std::abs(l[n])));
    // the coefficients from the table against those from the derivatives, between the nodes
    for (const double q : {0.0, 0.00012207, 0.24999, 0.25, 0.5, 0.731, 0.99999}) {
        const auto t = zeta_detail::riemann_siegel_coefficients(q), e = zeta_detail::riemann_siegel_coefficients_exact(q);
        for (int k = 0; k < 5; ++k) CHECK(std::abs(t[k] - e[k]) < 1e-12);
    }
    // the second route: the proved Z, whose interval the formula's value lies in up to its own error
    for (const double t : {1000.0, 5000.0, 7005.08, 1e5}) {
        const I z = hardy_z(I(t));
        INFO(t);
        CHECK(std::abs(riemann_siegel_z(t) - z.mid()) < 1e-10 + z.width());
    }
    // the case it must fail: Z at a neighbouring height
    CHECK(std::abs(riemann_siegel_z(1000.01) - hardy_z(I(1000.0)).mid()) > 1e-3);
    // Gabcke's bound on the cut holds the difference from the proved Z, and is no bound below 200
    for (const double t : {250.0, 1000.0, 7005.08, 1e5}) {
        const I z = hardy_z(I(t));
        INFO(t);
        CHECK(std::abs(riemann_siegel_z(t) - z.mid()) <= riemann_siegel_truncation(t) + z.width());
    }
    CHECK(riemann_siegel_truncation(1e6) < 1e-10);
    CHECK(std::isinf(riemann_siegel_truncation(150.0)));
}

TEST_CASE("zeta' and zeta'' by dual numbers: their known values, and the difference quotient", "[zeta]") {
    const auto half = zeta_jet({0.5, 0});
    CHECK(std::abs(half.first.re - -3.92264613920915172) < 1e-10);   // zeta'(1/2)
    const auto two = zeta_jet({2, 0});
    CHECK(std::abs(two.first.re - -0.93754825431584375) < 1e-11);    // zeta'(2)
    CHECK(std::abs(two.second.re - 1.98928023429890130) < 1e-10);    // zeta''(2)
    // the second route: a central difference of zeta in doubles, off the axis
    const Complex<double> s{0.7, 21.3};
    const double h = 1e-5;
    const auto q = (zeta(Complex<double>{s.re + h, s.im}) - zeta(Complex<double>{s.re - h, s.im})) / Complex<double>{2 * h, 0};
    const auto j = zeta_jet(s);
    CHECK(std::hypot(q.re - j.first.re, q.im - j.first.im) < 1e-8);
    // the case it must fail: the derivative is not the value
    CHECK(std::hypot(j.value.re - j.first.re, j.value.im - j.first.im) > 1e-3);
}

TEST_CASE("the zeros of zeta' below 50 are right of the critical line, the first where it is known", "[zeta]") {
    const auto zeros = zeta_prime_zeros({-1, 1}, {3, 50});
    REQUIRE(zeros.size() == 5);
    CHECK(std::abs(zeros.front().re - 2.46316186945) < 1e-8);
    CHECK(std::abs(zeros.front().im - 23.29832049) < 1e-7);
    for (const auto& z : zeros) {
        CHECK(z.re > 0.5);
        CHECK(std::hypot(zeta_jet(z).first.re, zeta_jet(z).first.im) < 1e-10);
    }
    // a start right of every zero leaves to the right, where zeta' is near -log 2 / 2^s
    const auto away = zeta_prime_newton({5, 10});
    CHECK_FALSE(away.settled);
    CHECK(away.escaped == 1);
}

TEST_CASE("the first moment of zeta on the line against Hardy and Littlewood's asymptotics", "[zeta]") {
    const double T = 50, gamma = std::numbers::egamma, pi = std::numbers::pi;
    const auto m1 = zeta_moment(1, T);
    const double main = T * std::log(T / (2 * pi)) + (2 * gamma - 1) * T;
    CHECK(std::abs(m1.value / main - 1) < 0.06);
    CHECK(m1.error < 1e-6);
    // M_2 is larger than M_1 squared over T (Cauchy-Schwarz), the second route's bound
    const auto m2 = zeta_moment(2, T);
    CHECK(m2.value * T >= m1.value * m1.value);
}
