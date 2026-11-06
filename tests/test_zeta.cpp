// The proved zeta: its values where they are known in closed form, two routes to one number that
// must agree, and the zeros on the critical line below 100 against those mpmath gives to 25 digits.

#include <catch2/catch_test_macros.hpp>
#include <rsc/spatium/zeta_request.hpp>
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
    const auto zeros = rsc::zeta_request(R"({"what": "zeros", "t0": 0, "t1": 30})");
    CHECK(zeros.find("zeros")->as_array().size() == 3);
    CHECK(zeros.find("complete")->as_bool());
    CHECK(rsc::zeta_request(R"({"what": "hardy", "t0": 5, "t1": 1})").find("refused") != nullptr);
}
