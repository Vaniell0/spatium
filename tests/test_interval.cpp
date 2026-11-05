#include <catch2/catch_test_macros.hpp>
#include <spatium/core/interval.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/algebra/polynomial.hpp>
#include <spatium/algebra/roots.hpp>

#include <cmath>
#include <format>
#include <string>

using namespace spatium;
using I = Interval<double>;
using R = Real<60>;

namespace {

// One formula, written once over any scalar: through every operation an interval has.
template<class T>
T formula(const T& x) {
    using std::sin, std::cos, std::exp, std::log, std::sqrt, std::pow, std::abs;
    return sin(x) * exp(x / T(3)) - sqrt(abs(x) + T(1)) / (x * x + T(2)) + log(x * x + T(1)) * cos(T(2) * x) +
           pow(x, T(3)) - pow(abs(x) + T(1), x / T(5));
}

// Whether the sixty-digit value is inside the interval, compared at sixty digits.
bool holds(const I& x, const R& truth) { return R(x.lo) <= truth && truth <= R(x.hi); }

} // namespace

TEST_CASE("an interval holds the value of what was written, against sixty digits", "[interval]") {
    // The decimals -4.000 ... 4.000 by thousandths, each read from its text into both.
    int held = 0, cases = 0, midpoint_exact = 0;
    for (int k = -4000; k <= 4000; k += 7) {
        const std::string text = std::format("{}.{:03}", k / 1000, std::abs(k % 1000));
        const std::string written = k < 0 && k > -1000 ? "-" + text : text;
        const I x(written.c_str());
        const R truth = formula(R(written));
        const I v = formula(x);
        ++cases;
        held += holds(v, truth);
        midpoint_exact += R(v.mid()) == truth;
        CHECK(v.width() < 1e-12 * (std::abs(v.mid()) + 1));
    }
    // a double's sine near its zeros: glibc reduces the argument exactly, so a few units hold it
    for (const char* x : {"3.141592653589793", "6.283185307179586", "-9.42477796076938", "31.41592653589793", "1e10"}) {
        INFO(x);
        using std::sin, std::cos;
        CHECK(holds(sin(I(x)), sin(R(x))));
        CHECK(holds(cos(I(x)), cos(R(x))));
    }
    CHECK(held == cases);
    // The check has teeth: the double the same formula gives is not the value it is held against.
    CHECK(midpoint_exact == 0);
}

TEST_CASE("fifty digits hold too, against a hundred and twenty, near a sine's zeros most of all", "[interval]") {
    using I50 = Interval<Real<50>>;
    using R120 = Real<120>;
    int held = 0, cases = 0;
    const auto check = [&](const std::string& written, auto&& f) {
        const I50 v = f(I50(written.c_str()));
        const R120 truth = f(R120(written));
        ++cases;
        const bool in = R120(v.lo) <= truth && truth <= R120(v.hi);
        held += in;
        if (!in) FAIL_CHECK(written);
    };
    // the decimals of pi and its multiples cut short, where sin is the digits cut off
    for (const char* x : {"3.141592653589793", "6.283185307179586", "3.14159265358979323846264338", "-9.42477796076938",
                          "1.5707963267948966", "31.41592653589793", "0.000001", "40", "-40"}) {
        using std::sin, std::cos, std::exp, std::log;
        check(x, [](const auto& t) { return sin(t); });
        check(x, [](const auto& t) { return cos(t); });
        check(x, [](const auto& t) { return exp(t); });
        check(x, [](const auto& t) { return log(abs(t)); });
    }
    for (int k = -3000; k <= 3000; k += 37) {
        const std::string text = std::format("{}{}.{:03}", k < 0 ? "-" : "", std::abs(k) / 1000, std::abs(k) % 1000);
        check(text, [](const auto& t) { return formula(t); });
    }
    CHECK(held == cases);
}

TEST_CASE("a decimal enters as written: exact when the double is it, else between its neighbours", "[interval]") {
    const I half("0.5"), tenth("0.1"), big("123456789"), small("-2.75e-3");
    CHECK(half.lo == half.hi);
    CHECK(big.lo == big.hi);
    CHECK(tenth.lo < tenth.hi);
    CHECK(holds(tenth, R("0.1")));
    CHECK(holds(small, R("-0.00275")));
    // a double read from its own decimal is the same double at both ends only when it is that decimal
    CHECK(I(0.1).lo == I(0.1).hi);
    CHECK_FALSE(holds(I(0.1), R("0.1")));
}

TEST_CASE("a comparison of overlapping intervals is undecided and says so", "[interval]") {
    const I a{1, 2}, b{3, 4}, c{1, 3}, d{2, 4};
    {
        Undecided w;
        CHECK(a < b);
        CHECK_FALSE(b < a);
        CHECK(a <= I{2, 5});
        CHECK(a != b);
        CHECK_FALSE(w.raised());
    }
    CHECK(less(c, d) == Truth::Undecided);
    CHECK(equal(I{2}, I{2}) == Truth::Yes);
    CHECK(equal(c, d) == Truth::Undecided);
    Undecided w;
    (void)(c < d);
    CHECK(w.raised());
    // a root reaching below the domain is a question part of the interval has no answer to
    Undecided root;
    CHECK(sqrt(I{-1, 4}).hi >= 2);
    CHECK(root.raised());
}

TEST_CASE("a power is not a product: an even power keeps its sign", "[interval]") {
    const I x{-1, 2};
    const I square = pow(x, I{2}), product = x * x;
    CHECK(square.lo == 0);
    CHECK(square.hi >= 4);
    CHECK(product.lo <= -2);
    const I cube = pow(x, I{3});
    CHECK(cube.lo <= -1);
    CHECK(cube.hi >= 8);
}

TEST_CASE("a sine reaches one where the interval holds a crest, and is tight where it does not", "[interval]") {
    const I s = sin(I{0, 4});
    CHECK(s.hi == 1);
    CHECK(s.lo < std::sin(4.0));
    const I t = sin(I{0.1, 0.2});
    CHECK(t.lo <= std::sin(0.1));
    CHECK(t.hi >= std::sin(0.2));
    CHECK(t.width() < 0.1);
    CHECK(cos(I{3, 3.3}).lo == -1);
    CHECK(sin(I{-1e300, 1e300}).hi == 1);
}

TEST_CASE("a divisor holding zero bounds nothing, and one that is zero gives nothing", "[interval]") {
    const I q = I{1} / I{-1, 1};
    CHECK(std::isinf(q.lo));
    CHECK(std::isinf(q.hi));
    CHECK((I{1} / I{0}).is_empty());
    const I r = I{1} / I{3};
    CHECK(holds(r, R(1) / R(3)));
}

TEST_CASE("an interval proves the digits its two ends round to alike", "[interval]") {
    const Proved third = proved(I{1} / I{3});
    CHECK(third.digits == 15);
    CHECK(third.text == "0.333333333333333");
    // one number is exact, written in its shortest digits
    const Proved exact = proved(I("0.5"));
    CHECK(exact.exact);
    CHECK(exact.text == "0.5");
    CHECK(exact.digits == 1);
    const Proved four = proved(I{2} + I{2});
    CHECK(four.exact);
    CHECK(four.text == "4");
    CHECK(proved(I("1.25") * I{8}).text == "10");
    CHECK(proved(sqrt(I{2.25})).text == "1.5");
    CHECK_FALSE(proved(I("0.1") * I{3}).exact);
    // a whole power of one number is a product: exact where it is, 10^20 among them
    CHECK(proved(pow(I{10}, I{20})).text == "100000000000000000000");
    CHECK(proved(pow(I{2}, I{-3})).text == "0.125");
    CHECK_FALSE(proved(pow(I{3}, I{40})).exact);
    CHECK(holds(pow(I{3}, I{40}), pow(R(3), 40)));
    // about zero there is no digit, only how far from zero it reaches
    const Proved zero = proved(sin(I("3.141592653589793")) - sin(I("3.141592653589793")));
    CHECK(zero.digits == 0);
    CHECK(zero.bound < 1e-15);
    // fifty digits deep: the same third over Real<50>
    const Proved deep = proved(Interval<Real<50>>{1} / Interval<Real<50>>{3});
    CHECK(deep.digits >= 48);
    CHECK(deep.text.starts_with("0.33333333333333333333333333333333333333333333333"));
}

TEST_CASE("a statement on a box: proved, refuted with a region, or undecided where it is", "[interval][prove_on]") {
    using B1 = Box<double, 1>;
    // x*x + 1 > 0 on [-10, 10]: true, but x*x of [-10, 10] is [-100, 100], so only cut pieces show it
    const auto positive = prove_on(B1{I{-10, 10}}, [](const B1& b) { return less(I{0}, b[0] * b[0] + I{1}); });
    CHECK(positive.verdict == Truth::Yes);
    CHECK(positive.evaluations > 1);
    // x != 0 on [-1, 1]: the piece holding zero is never excluded
    const auto nonzero = prove_on(B1{I{-1, 1}}, [](const B1& b) {
        return b[0].hi < 0 || 0 < b[0].lo ? Truth::Yes : Truth::Undecided;
    });
    CHECK(nonzero.verdict == Truth::Undecided);
    int holding_zero = 0;
    for (const auto& [b, t] : nonzero.pieces) {
        if (t == Truth::Undecided) holding_zero += b[0].contains(0);
        if (b[0].contains(0)) CHECK(t == Truth::Undecided);
    }
    CHECK(holding_zero >= 1);
    // x*x < 50 on [-10, 10] is false past the square root of 50: refuted, on pieces wholly past it
    const auto below = prove_on(B1{I{-10, 10}}, [](const B1& b) { return less(b[0] * b[0], I{50}); });
    CHECK(below.verdict == Truth::No);
    for (const auto& [b, t] : below.pieces)
        if (t == Truth::No) CHECK(std::min(b[0].lo * b[0].lo, b[0].hi * b[0].hi) >= 50 - 1e-12);
    // a comparison taken by guess is no answer, whatever it said
    const auto guessed = prove_on(B1{I{-1, 1}}, [](const B1& b) { return b[0] < I{0.5} ? Truth::Yes : Truth::No; },
                                  Budget{.evaluations = 64});
    for (const auto& [b, t] : guessed.pieces)
        if (b[0].contains(0.5) && b[0].lo < 0.5 && 0.5 < b[0].hi) CHECK(t == Truth::Undecided);
}

TEST_CASE("two dimensions: (x - y)^2 >= 0 written out is proved off the diagonal, undecided along it", "[interval][prove_on]") {
    using B2 = Box<double, 2>;
    // x*x + y*y - 2*x*y reads x and y twice each, so its interval is wider than the square's and a
    // strip about x = y stays undecided; it narrows as the budget grows, slowly (measured: undecided
    // area 1.47, 1.06, 0.73 of 4 at 250, 1000, 4000 evaluations)
    double area[3] = {};
    const std::size_t budgets[3] = {250, 1000, 4000};
    for (int n = 0; n < 3; ++n) {
        const auto square = prove_on(B2{I{-1, 1}, I{-1, 1}}, [](const B2& b) {
            const I& x = b[0];
            const I& y = b[1];
            return less_equal(I{0}, x * x + y * y - I{2} * x * y);
        }, Budget{.evaluations = budgets[n]});
        CHECK(square.verdict == Truth::Undecided);
        double total = 0;
        for (const auto& [b, t] : square.pieces) {
            const double a = b[0].width() * b[1].width();
            total += a;
            CHECK(t != Truth::No);
            if (t == Truth::Undecided) area[n] += a;
            // a piece the diagonal crosses holds points where the square is zero: never proved
            if (b[0].lo < b[1].hi && b[1].lo < b[0].hi) CHECK(t == Truth::Undecided);
        }
        CHECK(total == 4);
    }
    CHECK(area[0] > area[1]);
    CHECK(area[1] > area[2]);
    CHECK(area[2] < 0.75);
    // written as a square it is one evaluation
    const auto direct = prove_on(B2{I{-1, 1}, I{-1, 1}}, [](const B2& b) { return less_equal(I{0}, pow(b[0] - b[1], I{2})); });
    CHECK(direct.verdict == Truth::Yes);
    CHECK(direct.evaluations == 1);
}

namespace {
// The whole range of the doubles, as one box.
const I doubles{-std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
// A constant of whichever arithmetic f is asked in: an interval, or a dual number of intervals.
template<class X> X k(const X&, const I& c) { return X(c); }
} // namespace

TEST_CASE("every root on a box: each proved alone in its interval, against the closed form and sixty digits", "[interval][roots]") {
    int quadratics = 0;
    for (int a = 1; a <= 3; ++a)
        for (int b = -7; b <= 7; b += 2)
            for (int c = -9; c <= 3; c += 3) {
                if (b * b - 4 * a * c <= 0) continue;
                INFO(a << "x^2 + " << b << "x + " << c);
                const auto found = roots_on([&](const auto& x) { return k(x, a) * x * x + k(x, b) * x + k(x, c); }, I{-100, 100});
                REQUIRE(found.complete());
                REQUIRE(found.roots.size() == 2);
                // the roots at sixty digits, as q / a and c / q, so a root that is zero is zero
                const R q = (R(-b) - (b < 0 ? R(-1) : R(1)) * sqrt(R(b * b - 4 * a * c))) / R(2);
                const R one = q / R(a), other = R(c) / q;
                const R exact[2] = {one < other ? one : other, one < other ? other : one};
                // the formula's roots, in doubles, against the proved intervals: within a unit of them
                const auto closed = solve_quadratic(double(a), double(b), double(c));
                REQUIRE(closed.size() == 2);
                const double formula[2] = {std::min(closed[0].re, closed[1].re), std::max(closed[0].re, closed[1].re)};
                for (int r = 0; r < 2; ++r) {
                    const I& root = found.roots[r];
                    CHECK(holds(root, exact[r]));
                    CHECK((c == 0 && exact[r] == 0 ? root.width() < 1e-300 : proved(root).digits >= 14));
                    CHECK(std::nextafter(root.lo, -1e300) <= formula[r]);
                    CHECK(formula[r] <= std::nextafter(root.hi, 1e300));
                }
                ++quadratics;
            }
    CHECK(quadratics == 106);
    // a root on a cut: [-8, 8] is cut at 0, 4, 2, so 2 is the end of two pieces, and found once
    const auto cut = roots_on([](const auto& x) { return x * x - k(x, 4); }, I{-8, 8});
    REQUIRE(cut.complete());
    REQUIRE(cut.roots.size() == 2);
    CHECK(cut.roots[1].contains(2));
    CHECK(cut.roots[0].contains(-2));
    CHECK(proved(cut.roots[1]).digits == 16);
}

TEST_CASE("roots over the doubles' whole range, and none where there are none", "[interval][roots]") {
    using std::pow, std::sin, std::sqrt;
    // 2^x = 8: one root, 3, on the whole line
    const auto power = roots_on([](const auto& x) { return pow(k(x, 2), x) - k(x, 8); }, doubles);
    CHECK(power.complete());
    REQUIRE(power.roots.size() == 1);
    CHECK(power.roots[0].contains(3));
    CHECK(proved(power.roots[0]).digits >= 15);
    CHECK(power.evaluations < 2000);
    // x^2 + 1 has none: the proof is that every piece of the line was excluded
    const auto none = roots_on([](const auto& x) { return x * x + k(x, 1); }, doubles);
    CHECK(none.complete());
    CHECK(none.roots.empty());
    // its pieces, counted: with no root, every one is proved to hold none
    CHECK(none.pieces > 0);
    CHECK(none.excluded == none.pieces);
    // sin x = 1/2 on [0, 10]: pi/6, 5pi/6, 13pi/6, 17pi/6
    const auto wave = roots_on([](const auto& x) { return sin(x) - k(x, I("0.5")); }, I{0, 10});
    CHECK(wave.complete());
    REQUIRE(wave.roots.size() == 4);
    // the case it must fail: with roots, some pieces hold one and are not excluded
    CHECK(wave.excluded + wave.roots.size() <= wave.pieces);
    const R pi = acos(R(-1));
    const int sixths[4] = {1, 5, 13, 17};
    for (int r = 0; r < 4; ++r) {
        CHECK(holds(wave.roots[r], pi * sixths[r] / 6));
        CHECK(proved(wave.roots[r]).digits >= 14);
    }
    // a root at the edge of a square root's domain: the part of a piece below zero has no value and no root
    const auto edge = roots_on([](const auto& x) { return sqrt(x) - k(x, 1); }, I{-1, 4});
    CHECK(edge.complete());
    REQUIRE(edge.roots.size() == 1);
    CHECK(edge.roots[0].contains(1));
    CHECK(proved(edge.roots[0]).digits >= 14);
}

TEST_CASE("a kink is no root: abs(x - 2) = 1 has two, and the corner between them is excluded", "[interval][roots]") {
    using std::abs;
    // the derivative of abs guesses a branch at the corner; the value alone excludes it, and is asked first
    const auto kink = roots_on([](const auto& x) { return abs(x - k(x, 2)) - k(x, 1); }, I{0, 5});
    CHECK(kink.complete());
    REQUIRE(kink.roots.size() == 2);
    CHECK(kink.roots[0].contains(1));
    CHECK(kink.roots[1].contains(3));
}

TEST_CASE("a double root stays undecided, and says where it may be", "[interval][roots]") {
    const auto twice = roots_on([](const auto& x) { return x * x; }, I{-1, 1});
    CHECK_FALSE(twice.complete());
    CHECK(twice.roots.empty());
    REQUIRE(twice.undecided.size() == 1);
    CHECK(twice.undecided[0].contains(0));
    CHECK(twice.undecided[0].width() < 1e-8);
    // the same square moved off zero, beside a simple root: the simple one is proved, the double one not
    const auto both = roots_on([](const auto& x) { return (x - k(x, 1)) * (x - k(x, 1)) * (x + k(x, 2)); }, I{-5, 5});
    REQUIRE(both.roots.size() == 1);
    CHECK(both.roots[0].contains(-2));
    REQUIRE(both.undecided.size() == 1);
    CHECK(both.undecided[0].contains(1));
}
