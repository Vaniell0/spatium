#include <catch2/catch_test_macros.hpp>
#include <spatium/core/interval.hpp>
#include <spatium/core/precision.hpp>

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
