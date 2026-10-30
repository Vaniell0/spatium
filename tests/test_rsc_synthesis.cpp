// Synthesis of particular analytic solutions on a region, and piecewise ones.
//
// Each claim is held against a second route: a found expression against the function it was
// asked for, at points it was not fitted on; the field the synthesis writes against the same
// field lowered to plain data (io/field_pod.hpp), so what is found is an IR artifact and not a
// closure; the deflection of a light ray by a Schwarzschild hole, the reference, by quadrature
// against its closed form in elliptic integrals. Each check gets a case it must fail: near the
// critical impact parameter the deflection diverges like a logarithm, which the operations of
// the IR do not hold, and the pieces there stop at the narrowest width allowed.

#include <rsc/spatium/synthesis.hpp>
#include <spatium/algebra/quadrature.hpp>
#include <spatium/io/field_pod.hpp>
#include <spatium/io/field_text.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <functional>
#include <numbers>
#include <random>

using Catch::Matchers::WithinAbs;

namespace {

// The deflection of a null ray of impact parameter b by a hole of mass 1 (b_c = 3 sqrt 3).
double closest_approach(double b) {
    double r = b;
    for (int i = 0; i < 100; ++i) { // the largest root of r^3 - b^2 (r - 2)
        const double step = (r * r * r - b * b * (r - 2)) / (3 * r * r - b * b);
        r -= step;
        if (std::abs(step) < 1e-15 * r)
            break;
    }
    return r;
}

// Darwin's closed form, through the elliptic integrals of <cmath>.
double deflection(double b) {
    const double r0 = closest_approach(b), q = std::sqrt((r0 - 2) * (r0 + 6));
    const double k = std::sqrt((q - r0 + 6) / (2 * q));
    const double psi = std::asin(std::sqrt((q - r0 + 2) / (q - r0 + 6)));
    return -std::numbers::pi + 4 * std::sqrt(r0 / q) * (std::comp_ellint_1(k) - std::ellint_1(k, psi));
}

// 2 int_0^u0 du / sqrt(1/b^2 - u^2 + 2 u^3) - pi, with u = u0 - s^2 taking the root at the end
// away, and the radicand written as a product so that nothing cancels.
double deflection_by_quadrature(double b) {
    const double u0 = 1 / closest_approach(b);
    const auto f = [&](double s) {
        const double u = u0 - s * s;
        return 2 / std::sqrt((u0 + u) - 2 * (u0 * u0 + u0 * u + u * u));
    };
    return 2 * spatium::integrate(f, spatium::Finite<double>{0, std::sqrt(u0)}).value - std::numbers::pi;
}

const double critical = 3 * std::sqrt(3.0);

std::size_t accepted(const std::vector<rsc::Solution>& pieces) {
    std::size_t n = 0;
    for (const auto& p : pieces)
        n += p.accepted;
    return n;
}

} // namespace

TEST_CASE("synthesis: an expression the operations can write is found exactly", "[rsc][synthesis]") {
    const rsc::SynthesisSpec spec{.lo = 0.1, .hi = 1, .tolerance = 1e-10,
                                  .unary = {rsc::Unary::Sqrt, rsc::Unary::Sin}};
    const auto f = [](double x) { return 3 * x * x + 2 * std::sin(x); };
    const auto s = rsc::synthesize(f, spec);
    REQUIRE(s.accepted);
    CHECK(s.terms == 2);
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> u(0.1, 1);
    const auto pod = spatium::io::build::lower(s.expr);
    REQUIRE(pod.has_value());
    for (int i = 0; i < 200; ++i) {
        const double x = u(rng);
        CHECK_THAT(s.expr(x, 0.0), WithinAbs(f(x), 1e-10 * f(x)));
        CHECK(spatium::io::build::interpret(*pod, x, 0.0) == s.expr(x, 0.0)); // the lowered form
    }
}

TEST_CASE("synthesis: a function the operations cannot write takes more pieces the tighter it is asked", "[rsc][synthesis]") {
    // log(1 + x) with no logarithm among the operations
    const auto f = [](double x) { return std::log1p(x); };
    std::size_t last = 0;
    for (double tol : {1e-3, 1e-5, 1e-7}) {
        INFO("tolerance " << tol);
        const rsc::SynthesisSpec spec{.lo = 0, .hi = 1, .tolerance = tol, .relative = false, .max_terms = 3};
        const auto pieces = rsc::synthesize_piecewise(f, spec, 1e-3);
        CHECK(accepted(pieces) == pieces.size());
        CHECK(pieces.size() >= last);
        last = pieces.size();
        // the pieces cover [0, 1] in order and each holds at points it was not certified on
        CHECK(pieces.front().lo == 0.0);
        CHECK(pieces.back().hi == 1.0);
        std::mt19937 rng(7);
        for (std::size_t i = 0; i < pieces.size(); ++i) {
            if (i > 0)
                CHECK(pieces[i].lo == pieces[i - 1].hi);
            std::uniform_real_distribution<double> u(pieces[i].lo, pieces[i].hi);
            for (int k = 0; k < 20; ++k) {
                const double x = u(rng);
                CHECK(std::abs(pieces[i].expr(x, 0.0) - f(x)) <= 2 * tol);
            }
        }
    }
    CHECK(last > 4); // measured: 1, 4 and 16 pieces
}

TEST_CASE("synthesis: the deflection of light, its reference held two ways", "[rsc][synthesis][relativity]") {
    for (double b : {critical * (1 + 1e-4), critical * 1.1, 7.0, 30.0, 1000.0})
        CHECK_THAT(deflection_by_quadrature(b), WithinAbs(deflection(b), 1e-12 * deflection(b) + 1e-14));
    // the weak field: 4 / b + 15 pi / (4 b^2) + 128 / (3 b^3) + O(b^-4); at b = 1e4 the closed
    // form subtracts two elliptic integrals that agree to four digits, so it holds about 1e-13
    CHECK_THAT(deflection(1e4), WithinAbs(4e-4 + 15 * std::numbers::pi / 4e8 + 128 / 3e12, 1e-12));
}

TEST_CASE("synthesis: the deflection of light in pieces of x = 1/b", "[rsc][synthesis][relativity]") {
    const auto alpha = [](double x) { return deflection(1 / x); };
    const rsc::SynthesisSpec spec{.lo = 1e-4, .hi = (1 / critical) * (1 - 1e-4), .tolerance = 1e-2, .max_terms = 4};
    const auto pieces = rsc::synthesize_piecewise(alpha, spec, 1e-5);
    INFO(pieces.size() << " pieces");
    CHECK(accepted(pieces) == pieces.size());
    // the pieces crowd towards the critical parameter, where the deflection diverges: more than
    // half of them lie within b < 6 (measured: 11 pieces, the first from b = 20.8 to 10^4)
    std::size_t near = 0;
    for (const auto& p : pieces)
        near += 1 / p.lo < 6;
    CHECK(2 * near > pieces.size());
    CHECK(1 / pieces.front().hi < 21);
}

TEST_CASE("synthesis: next to the critical parameter the operations run out", "[rsc][synthesis][relativity]") {
    // The deflection diverges like -log(b / b_c - 1) there and the logarithm is not among the
    // operations: each decade nearer b_c costs more pieces at the same tolerance, and an end
    // nearer than the narrowest piece allowed is not reached at all.
    const auto alpha = [](double x) { return deflection(1 / x); };
    const double xc = 1 / critical;
    std::size_t last = 0;
    for (double gap : {1e-2, 1e-3, 1e-4}) {
        INFO("gap " << gap);
        const auto pieces = rsc::synthesize_piecewise(
            alpha, {.lo = xc * (1 - 3e-2), .hi = xc * (1 - gap), .tolerance = 1e-4, .max_terms = 4}, 1e-9);
        CHECK(accepted(pieces) == pieces.size());
        CHECK(pieces.size() > last);
        last = pieces.size();
    }
    // 2e-8 in x from the divergence, with no piece narrower than 1e-5: the last one fails
    const auto out = rsc::synthesize_piecewise(
        alpha, {.lo = xc * (1 - 3e-2), .hi = xc * (1 - 1e-7), .tolerance = 1e-4, .max_terms = 4}, 1e-5);
    CHECK(!out.back().accepted);
    CHECK(out.back().hi == xc * (1 - 1e-7));
}

TEST_CASE("synthesis: with the logarithm among the operations the critical neighbourhood is one piece", "[rsc][synthesis][relativity]") {
    // The job's own variable there is y = b / b_c - 1 (b_c is a constant of the job, not a hint),
    // and the deflection is -log(y) and a constant to leading order. Measured from y = 1e-6 to
    // 1e-2: at 1e-4 one piece of 19 nodes against five of 106, with sets of three terms searched
    // whole (eight of 120 before them); before them, at 1e-6 four of 75 against 24 of 537, and at
    // 1e-8 six of 161 against 62 of 1555.
    const auto alpha = [](double y) { return deflection(critical * (1 + y)); };
    const rsc::SynthesisSpec without{.lo = 1e-6, .hi = 1e-2, .tolerance = 1e-4, .max_terms = 4};
    rsc::SynthesisSpec with = without;
    with.unary = {rsc::Unary::Sqrt, rsc::Unary::Log};
    const auto a = rsc::synthesize_piecewise(alpha, without, 1e-9);
    const auto b = rsc::synthesize_piecewise(alpha, with, 1e-9);
    CHECK(accepted(a) == a.size());
    CHECK(accepted(b) == b.size());
    CHECK(b.size() == 1);
    CHECK(a.size() >= 4);
    std::size_t na = 0, nb = 0;
    for (const auto& p : a)
        na += p.nodes;
    for (const auto& p : b)
        nb += p.nodes;
    CHECK(na > 5 * nb);
}

TEST_CASE("synthesis: the deflection of light at 1e-4 and 1e-6", "[.deep][rsc][synthesis][relativity]") {
    // measured: 19 pieces, all accepted, at 1e-4; 59 pieces, the last one not, at 1e-6
    const auto alpha = [](double x) { return deflection(1 / x); };
    const double lo = 1e-4, hi = (1 / critical) * (1 - 1e-4);
    const auto mid = rsc::synthesize_piecewise(alpha, {.lo = lo, .hi = hi, .tolerance = 1e-4, .max_terms = 4}, 1e-5);
    CHECK(accepted(mid) == mid.size());
    const auto fine = rsc::synthesize_piecewise(alpha, {.lo = lo, .hi = hi, .tolerance = 1e-6, .max_terms = 4}, 1e-5);
    CHECK(!fine.back().accepted);
    CHECK(fine.size() > mid.size());
}

TEST_CASE("synthesis: an answer that x itself writes is written in x, and reads as text", "[rsc][synthesis]") {
    // On [0, 10] the interval's variable t = (x - 5) * 0.2 fits a line as well as x does; at equal
    // error the form of fewer IR nodes is kept, so the answer reads in x.
    const auto s = rsc::synthesize([](double x) { return 1 + x / 6; }, {.lo = 0, .hi = 10, .tolerance = 1e-12, .relative = false});
    REQUIRE(s.accepted);
    const auto text = spatium::io::build::to_text(s.expr, "p");
    REQUIRE(text.has_value());
    INFO(*text);
    CHECK(text->find("p") != std::string::npos);
    CHECK(text->find("- 5") == std::string::npos); // no t in it
    // the text of a known field, and a field that has no text
    using F = spatium::io::build::Field<double>;
    CHECK(*spatium::io::build::to_text(sqrt(F::u()) + F(2.0) * log(F::v()), "x", "y") == "sqrt(x) + 2 * log(y)");
    CHECK_FALSE(spatium::io::build::to_text(F([](double, double) { return 0.0; })).has_value());
}

TEST_CASE("synthesis: a coefficient that is a small fraction is written as one, if the answer still holds",
          "[rsc][synthesis]") {
    using rsc::synthesis_detail::simplest_fraction;
    CHECK(simplest_fraction(0.16666666666666663) == std::pair{1LL, 6LL});
    CHECK(simplest_fraction(1.0000000000000002) == std::pair{1LL, 1LL});
    CHECK(simplest_fraction(-2.5) == std::pair{-5LL, 2LL});
    CHECK(simplest_fraction(0.333333333) == std::pair{1LL, 3LL});   // 1/3 to nine digits
    // the cases it must fail: a number that is no small fraction to nine digits
    CHECK_FALSE(simplest_fraction(std::sqrt(2.0)).has_value());
    CHECK_FALSE(simplest_fraction(0.3333).has_value());

    // the prism's volume, 1 + p/6: the fitted 0.16666666666666663 and 1.0000000000000002 become
    // 1/6 and 1, and the answer is certified again in that form
    const rsc::SynthesisSpec spec{.lo = 0, .hi = 10, .tolerance = 1e-12, .relative = false};
    const auto s = rsc::synthesize([](double x) { return 1 + x / 6; }, spec);
    REQUIRE(s.accepted);
    CHECK(s.error <= 1e-12);
    CHECK(*spatium::io::build::to_text(s.expr, "p") == "1 / 6 * p + 1");
    CHECK(*spatium::io::build::to_latex(s.expr, "p") == "\\frac{1}{6} p + 1");
    // and one whose coefficient is not a fraction keeps the number
    const auto r = rsc::synthesize([](double x) { return std::sqrt(2.0) * x + 1; }, spec);
    REQUIRE(r.accepted);
    CHECK(spatium::io::build::to_text(r.expr, "p")->find('/') == std::string::npos);
}

TEST_CASE("synthesis: a size begun is searched whole, and a term that only refits the noise is dropped",
          "[rsc][synthesis]") {
    // 2n / (log n)^2, the Hardy-Littlewood growth of Goldbach pairs, is one term of seven nodes,
    // n / (log(n) * log(n)); a budget cut inside size seven missed it, and once found by a swap the
    // term beside it stayed with a coefficient of 1e-17
    rsc::SynthesisSpec spec{.lo = 10, .hi = 1e6, .max_size = 7, .unary = {rsc::Unary::Sqrt, rsc::Unary::Log}};
    const auto reference = [](double n) { return 2 * n / (std::log(n) * std::log(n)); };
    const auto s = rsc::synthesize(reference, spec);
    REQUIRE(s.accepted);
    CHECK(s.terms == 1);
    CHECK(*spatium::io::build::to_text(s.expr, "n") == "2 * (n / log(n)^2)");

    // every expression of the sizes held is there, whatever the budget once the size is begun
    const auto x = rsc::synthesis_detail::grid(spec.lo, spec.hi, spec.samples);
    const auto count_of = [&](std::size_t budget, std::size_t size) {
        rsc::SynthesisSpec b = spec;
        b.max_candidates = budget;
        return std::ranges::count_if(rsc::enumerate(b, x), [&](const rsc::Term& t) { return t.size == size; });
    };
    CHECK(count_of(4000, 7) == count_of(1000000, 7));
    CHECK(count_of(10, 7) == 0);   // a size not begun is not there at all, not in part

    // the case it must fail: six nodes cannot write it, and no single term is accepted
    spec.max_size = 6;
    const auto six = rsc::synthesize(reference, spec);
    CHECK_FALSE((six.accepted && six.terms == 1));
}

TEST_CASE("synthesis: three short terms are tried together, not grown one at a time", "[rsc][synthesis]") {
    // the Thomas-Fermi solution at the origin, 1 + B x + 4/3 x^(3/2): three terms of six nodes in
    // all, which on a short interval are nearly collinear, so a third term picked greedily after
    // the best pair took whichever leaned toward the residual and the answer was not found
    const double B = -1.5880710226113753;
    const auto origin = [&](double x) { return 1 + B * x + 4.0 / 3 * x * std::sqrt(x); };
    for (double hi : {1e-3, 0.1, 1.0}) {
        INFO("on [0, " << hi << "]");
        rsc::SynthesisSpec spec{.lo = 0, .hi = hi, .tolerance = 1e-12};
        const auto s = rsc::synthesize(origin, spec);
        REQUIRE(s.accepted);
        CHECK(s.terms == 3);
        CHECK(spatium::io::build::to_text(s.expr, "x")->find("x * sqrt(x)") != std::string::npos);
        // the case it must fail: without the sets of three the short intervals are not answered
        if (hi < 1) {
            spec.exhaustive_triples = 0;
            CHECK_FALSE(rsc::synthesize(origin, spec).accepted);
        }
    }
}

TEST_CASE("synthesis: a sum of its own terms comes back, the inverse of evaluation", "[rsc][synthesis]") {
    using rsc::Expr;
    // Synthesis is the inverse of evaluating an expression: a sum of up to three terms of its own
    // language (x, 1, + - * /, sqrt, up to four nodes each) with coefficients no fraction writes,
    // evaluated on an interval, must come back accepted, in no more terms than it has. This is
    // the measure the search is held to: before the sets of three, 46 of 60 sums of three terms
    // came back on [0.5, 3] and 36 of 60 of four; with them all of three and 56 of four.
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> coef(-3, 3);
    const std::function<Expr(int)> grow = [&](int size) -> Expr {
        if (size == 1) return std::uniform_int_distribution<int>(0, 2)(rng) ? Expr::u() : Expr(1.0);
        if (size == 2) return sqrt(grow(1));
        if (std::uniform_int_distribution<int>(0, 4)(rng) == 0) return sqrt(grow(size - 1));
        const int a = std::uniform_int_distribution<int>(1, size - 2)(rng);
        const Expr l = grow(a), r = grow(size - 1 - a);
        switch (std::uniform_int_distribution<int>(0, 3)(rng)) {
            case 0: return l + r;
            case 1: return l - r;
            case 2: return l * r;
            default: return l / r;
        }
    };
    for (auto [lo, hi] : {std::pair{1.0, 1.01}, std::pair{0.5, 3.0}}) {
        for (int k = 1; k <= 3; ++k) {
            int tried = 0;
            while (tried < 12) {
                Expr e(0.0);
                for (int j = 0; j < k; ++j) {
                    const Expr t = grow(std::uniform_int_distribution<int>(1, 4)(rng));
                    e = j == 0 ? Expr(coef(rng) + 1e-3 * std::sqrt(2.0)) * t : e + Expr(coef(rng) + 1e-3 * std::sqrt(3.0)) * t;
                }
                bool usable = true;
                for (int i = 0; i <= 64 && usable; ++i) {
                    const double v = e(lo + (hi - lo) * i / 64, 0.0);
                    usable = std::isfinite(v) && std::abs(v) < 1e6 && std::abs(v) > 1e-6;
                }
                if (!usable) continue;
                ++tried;
                const auto s = rsc::synthesize([&](double x) { return e(x, 0.0); }, {.lo = lo, .hi = hi, .tolerance = 1e-10});
                INFO("[" << lo << ", " << hi << "] " << *spatium::io::build::to_text(e, "x"));
                CHECK(s.accepted);
                CHECK(s.terms <= std::size_t(k));
            }
        }
    }
}

TEST_CASE("synthesis: where no short expression answers, a Chebyshev polynomial does", "[rsc][synthesis]") {
    // the error function on [0, 3]: no closed form in the operations; the fallback interpolates it
    // at Chebyshev points, and the polynomial is a field certified on the same grid
    rsc::SynthesisSpec spec{.lo = 0, .hi = 3, .tolerance = 1e-10, .relative = false, .polynomial_degree = 32};
    const auto erf = rsc::chebyshev([](double x) { return std::erf(x); }, spec);
    REQUIRE(erf.accepted);
    CHECK(erf.form == rsc::Form::Chebyshev);
    CHECK(erf.terms <= 33);
    for (double x : {0.0, 0.3, 1.0, 2.2, 3.0}) CHECK_THAT(erf.expr(x, 0.0), WithinAbs(std::erf(x), 1e-10));
    // a field like any other: lowered to plain data, it gives the interpreter's bits
    const auto pod = spatium::io::build::lower(erf.expr);
    REQUIRE(pod.has_value());
    for (double x : {0.1, 1.7, 2.9})
        CHECK(spatium::io::build::interpret(*pod, x, 0.0) == erf.expr(x, 0.0));

    // pieces: the error function on [0.25, 3] took 46 pieces of short terms; with the fallback, few
    rsc::SynthesisSpec piecewise{.lo = 0.25, .hi = 3, .tolerance = 1e-8, .polynomial_degree = 32};
    const auto pieces = rsc::synthesize_piecewise([](double x) { return std::erf(x); }, piecewise, 1e-6);
    CHECK(pieces.size() <= 2);
    CHECK(std::ranges::all_of(pieces, [](const rsc::Solution& s) { return s.accepted; }));

    // the case it must fail: a kink, |x - 1| on [0, 2], is no polynomial of degree 32 to 1e-10
    rsc::SynthesisSpec kink{.lo = 0, .hi = 2, .tolerance = 1e-10, .relative = false, .polynomial_degree = 32};
    CHECK_FALSE(rsc::chebyshev([](double x) { return std::abs(x - 1); }, kink).accepted);
}
