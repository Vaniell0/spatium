// A law from noisy data, chosen by its description length, its coefficients read as named constants
// only where they could not have landed there by chance. Held on the forms of two laws of the zeta
// function with noise of a stated size: the count of zeros, (T / 2pi) log(T / 2pi) - T / 2pi + 7/8,
// and Hardy and Littlewood's mean square, T log T + (2 gamma - 1 - log 2pi) T. The case it must fail:
// the same data with a coefficient that is no such constant, which keeps its digits.

#include <catch2/catch_test_macros.hpp>
#include <spatium/answer/law.hpp>
#include <spatium/io/field_text.hpp>

#include <cmath>
#include <cstdint>
#include <numbers>
#include <random>
#include <string>
#include <vector>

namespace {

namespace sa = spatium::answer;

// noise of unit variance at each point, the same on every platform: a 64-bit linear congruential
// generator (Knuth's MMIX constants) through Box and Muller
struct Noise {
    std::uint64_t state = 0x9E3779B97F4A7C15ull;
    double uniform() {
        state = state * 6364136223846793005ull + 1442695040888963407ull;
        return (static_cast<double>(state >> 11) + 0.5) / 9007199254740992.0;
    }
    double operator()() {
        const double r = std::sqrt(-2 * std::log(uniform())), a = 2 * std::numbers::pi * uniform();
        return r * std::cos(a);
    }
};

struct Data {
    std::vector<double> x, y, sigma;
};

template <class F, class S>
Data sampled(F law, S noise) {
    Data d;
    Noise unit;
    for (double t = 1e3; t <= 5e6; t *= 1.15) {
        d.x.push_back(t);
        d.sigma.push_back(noise(t));
        d.y.push_back(law(t) + noise(t) * unit());
    }
    return d;
}

std::vector<std::string> terms_of(const sa::Law& l) {
    std::vector<std::string> out;
    for (const auto& t : l.terms) out.push_back(*spatium::io::build::to_text(t, "x"));
    std::ranges::sort(out);
    return out;
}

const double euler = std::numbers::egamma, log2pi = std::log(2 * std::numbers::pi), pi = std::numbers::pi;

} // namespace

TEST_CASE("law: the count of zeros, its terms found and two coefficients named", "[answer][law]") {
    const auto d = sampled([](double t) { return t * std::log(t) / (2 * pi) - (1 + log2pi) / (2 * pi) * t + 7.0 / 8; },
                           [](double) { return 0.05; });   // at 0.3, log(x) and 1 are told apart in some draws and not in others
    const sa::Law l = sa::synthesize_law(d.x, d.y, d.sigma);
    INFO(l.text);
    CHECK(terms_of(l) == std::vector<std::string>{"1", "x", "x * log(x)"});
    std::size_t named = 0;
    for (std::size_t j = 0; j < l.terms.size(); ++j) {
        const std::string t = *spatium::io::build::to_text(l.terms[j], "x");
        if (t == "x * log(x)") {
            REQUIRE(l.named[j]);
            CHECK(l.named[j]->text == "1/(2 pi)");
            ++named;
        }
        if (t == "x") {
            REQUIRE(l.named[j]);
            CHECK(l.named[j]->text == "-1/(2 pi) - log(2 pi)/(2 pi)");
            ++named;
        }
        if (l.named[j]) CHECK(l.named[j]->by_chance < 0.01);
    }
    CHECK(named == 2);
}

TEST_CASE("law: Hardy and Littlewood's mean square, gamma and log 2pi named", "[answer][law]") {
    const auto d = sampled([](double t) { return t * std::log(t) + (2 * euler - 1 - log2pi) * t; },
                           [](double t) { return std::pow(t, 0.25); });
    const sa::Law l = sa::synthesize_law(d.x, d.y, d.sigma);
    INFO(l.text);
    CHECK(terms_of(l) == std::vector<std::string>{"x", "x * log(x)"});
    for (std::size_t j = 0; j < l.terms.size(); ++j) {
        REQUIRE(l.named[j]);
        const std::string t = *spatium::io::build::to_text(l.terms[j], "x");
        CHECK(l.named[j]->text == (t == "x" ? "-1 + 2 gamma - log(2 pi)" : "1"));
    }
}

TEST_CASE("law: a coefficient that is no named constant keeps its digits", "[answer][law]") {
    // the case it must fail: -1.7 lies near many combinations, and none is cheap enough to be trusted
    const auto d = sampled([](double t) { return t * std::log(t) - 1.7 * t; }, [](double t) { return std::pow(t, 0.25); });
    const sa::Law l = sa::synthesize_law(d.x, d.y, d.sigma);
    INFO(l.text);
    for (std::size_t j = 0; j < l.terms.size(); ++j)
        if (*spatium::io::build::to_text(l.terms[j], "x") == "x") {
            CHECK_FALSE(l.named[j]);
            CHECK(std::abs(l.coefficients[j] + 1.7) < 3 * l.errors[j] + 1e-9);
        }
}

TEST_CASE("law: no data, data of unequal lengths or a point with no noise gives no law", "[answer][law]") {
    const std::vector<double> none, six{1, 2, 3, 4, 5, 6};
    const std::vector<sa::Law> laws{
        sa::synthesize_law(none, none, none),
        sa::synthesize_law(six, std::vector<double>(5, 1.0), std::vector<double>(6, 1.0)),   // one value short
        sa::synthesize_law(six, six, std::vector<double>(3, 1.0)),                           // too few noises
        sa::synthesize_law(six, six, std::vector<double>{1, 1, 0, 1, 1, 1}),                 // a point of noise 0 weighs infinitely
        sa::synthesize_law(six, six, std::vector<double>{1, 1, -1, 1, 1, 1}),
    };
    for (const sa::Law& l : laws) {
        CHECK(l.terms.empty());
        CHECK(std::isinf(l.length));
        CHECK(l.text.empty());
    }
}

TEST_CASE("law: powers the caller names are terms the four rules cannot build", "[answer][law]") {
    // two terms of a phase law, 3 x^(-5/3) + 0.5 x^(-2/3): no sum of 1, x, +, -, *, /, log and sqrt in five nodes is a
    // power with a third in its exponent, so without the powers the law is not found
    const auto truth = [](double t) { return 3 * std::pow(t, -5.0 / 3) + 0.5 * std::pow(t, -2.0 / 3); };
    Data d;
    Noise unit;
    for (double t = 20; t <= 400; t *= 1.04) {
        d.x.push_back(t);
        d.sigma.push_back(1e-5);
        d.y.push_back(truth(t) + 1e-5 * unit());
    }
    sa::LawSpec with;
    with.max_terms = 2;
    with.powers = {-5.0 / 3, -4.0 / 3, -1.0, -2.0 / 3, -1.0 / 3};
    const sa::Law found = sa::synthesize_law(d.x, d.y, d.sigma, with);
    INFO(found.text);
    REQUIRE(found.terms.size() == 2);
    CHECK(found.rms < 1.5);   // the residual in units of the noise
    for (double t : {25.0, 100.0, 380.0}) CHECK(std::abs(found.expr(t, 0.0) - truth(t)) < 1e-3 * truth(t));
    for (std::size_t j = 0; j < found.terms.size(); ++j) CHECK(std::abs(found.coefficients[j]) > 5 * found.errors[j]);

    sa::LawSpec without;
    without.max_terms = 2;
    const sa::Law plain = sa::synthesize_law(d.x, d.y, d.sigma, without);
    INFO(plain.text);
    CHECK(plain.rms > 5);                  // the case it must fail: no two terms of the four rules reach the noise
    CHECK(plain.length > found.length);
}

TEST_CASE("law: more than three terms are searched when the caller allows them", "[answer][law]") {
    // four terms, each a named power or x: a search that stops at triples leaves a residual far above the noise
    const auto truth = [](double t) { return 3 * std::pow(t, -5.0 / 3) + 2 * std::pow(t, -1.0) - 1.5 * std::pow(t, -1.0 / 3) + 0.01 * t; };
    Data d;
    Noise unit;
    for (double t = 5; t <= 400; t *= 1.04) {
        d.x.push_back(t);
        d.sigma.push_back(1e-6);
        d.y.push_back(truth(t) + 1e-6 * unit());
    }
    sa::LawSpec spec;
    spec.max_terms = 4;
    spec.powers = {-5.0 / 3, -1.0, -1.0 / 3};
    const sa::Law found = sa::synthesize_law(d.x, d.y, d.sigma, spec);
    INFO(found.text);
    REQUIRE(found.terms.size() == 4);
    CHECK(found.rms < 1.5);
    for (double t : {6.0, 50.0, 390.0}) CHECK(std::abs(found.expr(t, 0.0) - truth(t)) < 1e-4 * std::abs(truth(t)));

    spec.max_terms = 3;
    const sa::Law three = sa::synthesize_law(d.x, d.y, d.sigma, spec);
    CHECK(three.terms.size() <= 3);
    CHECK(three.rms > 5);
}

TEST_CASE("law: correlated noise given as a covariance is whitened, and a bad covariance gives no law", "[answer][law]") {
    // 2 + 0.5 log x over a decade with noise that is smooth from point to point (AR(1), rho 0.9): read as independent,
    // the residual looks too small and a flexible form wins; with the covariance the form is the true one
    const int n = 24;
    const double rho = 0.9, sg = 0.03;
    std::vector<double> cov(n * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) cov[i * n + j] = sg * sg * std::pow(rho, std::abs(i - j));
    int plain = 0, whitened = 0;
    double rms = 0;
    const int trials = 40;
    for (int s = 0; s < trials; ++s) {
        std::mt19937 g(100 + s);
        std::normal_distribution<double> unit(0, 1);
        std::vector<double> x, y;
        double e = sg * unit(g);
        for (int i = 0; i < n; ++i) {
            x.push_back(10 * std::pow(10.0, i / double(n - 1)));
            if (i) e = rho * e + sg * std::sqrt(1 - rho * rho) * unit(g);
            y.push_back(2.0 + 0.5 * std::log(x.back()) + e);
        }
        const auto is_truth = [](const sa::Law& l) { return l.terms.size() == 2 && l.text.find("log") != std::string::npos; };
        sa::LawSpec spec;
        spec.max_terms = 3;
        plain += is_truth(sa::synthesize_law(x, y, std::vector<double>(n, sg), spec));
        spec.covariance = cov;
        const sa::Law found = sa::synthesize_law(x, y, {}, spec);
        whitened += is_truth(found);
        rms += found.rms;
    }
    INFO("independent reading " << plain << " of " << trials << ", covariance " << whitened << " of " << trials);
    CHECK(whitened >= 30);
    CHECK(whitened > plain + 10);
    CHECK(std::abs(rms / trials - 1.0) < 0.25);   // the whitened residual is the noise, in units of its own sigma

    // a diagonal covariance is the same as the sigma
    std::vector<double> x{1, 2, 3, 4, 5, 6, 7, 8}, y, sig(8, 0.1), diag(64, 0.0);
    for (int i = 0; i < 8; ++i) { y.push_back(1.0 + 2.0 * x[i] + 0.1 * ((i % 3) - 1)); diag[i * 8 + i] = 0.01; }
    sa::LawSpec a;
    sa::LawSpec b;
    b.covariance = diag;
    const sa::Law la = sa::synthesize_law(x, y, sig, a), lb = sa::synthesize_law(x, y, {}, b);
    CHECK(la.text == lb.text);
    CHECK(std::abs(la.length - lb.length) < 1e-9);

    // not n x n, not symmetric, not positive definite
    auto wrong = diag; wrong.pop_back();
    auto asym = diag; asym[1] = 0.005;
    auto notpd = diag; notpd[1] = notpd[8] = 0.02;   // |rho| > 1
    for (const auto* c : {&wrong, &asym, &notpd}) {
        sa::LawSpec bad;
        bad.covariance = *c;
        const sa::Law l = sa::synthesize_law(x, y, sig, bad);
        CHECK(l.terms.empty());
    }
}
