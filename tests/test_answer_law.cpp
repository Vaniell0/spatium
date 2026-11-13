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
