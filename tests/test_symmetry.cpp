// The registry of pairs of paths (tests/symmetry/pairs.hpp): every answer the library
// computes two ways, held to the tolerance its paths can promise, on inputs drawn by
// measure.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "symmetry/pairs.hpp"

using Catch::Matchers::WithinAbs;

TEST_CASE("every registered pair of paths agrees within what its paths can promise", "[symmetry]") {
    for (const auto& entry : symmetry::all_pairs()) {
        const auto r = entry.run();
        INFO(r.name << ": " << r.samples << " inputs, " << r.failures << " failures, worst "
                    << r.worst_ratio << " of the bound, first failure at input " << r.first_failure);
        CHECK(r.failures == 0);
        CHECK(r.samples >= entry.min_samples);
        CHECK(r.unmatched_measure() <= 3.0 / static_cast<double>(entry.min_samples) + 1e-12);
    }
}

TEST_CASE("a pair that disagrees is reported, with where and how much", "[symmetry]") {
    // Every verifier gets a case it must fail. Here the route is the reference plus a
    // fault that is above the bound on one input in ten.
    const auto r = symmetry::check<double>(
        "a route with a fault", symmetry::Kind::Function, 1000,
        [](std::size_t i) { return symmetry::halton(i, 0); },
        [](double x) { return x * x; },
        [](double x) { return x * x + (x > 0.9 ? 1e-3 : 0.0); },
        [](double a, double b) { return std::abs(a - b); },
        [](double) { return 1e-6; });
    CHECK(r.failures > 0);
    CHECK(r.first_failure > 0);
    CHECK(r.worst_ratio > 100.0);
    CHECK_THAT(r.unmatched_measure(), WithinAbs(0.1, 0.02));

    // A distance that is not finite is a failure, never a pass.
    const auto n = symmetry::check<double>(
        "a route that is NaN", symmetry::Kind::Relation, 10,
        [](std::size_t i) { return symmetry::halton(i, 0); },
        [](double x) { return x; },
        [](double) { return std::nan(""); },
        [](double a, double b) { return std::abs(a - b); },
        [](double) { return 1.0; });
    CHECK(n.failures == 10);
}

TEST_CASE("a pair with no failure bounds the failing measure by the rule of three", "[symmetry]") {
    const auto r = symmetry::check<double>(
        "an exact pair", symmetry::Kind::Function, 3000,
        [](std::size_t i) { return symmetry::halton(i, 0); },
        [](double x) { return x; }, [](double x) { return x; },
        [](double a, double b) { return std::abs(a - b); }, [](double) { return 1e-12; });
    CHECK(r.failures == 0);
    CHECK_THAT(r.unmatched_measure(), WithinAbs(1e-3, 1e-12));
}

TEST_CASE("ray_torus on the ray that used to be 4e-4 off", "[symmetry]") {
    // Input 5037 of the registry's sequence: both hits were off by 4.3e-4 against fifty
    // digits, with the other roots a gap of 3.8 away, before the roots were polished.
    const auto s = symmetry::torus_ray(5037);
    const auto d = symmetry::torus_hits_double(s);
    const auto r = symmetry::torus_hits_real50(s);
    REQUIRE(d.size() == 2);
    REQUIRE(r.size() == 2);
    CHECK_THAT(d[0], WithinAbs(r[0], 1e-9));
    CHECK_THAT(d[1], WithinAbs(r[1], 1e-9));
}
