// RSC labelling a continuous-collision choice: the one-sided rule has to
// reach the labels, or a model is trained to prefer answers that are late.
#include <catch2/catch_test_macros.hpp>

#include <comparison_task.hpp>

#include <cmath>
#include <random>

namespace {

// Three candidates for a time of impact, cheapest first: one that answers a
// hair late, one that answers far too early to be of use, one that answers
// a hair early. Only the last is admissible.
struct Toi { double toi; };

rsc::ComparisonTaskGenerator<double, Toi> make_generator(bool one_sided) {
    std::vector<rsc::Candidate<double, Toi>> c{
        {"late by 1e-6", [](const double& t) { return Toi{t + 1e-6}; }},
        {"early by 0.05", [](const double& t) { return Toi{t - 0.05}; }},
        {"early by 1e-5", [](const double& t) { return Toi{t - 1e-5}; }},
    };
    rsc::ComparisonTaskGenerator<double, Toi> g(
        std::move(c), [](std::mt19937_64& rng) { return std::uniform_real_distribution<double>(0.2, 0.8)(rng); },
        [](const double& t) { return std::vector<double>{t}; }, [](const double& t) { return Toi{t}; },
        [](const Toi& a, const Toi& b) { return std::abs(a.toi - b.toi); }, 1e-4, 3);
    if (one_sided)
        g.with_gate([](const Toi& cand, const Toi& ref) {
            // Never late; early, but by no more than it takes to stay useful.
            return cand.toi <= ref.toi && cand.toi >= ref.toi - 1e-3;
        });
    return g;
}

}  // namespace

// The symmetric gate took the cheapest candidate within 1e-4 of the
// reference either way, and the cheapest is late: every label it made
// taught a model to prefer the one answer continuous collision must never
// give.
TEST_CASE("A continuous-collision label is never a late answer", "[rsc][ccd]") {
    auto g = make_generator(/*one_sided=*/true);
    for (int i = 0; i < 50; ++i) {
        const auto task = g.sample();
        INFO("sample " << i << " labelled " << g.candidates()[task.op_index].name);
        CHECK(task.op_index == 2);
    }
}

// And without a gate the template behaves exactly as the seven domains
// built on it rely on.
TEST_CASE("Without a gate the comparison is the symmetric one it always was", "[rsc][ccd]") {
    auto g = make_generator(/*one_sided=*/false);
    CHECK(g.sample().op_index == 0);
}

// ── The primitives a chain is made of ───────────────────────────
#include <ccd_ops.hpp>

namespace {

using rsc::ccd::V3;

std::shared_ptr<const spatium::ParametricSurface<double>> sphere_chart(double r) {
    return std::make_shared<spatium::ParametricSurface<double>>(
        [r](double u, double v) { return V3{r * std::sin(v) * std::cos(u), r * std::sin(v) * std::sin(u), r * std::cos(v)}; },
        spatium::ParametricSurface<double>::Domain{0.0, 2 * std::numbers::pi, 0.0, std::numbers::pi}, true, false);
}

// The true first contact of a point with a sphere of radius r, if any.
double truth(const V3& p0, const V3& disp, double r) {
    const double len = disp.norm();
    const spatium::geometry::Ray<3, double> ray{p0, V3{disp * (1.0 / len)}};
    double first = 2.0;
    for (const auto& h : spatium::geometry::ray_quadric(ray, spatium::geometry::Quadric<double>::sphere(r)))
        if (h.t >= 0 && h.t <= len) first = std::min(first, h.t / len);
    return first;
}

}  // namespace

// Every chain, on queries whose truth is known: never later than the
// truth, never a miss where there is contact -- whatever the primitives,
// in whatever order, since each keeps the rule on its own.
TEST_CASE("A chain of CCD primitives is never late", "[rsc][ccd]") {
    using namespace rsc::ccd;
    const double r = 1.0;
    Obstacle exact{Shape::Sphere, r, {}, sphere_chart(r)};
    Obstacle chart{Shape::Chart, r, {}, sphere_chart(r)};
    chart.lipschitz = r;
    const std::vector<Chain> chains{
        {{Op::CF}}, {{Op::A, 16}}, {{Op::A, 8}, {Op::S, 1 << 14}}, {{Op::P, 0.25}, {Op::P, 0.25}, {Op::A, 8}, {Op::S, 1 << 14}},
        {{Op::S, 1 << 16}}, {{Op::P, 0.5}}, {}};
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    int contacts = 0;
    for (int i = 0; i < 300; ++i) {
        V3 a{uni(rng), uni(rng), uni(rng)}, b{uni(rng), uni(rng), uni(rng)};
        a = V3{a * 3.0};
        b = V3{b * 3.0};
        if (a.norm() < 1.05) continue;
        const V3 disp{b - a};
        const double t_true = truth(a, disp, r);
        const bool hit = t_true <= 1.0;
        contacts += hit;
        for (const Obstacle* ob : {&exact, &chart})
            for (const auto& c : chains) {
                if (c.size() == 1 && c[0].op == Op::CF && ob->shape == Shape::Chart) continue;
                const auto ans = run(Query{a, disp, 0.0, ob}, c);
                INFO("query " << i << " chain '" << name(c) << "' on " << (ob == &exact ? "sphere" : "chart")
                              << ": truth " << t_true << ", answer hit=" << ans.hit << " toi=" << ans.toi);
                if (hit) {
                    REQUIRE(ans.hit);
                    REQUIRE(ans.toi <= t_true + 1e-9);
                }
            }
    }
    CHECK(contacts > 30);
}

// The closed form decides and is exact up to the tolerance it answers
// before; a chain that ends undecided says contact where it stopped.
TEST_CASE("CF is exact, and an undecided chain answers where it stopped", "[rsc][ccd]") {
    using namespace rsc::ccd;
    Obstacle ob{Shape::Sphere, 1.0, {}, sphere_chart(1.0)};
    const auto cf = run(Query{V3{0, 0, 3}, V3{0, 0, -6}, 0.0, &ob}, {{Op::CF}});
    REQUIRE(cf.hit);
    CHECK(cf.toi <= 1.0 / 3.0);
    CHECK(cf.toi > 1.0 / 3.0 - 1e-7);
    CHECK(cf.cost == 1);
    const auto nothing = run(Query{V3{0, 0, 3}, V3{0, 0, -6}, 0.0, &ob}, {});
    CHECK(nothing.hit);
    CHECK(nothing.toi == 0.0);
}
