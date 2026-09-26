// ccd_synthesize -- find, for each class of continuous-collision query, the
// cheapest chain of primitives (rsc/include/ccd_ops.hpp) that answers every
// query of the class right.
//
// Right is the one-sided rule plus usefulness: never later than the truth,
// never a miss where there is contact, never a contact where there is
// none, and at most 1e-3 of the step early. Without the last two "answer
// contact at once" would pass and cost nothing. The truth is exact: every
// obstacle here has one -- a sphere's quadric, a torus's quartic -- also
// when it is presented to the chain as a chart, so the chart's primitives
// are judged against an answer that does not come from themselves.
//
// A class is an obstacle and a way of moving past it. The result is a
// table: the chain each class wants, what it costs, and what the general
// chains -- the hybrid rigid_contact.hpp's first_contact runs, and the
// search alone -- cost on the same queries.

#include <ccd_ops.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <format>
#include <map>
#include <mutex>
#include <print>
#include <random>
#include <string>
#include <thread>
#include <vector>

using namespace rsc::ccd;
namespace geo = spatium::geometry;

namespace {

enum class Motion { HeadOn, Graze, Miss, Through };
const char* name(Motion m) {
    return m == Motion::HeadOn ? "head-on" : m == Motion::Graze ? "graze" : m == Motion::Miss ? "miss" : "through";
}

struct Case {
    V3 p0, disp;
    double truth;   // fraction of the step at first contact, > 1 for none
};

std::shared_ptr<const spatium::ParametricSurface<double>> sphere_chart(double r) {
    return std::make_shared<spatium::ParametricSurface<double>>(
        [r](double u, double v) { return V3{r * std::sin(v) * std::cos(u), r * std::sin(v) * std::sin(u), r * std::cos(v)}; },
        spatium::ParametricSurface<double>::Domain{0.0, 2 * std::numbers::pi, 0.0, std::numbers::pi}, true, false);
}
std::shared_ptr<const spatium::ParametricSurface<double>> torus_chart(double R, double r) {
    return std::make_shared<spatium::ParametricSurface<double>>(
        [R, r](double u, double v) {
            return V3{(R + r * std::cos(v)) * std::cos(u), (R + r * std::cos(v)) * std::sin(u), r * std::sin(v)};
        },
        spatium::ParametricSurface<double>::Domain{0.0, 2 * std::numbers::pi, 0.0, 2 * std::numbers::pi}, true, true);
}

double truth_of(const Obstacle& ob, const V3& p0, const V3& disp) {
    const double len = disp.norm();
    const geo::Ray<3, double> ray{p0, V3{disp * (1.0 / len)}};
    double first = 2.0;
    const bool torus = ob.torus.minor_radius > 0 && (ob.shape == Shape::Torus || ob.radius < 0);
    if (!torus) {
        for (const auto& h : geo::ray_quadric(ray, geo::Quadric<double>::sphere(std::abs(ob.radius))))
            if (h.t >= 0 && h.t <= len) first = std::min(first, h.t / len);
    } else {
        for (const auto& h : geo::ray_torus(ray, ob.torus))
            if (h.t >= 0 && h.t <= len) first = std::min(first, h.t / len);
    }
    return first;
}

// A surface point and its outward normal, for aiming a query at it.
std::pair<V3, V3> surface_point(const Obstacle& ob, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> ang(0.0, 2 * std::numbers::pi);
    const bool torus = ob.radius < 0 || ob.shape == Shape::Torus;
    if (!torus) {
        const double u = ang(rng), v = std::acos(std::uniform_real_distribution<double>(-1, 1)(rng));
        const V3 n{std::sin(v) * std::cos(u), std::sin(v) * std::sin(u), std::cos(v)};
        return {V3{n * std::abs(ob.radius)}, n};
    }
    const double u = ang(rng), v = ang(rng), R = ob.torus.major_radius, r = ob.torus.minor_radius;
    const V3 n{std::cos(v) * std::cos(u), std::cos(v) * std::sin(u), std::sin(v)};
    return {V3{V3{std::cos(u), std::sin(u), 0.0} * R + n * r}, n};
}

std::vector<Case> sample(const Obstacle& ob, Motion m, int count, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    std::vector<Case> out;
    while (static_cast<int>(out.size()) < count) {
        const auto [foot, n] = surface_point(ob, rng);
        V3 tan{uni(rng), uni(rng), uni(rng)};
        tan = V3{tan - n * tan.dot(n)};
        tan = V3{tan * (1.0 / tan.norm())};
        V3 p0, disp;
        if (m == Motion::HeadOn) {           // straight at the surface from a few units out
            p0 = V3{foot + n * (1.0 + 2.0 * (uni(rng) + 1.0))};
            disp = V3{V3{foot - p0} * 1.5};
        } else if (m == Motion::Graze) {     // tangent, a hair inside
            const double depth = std::pow(10.0, -2.0 - 8.0 * (uni(rng) + 1.0) / 2.0);
            const V3 through{foot - n * depth};
            p0 = V3{through - tan * 2.0};
            disp = V3{tan * 4.0};
        } else if (m == Motion::Miss) {      // passes by at a distance
            const V3 by{foot + n * (0.05 + 0.5 * (uni(rng) + 1.0))};
            p0 = V3{by - tan * 2.0};
            disp = V3{tan * 4.0};
        } else {                             // fast, a long step clean through the middle
            p0 = V3{foot + n * 6.0};
            disp = V3{n * -12.0};
        }
        double d;
        if (ob.shape == Shape::Torus || ob.radius < 0) d = spatium::physics::mechanics::point_to(p0, ob.torus).distance;
        else d = std::abs(p0.norm() - std::abs(ob.radius));
        if (d < 1e-3) continue;              // not starting on the surface
        out.push_back({p0, disp, truth_of(ob, p0, disp)});
    }
    return out;
}

bool admissible(const Answer& a, double truth) {
    if (truth > 1.0) return !a.hit;
    return a.hit && a.toi <= truth + 1e-9 && a.toi >= truth - 1e-3;
}

std::vector<Chain> chains() {
    const std::vector<Step> prims{{Op::CF},       {Op::A, 4},       {Op::A, 32},   {Op::S, 1 << 12},
                                  {Op::S, 1 << 16}, {Op::S, 1 << 20}, {Op::P, 0.10}, {Op::P, 0.30},
                                  {Op::P, 0.60}};
    std::vector<Chain> out;
    for (const auto& a : prims) {
        out.push_back({a});
        for (const auto& b : prims) {
            out.push_back({a, b});
            for (const auto& c : prims) out.push_back({a, b, c});
        }
    }
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    const int per_class = argc > 1 ? std::atoi(argv[1]) : 40;
    // Four obstacles, each with a truth: radius < 0 marks "a torus, shown
    // to the chain as a chart", so truth_of still knows what it is.
    Obstacle sphere{Shape::Sphere, 1.0, {}, sphere_chart(1.0)};
    Obstacle torus{Shape::Torus, 0.0, spatium::geometry::Torus<double>{V3{}, V3{0, 0, 1}, 1.0, 0.3}, torus_chart(1.0, 0.3)};
    Obstacle sphere_as_chart{Shape::Chart, 1.0, {}, sphere_chart(1.0)};
    sphere_as_chart.lipschitz = 1.0;
    Obstacle torus_as_chart{Shape::Chart, -1.0, spatium::geometry::Torus<double>{V3{}, V3{0, 0, 1}, 1.0, 0.3}, torus_chart(1.0, 0.3)};
    torus_as_chart.lipschitz = 1.3;
    const std::vector<std::pair<const char*, const Obstacle*>> obstacles{
        {"sphere", &sphere}, {"torus", &torus}, {"sphere as chart", &sphere_as_chart}, {"torus as chart", &torus_as_chart}};
    const auto all = chains();
    const Chain hybrid{{Op::A, 8}, {Op::S, 1 << 16}}, search{{Op::S, 1 << 16}};

    struct Row { std::string obstacle, motion, best; double best_cost, hybrid_cost, search_cost; int hybrid_ok, search_ok, admissible_chains; };
    std::vector<Row> rows;
    std::mutex mu;
    std::vector<std::thread> pool;
    for (std::size_t oi = 0; oi < obstacles.size(); ++oi)
        for (Motion m : {Motion::HeadOn, Motion::Graze, Motion::Miss, Motion::Through})
            pool.emplace_back([&, oi, m] {
                const auto& [oname, ob] = obstacles[oi];
                const auto cases = sample(*ob, m, per_class, 1000 * oi + static_cast<int>(m));
                auto eval = [&](const Chain& c, int& ok) {
                    double cost = 0;
                    ok = 0;
                    for (const auto& k : cases) {
                        const auto a = run(Query{k.p0, k.disp, 0.0, ob}, c);
                        cost += static_cast<double>(a.cost);
                        ok += admissible(a, k.truth);
                    }
                    return cost / cases.size();
                };
                Row row{oname, name(m), "(none)", 0, 0, 0, 0, 0, 0};
                double best = std::numeric_limits<double>::infinity();
                std::size_t best_len = 99;
                for (const auto& c : all) {
                    int ok;
                    const double cost = eval(c, ok);
                    if (ok != static_cast<int>(cases.size())) continue;
                    ++row.admissible_chains;
                    // Ties to the shorter chain: a CF on a chart does nothing
                    // and costs nothing, and should not be kept for it.
                    if (cost < best || (cost == best && c.size() < best_len)) {
                        best = cost;
                        best_len = c.size();
                        row.best = rsc::ccd::name(c);
                        row.best_cost = cost;
                    }
                }
                row.hybrid_cost = eval(hybrid, row.hybrid_ok);
                row.search_cost = eval(search, row.search_ok);
                std::lock_guard lk(mu);
                rows.push_back(row);
            });
    for (auto& t : pool) t.join();
    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
        return std::tie(a.obstacle, a.motion) < std::tie(b.obstacle, b.motion);
    });
    std::println("{} queries a class, {} chains; admissible = right on every query of the class", per_class, all.size());
    std::println("{:<16} {:<8} | {:<14} {:>10} | {:>18} | {:>18} | {:>6}", "obstacle", "motion", "best chain", "cost",
                 "hybrid A8 S16", "search S16", "admis.");
    for (const auto& r : rows)
        std::println("{:<16} {:<8} | {:<14} {:>10.1f} | {:>9.1f} ({:>2}/{:<2}) | {:>9.1f} ({:>2}/{:<2}) | {:>6}", r.obstacle,
                     r.motion, r.best, r.best_cost, r.hybrid_cost, r.hybrid_ok, per_class, r.search_cost, r.search_ok,
                     per_class, r.admissible_chains);
    return 0;
}
