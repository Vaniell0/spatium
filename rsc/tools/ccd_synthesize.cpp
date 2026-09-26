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
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdio>
#include <format>
#include <map>
#include <mutex>
#include <print>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <vector>

using namespace rsc::ccd;
namespace geo = spatium::geometry;
namespace mech = spatium::physics::mechanics;

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

    // A class's queries run as one step of a solver would: every query of
    // the class against one cell tree, fresh for each chain, so a chain
    // pays for the cells it makes and for walking the ones already made,
    // and no chain inherits another's.
    struct Class { std::string obstacle, motion; const Obstacle* ob; std::vector<Case> cases; };
    std::vector<Class> classes;
    for (std::size_t oi = 0; oi < obstacles.size(); ++oi)
        for (Motion m : {Motion::HeadOn, Motion::Graze, Motion::Miss, Motion::Through})
            classes.push_back({obstacles[oi].first, name(m), obstacles[oi].second,
                               sample(*obstacles[oi].second, m, per_class, 1000 * oi + static_cast<int>(m))});
    // What a chain did on a class, as three fingerprints: its answers bit
    // for bit, its answers and what each cost, and its answers with the
    // time rounded to 1e-9 of the step. Two chains with one fingerprint
    // are one algorithm as far as the class can tell -- the measure of how
    // much of the chain space is really there, as mesh_collisions.cpp
    // measures it for meshes. Bits flatter the space, rounding flatters
    // its collapse; the gap between them is part of the answer.
    struct Print { std::uint64_t exact = 1469598103934665603u, costed = 1469598103934665603u, rounded = 1469598103934665603u; };
    auto mix = [](std::uint64_t& h, std::uint64_t v) { h = (h ^ v) * 1099511628211u; };
    auto eval = [&](const Class& k, const Chain& c, int& ok, Print* print = nullptr) {
        mech::ChartCellTree<double> tree(k.ob->bound());
        double cost = 0;
        ok = 0;
        for (const auto& q : k.cases) {
            const auto a = run(Query{q.p0, q.disp, 0.0, k.ob, &tree}, c);
            cost += static_cast<double>(a.cost);
            ok += admissible(a, q.truth);
            if (print) {
                const std::uint64_t bits = std::bit_cast<std::uint64_t>(a.toi) ^ (a.hit ? 1u : 0u);
                mix(print->exact, bits);
                mix(print->costed, bits);
                mix(print->costed, a.cost);
                mix(print->rounded, static_cast<std::uint64_t>(std::llround(a.toi * 1e9)) * 2 + (a.hit ? 1u : 0u));
            }
        }
        return cost / static_cast<double>(k.cases.size());
    };

    struct Result { double cost; int ok; Print print; };
    std::vector<std::vector<Result>> results(classes.size(), std::vector<Result>(all.size()));
    std::vector<std::atomic<std::size_t>> left(classes.size());
    for (auto& l : left) l = all.size();
    std::atomic<std::size_t> next{0};
    std::mutex mu;
    std::println("{} queries a class, {} chains; admissible = right on every query of the class", per_class, all.size());
    std::println("{:<16} {:<8} | {:<14} {:>10} | {:>18} | {:>18} | {:>6}", "obstacle", "motion", "best chain", "cost",
                 "hybrid A8 S16", "search S16", "admis.");
    std::println("{:>104} distinct: exact +cost 1e-9", "");
    std::fflush(stdout);
    auto report = [&](std::size_t ci) {
        const auto& k = classes[ci];
        std::string best = "(none)";
        double best_cost = 0, lowest = std::numeric_limits<double>::infinity();
        std::size_t best_len = 99;
        int admissible_chains = 0;
        std::set<std::uint64_t> exact, costed, rounded;
        for (std::size_t j = 0; j < all.size(); ++j) {
            const auto& r = results[ci][j];
            exact.insert(r.print.exact);
            costed.insert(r.print.costed);
            rounded.insert(r.print.rounded);
            if (r.ok != static_cast<int>(k.cases.size())) continue;
            ++admissible_chains;
            // Ties to the shorter chain: a CF on a chart does nothing and
            // costs nothing, and should not be kept for it.
            if (r.cost < lowest || (r.cost == lowest && all[j].size() < best_len)) {
                lowest = r.cost;
                best_len = all[j].size();
                best = rsc::ccd::name(all[j]);
                best_cost = r.cost;
            }
        }
        int hok, sok;
        const double hc = eval(k, hybrid, hok), sc = eval(k, search, sok);
        std::lock_guard lk(mu);
        std::println("{:<16} {:<8} | {:<14} {:>10.1f} | {:>9.1f} ({:>2}/{:<2}) | {:>9.1f} ({:>2}/{:<2}) | {:>6} | {:>5} {:>5} {:>5}",
                     k.obstacle, k.motion, best, best_cost, hc, hok, per_class, sc, sok, per_class, admissible_chains,
                     exact.size(), costed.size(), rounded.size());
        std::fflush(stdout);
    };
    std::vector<std::thread> pool;
    const unsigned workers = std::max(1u, std::thread::hardware_concurrency());
    for (unsigned w = 0; w < workers; ++w)
        pool.emplace_back([&] {
            for (std::size_t i; (i = next++) < classes.size() * all.size();) {
                const std::size_t ci = i / all.size(), j = i % all.size();
                int ok;
                Print print;
                const double cost = eval(classes[ci], all[j], ok, &print);
                results[ci][j] = {cost, ok, print};
                if (--left[ci] == 0) report(ci);
            }
        });
    for (auto& t : pool) t.join();
    return 0;
}
