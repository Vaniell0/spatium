// A scatter scene whose count is changed a step at a time, laid out by
// updating the previous scene's instance tree instead of building it.
//
// Each step is what an editor does when a "count" slider moves: the Trace
// is built again (building is nothing -- the work is in the next stages),
// cooked, and laid out. The stage that grows with the number of copies and
// has no memory is the instance tree; `render::lay_out(previous, ...)`
// keeps the copies that did not change and removes or inserts the rest.
//
// Every step is also laid out from scratch, for two reasons: the time of
// that is what the update is compared to, and the two scenes must give the
// same hit on every ray, which is checked and reported.
//
//   scatter_regen [--count N] [--steps +1,+10,-3,...] [--seed S] [--rays R]
//
// A step is a signed change of the count, or `seed` to change the seed
// (every copy then moves, so the tree is built again).

#include <spatium/io/build.hpp>
#include <spatium/render/cooked_scene.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <print>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace spatium;
namespace bd = spatium::io::build;

namespace {

using clk = std::chrono::steady_clock;
template<class F>
double ms(F&& f) {
    const auto t0 = clk::now();
    f();
    return std::chrono::duration<double, std::milli>(clk::now() - t0).count();
}

// The Trace for a count and a seed; ready to cook.
struct Scene {
    bd::Trace<double> trace;
    std::size_t root = 0;
    Scene(std::size_t count, std::uint32_t seed) {
        auto torus = trace.torus(2.0, 1.0, 240, 120);
        auto item = trace.cylinder(0.011, 0.062, 8, 2);
        auto sc = trace.scatter(item, torus, count, seed, 0.45, bd::SeatAxis::X);
        root = trace.compose({torus, sc}).index;
    }
};

std::vector<geometry::Ray<3, double>> rays(std::size_t n) {
    std::mt19937_64 g(7);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<geometry::Ray<3, double>> out;
    for (std::size_t i = 0; i < n; ++i) {
        const double a = u(g) * 6.2831853, r = 5.0 + u(g), b = u(g) * 6.2831853;
        const Vec<double, 3> o{r * std::cos(a), r * std::sin(a), (u(g) - 0.5) * 4.0};
        const Vec<double, 3> target{2.0 * std::cos(b), 2.0 * std::sin(b), (u(g) - 0.5) * 1.5};
        out.push_back(*geometry::Ray<3, double>::from(o, Vec<double, 3>{target - o}));
    }
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    std::size_t count = 20000, n_rays = 4000;
    std::uint32_t seed = 42;
    std::string steps = "+1,+10,+100,-1,-10,-100,seed,+1";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--count" && i + 1 < argc) count = std::strtoull(argv[++i], nullptr, 10);
        else if (a == "--steps" && i + 1 < argc) steps = argv[++i];
        else if (a == "--seed" && i + 1 < argc) seed = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        else if (a == "--rays" && i + 1 < argc) n_rays = std::strtoull(argv[++i], nullptr, 10);
        else {
            std::print("scatter_regen [--count N] [--steps +1,+10,-3,seed,...] [--seed S] [--rays R]\n");
            return a == "--help" ? 0 : 1;
        }
    }

    render::StaticCache<double> cache;   // the torus is the same at every step
    const auto probe = rays(n_rays);

    std::println("{:>8} {:>8} {:>7} {:>9} {:>9} {:>9} {:>9}  {}", "count", "step", "kept", "removed+new",
                 "cook ms", "update ms", "fresh ms", "rays");

    render::CookedScene<double> scene;
    {
        Scene s(count, seed);
        scene = render::lay_out(bd::cook(s.trace, s.root, 0.0), {}, nullptr, &cache);
    }

    std::istringstream in(steps);
    std::string step;
    int failures = 0;
    while (std::getline(in, step, ',')) {
        if (step == "seed") ++seed;
        else count = static_cast<std::size_t>(static_cast<long long>(count) + std::atoll(step.c_str()));

        Scene s(count, seed);   // the control moved: build the Trace again
        bd::Cooked<double> cooked;
        const double cook_ms = ms([&] { cooked = bd::cook(s.trace, s.root, 0.0); });

        render::InstanceUpdate u;
        render::CookedScene<double> next;
        const double update_ms = ms([&] {
            next = render::lay_out(std::move(scene), cooked, {}, nullptr, &cache, &u);
        });
        render::CookedScene<double> fresh;
        const double fresh_ms = ms([&] { fresh = render::lay_out(cooked, {}, nullptr, &cache); });

        // Same hit, or same miss, on every ray; the same copy by position.
        std::size_t hits = 0, bad = 0;
        for (const auto& r : probe) {
            const auto a = next.instances.ray_cast(r), b = fresh.instances.ray_cast(r);
            if (a.has_value() != b.has_value()) { ++bad; continue; }
            if (!a) continue;
            ++hits;
            const auto &ia = next.instances.shapes()[a->index], &ib = fresh.instances.shapes()[b->index];
            if (a->t != b->t || std::memcmp(&ia.translation, &ib.translation, sizeof(ia.translation)) != 0) ++bad;
        }
        failures += bad > 0;
        std::println("{:>8} {:>8} {:>7} {:>5}+{:<5} {:>9.1f} {:>9.1f} {:>9.1f}  {}/{} hit, {}", count, step,
                     u.kept, u.removed, u.inserted, cook_ms, update_ms, fresh_ms, hits, probe.size(),
                     bad ? "DISAGREE" : (u.rebuilt ? "agree (rebuilt)" : "agree"));
        scene = std::move(next);
    }
    return failures ? 1 : 0;
}
