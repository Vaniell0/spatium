// blackhole_live -- black holes described in the scene language and traced
// on the device. One program for one hole or a pair, replacing the two
// demos that each hard-coded a metric.
//
// This first form has one mode:
//
//   --check [--rays N] [--steps S]
//       For a single Kerr hole and for a superposed binary, sends N light
//       rays from a camera S RK4 steps through the scene's metric three
//       ways -- on the device (the metric lowered to GLSL with its
//       derivatives, physics/relativity/metric_glsl.hpp), on the host in
//       float (geodesic.hpp itself, on Dual<float>), and on the host in
//       double -- and reports how far the three end points are apart, and
//       what the device took.

#include <spatium/physics/relativity/geodesic.hpp>
#include <spatium/physics/relativity/metric_glsl.hpp>
#include <spatium/physics/relativity/spacetime_scene.hpp>
#include <spatium/viewer/compute.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <print>
#include <random>
#include <string>
#include <vector>

namespace rel = spatium::physics::relativity;
namespace vc = spatium::viewer::compute;
using spatium::Vec;

namespace {

struct Push {
    std::uint32_t count, steps;
    float dl, pad;
};

const char* kCheckMain = R"GLSL(
layout(std430, binding = 0) buffer States { vec4 s[]; };   // x, u per ray
layout(push_constant) uniform P { uint count; uint steps; float dl; float pad; } pc;

void main() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pc.count) return;
    vec4 x = s[2u * i], u = s[2u * i + 1u];
    for (uint k = 0u; k < pc.steps; ++k) geodesic_step(x, u, pc.dl);
    s[2u * i] = x;
    s[2u * i + 1u] = u;
}
)GLSL";

// The future-pointing null 4-velocity at x with spatial part n: the root
// of g(u, u) = 0 in u^0 that is positive.
template<typename Metric>
Vec<double, 4> null_velocity(const Metric& m, const Vec<double, 4>& x, const Vec<double, 3>& n) {
    const auto g = m(x);
    const double A = g(0, 0);
    double B = 0, C = 0;
    for (std::size_t i = 0; i < 3; ++i) {
        B += 2 * g(0, i + 1) * n[i];
        for (std::size_t j = 0; j < 3; ++j) C += g(i + 1, j + 1) * n[i] * n[j];
    }
    const double disc = std::sqrt(B * B - 4 * A * C);
    const double r1 = (-B + disc) / (2 * A), r2 = (-B - disc) / (2 * A);
    return {std::max(r1, r2), n[0], n[1], n[2]};
}

int check_scene(vc::Context& ctx, const char* label, const rel::SpacetimeScene<double>& scene, int rays,
                int steps) {
    const auto metric = scene.metric();
    if (!metric) {
        std::println(stderr, "{}: {}", label, metric.error().message);
        return 1;
    }
    const std::string glsl = std::string("#version 450\nlayout(local_size_x = 64) in;\n") +
                             rel::emit_metric_glsl(*metric, "metric") + rel::geodesic_glsl("metric") + kCheckMain;

    // Rays from a camera 30 M off, in a cone around the direction to the
    // origin wide enough that some pass close to a hole.
    std::mt19937_64 rng(3);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    const Vec<double, 4> cam{0.0, 0.0, -30.0, 6.0};
    std::vector<Vec<double, 8>> start;
    for (int k = 0; k < rays; ++k) {
        Vec<double, 3> n{-cam[1] / 30.0 + 0.35 * uni(rng), -cam[2] / 30.0 + 0.35 * uni(rng), -cam[3] / 30.0 + 0.35 * uni(rng)};
        n = Vec<double, 3>{n * (1.0 / n.norm())};
        const auto u = null_velocity(*metric, cam, n);
        start.push_back({cam[0], cam[1], cam[2], cam[3], u[0], u[1], u[2], u[3]});
    }
    const double dl = 0.2;

    // Device.
    std::vector<float> states(8 * static_cast<std::size_t>(rays));
    for (int k = 0; k < rays; ++k)
        for (int c = 0; c < 8; ++c) states[8 * k + c] = static_cast<float>(start[k][c]);
    auto buf = vc::Buffer::from(ctx, std::span<const float>(states));
    vc::Kernel kernel(ctx, glsl.c_str(), "blackhole_check.comp", 1, sizeof(Push));
    vc::Buffer* bufs[] = {&buf};
    kernel.bind(bufs);
    const Push push{static_cast<std::uint32_t>(rays), static_cast<std::uint32_t>(steps), static_cast<float>(dl), 0.0f};
    const double ms = ctx.run([&](VkCommandBuffer cmd) {
        kernel.dispatch(cmd, &push, static_cast<std::uint32_t>((rays + 63) / 64), 1);
    });
    const auto* dev = static_cast<const float*>(buf.data());

    // Host, float and double, the same steps. A ray that comes within 3 M
    // of a hole is set aside: inside the horizon it heads for the
    // singularity, where float and double part ways on the host just as
    // much, and what it says is about the spacetime, not the device.
    auto closest_hole = [&](const Vec<double, 8>& st) {
        double best = 1e300;
        for (const auto& h : scene.holes()) {
            const double t = st[0];
            const double dx = st[1] - h.x(t, 0.0), dy = st[2] - h.y(t, 0.0), dz = st[3] - h.z(t, 0.0);
            best = std::min(best, std::sqrt(dx * dx + dy * dy + dz * dz));
        }
        return best;
    };
    double worst_df = 0, worst_fd = 0, worst_dd = 0;
    std::vector<double> dev_err;
    int finite = 0, near = 0;
    for (int k = 0; k < rays; ++k) {
        Vec<double, 8> sd = start[k];
        Vec<float, 8> sf{};
        for (int c = 0; c < 8; ++c) sf[c] = static_cast<float>(start[k][c]);
        double closest = closest_hole(sd);
        for (int s = 0; s < steps; ++s) {
            sd = rel::geodesic_step(*metric, sd, dl);
            sf = rel::geodesic_step(*metric, sf, static_cast<float>(dl));
            closest = std::min(closest, closest_hole(sd));
        }
        if (!(closest > 3.0)) { ++near; continue; }
        bool ok = true;
        for (int c = 0; c < 4; ++c) ok = ok && std::isfinite(sd[c]) && std::isfinite(sf[c]) && std::isfinite(dev[8 * k + c]);
        if (!ok) continue;
        ++finite;
        // Distance between end points, relative to the path's length.
        const double len = dl * steps;
        double df = 0, fd = 0, dd = 0;
        for (int c = 1; c < 4; ++c) {
            df = std::max(df, static_cast<double>(std::abs(dev[8 * k + c] - sf[c])));
            fd = std::max(fd, std::abs(static_cast<double>(sf[c]) - sd[c]));
            dd = std::max(dd, std::abs(dev[8 * k + c] - sd[c]));
        }
        worst_df = std::max(worst_df, df / len);
        worst_fd = std::max(worst_fd, fd / len);
        worst_dd = std::max(worst_dd, dd / len);
        dev_err.push_back(dd / len);
    }
    std::sort(dev_err.begin(), dev_err.end());
    const double median = dev_err.empty() ? 0.0 : dev_err[dev_err.size() / 2];
    std::println("{}: {} rays x {} steps; {} came within 3 M of a hole and are set aside; for the other {},"
                 " end points apart relative to path length:", label, rays, steps, near, finite);
    std::println("  device vs host float  worst {:.2e}", worst_df);
    std::println("  host float vs double  worst {:.2e}", worst_fd);
    std::println("  device vs double      worst {:.2e}, median {:.2e}", worst_dd, median);
    std::println("  device time {:.2f} ms ({:.1f} ns per ray-step); metric pool {} ops", ms,
                 1e6 * ms / (double(rays) * steps), metric->pool().size());
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    bool check = false;
    int rays = 4096, steps = 300;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--check") check = true;
        else if (a == "--rays" && i + 1 < argc) rays = std::atoi(argv[++i]);
        else if (a == "--steps" && i + 1 < argc) steps = std::atoi(argv[++i]);
        else {
            std::println(stderr, "usage: blackhole_live --check [--rays N] [--steps S]");
            return 1;
        }
    }
    if (!check) {
        std::println(stderr, "blackhole_live: only --check exists so far");
        return 1;
    }
    vc::Context ctx("blackhole_live");

    rel::SpacetimeScene<double> one;
    one.hole(1.0, 0.9);
    rel::SpacetimeScene<double> pair;
    pair.binary(0.5, 0.5, 12.0, /*inspiral=*/false);

    int rc = check_scene(ctx, "one Kerr hole, spin 0.9", one, rays, steps);
    rc |= check_scene(ctx, "binary, 12 M apart", pair, rays, steps);
    return rc;
}
