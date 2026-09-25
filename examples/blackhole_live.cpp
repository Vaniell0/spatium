// blackhole_live -- black holes described in the scene language and traced
// on the device. One program for one hole or a pair, replacing the two
// demos that each hard-coded a metric.
//
// Modes, all headless so far:
//
//   --frame PATH [--scene one|binary] [--spin a] [--t T] [--width W]
//       [--height H] [--steps N] [--no-disk]
//       Renders one frame of the scene to PNG and reports the device time.
//
//   --shadow-check
//       Renders a Schwarzschild hole with no disk and measures its shadow's
//       angular radius from the pixels a horizon ended, against the static
//       observer's sin(alpha) = 3 sqrt(3) M / D sqrt(1 - 2M/D).
//
//   --check [--rays N] [--steps S]
//       For a single Kerr hole and for a superposed binary, sends N light
//       rays from a camera S RK4 steps through the scene's metric three
//       ways -- on the device (the metric lowered to GLSL with its
//       derivatives, physics/relativity/metric_glsl.hpp), on the host in
//       float (geodesic.hpp itself, on Dual<float>), and on the host in
//       double -- and reports how far the three end points are apart, and
//       what the device took.

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <spatium/physics/relativity/geodesic.hpp>
#include <spatium/render/blackhole_glsl.hpp>
#include <spatium/render/write_image.hpp>
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

struct Frame {
    std::vector<std::uint32_t> pixels;
    double ms = 0;
};

// One frame of `scene` at coordinate time t, seen by an observer at rest.
Frame render_frame(vc::Context& ctx, const rel::SpacetimeScene<double>& scene, int W, int H, double t,
                   int max_steps, bool disk, double exposure) {
    const auto src = spatium::render::blackhole_shader(scene);
    if (!src) {
        std::println(stderr, "shader: {}", src.error().message);
        std::exit(1);
    }
    // The tetrad in the coordinates the shader traces in.
    const auto metric = *scene.metric(rel::SpacetimeScene<double>::Form::outgoing);
    const auto cam = scene.camera();
    const double deg = std::numbers::pi / 180.0;
    const double ca = std::cos(cam.azimuth_deg * deg), sa = std::sin(cam.azimuth_deg * deg);
    const double ce = std::cos(cam.elevation_deg * deg), se = std::sin(cam.elevation_deg * deg);
    const Vec<double, 4> at{t, cam.distance * ce * ca, cam.distance * ce * sa, cam.distance * se};
    const auto tet = spatium::render::observer_tetrad(metric, at, Vec<double, 3>{0.0, 0.0, 1.0});

    spatium::render::BlackholePush push{};
    for (int c = 0; c < 4; ++c) {
        push.e0[c] = static_cast<float>(tet[0][c]);
        push.e1[c] = static_cast<float>(tet[1][c]);
        push.e2[c] = static_cast<float>(tet[2][c]);
        push.e3[c] = static_cast<float>(tet[3][c]);
        push.cam[c] = static_cast<float>(at[c]);
    }
    push.view[0] = static_cast<float>(W);
    push.view[1] = static_cast<float>(H);
    push.view[2] = static_cast<float>(std::tan(0.5 * cam.fov_deg * deg));
    push.view[3] = static_cast<float>(exposure);
    const auto& d = scene.disk();
    const double M = scene.total_mass();
    push.disk[0] = static_cast<float>(d.inner * M);
    push.disk[1] = static_cast<float>(d.outer * M);
    push.disk[2] = 0.018f;
    push.disk[3] = static_cast<float>(d.temperature);
    push.flags[0] = 0;
    push.flags[1] = static_cast<std::uint32_t>(max_steps);
    push.flags[2] = disk && d.on ? 1u : 0u;
    push.flags[3] = scene.sky().seed;

    const auto table = spatium::render::blackbody_table();
    vc::Buffer pixels(ctx, static_cast<std::size_t>(W) * H * 4);
    auto bb = vc::Buffer::from(ctx, std::span<const float>(table));
    const std::vector<std::uint32_t> none{0};
    auto perm = vc::Buffer::from(ctx, std::span<const std::uint32_t>(none));
    auto pts = vc::Buffer::from(ctx, std::span<const std::uint32_t>(none));
    vc::Kernel kernel(ctx, src->c_str(), "blackhole.comp", 4, sizeof(push));
    vc::Buffer* bufs[] = {&pixels, &bb, &perm, &pts};
    kernel.bind(bufs);
    Frame f;
    f.ms = ctx.run([&](VkCommandBuffer cmd) {
        kernel.dispatch(cmd, &push, static_cast<std::uint32_t>((W + 7) / 8), static_cast<std::uint32_t>((H + 7) / 8));
    });
    const auto* px = static_cast<const std::uint32_t*>(pixels.data());
    f.pixels.assign(px, px + static_cast<std::size_t>(W) * H);
    return f;
}

bool write_frame(const std::string& path, const Frame& f, int W, int H) {
    std::vector<std::uint8_t> rgb(static_cast<std::size_t>(W) * H * 3);
    for (std::size_t i = 0; i < f.pixels.size(); ++i)
        for (int c = 0; c < 3; ++c) rgb[3 * i + c] = static_cast<std::uint8_t>((f.pixels[i] >> (8 * c)) & 0xffu);
    return spatium::render::write_png_rgb(path, W, H, rgb);
}

}  // namespace

int main(int argc, char** argv) {
    std::string mode, out;
    int rays = 4096, steps = 300, W = 640, H = 360, max_steps = 1500;
    std::string which = "one";
    double spin = 0.9, t = 0.0, exposure = 1.0;
    bool disk = true;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&] { return std::string(i + 1 < argc ? argv[++i] : ""); };
        if (a == "--check") mode = "check";
        else if (a == "--shadow-check") mode = "shadow";
        else if (a == "--frame") { mode = "frame"; out = next(); }
        else if (a == "--rays") rays = std::stoi(next());
        else if (a == "--steps") steps = max_steps = std::stoi(next());
        else if (a == "--width") W = std::stoi(next());
        else if (a == "--height") H = std::stoi(next());
        else if (a == "--scene") which = next();
        else if (a == "--spin") spin = std::stod(next());
        else if (a == "--t") t = std::stod(next());
        else if (a == "--exposure") exposure = std::stod(next());
        else if (a == "--no-disk") disk = false;
        else {
            std::println(stderr, "usage: blackhole_live --check | --shadow-check | --frame PATH [options]");
            return 1;
        }
    }
    vc::Context ctx("blackhole_live");

    if (mode == "check") {
        rel::SpacetimeScene<double> one;
        one.hole(1.0, 0.9);
        rel::SpacetimeScene<double> pair;
        pair.binary(0.5, 0.5, 12.0, /*inspiral=*/false);
        int rc = check_scene(ctx, "one Kerr hole, spin 0.9", one, rays, steps);
        rc |= check_scene(ctx, "binary, 12 M apart", pair, rays, steps);
        return rc;
    }

    if (mode == "shadow") {
        // A Schwarzschild hole, no disk, observer at rest at D: the shadow's
        // edge is where rays stop ending at the horizon.
        const double D = 30.0;
        rel::SpacetimeScene<double> scene;
        scene.hole(1.0, 0.0);
        scene.camera() = {.distance = D, .azimuth_deg = 0.0, .elevation_deg = 0.0, .fov_deg = 30.0};
        const int w = 2000, h = 1125;
        const auto f = render_frame(ctx, scene, w, h, 0.0, 4000, false, 1.0);
        int first = -1, last = -1;
        for (int x = 0; x < w; ++x)
            if ((f.pixels[static_cast<std::size_t>(h / 2) * w + x] >> 24) == 0u) {
                if (first < 0) first = x;
                last = x;
            }
        const double tan_half = std::tan(0.5 * 30.0 * std::numbers::pi / 180.0), aspect = double(w) / h;
        auto angle = [&](double px) { return std::atan(std::abs((px + 0.5) / w * 2.0 - 1.0) * tan_half * aspect); };
        const double measured = 0.5 * (angle(first - 0.5) + angle(last + 0.5));
        const double expected = std::asin(3.0 * std::sqrt(3.0) / D * std::sqrt(1.0 - 2.0 / D));
        // The angle one pixel spans at the shadow's edge.
        const double pixel = angle(last + 1.0) - angle(last);
        std::println("shadow of a Schwarzschild hole seen from D = {} M: measured {:.5f} rad, expected {:.5f} rad,"
                     " apart {:.2e} (a pixel is {:.2e}); {:.1f} ms",
                     D, measured, expected, std::abs(measured - expected), pixel, f.ms);
        return std::abs(measured - expected) < 2.0 * pixel ? 0 : 1;
    }

    if (mode == "frame") {
        rel::SpacetimeScene<double> scene;
        if (which == "binary") scene.binary(0.5, 0.5, 14.0, /*inspiral=*/false, 10.0, spin * 0.5, spin * 0.5);
        else scene.hole(1.0, spin);
        const auto f = render_frame(ctx, scene, W, H, t, max_steps, disk, exposure);
        if (!write_frame(out, f, W, H)) {
            std::println(stderr, "cannot write {}", out);
            return 1;
        }
        std::println("{}: {}x{}, {} steps at most, {:.1f} ms", out, W, H, max_steps, f.ms);
        return 0;
    }
    std::println(stderr, "usage: blackhole_live --check | --shadow-check | --frame PATH [options]");
    return 1;
}
