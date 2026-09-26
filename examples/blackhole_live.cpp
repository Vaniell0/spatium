// blackhole_live -- black holes described in the scene language and traced
// on the device. One program for one hole or a pair, replacing the two
// demos that each hard-coded a metric.
//
// Modes:
//
//   --live [--frames N] [--screenshot PATH]
//       A window: the scene rendered at a chosen resolution (144p by
//       default, the whole frame each time) and scaled to the window, with
//       the few settings that change what is seen -- one hole or a pair,
//       spin, the pair's masses, separation and inspiral, the disk and its
//       colour, the camera's angles, exposure, time -- and a button that
//       saves the current view rendered afresh at 1920x1080.
//
//   --video DIR [--frames N] [--fps F] [--orbit DEG] [--scene one|binary]
//       [--spin a] [--width W] [--height H] [--inspiral]
//       A sequence DIR/frame_%04d.png, time advancing with the frames and
//       the camera turning DEG degrees of azimuth over the whole sequence.
//       Vulkan, not CUDA: it runs on any device with a Vulkan driver, a
//       rented NVIDIA card included. Assemble with
//         ffmpeg -framerate F -i DIR/frame_%04d.png -pix_fmt yuv420p out.mp4
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
//   --dust-check
//       2048 dust particles on circular orbits about a Kerr hole (spin
//       0.9), moved on the device to coordinate time 200 M, against the
//       host integrating the same geodesics in double with a fine step.
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

#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>
#if SPATIUM_HAS_IMGUI
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_vulkan.h>
#endif
#include <chrono>
#include <ctime>
#include <filesystem>
#include <memory>

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

// ── What a person can set ────────────────────────────────────────
//
// Few on purpose: what changes the picture, nothing that tunes the solver.
struct Settings {
    int scene = 0;                 // 0 one hole, 1 a pair
    float spin = 0.9f;             // per hole
    float mass_ratio = 0.8f;       // the pair: m2 / m1
    float separation = 14.0f;      // the pair, in units of the total mass
    bool inspiral = false;
    bool disk = true;
    float temperature = 1.0f;
    float disk_outer = 30.0f;
    float distance = 60.0f, azimuth = 30.0f, elevation = 10.0f, fov = 50.0f;
    float exposure = 1.6f;
    int max_steps = 1500;
    int dust = 1;                  // index into kDustCounts
};

constexpr std::uint32_t kDustCounts[] = {0, 100000, 300000, 1000000};

rel::SpacetimeScene<double> make_scene(const Settings& st) {
    rel::SpacetimeScene<double> scene;
    if (st.scene == 1) {
        const double m1 = 1.0 / (1.0 + st.mass_ratio), m2 = 1.0 - m1;
        scene.binary(m1, m2, st.separation, st.inspiral, 6.0, st.spin * m1, st.spin * m2);
    } else {
        scene.hole(1.0, st.spin);
    }
    scene.disk().on = st.disk;
    scene.disk().temperature = st.temperature;
    scene.disk().outer = st.disk_outer;
    scene.dust().count = kDustCounts[st.dust];
    scene.dust().outer = std::max(40.0, 1.3 * st.disk_outer);
    scene.camera() = {.distance = st.distance, .azimuth_deg = st.azimuth, .elevation_deg = st.elevation,
                      .fov_deg = st.fov};
    return scene;
}

// The shader for one scene, its buffers, and one frame's push constants.
struct Renderer {
    vc::Context& ctx;
    rel::SpacetimeScene<double> scene;
    std::unique_ptr<vc::Kernel> kernel;
    std::unique_ptr<vc::Buffer> pixels, bb, perm, pts;
    std::uint32_t W = 0, H = 0;
    // Dust: particles moved on the device each frame, counted into a grid
    // the rays read.
    std::unique_ptr<vc::Kernel> move, deposit;
    std::unique_ptr<vc::Buffer> particles, grid;
    spatium::render::DustGrid grid_shape;
    std::uint32_t dust_count = 0;
    std::uint32_t frame_no = 0;

    Renderer(vc::Context& c, rel::SpacetimeScene<double> sc, std::uint32_t w, std::uint32_t h) : ctx(c), scene(std::move(sc)) {
        const auto src = spatium::render::blackhole_shader(scene);
        if (!src) {
            std::println(stderr, "shader: {}", src.error().message);
            std::exit(1);
        }
        const auto table = spatium::render::blackbody_table();
        bb = std::make_unique<vc::Buffer>(vc::Buffer::from(ctx, std::span<const float>(table)));
        const std::vector<std::uint32_t> none{0};
        perm = std::make_unique<vc::Buffer>(vc::Buffer::from(ctx, std::span<const std::uint32_t>(none)));
        pts = std::make_unique<vc::Buffer>(vc::Buffer::from(ctx, std::span<const std::uint32_t>(none)));
        kernel = std::make_unique<vc::Kernel>(ctx, src->c_str(), "blackhole.comp", 5,
                                              static_cast<std::uint32_t>(sizeof(spatium::render::BlackholePush)));
        grid_shape = spatium::render::dust_grid(scene);
        dust_count = scene.dust().count;
        grid = std::make_unique<vc::Buffer>(ctx, std::size_t{grid_shape.cells()} * 4);
        std::memset(grid->data(), 0, grid->size());
        if (dust_count > 0) {
            const auto msrc = spatium::render::dust_move_shader(scene);
            if (!msrc) {
                std::println(stderr, "dust shader: {}", msrc.error().message);
                std::exit(1);
            }
            move = std::make_unique<vc::Kernel>(ctx, msrc->c_str(), "dust_move.comp", 4, 16);
            const auto dsrc = spatium::render::dust_deposit_shader(grid_shape);
            deposit = std::make_unique<vc::Kernel>(ctx, dsrc.c_str(), "dust_deposit.comp", 2, 16);
            reset_dust();
            vc::Buffer* mb[] = {particles.get(), bb.get(), perm.get(), pts.get()};
            move->bind(mb);
            vc::Buffer* db[] = {particles.get(), grid.get()};
            deposit->bind(db);
        }
        resize(w, h);
    }
    void resize(std::uint32_t w, std::uint32_t h) {
        W = w;
        H = h;
        pixels = std::make_unique<vc::Buffer>(ctx, std::size_t{W} * H * 4);
        vc::Buffer* bufs[] = {pixels.get(), bb.get(), perm.get(), pts.get(), grid.get()};
        kernel->bind(bufs);
    }
    // Back to the first state: circular orbits at t = 0.
    void reset_dust() {
        if (dust_count == 0) return;
        const auto init = spatium::render::dust_initial(scene);
        if (!particles) particles = std::make_unique<vc::Buffer>(ctx, init.size() * sizeof(float));
        std::memcpy(particles->data(), init.data(), init.size() * sizeof(float));
    }
    spatium::render::BlackholePush push(double t, double exposure, int max_steps, bool bgr) const {
        const auto metric = *scene.metric(rel::SpacetimeScene<double>::Form::outgoing);
        const auto& cam = scene.camera();
        const double deg = std::numbers::pi / 180.0;
        const double ca = std::cos(cam.azimuth_deg * deg), sa = std::sin(cam.azimuth_deg * deg);
        const double ce = std::cos(cam.elevation_deg * deg), se = std::sin(cam.elevation_deg * deg);
        const Vec<double, 4> at{t, cam.distance * ce * ca, cam.distance * ce * sa, cam.distance * se};
        const auto tet = spatium::render::observer_tetrad(metric, at, Vec<double, 3>{0.0, 0.0, 1.0});
        spatium::render::BlackholePush p{};
        for (int c = 0; c < 4; ++c) {
            p.e0[c] = static_cast<float>(tet[0][c]);
            p.e1[c] = static_cast<float>(tet[1][c]);
            p.e2[c] = static_cast<float>(tet[2][c]);
            p.e3[c] = static_cast<float>(tet[3][c]);
            p.cam[c] = static_cast<float>(at[c]);
        }
        p.view[0] = static_cast<float>(W);
        p.view[1] = static_cast<float>(H);
        p.view[2] = static_cast<float>(std::tan(0.5 * cam.fov_deg * deg));
        p.view[3] = static_cast<float>(exposure);
        const auto& d = scene.disk();
        const double M = scene.total_mass();
        p.disk[0] = static_cast<float>(d.inner * M);
        p.disk[1] = static_cast<float>(d.outer * M);
        p.disk[2] = 0.018f;
        p.disk[3] = static_cast<float>(d.temperature);
        p.flags[0] = bgr ? 1u : 0u;
        p.flags[1] = static_cast<std::uint32_t>(max_steps);
        p.flags[2] = d.on ? 1u : 0u;
        p.flags[3] = scene.sky().seed;
        return p;
    }
    static void barrier(VkCommandBuffer cmd, VkPipelineStageFlags from, VkAccessFlags src) {
        VkMemoryBarrier mb{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
        mb.srcAccessMask = src;
        mb.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
        vkCmdPipelineBarrier(cmd, from, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &mb, 0, nullptr, 0, nullptr);
    }
    // The dust moved to the frame's time and counted, then the frame.
    void record(VkCommandBuffer cmd, const spatium::render::BlackholePush& p) {
        if (dust_count > 0) {
            struct { float target; std::uint32_t count, frame; float outer; } mp{
                p.cam[0], dust_count, frame_no++, static_cast<float>(scene.dust().outer * scene.total_mass())};
            move->dispatch(cmd, &mp, (dust_count + 63) / 64, 1);
            barrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_SHADER_WRITE_BIT);
            vkCmdFillBuffer(cmd, grid->handle(), 0, VK_WHOLE_SIZE, 0u);
            barrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_ACCESS_TRANSFER_WRITE_BIT);
            const std::uint32_t dp[4] = {dust_count, 0, 0, 0};
            deposit->dispatch(cmd, dp, (dust_count + 63) / 64, 1);
            barrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_SHADER_WRITE_BIT);
        }
        kernel->dispatch(cmd, &p, (W + 7) / 8, (H + 7) / 8);
    }
    double render(const spatium::render::BlackholePush& p) {
        return ctx.run([&](VkCommandBuffer cmd) { record(cmd, p); });
    }
    std::vector<std::uint32_t> read() const {
        const auto* px = static_cast<const std::uint32_t*>(pixels->data());
        return {px, px + std::size_t{W} * H};
    }
};

bool write_png(const std::string& path, const std::vector<std::uint32_t>& px, std::uint32_t W, std::uint32_t H,
               bool bgr) {
    std::vector<std::uint8_t> rgb(std::size_t{W} * H * 3);
    for (std::size_t i = 0; i < px.size(); ++i)
        for (int c = 0; c < 3; ++c) {
            const int src = bgr ? 2 - c : c;
            rgb[3 * i + c] = static_cast<std::uint8_t>((px[i] >> (8 * src)) & 0xffu);
        }
    return spatium::render::write_png_rgb(path, static_cast<int>(W), static_cast<int>(H), rgb);
}

// The low-resolution frame scaled up to the window, bilinearly.
const char* kUpscaleGlsl = R"GLSL(
#version 450
layout(local_size_x = 8, local_size_y = 8) in;
layout(std430, binding = 0) readonly buffer Src { uint src[]; };
layout(std430, binding = 1) writeonly buffer Dst { uint dst[]; };
layout(push_constant) uniform P { uvec4 size; } pc;   // source w h, target w h

vec4 at(ivec2 p) {
    p = clamp(p, ivec2(0), ivec2(pc.size.xy) - 1);
    return unpackUnorm4x8(src[uint(p.y) * pc.size.x + uint(p.x)]);
}
void main() {
    uvec2 d = gl_GlobalInvocationID.xy;
    if (d.x >= pc.size.z || d.y >= pc.size.w) return;
    vec2 uv = (vec2(d) + 0.5) / vec2(pc.size.zw) * vec2(pc.size.xy) - 0.5;
    ivec2 i = ivec2(floor(uv));
    vec2 f = uv - vec2(i);
    vec4 c = mix(mix(at(i), at(i + ivec2(1, 0)), f.x), mix(at(i + ivec2(0, 1)), at(i + ivec2(1, 1)), f.x), f.y);
    c.a = 1.0;
    dst[d.y * pc.size.z + d.x] = packUnorm4x8(c);
}
)GLSL";

struct Preset { const char* name; std::uint32_t h; };
// "auto" holds a frame rate by moving the render height; "window" is the
// window's own resolution.
constexpr std::uint32_t kAuto = 1;
constexpr Preset kPresets[] = {{"auto", kAuto}, {"144p", 144}, {"240p", 240}, {"360p", 360}, {"540p", 540},
                               {"720p", 720}, {"window", 0}};

std::string timestamp_name() {
    const std::time_t now = std::time(nullptr);
    char buf[64];
    std::strftime(buf, sizeof buf, "blackhole_%Y%m%d_%H%M%S.png", std::localtime(&now));
    return buf;
}

int run_live(int max_frames, const std::string& screenshot) {
    if (!glfwInit()) {
        std::println(stderr, "blackhole_live: glfwInit failed");
        return 1;
    }
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    GLFWwindow* window = glfwCreateWindow(1280, 720, "blackhole_live", nullptr, nullptr);
    if (!window) {
        std::println(stderr, "blackhole_live: could not open a window");
        glfwTerminate();
        return 1;
    }
    {
        vc::Context ctx("blackhole_live", window);
        vc::Presenter present(ctx, window);
        Settings st;
        int preset = 0;
        float auto_h = 144.0f, target_fps = 30.0f;
        bool fullscreen = false;
        int windowed[4] = {0, 0, 1280, 720};   // x, y, w, h before going full screen
        auto set_fullscreen = [&](bool on) {
            if (on == fullscreen) return;
            if (on) {
                glfwGetWindowPos(window, &windowed[0], &windowed[1]);
                glfwGetWindowSize(window, &windowed[2], &windowed[3]);
                GLFWmonitor* mon = glfwGetPrimaryMonitor();
                const GLFWvidmode* mode = glfwGetVideoMode(mon);
                glfwSetWindowMonitor(window, mon, 0, 0, mode->width, mode->height, mode->refreshRate);
            } else {
                glfwSetWindowMonitor(window, nullptr, windowed[0], windowed[1], windowed[2], windowed[3], 0);
            }
            fullscreen = on;
        };
        bool f11_was = false;
        std::uint32_t WW = present.width(), WH = present.height();
        auto render_size = [&] {
            const std::uint32_t want = kPresets[preset].h == kAuto ? static_cast<std::uint32_t>(auto_h)
                                                                     : kPresets[preset].h;
            const std::uint32_t h = want == 0 ? WH : std::min(want, WH);
            const std::uint32_t w = std::max(1u, static_cast<std::uint32_t>(std::lround(double(h) * WW / WH)));
            return std::pair{w, h};
        };
        auto [rw, rh] = render_size();
        auto renderer = std::make_unique<Renderer>(ctx, make_scene(st), rw, rh);
        auto window_px = std::make_unique<vc::Buffer>(ctx, std::size_t{WW} * WH * 4);
        vc::Kernel upscale(ctx, kUpscaleGlsl, "upscale.comp", 2, 16);
        auto rebind = [&] {
            vc::Buffer* b[] = {renderer->pixels.get(), window_px.get()};
            upscale.bind(b);
        };
        rebind();

        double t = 0.0, gpu_ms = 0.0, fps = 0.0;
        float rate = 1.0f;
        bool playing = true;
        std::string saved;
        auto t_prev = std::chrono::steady_clock::now();
        std::println("device: {}", ctx.device_name());

        for (int frame = 0; !glfwWindowShouldClose(window); ++frame) {
            if (max_frames > 0 && frame >= max_frames) break;
            glfwPollEvents();
            const auto now = std::chrono::steady_clock::now();
            const double dt = std::chrono::duration<double>(now - t_prev).count();
            t_prev = now;
            fps = fps == 0 ? 1.0 / std::max(dt, 1e-6) : 0.9 * fps + 0.1 / std::max(dt, 1e-6);
            if (playing) t += dt * rate * 20.0;   // 20 M of coordinate time a second at rate 1

            bool rebuild = false, resize = false, save = false;
#if SPATIUM_HAS_IMGUI
            ImGui_ImplVulkan_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
            ImGui::Begin("black holes");
            ImGui::Text("%s", ctx.device_name().c_str());
            ImGui::Text("%.0f fps, render %ux%u, %.1f ms", fps, renderer->W, renderer->H, gpu_ms);
            auto changed = [&](bool edited) { if (edited && !ImGui::IsItemActive()) rebuild = true;
                                              if (ImGui::IsItemDeactivatedAfterEdit()) rebuild = true; };
            ImGui::SeparatorText("holes");
            changed(ImGui::RadioButton("one", &st.scene, 0));
            ImGui::SameLine();
            changed(ImGui::RadioButton("pair", &st.scene, 1));
            changed(ImGui::SliderFloat("spin", &st.spin, 0.0f, 0.99f, "%.2f"));
            if (st.scene == 1) {
                changed(ImGui::SliderFloat("mass ratio", &st.mass_ratio, 0.1f, 1.0f, "%.2f"));
                changed(ImGui::SliderFloat("separation", &st.separation, 6.0f, 40.0f, "%.1f M"));
                changed(ImGui::Checkbox("inspiral", &st.inspiral));
            }
            ImGui::SeparatorText("disk");
            // The disk is push constants: no rebuild.
            ImGui::Checkbox("disk", &st.disk);
            ImGui::SliderFloat("temperature", &st.temperature, 0.3f, 3.0f, "%.2f");
            ImGui::SliderFloat("outer edge", &st.disk_outer, 10.0f, 60.0f, "%.0f M");
            ImGui::SeparatorText("dust");
            {
                const char* names[] = {"none", "100k", "300k", "1M"};
                for (int i = 0; i < 4; ++i) {
                    if (i) ImGui::SameLine();
                    if (ImGui::RadioButton(names[i], &st.dust, i)) rebuild = true;
                }
            }
            ImGui::SeparatorText("camera");
            bool cam = false;
            cam |= ImGui::SliderFloat("distance", &st.distance, 12.0f, 200.0f, "%.0f M");
            cam |= ImGui::SliderFloat("azimuth", &st.azimuth, -180.0f, 180.0f, "%.0f deg");
            cam |= ImGui::SliderFloat("elevation", &st.elevation, -89.0f, 89.0f, "%.0f deg");
            cam |= ImGui::SliderFloat("field of view", &st.fov, 10.0f, 120.0f, "%.0f deg");
            ImGui::SliderFloat("exposure", &st.exposure, 0.1f, 8.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
            ImGui::SeparatorText("time");
            ImGui::Checkbox("play", &playing);
            ImGui::SameLine();
            ImGui::SliderFloat("rate", &rate, 0.05f, 5.0f, "%.2fx", ImGuiSliderFlags_Logarithmic);
            ImGui::Text("t = %.0f M", t);
            if (ImGui::Button("t = 0")) { t = 0.0; renderer->reset_dust(); }
            ImGui::SeparatorText("render");
            for (int i = 0; i < static_cast<int>(std::size(kPresets)); ++i) {
                if (i) ImGui::SameLine();
                if (ImGui::RadioButton(kPresets[i].name, &preset, i)) resize = true;
            }
            if (kPresets[preset].h == kAuto)
                ImGui::SliderFloat("target fps", &target_fps, 10.0f, 60.0f, "%.0f");
            ImGui::SeparatorText("window");
            ImGui::Text("%ux%u", WW, WH);
            if (!fullscreen) {
                if (ImGui::Button("1280x720")) glfwSetWindowSize(window, 1280, 720);
                ImGui::SameLine();
                if (ImGui::Button("1600x900")) glfwSetWindowSize(window, 1600, 900);
                ImGui::SameLine();
                if (ImGui::Button("1920x1080")) glfwSetWindowSize(window, 1920, 1080);
            }
            if (ImGui::Button(fullscreen ? "leave full screen (F11)" : "full screen (F11)")) set_fullscreen(!fullscreen);
            if (ImGui::Button("save frame (1920x1080)")) save = true;
            if (!saved.empty()) ImGui::Text("saved %s", saved.c_str());
            ImGui::End();
            ImGui::Render();
            // Camera settings are read from the scene each frame; no rebuild.
            if (cam) renderer->scene.camera() = make_scene(st).camera();
#endif
            renderer->scene.disk().on = st.disk;
            renderer->scene.disk().temperature = st.temperature;
            renderer->scene.disk().outer = st.disk_outer;
            {
                const bool f11 = glfwGetKey(window, GLFW_KEY_F11) == GLFW_PRESS;
                if (f11 && !f11_was) set_fullscreen(!fullscreen);
                f11_was = f11;
            }
            // Auto: every half second, move the render height toward what the
            // frame budget allows -- down quickly when over it, up slowly
            // when well under -- and resize only on a change worth it.
            if (kPresets[preset].h == kAuto && frame % 30 == 29 && gpu_ms > 0.0) {
                const double budget = 1000.0 / target_fps;
                float h = auto_h;
                if (gpu_ms > 0.9 * budget) h *= 0.8f;
                else if (gpu_ms < 0.5 * budget) h *= 1.15f;
                h = std::clamp(h, 72.0f, static_cast<float>(WH));
                if (std::abs(h - auto_h) > 0.08f * auto_h) {
                    auto_h = h;
                    resize = true;
                }
            }
            if (rebuild) {
                vkDeviceWaitIdle(ctx.device());
                auto [w, h] = render_size();
                renderer = std::make_unique<Renderer>(ctx, make_scene(st), w, h);
                rebind();
            }
            if (resize) {
                vkDeviceWaitIdle(ctx.device());
                auto [w, h] = render_size();
                renderer->resize(w, h);
                rebind();
            }
            if (save || (!screenshot.empty() && max_frames > 0 && frame == max_frames - 1)) {
                vkDeviceWaitIdle(ctx.device());
                Renderer big(ctx, make_scene(st), 1920, 1080);
                // The dust as it is now, not as it started.
                if (big.dust_count > 0 && big.dust_count == renderer->dust_count)
                    std::memcpy(big.particles->data(), renderer->particles->data(), big.particles->size());
                big.render(big.push(t, st.exposure, std::max(st.max_steps, 3000), false));
                saved = screenshot.empty() ? timestamp_name() : screenshot;
                write_png(saved, big.read(), 1920, 1080, false);
                std::println("saved {}", saved);
            }

            const auto p = renderer->push(t, st.exposure, st.max_steps, present.bgra());
            const std::uint32_t up[4] = {renderer->W, renderer->H, WW, WH};
            const bool ok = present.frame(*window_px, [&](VkCommandBuffer cmd) {
                renderer->record(cmd, p);
                VkMemoryBarrier mb{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
                mb.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
                mb.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
                vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                                     0, 1, &mb, 0, nullptr, 0, nullptr);
                upscale.dispatch(cmd, up, (WW + 7) / 8, (WH + 7) / 8);
            }, gpu_ms);
            if (!ok && (present.width() != WW || present.height() != WH)) {
                vkDeviceWaitIdle(ctx.device());
                WW = present.width();
                WH = present.height();
                window_px = std::make_unique<vc::Buffer>(ctx, std::size_t{WW} * WH * 4);
                auto [w, h] = render_size();
                renderer->resize(w, h);
                rebind();
            }
        }
        vkDeviceWaitIdle(ctx.device());
        std::println("{:.0f} fps at the end, render {}x{} in {:.1f} ms", fps, renderer->W, renderer->H, gpu_ms);
    }
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    std::string mode, out, screenshot;
    double fps_video = 30.0, orbit = 90.0;
    bool inspiral = false;
    int rays = 4096, steps = 300, W = 640, H = 360, max_steps = 1500, frames = 0;
    std::string which = "one";
    double spin = 0.9, t = 0.0, exposure = 1.6;
    bool disk = true;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&] { return std::string(i + 1 < argc ? argv[++i] : ""); };
        if (a == "--check") mode = "check";
        else if (a == "--shadow-check") mode = "shadow";
        else if (a == "--dust-check") mode = "dust";
        else if (a == "--live") mode = "live";
        else if (a == "--frame") { mode = "frame"; out = next(); }
        else if (a == "--video") { mode = "video"; out = next(); }
        else if (a == "--fps") fps_video = std::stod(next());
        else if (a == "--orbit") orbit = std::stod(next());
        else if (a == "--inspiral") inspiral = true;
        else if (a == "--frames") frames = std::stoi(next());
        else if (a == "--screenshot") screenshot = next();
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
            std::println(stderr, "usage: blackhole_live --live | --check | --shadow-check | --frame PATH [options]");
            return 1;
        }
    }
    if (mode == "live") return run_live(frames, screenshot);
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

    if (mode == "dust") {
        Settings st;
        st.dust = 0;
        rel::SpacetimeScene<double> scene = make_scene(st);
        scene.dust().count = 2048;
        scene.dust().outer = 40.0;
        Renderer r(ctx, scene, 16, 16);
        const auto init = spatium::render::dust_initial(r.scene);
        const double T = 200.0;
        // The device moves at most 32 steps a dispatch; dispatch until every
        // particle has reached T.
        double ms = 0;
        for (int k = 0; k < 64; ++k) {
            struct { float target; std::uint32_t count, frame; float outer; } mp{
                static_cast<float>(T), r.dust_count, 0u, static_cast<float>(r.scene.dust().outer * r.scene.total_mass())};
            ms += ctx.run([&](VkCommandBuffer cmd) { r.move->dispatch(cmd, &mp, (r.dust_count + 63) / 64, 1); });
        }
        const auto* dev = static_cast<const float*>(r.particles->data());
        const auto metric = *r.scene.metric(rel::SpacetimeScene<double>::Form::outgoing);
        std::vector<double> err;
        double worst = 0;
        int reborn = 0;
        for (std::uint32_t i = 0; i < r.dust_count; ++i) {
            Vec<double, 8> st8{};
            for (int c = 0; c < 8; ++c) st8[c] = init[8 * i + c];
            const double R0 = std::hypot(st8[1], st8[2]);
            Vec<double, 8> prev = st8;
            while (st8[0] < T) {
                prev = st8;
                st8 = rel::geodesic_step(metric, st8, 0.02 * std::max(1.0, R0 / 4.0));
            }
            // Linear in coordinate time between the last two steps.
            const double a = (T - prev[0]) / (st8[0] - prev[0]);
            Vec<double, 3> host{};
            for (int c = 0; c < 3; ++c) host[c] = prev[1 + c] + a * (st8[1 + c] - prev[1 + c]);
            const Vec<double, 3> d{dev[8 * i + 1], dev[8 * i + 2], dev[8 * i + 3]};
            if (std::abs(dev[8 * i + 0] - T) > 1.0) { ++reborn; continue; }
            const double e = Vec<double, 3>{d - host}.norm() / R0;
            err.push_back(e);
            worst = std::max(worst, e);
        }
        std::sort(err.begin(), err.end());
        std::println("dust: {} particles on circular orbits about a Kerr hole (spin 0.9), moved to t = {} M; "
                     "{} reborn; device against host double, position apart relative to the orbit's radius: "
                     "worst {:.2e}, median {:.2e}; device {:.1f} ms",
                     r.dust_count, T, reborn, worst, err.empty() ? 0.0 : err[err.size() / 2], ms);
        return worst < 1e-2 ? 0 : 1;
    }

    if (mode == "shadow") {
        // A Schwarzschild hole, no disk, observer at rest at D: the shadow's
        // edge is where rays stop ending at the horizon.
        const double D = 30.0;
        Settings st;
        st.spin = 0.0f;
        st.disk = false;
        st.distance = static_cast<float>(D);
        st.azimuth = st.elevation = 0.0f;
        st.fov = 30.0f;
        const std::uint32_t w = 2000, h = 1125;
        Renderer r(ctx, make_scene(st), w, h);
        const double ms = r.render(r.push(0.0, 1.0, 4000, false));
        const auto px = r.read();
        int first = -1, last = -1;
        for (std::uint32_t x = 0; x < w; ++x)
            if ((px[std::size_t{h / 2} * w + x] >> 24) == 0u) {
                if (first < 0) first = static_cast<int>(x);
                last = static_cast<int>(x);
            }
        const double tan_half = std::tan(0.5 * 30.0 * std::numbers::pi / 180.0), aspect = double(w) / h;
        auto angle = [&](double p) { return std::atan(std::abs((p + 0.5) / w * 2.0 - 1.0) * tan_half * aspect); };
        const double measured = 0.5 * (angle(first - 0.5) + angle(last + 0.5));
        const double expected = std::asin(3.0 * std::sqrt(3.0) / D * std::sqrt(1.0 - 2.0 / D));
        const double pixel = angle(last + 1.0) - angle(last);
        std::println("shadow of a Schwarzschild hole seen from D = {} M: measured {:.5f} rad, expected {:.5f} rad,"
                     " apart {:.2e} (a pixel is {:.2e}); {:.1f} ms",
                     D, measured, expected, std::abs(measured - expected), pixel, ms);
        return std::abs(measured - expected) < 2.0 * pixel ? 0 : 1;
    }

    if (mode == "video") {
        Settings st;
        st.scene = which == "binary" ? 1 : 0;
        st.spin = static_cast<float>(which == "binary" ? spin * 0.5 : spin);
        st.disk = disk;
        st.inspiral = inspiral;
        const int n = frames > 0 ? frames : 300;
        std::filesystem::create_directories(out);
        const float azimuth0 = st.azimuth;
        double total_ms = 0;
        // One shader for the whole sequence: the camera turning and time
        // advancing are push constants.
        Renderer r(ctx, make_scene(st), static_cast<std::uint32_t>(W), static_cast<std::uint32_t>(H));
        for (int f = 0; f < n; ++f) {
            r.scene.camera().azimuth_deg = azimuth0 + orbit * f / std::max(1, n - 1);
            // 20 M of coordinate time a second of video, as the window's rate 1.
            const double at = t + 20.0 * f / fps_video;
            total_ms += r.render(r.push(at, exposure, std::max(max_steps, 3000), false));
            char name[64];
            std::snprintf(name, sizeof name, "frame_%04d.png", f);
            write_png((std::filesystem::path(out) / name).string(), r.read(), r.W, r.H, false);
            if (f % 10 == 0) std::println("frame {}/{}", f, n);
        }
        std::println("{} frames {}x{} in {}, {:.0f} ms of device time a frame; "
                     "ffmpeg -framerate {} -i {}/frame_%04d.png -pix_fmt yuv420p out.mp4",
                     n, W, H, out, total_ms / n, fps_video, out);
        return 0;
    }

    if (mode == "frame") {
        Settings st;
        st.scene = which == "binary" ? 1 : 0;
        st.spin = static_cast<float>(which == "binary" ? spin * 0.5 : spin);
        st.disk = disk;
        Renderer r(ctx, make_scene(st), static_cast<std::uint32_t>(W), static_cast<std::uint32_t>(H));
        const double ms = r.render(r.push(t, exposure, max_steps, false));
        if (!write_png(out, r.read(), r.W, r.H, false)) {
            std::println(stderr, "cannot write {}", out);
            return 1;
        }
        std::println("{}: {}x{}, {} steps at most, {:.1f} ms", out, W, H, max_steps, ms);
        return 0;
    }
    std::println(stderr, "usage: blackhole_live --live | --check | --shadow-check | --frame PATH [options]");
    return 1;
}
