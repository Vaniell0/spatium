// donut_live -- the donut scene traced by a Vulkan compute shader.
//
// The first thing the device does in this repository that is not CUDA,
// and the thing it has to earn before a window is worth opening: the same
// picture as the host. So the default run is headless and prints three
// comparisons, each isolating one layer:
//
//   fp64 trees vs the fp32 host trace   -- what dropping to fp32 costs
//   the fp32 host trace vs the shader   -- whether the shader is the spec
//   frame time, shader against host     -- what the device buys
//
//   donut_live [--t seconds] [--photo PATH] [--runs N]
//
// Shading is primary rays with the demo's key and fill lights -- no
// shadows, highlights, reflections or see-through dust yet -- so this is
// a check of the traversal and not the finished picture; `donut_demo
// --photo` remains the render.
#include "donut_scene.hpp"

#include <spatium/render/camera.hpp>
#include <spatium/render/cooked_scene.hpp>
#include <spatium/render/gpu_scene.hpp>
#include <spatium/render/gpu_trace_glsl.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/write_image.hpp>
#include <spatium/viewer/compute.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <print>
#include <string>
#include <string_view>
#include <vector>

using namespace spatium;
namespace bd = spatium::io::build;
namespace gpu = spatium::render::gpu;
namespace vc = spatium::viewer::compute;

namespace {

// The push-constant block, laid out as the shader's `Push`.
struct Push {
    float cam_pos[4], fwd[4], right[4], up[4], key[4], fill[4], background[4];
    std::uint32_t size[4];
};
static_assert(sizeof(Push) == 128);

// The shader's primary ray, in the shader's arithmetic. The host's own
// `camera_pixel_dir` works in fp64, and a comparison against it would
// charge the camera's rounding to the traversal.
gpu::Ray32 pixel_ray(const Push& pc, int x, int y) {
    const float W = static_cast<float>(pc.size[0]), H = static_cast<float>(pc.size[1]);
    const float nx = (2.0f * (static_cast<float>(x) + 0.5f) / W - 1.0f) * pc.fwd[3] * pc.cam_pos[3];
    const float ny = (1.0f - 2.0f * (static_cast<float>(y) + 0.5f) / H) * pc.cam_pos[3];
    gpu::V3 d = gpu::v3(pc.fwd) + gpu::v3(pc.right) * nx + gpu::v3(pc.up) * ny;
    d = d * (1.0f / std::sqrt(gpu::dot(d, d)));
    return {gpu::v3(pc.cam_pos), d};
}

std::uint32_t pack_rgba(gpu::V3 c) {
    auto q = [](float v) { return static_cast<std::uint32_t>(std::clamp(v, 0.0f, 1.0f) * 255.0f + 0.5f); };
    return q(c.x) | (q(c.y) << 8) | (q(c.z) << 16) | (255u << 24);
}

std::uint32_t hit_id(const gpu::Hit& h) {
    const std::uint32_t kind = h.kind == gpu::HitKind::None ? 0u
                             : h.kind == gpu::HitKind::Triangle ? 1u : 2u;
    return (kind << 30) | (kind ? h.index : 0u);
}

double ms_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

int main(int argc, char** argv) {
    double t = 1.5;
    std::string photo;
    int runs = 5;
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--t" && i + 1 < argc) { t = std::atof(argv[++i]); continue; }
        if (a == "--photo" && i + 1 < argc) { photo = argv[++i]; continue; }
        if (a == "--runs" && i + 1 < argc) { runs = std::max(1, std::atoi(argv[++i])); continue; }
        std::print("donut_live [--t seconds] [--photo PATH] [--runs N]\n");
        return a == "--help" ? 0 : 1;
    }

    auto boom = donut::load_boom_points("examples/data/boom_points.txt");
    if (boom.empty()) boom = donut::load_boom_points("data/boom_points.txt");
    if (boom.empty()) {
        std::println(stderr, "donut_live: examples/data/boom_points.txt not found -- run from the "
                             "repository root or from examples/");
        return 1;
    }
    bd::Trace<double> scene;
    const auto root = donut::build_scene(scene, 11000, 35200, boom);

    auto t0 = std::chrono::steady_clock::now();
    const auto cooked = bd::cook(scene, root, t);
    const double ms_cook = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const auto laid = render::lay_out(cooked);
    const double ms_lay = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const auto packed = gpu::pack(laid);
    const double ms_pack = ms_since(t0);
    std::println("t={}: cook {:.1f} ms, lay_out {:.1f} ms, pack {:.1f} ms -- {} triangles, {} "
                 "instances, {:.1f} MB on the device",
                 t, ms_cook, ms_lay, ms_pack, packed.triangles.size(), packed.instances.size(),
                 static_cast<double>(packed.bytes()) / 1e6);

    constexpr int W = 960, H = 720;
    const auto cam = donut::hero_camera();
    const auto basis = render::make_camera_basis(cam);
    const Vec<double, 3> key = Vec<double, 3>{Vec<double, 3>{0.85, 0.45, 0.55}.normalized()};
    const Vec<double, 3> fill = Vec<double, 3>{Vec<double, 3>{0.62, -0.70, 0.35}.normalized()};
    Push pc{};
    gpu::detail::put3(pc.cam_pos, cam.position, static_cast<float>(basis.tan_half));
    gpu::detail::put3(pc.fwd, basis.fwd, static_cast<float>(W) / static_cast<float>(H));
    gpu::detail::put3(pc.right, basis.right);
    gpu::detail::put3(pc.up, basis.up);
    gpu::detail::put3(pc.key, key, 0.38f);
    gpu::detail::put3(pc.fill, fill);
    gpu::detail::put3(pc.background, Vec<double, 3>{0.04, 0.04, 0.05});
    pc.size[0] = W;
    pc.size[1] = H;
    pc.size[2] = static_cast<std::uint32_t>(packed.tri_nodes.size());
    pc.size[3] = static_cast<std::uint32_t>(packed.inst_nodes.size());
    const gpu::Lights lights{gpu::v3(pc.key), gpu::v3(pc.fill), pc.key[3], gpu::v3(pc.background)};

    // ── Host: the fp32 specification, and the fp64 trees beside it ──
    const std::size_t npix = static_cast<std::size_t>(W) * H;
    std::vector<std::uint32_t> host_rgba(npix), host_ids(npix);
    std::vector<char> fp64_agrees(npix, 0);
    t0 = std::chrono::steady_clock::now();
    render::parallel_for_rows(H, [&](int y) {
        for (int x = 0; x < W; ++x) {
            const auto r = pixel_ray(pc, x, y);
            const auto h = gpu::trace(packed, r);
            const std::size_t i = static_cast<std::size_t>(y) * W + x;
            host_rgba[i] = pack_rgba(gpu::shade(packed, r, h, lights));
            host_ids[i] = hit_id(h);
        }
    });
    const double ms_host = ms_since(t0);

    render::parallel_for_rows(H, [&](int y) {
        for (int x = 0; x < W; ++x) {
            const auto r32 = pixel_ray(pc, x, y);
            geometry::Ray<3, double> r{cam.position,
                                       Vec<double, 3>{r32.d.x, r32.d.y, r32.d.z}};
            const auto a = laid.triangles.ray_cast(r);
            const auto b = laid.instances.ray_cast(r);
            const std::size_t i = static_cast<std::size_t>(y) * W + x;
            const std::uint32_t id = host_ids[i];
            const std::uint32_t kind = id >> 30, idx = id & 0x3fffffffu;
            if (b && (!a || b->t < a->t))
                fp64_agrees[i] = kind == 2 && packed.instance_source[idx] == b->index;
            else if (a)
                fp64_agrees[i] = kind == 1 && packed.triangle_source[idx] == a->index;
            else
                fp64_agrees[i] = kind == 0;
        }
    });
    const auto fp64_diff = static_cast<std::size_t>(std::count(fp64_agrees.begin(), fp64_agrees.end(), 0));

    // ── Device ──
    vc::Context ctx("donut_live");
    auto b_tri_nodes = vc::Buffer::from(ctx, std::span<const gpu::Node>(packed.tri_nodes));
    auto b_tris = vc::Buffer::from(ctx, std::span<const gpu::Triangle>(packed.triangles));
    auto b_inst_nodes = vc::Buffer::from(ctx, std::span<const gpu::Node>(packed.inst_nodes));
    auto b_quads = vc::Buffer::from(ctx, std::span<const gpu::Quadric>(packed.quadrics));
    auto b_insts = vc::Buffer::from(ctx, std::span<const gpu::Instance>(packed.instances));
    vc::Buffer b_pixels(ctx, npix * sizeof(std::uint32_t));
    vc::Buffer b_ids(ctx, npix * sizeof(std::uint32_t));
    vc::Kernel kernel(ctx, gpu::kTraceGlsl, "gpu_trace.comp", 7, sizeof(Push));
    vc::Buffer* bufs[] = {&b_tri_nodes, &b_tris, &b_inst_nodes, &b_quads, &b_insts, &b_pixels, &b_ids};
    kernel.bind(bufs);

    std::vector<double> gpu_ms;
    for (int k = 0; k < runs; ++k)
        gpu_ms.push_back(ctx.run([&](VkCommandBuffer cmd) {
            kernel.dispatch(cmd, &pc, (W + 7) / 8, (H + 7) / 8);
        }));
    std::sort(gpu_ms.begin(), gpu_ms.end());

    const auto* dev_rgba = static_cast<const std::uint32_t*>(b_pixels.data());
    const auto* dev_ids = static_cast<const std::uint32_t*>(b_ids.data());
    std::size_t id_diff = 0, px_diff = 0;
    int worst = 0;
    for (std::size_t i = 0; i < npix; ++i) {
        if (dev_ids[i] != host_ids[i]) ++id_diff;
        if (dev_rgba[i] != host_rgba[i]) {
            ++px_diff;
            for (int c = 0; c < 3; ++c) {
                const int a = static_cast<int>((dev_rgba[i] >> (8 * c)) & 0xffu);
                const int b = static_cast<int>((host_rgba[i] >> (8 * c)) & 0xffu);
                worst = std::max(worst, std::abs(a - b));
            }
        }
    }

    std::println("device: {}", ctx.device_name());
    std::println("fp64 trees vs fp32 host: {} of {} pixels hit something else ({:.3f}%)",
                 fp64_diff, npix, 100.0 * static_cast<double>(fp64_diff) / static_cast<double>(npix));
    std::println("fp32 host vs shader:     {} pixels hit something else, {} differ in colour "
                 "(worst channel {} of 255)",
                 id_diff, px_diff, worst);
    std::println("frame, primary rays only: shader {:.2f} ms (median of {}), host {:.1f} ms on "
                 "all cores -- {:.0f}x",
                 gpu_ms[gpu_ms.size() / 2], runs, ms_host, ms_host / gpu_ms[gpu_ms.size() / 2]);

    if (!photo.empty()) {
        std::vector<std::uint8_t> rgb(npix * 3);
        for (std::size_t i = 0; i < npix; ++i)
            for (int c = 0; c < 3; ++c)
                rgb[i * 3 + c] = static_cast<std::uint8_t>((dev_rgba[i] >> (8 * c)) & 0xffu);
        render::write_png_rgb(photo, W, H, rgb);
        std::println("  -> {}", photo);
    }
}
