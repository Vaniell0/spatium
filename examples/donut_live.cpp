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
// And with --live, the same shader in a window: fly with WASD, Q/E for
// down and up, Shift to go faster, right mouse to look; the `t` slider
// re-cooks the scene on the host when it is let go.
//
//   donut_live --live [--t seconds] [--frames N] [--screenshot PATH]
//   either one: [--camera px py pz tx ty tz]
//
// Shading is primary rays with the demo's key and fill lights -- no
// shadows, highlights, reflections or see-through dust yet -- so this is
// a check of the traversal and not the finished picture; `donut_demo
// --photo` remains the render.
#include "donut_scene.hpp"

#include <spatium/render/camera.hpp>
#include <spatium/render/cooked_scene.hpp>
#include <spatium/render/gpu_instances.hpp>
#include <spatium/render/gpu_scene.hpp>
#include <spatium/render/lbvh.hpp>
#include <spatium/render/gpu_trace_glsl.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/write_image.hpp>
#include <spatium/viewer/compute.hpp>

#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>
#if SPATIUM_HAS_IMGUI
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_vulkan.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <numbers>
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


// ── Live ────────────────────────────────────────────────────────────

// The scene's arrays on the device, rebuilt as a unit when `t` changes.
struct SceneBuffers {
    vc::Buffer tri_nodes, tris, inst_nodes, quads, insts;
    SceneBuffers(vc::Context& ctx, const gpu::Scene& p)
        : tri_nodes(vc::Buffer::from(ctx, std::span<const gpu::Node>(p.tri_nodes))),
          tris(vc::Buffer::from(ctx, std::span<const gpu::Triangle>(p.triangles))),
          inst_nodes(vc::Buffer::from(ctx, std::span<const gpu::LNode>(p.inst_nodes))),
          quads(vc::Buffer::from(ctx, std::span<const gpu::Quadric>(p.quadrics))),
          insts(vc::Buffer::from(ctx, std::span<const gpu::Instance>(p.instances))) {}
};

gpu::Scene scene_at(const bd::Trace<double>& scene, std::size_t root, double t, double& ms) {
    const auto t0 = std::chrono::steady_clock::now();
    auto packed = gpu::pack(render::lay_out(bd::cook(scene, root, t)));
    ms = ms_since(t0);
    return packed;
}

// A free-flying camera: a position and two angles, z up. Kept as angles
// rather than as a basis so that looking around never rolls the horizon.
struct FlyCamera {
    Vec<double, 3> pos;
    double yaw = 0, pitch = 0;
    double fov_deg = 38;

    static FlyCamera from(const render::Camera<double>& c) {
        const Vec<double, 3> d = Vec<double, 3>{(c.target - c.position).normalized()};
        return {c.position, std::atan2(d[1], d[0]), std::asin(std::clamp(d[2], -1.0, 1.0)), c.fov_deg};
    }
    Vec<double, 3> forward() const {
        return {std::cos(pitch) * std::cos(yaw), std::cos(pitch) * std::sin(yaw), std::sin(pitch)};
    }
    render::Camera<double> camera() const {
        return {.position = pos, .target = Vec<double, 3>{pos + forward()}, .up = {0, 0, 1},
                .fov_deg = fov_deg};
    }
};

void fill_push(Push& pc, const render::Camera<double>& cam, std::uint32_t W, std::uint32_t H,
               bool bgra, const gpu::Scene& packed) {
    const auto basis = render::make_camera_basis(cam);
    const Vec<double, 3> key = Vec<double, 3>{Vec<double, 3>{0.85, 0.45, 0.55}.normalized()};
    const Vec<double, 3> fill = Vec<double, 3>{Vec<double, 3>{0.62, -0.70, 0.35}.normalized()};
    gpu::detail::put3(pc.cam_pos, cam.position, static_cast<float>(basis.tan_half));
    gpu::detail::put3(pc.fwd, basis.fwd, static_cast<float>(W) / static_cast<float>(H));
    gpu::detail::put3(pc.right, basis.right);
    gpu::detail::put3(pc.up, basis.up);
    gpu::detail::put3(pc.key, key, 0.38f);
    gpu::detail::put3(pc.fill, fill);
    gpu::detail::put3(pc.background, Vec<double, 3>{0.04, 0.04, 0.05}, bgra ? 1.0f : 0.0f);
    pc.size[0] = W;
    pc.size[1] = H;
    pc.size[2] = static_cast<std::uint32_t>(packed.tri_nodes.size());
    pc.size[3] = static_cast<std::uint32_t>(packed.inst_nodes.size());
}

int run_live(const bd::Trace<double>& scene, std::size_t root, double t, gpu::Scene packed,
             int max_frames, const std::string& screenshot, const render::Camera<double>& start) {
    if (!glfwInit()) {
        std::println(stderr, "donut_live: glfwInit failed");
        return 1;
    }
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    GLFWwindow* window = glfwCreateWindow(1280, 800, "donut_live", nullptr, nullptr);
    if (!window) {
        std::println(stderr, "donut_live: could not open a window");
        glfwTerminate();
        return 1;
    }
    int status = 0;
    {
        vc::Context ctx("donut_live", window);
        vc::Presenter present(ctx, window);
        auto sb = std::make_unique<SceneBuffers>(ctx, packed);
        std::uint32_t W = present.width(), H = present.height();
        auto pixels = std::make_unique<vc::Buffer>(ctx, std::size_t{W} * H * 4);
        auto ids = std::make_unique<vc::Buffer>(ctx, std::size_t{W} * H * 4);
        vc::Kernel kernel(ctx, gpu::kTraceGlsl, "gpu_trace.comp", 7, sizeof(Push));
        auto rebind = [&] {
            vc::Buffer* bufs[] = {&sb->tri_nodes, &sb->tris, &sb->inst_nodes, &sb->quads,
                                  &sb->insts, pixels.get(), ids.get()};
            kernel.bind(bufs);
        };
        rebind();

        FlyCamera cam = FlyCamera::from(start);
        const FlyCamera home = cam;
        float speed = 2.0f;
        float t_ui = static_cast<float>(t);
        double gpu_ms = 0, recook_ms = 0, fps = 0;
        double last_x = 0, last_y = 0;
        bool dragging = false;
        auto t_prev = std::chrono::steady_clock::now();
        std::println("device: {} -- WASD to fly, Q/E down/up, Shift faster, "
                     "right mouse to look", ctx.device_name());

        for (int frame = 0; !glfwWindowShouldClose(window); ++frame) {
            if (max_frames > 0 && frame >= max_frames) break;
            glfwPollEvents();
            const auto now = std::chrono::steady_clock::now();
            const double dt = std::chrono::duration<double>(now - t_prev).count();
            t_prev = now;
            fps = fps == 0 ? 1.0 / std::max(dt, 1e-6) : 0.9 * fps + 0.1 / std::max(dt, 1e-6);

            bool ui_mouse = false, ui_keys = false;
#if SPATIUM_HAS_IMGUI
            ImGui_ImplVulkan_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            ui_mouse = ImGui::GetIO().WantCaptureMouse;
            ui_keys = ImGui::GetIO().WantCaptureKeyboard;
            ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
            ImGui::Begin("donut_live");
            ImGui::Text("%s", ctx.device_name().c_str());
            ImGui::Text("trace %.2f ms, %.0f fps, %ux%u", gpu_ms, fps, W, H);
            ImGui::Text("%zu triangles, %zu instances", packed.triangles.size(),
                        packed.instances.size());
            ImGui::SliderFloat("t", &t_ui, 0.0f, 4.0f, "%.2f s");
            // Re-cooking costs a few hundred milliseconds on the host, so it
            // happens when the slider is let go rather than on every tick.
            if (ImGui::IsItemDeactivatedAfterEdit()) {
                packed = scene_at(scene, root, t_ui, recook_ms);
                vkDeviceWaitIdle(ctx.device());
                sb = std::make_unique<SceneBuffers>(ctx, packed);
                rebind();
            }
            if (recook_ms > 0) ImGui::Text("last re-cook %.0f ms (host)", recook_ms);
            ImGui::SliderFloat("speed", &speed, 0.1f, 20.0f, "%.1f /s", ImGuiSliderFlags_Logarithmic);
            if (ImGui::Button("reset camera")) cam = home;
            ImGui::Text("camera %.2f %.2f %.2f", cam.pos[0], cam.pos[1], cam.pos[2]);
            ImGui::End();
            ImGui::Render();
#endif
            double mx = 0, my = 0;
            glfwGetCursorPos(window, &mx, &my);
            const bool look = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
            if (look && !ui_mouse) {
                if (dragging) {
                    cam.yaw -= (mx - last_x) * 0.003;
                    cam.pitch = std::clamp(cam.pitch - (my - last_y) * 0.003, -1.55, 1.55);
                }
                dragging = true;
            } else {
                dragging = false;
            }
            last_x = mx;
            last_y = my;
            if (!ui_keys) {
                auto down = [&](int k) { return glfwGetKey(window, k) == GLFW_PRESS; };
                const Vec<double, 3> f = cam.forward();
                const Vec<double, 3> r = Vec<double, 3>{f.cross(Vec<double, 3>{0, 0, 1}).normalized()};
                const double step = speed * dt * (down(GLFW_KEY_LEFT_SHIFT) ? 4.0 : 1.0);
                Vec<double, 3> move{};
                if (down(GLFW_KEY_W)) move = Vec<double, 3>{move + f};
                if (down(GLFW_KEY_S)) move = Vec<double, 3>{move - f};
                if (down(GLFW_KEY_D)) move = Vec<double, 3>{move + r};
                if (down(GLFW_KEY_A)) move = Vec<double, 3>{move - r};
                if (down(GLFW_KEY_E)) move = Vec<double, 3>{move + Vec<double, 3>{0, 0, 1}};
                if (down(GLFW_KEY_Q)) move = Vec<double, 3>{move - Vec<double, 3>{0, 0, 1}};
                cam.pos = Vec<double, 3>{cam.pos + move * step};
                if (down(GLFW_KEY_ESCAPE)) glfwSetWindowShouldClose(window, 1);
            }

            Push pc{};
            fill_push(pc, cam.camera(), W, H, present.bgra(), packed);
            const bool ok = present.frame(*pixels, [&](VkCommandBuffer cmd) {
                kernel.dispatch(cmd, &pc, (W + 7) / 8, (H + 7) / 8);
            }, gpu_ms);
            if (!ok && (present.width() != W || present.height() != H)) {
                W = present.width();
                H = present.height();
                pixels = std::make_unique<vc::Buffer>(ctx, std::size_t{W} * H * 4);
                ids = std::make_unique<vc::Buffer>(ctx, std::size_t{W} * H * 4);
                rebind();
            }
        }

        if (!screenshot.empty()) {
            vkDeviceWaitIdle(ctx.device());
            const auto* px = static_cast<const std::uint32_t*>(pixels->data());
            std::vector<std::uint8_t> rgb(std::size_t{W} * H * 3);
            for (std::size_t i = 0; i < std::size_t{W} * H; ++i)
                for (int c = 0; c < 3; ++c) {
                    const int src = present.bgra() ? 2 - c : c;
                    rgb[i * 3 + c] = static_cast<std::uint8_t>((px[i] >> (8 * src)) & 0xffu);
                }
            render::write_png_rgb(screenshot, static_cast<int>(W), static_cast<int>(H), rgb);
            std::println("  -> {} (trace {:.2f} ms at {}x{})", screenshot, gpu_ms, W, H);
        }
    }
    glfwDestroyWindow(window);
    glfwTerminate();
    return status;
}


// ── The dust moved on the device ────────────────────────────────────

// Every Scatter whose objects are instanced exact forms, moved by a
// generated kernel and compared with what cook() put in the same slots.
// The comparison is within a tolerance and says what it found: the device
// is fp32 and its compiler may reorder arithmetic the host keeps in order.
int run_dust_check(const bd::Trace<double>& scene, std::size_t root, double t, int runs) {
    const auto cooked = bd::cook(scene, root, t);
    vc::Context ctx("donut_live");
    std::println("device: {}", ctx.device_name());

    for (std::size_t node = 0; node < scene.size(); ++node) {
        if (scene.node(node).kind != bd::Kind::Scatter) continue;
        std::vector<const bd::Object<double>*> objs;
        for (const auto& o : cooked.objects())
            if (o.source_node == node) objs.push_back(&o);
        if (objs.empty() || !objs.front()->instanceable) continue;
        const auto shape = objs.front()->shape;
        if (!cooked.shapes()[shape].exact.has_value()) continue;

        auto kernel_src = gpu::make_instance_kernel(
            scene, node, Vec<double, 3>{cooked.shapes()[shape].geometry.centroid()},
            static_cast<std::uint32_t>(shape));
        if (!kernel_src) {
            std::println("node {}: not movable on the device -- {}", node, kernel_src.error().message);
            continue;
        }
        auto& ks = *kernel_src;
        const std::size_t count = ks.sites.size();
        auto b_perm = vc::Buffer::from(ctx, std::span<const std::uint32_t>(ks.module.perm));
        auto b_points = vc::Buffer::from(ctx, std::span<const float>(ks.module.points));
        auto b_sites = vc::Buffer::from(ctx, std::span<const gpu::Site>(ks.sites));
        vc::Buffer b_out(ctx, count * sizeof(gpu::Instance));
        vc::Kernel kernel(ctx, ks.source.c_str(), "instances.comp", 4, sizeof(gpu::InstancePush));
        vc::Buffer* bufs[] = {&b_perm, &b_points, &b_sites, &b_out};
        kernel.bind(bufs);
        auto push = ks.push;
        push.rest_centroid_t[3] = static_cast<float>(t);
        std::vector<double> ms;
        for (int r = 0; r < runs; ++r)
            ms.push_back(ctx.run([&](VkCommandBuffer cmd) {
                kernel.dispatch(cmd, &push, static_cast<std::uint32_t>((count + 63) / 64), 1);
            }));
        std::sort(ms.begin(), ms.end());

        const auto* dev = static_cast<const gpu::Instance*>(b_out.data());
        double worst_t = 0, worst_r = 0, worst_s = 0, worst_c = 0, worst_e = 0;
        std::size_t over = 0;
        const std::size_t n_check = std::min(count, objs.size());
        for (std::size_t i = 0; i < n_check; ++i) {
            const auto& o = *objs[i];
            const auto& g = dev[i];
            const auto R = o.rotation_q.to_matrix();
            const float* rows[3] = {g.r0, g.r1, g.r2};
            double et = 0, er = 0;
            for (std::size_t a = 0; a < 3; ++a) {
                et = std::max(et, std::abs(rows[a][3] - o.translation[a]));
                for (std::size_t b = 0; b < 3; ++b)
                    er = std::max(er, std::abs(rows[a][b] - R(a, b)));
            }
            const double es = std::abs(g.scale_quadric[0] - o.scale);
            double ec = 0, ee = 0;
            for (std::size_t a = 0; a < 3; ++a) {
                ec = std::max(ec, std::abs(g.color_rough[a] - o.material.base_color[a]));
                ee = std::max(ee, std::abs(g.emissive_opacity[a] - o.material.emissive[a]));
            }
            worst_t = std::max(worst_t, et);
            worst_r = std::max(worst_r, er);
            worst_s = std::max(worst_s, es);
            worst_c = std::max(worst_c, ec);
            worst_e = std::max(worst_e, ee);
            if (et > 1e-3 || er > 1e-3 || es > 1e-4 || ec > 1e-3 || ee > 1e-3) ++over;
        }
        // The tree over what the device just wrote, built on the host from
        // the mapped buffer. Timed in two halves, because on shared memory
        // the read can cost as much as the build if the mapping is uncached.
        if (count >= 10000) {
            const auto packed = gpu::pack(render::lay_out(cooked));
            auto t0 = std::chrono::steady_clock::now();
            std::vector<gpu::Instance> copy(dev, dev + count);
            const double ms_read = ms_since(t0);
            t0 = std::chrono::steady_clock::now();
            const auto tree = gpu::build_lbvh(copy, packed.quadrics);
            const double ms_build = ms_since(t0);
            std::println("  LBVH on the host: read {:.1f} ms ({:.0f} MB), build {:.1f} ms -- {} leaves, "
                         "{} culled, {} nodes",
                         ms_read, count * sizeof(gpu::Instance) / 1e6, ms_build, tree.leaves,
                         tree.culled, tree.nodes.size());
            std::println("    boxes {:.1f}, keys {:.1f}, sort {:.1f}, hierarchy {:.1f}, bounds {:.1f} ms",
                         tree.ms[0], tree.ms[1], tree.ms[2], tree.ms[3], tree.ms[4]);
        }
        std::println("node {}: {} instances moved in {:.2f} ms (median of {}), {} lines of generated GLSL",
                     node, count, ms[ms.size() / 2], runs,
                     std::count(ks.source.begin(), ks.source.end(), '\n'));
        std::println("  against cook() at t={}: worst |translation| {:.2e}, |rotation| {:.2e}, "
                     "|scale| {:.2e}, |colour| {:.2e}, |glow| {:.2e}; {} of {} past 1e-3",
                     t, worst_t, worst_r, worst_s, worst_c, worst_e, over, n_check);
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    double t = 1.5;
    std::string photo;
    int runs = 5;
    bool live = false;
    bool dust_check = false;
    std::size_t dust = 35200;
    int frames = 0;
    std::string screenshot;
    render::Camera<double> start = donut::hero_camera();
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--t" && i + 1 < argc) { t = std::atof(argv[++i]); continue; }
        if (a == "--photo" && i + 1 < argc) { photo = argv[++i]; continue; }
        if (a == "--runs" && i + 1 < argc) { runs = std::max(1, std::atoi(argv[++i])); continue; }
        if (a == "--live") { live = true; continue; }
        if (a == "--dust-check") { dust_check = true; continue; }
        if (a == "--dust" && i + 1 < argc) { dust = static_cast<std::size_t>(std::atol(argv[++i])); continue; }
        if (a == "--camera" && i + 6 < argc) {
            for (std::size_t k = 0; k < 3; ++k) start.position[k] = std::atof(argv[++i]);
            for (std::size_t k = 0; k < 3; ++k) start.target[k] = std::atof(argv[++i]);
            continue;
        }
        if (a == "--frames" && i + 1 < argc) { frames = std::atoi(argv[++i]); continue; }
        if (a == "--screenshot" && i + 1 < argc) { screenshot = argv[++i]; continue; }
        std::print("donut_live [--t seconds] [--photo PATH] [--runs N]\n"
                   "donut_live --live [--t seconds] [--frames N] [--screenshot PATH]\n"
                   "  --camera px py pz tx ty tz   start from this position, looking at t\n"
                   "donut_live --dust-check [--t seconds] [--dust N] [--runs N]\n");
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
    const auto root = donut::build_scene(scene, 11000, dust, boom);
    if (dust_check) return run_dust_check(scene, root, t, runs);

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

    if (live) return run_live(scene, root, t, packed, frames, screenshot, start);

    constexpr int W = 960, H = 720;
    const auto cam = start;
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
    auto b_inst_nodes = vc::Buffer::from(ctx, std::span<const gpu::LNode>(packed.inst_nodes));
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
