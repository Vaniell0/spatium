// cloth_live -- an XPBD sheet dropped on an analytic obstacle, the scene
// cloth_visual_demo.cpp drew before April's removal, back with its scene
// kept apart (cloth_scene.hpp) and its failure measured.
//
//   cloth_live [--scene sphere|torus|klein] [--contact projection|triangles]
//       The window: the sheet falling, the obstacle drawn as a mesh, and
//       live stiffness, bending, contact band and gravity. Space pauses,
//       R resets, 1/2/3 pick the obstacle.
//
//   cloth_live --measure [--seconds S]
//       For each obstacle and contact method, a stiff and a soft sheet run
//       S simulated seconds, and the table says what held: the deepest
//       penetration and how many vertices ended inside, both by the exact
//       distance to the analytic obstacle, the worst stretch of a
//       structural edge, whether anything exploded, and the cost.

#include "cloth_scene.hpp"

#include <spatium/viewer/app.hpp>
#if SPATIUM_HAS_IMGUI
#  include <imgui.h>
#endif

#include <chrono>
#include <cstdlib>
#include <print>
#include <string>

using namespace cloth;
using Clock = std::chrono::steady_clock;

namespace {

mesh::Mesh<Euclidean<3>> render_mesh(const Sheet& s) {
    mesh::Mesh<Euclidean<3>> m;
    for (const auto& p : s.parts) m.vertices.push_back(p.x);
    for (const auto& f : s.faces) m.faces.push_back({f[0], f[1], f[2]});
    return m;
}

Metrics run(const Config& c, const Obstacles& ob, Obstacle o, Contact mode, double seconds) {
    Sheet s = make_sheet(c);
    TriangleObstacle tri;
    if (mode == Contact::Triangles) tri = TriangleObstacle::build(ob.tessellation(o));
    Metrics m;
    const int steps = static_cast<int>(seconds / c.substep_dt);
    const auto t0 = Clock::now();
    for (int k = 0; k < steps && !m.exploded; ++k) {
        substep(s, c, ob, o, mode, &tri, &m);
        measure_into(m, s, c, ob, o);
    }
    m.ms_per_substep = std::chrono::duration<double, std::milli>(Clock::now() - t0).count() / steps;
    for (const auto& p : s.parts) {
        double d;
        if (exact_distance(ob, o, p.x, d) && d < -c.dhat) ++m.inside;
    }
    return m;
}

int measure(double seconds) {
    const Obstacles ob;
    std::println("{:<7} {:<11} {:<11} | {:>11} {:>11} {:>9} {:>7} {:>9} {:>9} | {:>9}", "obstacle", "contact",
                 "sheet", "before", "after", "tunnelled", "inside", "stretch", "exploded", "substep");
    for (Obstacle o : {Obstacle::Sphere, Obstacle::Torus, Obstacle::Klein})
        for (Contact mode : {Contact::Projection, Contact::Triangles})
            for (int grid : {19, 31})
                for (double compl_ : {5e-6, 1e-8}) {
                    Config c;
                    c.grid = grid;
                    c.spacing = 1.35 / (grid - 1);
                    c.struct_compl = compl_;
                    const auto m = run(c, ob, o, mode, seconds);
                    const bool measured = o != Obstacle::Klein;
                    auto dist = [&](double v) { return measured ? std::format("{:.4f} m", v) : std::string("--"); };
                    std::println("{:<7} {:<11} {:>2}x{:<2} {:<5} | {:>11} {:>11} {:>9} {:>7} {:>8.1f}% {:>9} | {:>7.3f}ms",
                                 name(o), name(mode), grid, grid, compl_ < 1e-7 ? "rigid" : "stiff",
                                 dist(m.worst_before), dist(m.worst_penetration),
                                 measured ? std::to_string(m.tunnelled) : std::string("--"),
                                 measured ? std::to_string(m.inside) : std::string("--"), 100.0 * m.max_stretch,
                                 m.exploded ? "yes" : "no", m.ms_per_substep);
                }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    Obstacle scene = Obstacle::Sphere;
    Contact mode = Contact::Projection;
    bool do_measure = false;
    double seconds = 2.0;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&] { return std::string(i + 1 < argc ? argv[++i] : ""); };
        if (a == "--measure") do_measure = true;
        else if (a == "--seconds") seconds = std::stod(next());
        else if (a == "--scene") {
            const auto v = next();
            scene = v == "torus" ? Obstacle::Torus : v == "klein" ? Obstacle::Klein : Obstacle::Sphere;
        } else if (a == "--contact") mode = next() == "triangles" ? Contact::Triangles : Contact::Projection;
        else {
            std::println(stderr, "usage: cloth_live [--scene sphere|torus|klein] [--contact projection|triangles] | --measure [--seconds S]");
            return 1;
        }
    }
    if (do_measure) return measure(seconds);

    const Obstacles ob;
    Config cfg;
    Sheet sheet = make_sheet(cfg);
    TriangleObstacle tri = TriangleObstacle::build(ob.tessellation(scene));

    viewer::App app("Spatium -- cloth on an analytic obstacle", 1280, 720);
    app.add_mesh(render_mesh(sheet), Vec4f{0.85f, 0.55f, 0.30f, 1.0f});
    for (Obstacle o : {Obstacle::Sphere, Obstacle::Torus, Obstacle::Klein})
        app.add_mesh(mesh::parametric_mesh(ob.chart(o), 64, 32), Vec4f{0.30f, 0.55f, 0.85f, 1.0f});
    auto show = [&] {
        for (int k = 0; k < 3; ++k) app.set_mesh_visible(static_cast<std::size_t>(1 + k), static_cast<int>(scene) == k);
    };
    show();
    app.fit_camera(2.5f);

    bool paused = false;
    Metrics live{};
    double rolling_ms = 0;
    auto reset = [&] {
        sheet = make_sheet(cfg);
        tri = TriangleObstacle::build(ob.tessellation(scene));
        live = {};
        app.update_mesh_vertices(0, viewer::mesh_to_render_data(render_mesh(sheet)));
        show();
    };
    app.set_key_callback([&](int key, int action, int) {
        if (action == 0) return;
        if (key == ' ') paused = !paused;
        if (key == 'R' || key == 'r') reset();
        if (key == '1') { scene = Obstacle::Sphere; reset(); }
        if (key == '2') { scene = Obstacle::Torus; reset(); }
        if (key == '3') { scene = Obstacle::Klein; reset(); }
    });
    auto last = Clock::now();
    app.set_frame_callback([&] {
        const auto now = Clock::now();
        const double dt = std::min(1.0 / 30.0, std::chrono::duration<double>(now - last).count());
        last = now;
        if (!paused) {
            const auto t0 = Clock::now();
            const int n = std::clamp(static_cast<int>(dt / cfg.substep_dt), 1, 16);
            for (int k = 0; k < n; ++k) {
                substep(sheet, cfg, ob, scene, mode, &tri);
                measure_into(live, sheet, cfg, ob, scene);
            }
            rolling_ms = 0.9 * rolling_ms + 0.1 * std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
        }
        app.update_mesh_vertices(0, viewer::mesh_to_render_data(render_mesh(sheet)));
    });
#if SPATIUM_HAS_IMGUI
    app.enable_imgui();
    app.set_gui_callback([&] {
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::Begin("cloth");
        int si = static_cast<int>(scene), ci = static_cast<int>(mode);
        const char* scenes[] = {"sphere", "torus", "klein"};
        const char* modes[] = {"projection (analytic)", "triangles (tessellated)"};
        if (ImGui::Combo("obstacle", &si, scenes, 3)) { scene = static_cast<Obstacle>(si); reset(); }
        if (ImGui::Combo("contact", &ci, modes, 2)) { mode = static_cast<Contact>(ci); reset(); }
        ImGui::Text("%zu particles, sim %.2f ms a frame", sheet.parts.size(), rolling_ms);
        ImGui::Text("deepest penetration %.4f m, worst stretch %.1f%%%s", live.worst_penetration,
                    100.0 * live.max_stretch, live.exploded ? ", EXPLODED" : "");
        float st = static_cast<float>(cfg.struct_compl), bd = static_cast<float>(cfg.bend_compl);
        float dh = static_cast<float>(cfg.dhat), gr = static_cast<float>(cfg.gravity);
        bool edit = false;
        edit |= ImGui::SliderFloat("stretch compliance", &st, 1e-7f, 1e-1f, "%.2e", ImGuiSliderFlags_Logarithmic);
        edit |= ImGui::SliderFloat("bend compliance", &bd, 1e-7f, 1.0f, "%.2e", ImGuiSliderFlags_Logarithmic);
        edit |= ImGui::SliderFloat("contact band", &dh, 0.005f, 0.2f, "%.3f m");
        edit |= ImGui::SliderFloat("gravity", &gr, -25.0f, 0.0f, "%.2f");
        if (edit) {
            cfg.struct_compl = st;
            cfg.bend_compl = bd;
            cfg.dhat = dh;
            cfg.gravity = gr;
            for (std::size_t k = 0; k < sheet.cons.size(); ++k)
                sheet.cons[k].compliance = k < sheet.n_struct ? cfg.struct_compl : cfg.bend_compl;
        }
        if (ImGui::Button(paused ? "resume" : "pause")) paused = !paused;
        ImGui::SameLine();
        if (ImGui::Button("reset")) reset();
        ImGui::End();
    });
#endif
    app.run();
    return 0;
}
