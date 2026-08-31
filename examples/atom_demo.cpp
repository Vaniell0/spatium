// Atom visualization: rejection-sampled orbital point clouds + Bohr-model
// electron shells, by atomic number or element symbol (1-118). Vulkan
// viewer by default (--no-viewer for console-only), optional SVG export
// of the orbital cross-section and periodic table.

#include "io_helpers.hpp"

#include <spatium/spatium.hpp>
#include <spatium/viewer/app.hpp>
#include <spatium/physics/atomic/atom_model.hpp>
#include <spatium/physics/atomic/bohr_model.hpp>
#include <spatium/physics/atomic/atom_svg.hpp>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <filesystem>
#include <memory>
#include <mutex>
#include <print>
#include <string>
#include <thread>

#if SPATIUM_HAS_IMGUI
#include <imgui.h>
#endif

using namespace spatium;
using namespace spatium::viewer;
using namespace spatium::physics::atomic;

static void print_usage() {
    std::println("Usage: atom_demo [Z|symbol] [--points N] [--svg] [--no-viewer]"
                 " [--force] [--out-prefix PFX]");
    std::println("  Z or symbol: atomic number (1-118) or element symbol (e.g. Fe)");
    std::println("  --points N: points per orbital (default: 100000)");
    std::println("  --svg: export orbital cross-section and periodic table SVGs");
    std::println("  --no-viewer: console output only, no Vulkan window");
    std::println("  --force: overwrite existing SVG / screenshot files");
    std::println("  --out-prefix PFX: prefix prepended to every output filename");
}

static void append_superscript(std::string& s, int n) {
    static const char* super[] = {
        "\u2070", "\u00B9", "\u00B2", "\u00B3", "\u2074",
        "\u2075", "\u2076", "\u2077", "\u2078", "\u2079"
    };
    if (n < 10) {
        s += super[n];
    } else {
        s += super[n / 10];
        s += super[n % 10];
    }
}

static std::string electron_config_string(const Element& elem, bool abbreviated = true) {
    struct Core { const char* sym; int Z; int subs; };
    static const Core cores[] = {
        {"Rn", 86, 15}, {"Xe", 54, 11}, {"Kr", 36, 8},
        {"Ar", 18, 5},  {"Ne", 10, 3},  {"He", 2, 1},
    };

    auto config = elem.electron_config();
    int skip = 0;
    std::string result;

    if (abbreviated) {
        int total_e = elem.total_electrons();
        for (auto& c : cores) {
            if (total_e > c.Z && static_cast<int>(config.size()) > c.subs) {
                int core_sum = 0;
                for (int i = 0; i < c.subs; ++i)
                    core_sum += config[i].electrons;
                if (core_sum == c.Z) {
                    result = std::format("[{}] ", c.sym);
                    skip = c.subs;
                    break;
                }
            }
        }
    }

    for (int i = skip; i < static_cast<int>(config.size()); ++i) {
        auto& sub = config[i];
        if (i > skip) result += ' ';
        result += std::to_string(sub.n);
        result += subshell_letter(sub.l);
        append_superscript(result, sub.electrons);
    }
    return result;
}

static PointCloudData cloud_to_points(const OrbitalPointCloud<double>& cloud, bool positive) {
    PointCloudData pcd;
    for (std::size_t i = 0; i < cloud.positions.size(); ++i) {
        if (cloud.positive_lobe[i] != positive) continue;
        pcd.positions.push_back(static_cast<float>(cloud.positions[i][0]));
        pcd.positions.push_back(static_cast<float>(cloud.positions[i][1]));
        pcd.positions.push_back(static_cast<float>(cloud.positions[i][2]));
        pcd.point_count++;
    }
    return pcd;
}

// Per-orbital tracking for shell toggle + alpha fade-in
struct OrbitalEntry {
    int n;
    std::size_t pc_pos;
    std::size_t pc_neg;
    Vec4f target_pos;
    Vec4f target_neg;
    double spawn_time{-1.0};   // seconds since anim_start; -1 until first frame
    bool fade_done{false};
};

// Bohr model mesh indices
struct BohrEntry {
    std::size_t orbit_mesh_start;   // first orbit mesh index
    std::size_t orbit_mesh_count;
    std::size_t electron_mesh_start; // first electron mesh index
    std::size_t electron_mesh_count;
};

// View modes
enum ViewMode { VIEW_ORBITAL = 0, VIEW_BOHR = 1, VIEW_BOTH = 2 };

// Fast path: nucleus placeholder + Bohr model. No orbital sampling (cheap for any Z).
// Called synchronously on every element switch so the window never hangs.
static void rebuild_fast(App& app, int Z,
                         std::vector<OrbitalEntry>& orbital_entries,
                         BohrEntry& bohr_entry,
                         BohrModel<double>& bohr_model,
                         int& max_shell, int view_mode, double anim_t) {
    app.clear_meshes();
    app.clear_point_clouds();
    orbital_entries.clear();

    const auto& elem = element(Z);
    int max_n = 1;
    for (auto& sub : elem.electron_config())
        if (sub.n > max_n) max_n = sub.n;
    max_shell = max_n;

    double outer_bound = static_cast<double>(max_n * (max_n + 1)) + 3.0;
    double nucleus_r = std::max(outer_bound * 0.006, 0.04);
    auto nuc = mesh::uv_sphere_mesh<double>(8, 4, nucleus_r);
    app.add_mesh(nuc, Vec4f{0.8f, 0.8f, 0.8f, 1.0f});

    bohr_model = BohrModel<double>::build(Z);
    bohr_model.update_electrons(anim_t);  // prime spawn_time so initial meshes are scale~0

    bohr_entry.orbit_mesh_start = app.mesh_count();
    for (auto& orbit : bohr_model.orbits)
        app.add_mesh(orbit.ring, orbit.color);
    bohr_entry.orbit_mesh_count = bohr_model.orbits.size();

    bohr_entry.electron_mesh_start = app.mesh_count();
    for (auto& electron : bohr_model.electrons)
        app.add_mesh(electron.sphere, electron.color);
    bohr_entry.electron_mesh_count = bohr_model.electrons.size();

    float max_bound = std::max(0.3f, static_cast<float>(outer_bound));
    app.fit_camera(max_bound);

    bool show_bohr = (view_mode == VIEW_BOHR || view_mode == VIEW_BOTH);
    for (std::size_t i = 0; i < bohr_entry.orbit_mesh_count; ++i)
        app.set_mesh_visible(bohr_entry.orbit_mesh_start + i, show_bohr);
    for (std::size_t i = 0; i < bohr_entry.electron_mesh_count; ++i)
        app.set_mesh_visible(bohr_entry.electron_mesh_start + i, show_bohr);
}

// One orbital coming off the producer queue.
struct OrbitalResult {
    int gen;
    int n, l, m;
    OrbitalPointCloud<double> cloud;
    Vec4f color;
};

struct ProducerState {
    std::mutex mu;
    std::condition_variable cv;
    std::deque<int> jobs;              // element Z to build; latest wins
    std::deque<OrbitalResult> inbox;   // completed orbitals awaiting install
    std::atomic<int> generation{0};    // bumps on every new job; stale results are dropped
    bool quit{false};
};

// Install one streamed orbital as a point-cloud pair with alpha=0 (fades in via frame_callback).
static void install_one_orbital(App& app, const OrbitalResult& r,
                                std::vector<OrbitalEntry>& orbital_entries,
                                const bool shell_visible[10], int view_mode) {
    auto pos_cloud = cloud_to_points(r.cloud, true);
    auto neg_cloud = cloud_to_points(r.cloud, false);

    // Same vivid color for both lobes; the 0.5x neg-lobe variant washed everything to pastel.
    Vec4f color_pos = r.color;
    Vec4f color_neg = r.color;

    Vec4f cp_start = color_pos; cp_start[3] = 0.0f;
    Vec4f cn_start = color_neg; cn_start[3] = 0.0f;

    std::size_t idx_pos = app.point_cloud_count();
    app.add_point_cloud(std::move(pos_cloud), cp_start);
    std::size_t idx_neg = app.point_cloud_count();
    app.add_point_cloud(std::move(neg_cloud), cn_start);

    bool show_orbital = (view_mode == VIEW_ORBITAL || view_mode == VIEW_BOTH);
    bool vis = show_orbital && shell_visible[r.n];
    app.set_point_visible(idx_pos, vis);
    app.set_point_visible(idx_neg, vis);

    orbital_entries.push_back({r.n, idx_pos, idx_neg, color_pos, color_neg, -1.0, false});
}

// Background producer: pops jobs, samples orbitals one by one, pushes results.
// Aborts mid-atom as soon as a newer job arrives (generation bump).
static void producer_loop(std::shared_ptr<ProducerState> state, std::size_t num_points) {
    while (true) {
        int z = 0;
        int my_gen = 0;
        {
            std::unique_lock lk(state->mu);
            state->cv.wait(lk, [&]{ return state->quit || !state->jobs.empty(); });
            if (state->quit) return;
            z = state->jobs.back();
            state->jobs.clear();
            my_gen = state->generation.load();
        }

        const auto& e = element(z);
        auto cfg = e.electron_config();
        int total_subs = static_cast<int>(cfg.size());
        int sub_idx = 0;
        for (auto& sub : cfg) {
            if (state->generation.load() != my_gen) break;
            // One color per subshell (n,l) — all m variants share it so the full
            // orbital (e.g. 2p = three dumbbells) reads as one object.
            Vec4f sub_color = orbital_palette(sub_idx, total_subs);
            for (int m = -static_cast<int>(sub.l); m <= static_cast<int>(sub.l); ++m) {
                if (state->generation.load() != my_gen) break;

                auto cloud = sample_orbital_points<double>(sub.n, sub.l, m, num_points);
                if (cloud.positions.empty()) continue;
                if (state->generation.load() != my_gen) break;

                OrbitalResult r{my_gen, sub.n, static_cast<int>(sub.l), m,
                                std::move(cloud), sub_color};
                {
                    std::scoped_lock lk(state->mu);
                    if (state->generation.load() != my_gen) break;
                    state->inbox.push_back(std::move(r));
                }
            }
            ++sub_idx;
        }
    }
}

int main(int argc, char* argv[]) {
    int Z = 6;  // Carbon default
    std::size_t num_points = 100000;
    bool do_svg = false;
    bool do_viewer = true;
    bool force = false;
    std::string out_prefix;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") { print_usage(); return 0; }
        else if (arg == "--svg") do_svg = true;
        else if (arg == "--no-viewer") do_viewer = false;
        else if (arg == "--force") force = true;
        else if (arg == "--points" && i + 1 < argc) {
            num_points = std::stoul(argv[++i]);
        } else if (arg == "--out-prefix" && i + 1 < argc) {
            out_prefix = argv[++i];
        } else {
            try { Z = std::stoi(arg); }
            catch (...) {
                const auto& e = element(arg);
                Z = e.Z;
            }
        }
    }

    if (Z < 1 || Z > ELEMENT_COUNT) {
        std::println(stderr, "Error: Z must be 1-{}", ELEMENT_COUNT);
        return 1;
    }

    const auto& elem = element(Z);
    std::println("Element: {} ({}) Z={} mass={:.3f}", elem.name, elem.symbol, elem.Z, elem.atomic_mass);
    std::println("Config:  {}", electron_config_string(elem));

    std::println("\nBuilding atom model ({} points/orbital)...", num_points);
    auto t0 = std::chrono::steady_clock::now();
    auto model = AtomModel<>::build(Z, num_points);
    auto t1 = std::chrono::steady_clock::now();
    double build_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::println("Nucleus: {} vertices, {} faces",
        model.nucleus.vertex_count(), model.nucleus.face_count());

    std::size_t total_points = 0;
    for (auto& orb : model.orbitals) {
        std::size_t pos = 0, neg = 0;
        for (bool b : orb.cloud.positive_lobe) {
            if (b) ++pos; else ++neg;
        }
        total_points += orb.cloud.positions.size();
        std::println("  {}{}  m={:+d}: {:6} pts ({} +, {} -)",
            orb.n, subshell_letter(orb.l), orb.m,
            orb.cloud.positions.size(), pos, neg);
    }
    std::println("\n{} orbitals, {} total points, {:.0f} ms",
        model.orbitals.size(), total_points, build_ms);

    if (do_svg) {
        int max_n = 0;
        for (auto& orb : model.orbitals)
            if (orb.n > max_n) max_n = orb.n;

        for (auto& orb : model.orbitals) {
            if (orb.n != max_n) continue;
            auto svg = orbital_cross_section_svg(orb.n, orb.l, orb.m);
            std::string fname = out_prefix + std::format("{}_{}{}_{}.svg",
                elem.symbol, orb.n, subshell_letter(orb.l), orb.m);
            if (spatium::examples::confirm_overwrite(fname, force)) {
                svg.save(fname);
                std::println("  SVG: {}", fname);
            }
        }

        auto pt = periodic_table_svg();
        std::string pt_path = out_prefix + "periodic_table.svg";
        if (spatium::examples::confirm_overwrite(pt_path, force)) {
            pt.save(pt_path);
            std::println("  SVG: {}", pt_path);
        }
    }

    if (do_viewer) {
        std::string title = std::format("Spatium \u2014 {} ({}) Orbitals", elem.name, elem.symbol);
        App app(title, 1280, 720);

        std::vector<OrbitalEntry> orbital_entries;
        BohrEntry bohr_entry{};
        BohrModel<double> bohr_model;
        int max_shell = 0;
        int view_mode = VIEW_ORBITAL;

        bool shell_visible[10] = {true, true, true, true, true,
                                  true, true, true, true, true};

        auto anim_start = std::chrono::steady_clock::now();
        auto anim_now = [&] {
            return std::chrono::duration<double>(
                std::chrono::steady_clock::now() - anim_start).count();
        };

        // Persistent producer: streams orbital clouds one by one so points appear live.
        auto producer_state = std::make_shared<ProducerState>();
        std::thread producer(producer_loop, producer_state, num_points);

        auto kick_rebuild = [&](int z) {
            rebuild_fast(app, z, orbital_entries, bohr_entry, bohr_model, max_shell, view_mode,
                         anim_now());
            {
                std::scoped_lock lk(producer_state->mu);
                producer_state->generation.fetch_add(1);
                producer_state->inbox.clear();
                producer_state->jobs.clear();
                producer_state->jobs.push_back(z);
            }
            producer_state->cv.notify_one();
        };

        kick_rebuild(Z);

        // Apply view mode visibility helper
        auto apply_view_mode = [&]() {
            bool show_orbital = (view_mode == VIEW_ORBITAL || view_mode == VIEW_BOTH);
            bool show_bohr = (view_mode == VIEW_BOHR || view_mode == VIEW_BOTH);

            for (auto& oe : orbital_entries) {
                app.set_point_visible(oe.pc_pos, show_orbital && shell_visible[oe.n]);
                app.set_point_visible(oe.pc_neg, show_orbital && shell_visible[oe.n]);
            }
            for (std::size_t i = 0; i < bohr_entry.orbit_mesh_count; ++i)
                app.set_mesh_visible(bohr_entry.orbit_mesh_start + i, show_bohr);
            for (std::size_t i = 0; i < bohr_entry.electron_mesh_count; ++i)
                app.set_mesh_visible(bohr_entry.electron_mesh_start + i, show_bohr);
        };

#if SPATIUM_HAS_IMGUI
        app.enable_imgui();
        int current_z = Z;

        app.set_gui_callback([&]() {
            ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(300, 0), ImGuiCond_FirstUseEver);
            ImGui::Begin("Atom Explorer");

            // Element picker
            const auto& cur = element(current_z);
            std::string preview = std::format("{} {} (Z={})", cur.symbol, cur.name, cur.Z);
            if (ImGui::BeginCombo("Element", preview.c_str())) {
                for (int z = 1; z <= ELEMENT_COUNT; ++z) {
                    const auto& e = element(z);
                    bool selected = (z == current_z);
                    std::string label = std::format("{:3} {} {}", z, e.symbol, e.name);
                    if (ImGui::Selectable(label.c_str(), selected)) {
                        current_z = z;
                        for (auto& sv : shell_visible) sv = true;
                        kick_rebuild(current_z);
                    }
                    if (selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }

            // Electron config
            const auto& el = element(current_z);
            ImGui::TextWrapped("%s", electron_config_string(el).c_str());

            ImGui::Separator();

            // View mode
            ImGui::Text("View Mode");
            int prev_mode = view_mode;
            ImGui::RadioButton("Orbital", &view_mode, VIEW_ORBITAL); ImGui::SameLine();
            ImGui::RadioButton("Bohr", &view_mode, VIEW_BOHR); ImGui::SameLine();
            ImGui::RadioButton("Both", &view_mode, VIEW_BOTH);
            if (view_mode != prev_mode) apply_view_mode();

            ImGui::Separator();

            // Shell toggles (only relevant for orbital mode)
            if (view_mode == VIEW_ORBITAL || view_mode == VIEW_BOTH) {
                ImGui::Text("Shells");
                for (int n = 1; n <= max_shell; ++n) {
                    std::string label = std::format("n={}", n);
                    if (ImGui::Checkbox(label.c_str(), &shell_visible[n])) {
                        for (auto& oe : orbital_entries) {
                            if (oe.n == n) {
                                app.set_point_visible(oe.pc_pos, shell_visible[n]);
                                app.set_point_visible(oe.pc_neg, shell_visible[n]);
                            }
                        }
                    }
                    if (n < max_shell) ImGui::SameLine();
                }
                ImGui::Separator();
            }

            // Point size
            ImGui::SliderFloat("Point Size", &app.point_size_, 0.5f, 10.0f, "%.1f");

            // Quarter cutaway (world +X/+Z quadrant, no camera binding)
            ImGui::Checkbox("Quarter cutaway (C)", &app.cutaway_enabled_);

            // Screenshot
            if (ImGui::Button("Screenshot")) {
                auto now = std::chrono::system_clock::now();
                auto time = std::chrono::system_clock::to_time_t(now);
                char buf[64];
                std::strftime(buf, sizeof(buf), "screenshot_%Y%m%d_%H%M%S.png", std::localtime(&time));
                std::string shot_path = out_prefix + buf;
                if (spatium::examples::confirm_overwrite(shot_path, force))
                    app.save_screenshot(shot_path);
            }

            ImGui::End();
        });
#endif

        // Bohr electron animation + streamed orbital install + alpha fade-in
        app.set_frame_callback([&]() {
            double t = anim_now();

            // Drain any orbitals the producer has finished since last frame.
            std::deque<OrbitalResult> drained;
            {
                std::scoped_lock lk(producer_state->mu);
                drained.swap(producer_state->inbox);
            }
            int current_gen = producer_state->generation.load();
            for (auto& r : drained) {
                if (r.gen != current_gen) continue;  // stale from a superseded element
                install_one_orbital(app, r, orbital_entries, shell_visible, view_mode);
            }

            // Orbital alpha fade-in (smoothstep over 0.5s)
            constexpr double fade = 0.5;
            for (auto& oe : orbital_entries) {
                if (oe.fade_done) continue;
                if (oe.spawn_time < 0) oe.spawn_time = t;
                double u = std::clamp((t - oe.spawn_time) / fade, 0.0, 1.0);
                double s = u * u * (3.0 - 2.0 * u);
                Vec4f cp = oe.target_pos; cp[3] *= static_cast<float>(s);
                Vec4f cn = oe.target_neg; cn[3] *= static_cast<float>(s);
                app.set_point_color(oe.pc_pos, cp);
                app.set_point_color(oe.pc_neg, cn);
                if (u >= 1.0) oe.fade_done = true;
            }

            // Bohr electron motion + fade-in (update_electrons handles spawn scale internally).
            if (view_mode == VIEW_ORBITAL) return;
            bohr_model.update_electrons(t);
            for (std::size_t i = 0; i < bohr_model.electrons.size(); ++i) {
                std::size_t mesh_idx = bohr_entry.electron_mesh_start + i;
                app.update_mesh(mesh_idx, mesh_to_render_data(bohr_model.electrons[i].sphere));
            }
        });

        // Shell toggle: keys 1-9 (works with or without ImGui)
        app.set_key_callback([&](int key, int /*action*/, int /*mods*/) {
            if (key >= '1' && key <= '9') {
                int n = key - '0';
                shell_visible[n] = !shell_visible[n];
                for (auto& oe : orbital_entries) {
                    if (oe.n == n) {
                        app.set_point_visible(oe.pc_pos, shell_visible[n]);
                        app.set_point_visible(oe.pc_neg, shell_visible[n]);
                    }
                }
            }
        });

        std::println("\nControls: drag=orbit, right-drag=pan, scroll=zoom");
        std::println("  W=wireframe, C=cutaway, P/O=point size +/-, 1-9=toggle shells, Q=quit");
        app.run();

        // Stop producer: bump generation so current sample aborts, signal quit, join.
        {
            std::scoped_lock lk(producer_state->mu);
            producer_state->generation.fetch_add(1);
            producer_state->quit = true;
        }
        producer_state->cv.notify_all();
        producer.join();
    }

    return 0;
}
