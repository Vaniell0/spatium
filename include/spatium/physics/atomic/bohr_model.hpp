#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/physics/atomic/atom_palette.hpp>
#  include <spatium/physics/elements.hpp>
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/mesh/primitives.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <algorithm>
#  include <cmath>
#  include <numbers>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::atomic {

// Bohr model: circular orbits at r_n with electrons on the ring.
// Each shell n gets n electrons distributed evenly (simplified).
template<Scalar T = double>
struct BohrModel {
    struct Orbit {
        int n;
        T radius;
        mesh::Mesh<Euclidean<3, T>> ring;  // thin torus mesh
        Vec4f color;
    };

    struct Electron {
        int n;
        T orbit_radius;
        T speed;        // radians per second
        T phase;        // initial angle
        T base_radius;  // full sphere radius at scale=1
        T spawn_time{-1};  // set on first update_electrons call
        mesh::Mesh<Euclidean<3, T>> sphere;
        Vec4f color;
    };

    std::vector<Orbit> orbits;
    std::vector<Electron> electrons;

    // Build Bohr model for element Z.
    // Orbit radius = n^2 * scale (Bohr radius units).
    static BohrModel build(int Z, T scale = T{1}) {
        BohrModel model;
        const auto& elem = element(Z);
        auto config = elem.electron_config();

        // Determine max n and electrons per shell
        int max_n = 0;
        int shell_electrons[8] = {};
        for (auto& sub : config) {
            if (sub.n > max_n) max_n = sub.n;
            shell_electrons[sub.n] += sub.electrons;
        }

        constexpr T pi = std::numbers::pi_v<T>;

        for (int n = 1; n <= max_n; ++n) {
            if (shell_electrons[n] == 0) continue;

            T orbit_r = static_cast<T>(n * n) * scale;
            // Per-orbit scaling with clamps: inner shells stay tiny (clear of the nucleus),
            // outer shells don't blow up into fat tubes.
            T tube_r = std::clamp(orbit_r * T{0.015}, T{0.02}, T{0.15});
            T electron_r = std::clamp(tube_r * T{1.8}, T{0.06}, T{0.3});

            // Build torus ring mesh (Euclidean<3>)
            constexpr std::size_t ring_segs = 64;
            constexpr std::size_t tube_segs = 8;
            mesh::Mesh<Euclidean<3, T>> ring;

            for (std::size_t i = 0; i < ring_segs; ++i) {
                T u = T{2} * pi * static_cast<T>(i) / static_cast<T>(ring_segs);
                T cu = std::cos(u), su = std::sin(u);
                for (std::size_t j = 0; j < tube_segs; ++j) {
                    T v = T{2} * pi * static_cast<T>(j) / static_cast<T>(tube_segs);
                    T cv = std::cos(v), sv = std::sin(v);
                    ring.vertices.push_back(Vec<T, 3>{
                        (orbit_r + tube_r * cv) * cu,
                        tube_r * sv,
                        (orbit_r + tube_r * cv) * su
                    });
                }
            }

            // Faces (connect adjacent ring segments)
            for (std::size_t i = 0; i < ring_segs; ++i) {
                std::size_t ni = (i + 1) % ring_segs;
                for (std::size_t j = 0; j < tube_segs; ++j) {
                    std::size_t nj = (j + 1) % tube_segs;
                    auto a = static_cast<uint32_t>(i * tube_segs + j);
                    auto b = static_cast<uint32_t>(ni * tube_segs + j);
                    auto c = static_cast<uint32_t>(ni * tube_segs + nj);
                    auto d = static_cast<uint32_t>(i * tube_segs + nj);
                    ring.faces.push_back({a, b, c});
                    ring.faces.push_back({a, c, d});
                }
            }

            // Ring color: one vivid hue per shell from the palette.
            Vec4f color = orbital_palette(n - 1, max_n, 0.7f);

            model.orbits.push_back({n, orbit_r, std::move(ring), color});

            // Place electrons evenly around the orbit (electron_r computed once above)
            int ne = shell_electrons[n];
            T base_speed = T{2} * pi / (static_cast<T>(n) * T{1.5});  // slower for outer shells

            for (int k = 0; k < ne; ++k) {
                T phase = T{2} * pi * static_cast<T>(k) / static_cast<T>(ne);
                auto sphere = mesh::uv_sphere_mesh<T>(8, 4, electron_r);

                // Position sphere at initial angle
                T cx = orbit_r * std::cos(phase);
                T cz = orbit_r * std::sin(phase);
                for (auto& vert : sphere.vertices) {
                    vert[0] += cx;
                    vert[2] += cz;
                }

                model.electrons.push_back({n, orbit_r, base_speed, phase, electron_r, T{-1},
                                           std::move(sphere), {1.0f, 1.0f, 0.2f, 1.0f}});
            }
        }

        return model;
    }

    // Update electron positions for animation time t.
    // New electrons fade in (radius 0 -> base_radius) over fade_duration seconds.
    void update_electrons(T t, T fade_duration = T{0.45}) {
        for (auto& e : electrons) {
            if (e.spawn_time < T{0}) e.spawn_time = t;
            T age = t - e.spawn_time;
            T u = std::clamp(age / fade_duration, T{0}, T{1});
            T scale = u * u * (T{3} - T{2} * u);  // smoothstep
            T radius = e.base_radius * scale;
            if (radius <= T{0}) radius = e.base_radius * T{1e-3};  // keep mesh non-degenerate

            T angle = e.speed * t + e.phase;
            T cx = e.orbit_radius * std::cos(angle);
            T cz = e.orbit_radius * std::sin(angle);

            auto fresh = mesh::uv_sphere_mesh<T>(8, 4, radius);
            for (auto& v : fresh.vertices) {
                v[0] += cx;
                v[2] += cz;
            }
            e.sphere = std::move(fresh);
        }
    }
};

} // namespace spatium::physics::atomic
