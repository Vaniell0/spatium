#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/physics/atomic/atom_palette.hpp>
#  include <spatium/physics/atomic/orbital.hpp>
#  include <spatium/physics/elements.hpp>
#  include <spatium/mesh/primitives.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <algorithm>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::atomic {

template<Scalar T = double>
struct AtomModel {
    int Z;

    struct OrbitalInfo {
        int n, l, m;
        OrbitalPointCloud<T> cloud;
        Vec4f color;
    };

    mesh::Mesh<Euclidean<3, T>> nucleus;
    std::vector<OrbitalInfo> orbitals;

    // Default `num_points = 10000` keeps the call site cheap. Earlier
    // 100 000 was chosen for the standalone atom_demo screenshot but
    // is needlessly large for tests / smoke runs (≈ 100 MB of cloud
    // points across an Fe atom). Demos that want the dense version
    // pass an explicit larger argument.
    static AtomModel build(int Z, std::size_t num_points = 10000) {
        AtomModel model;
        model.Z = Z;

        // Build orbitals from electron config
        const auto& elem = element(Z);

        // Nucleus: proportional to outermost orbital extent
        int max_n = 1;
        for (auto& sub : elem.electron_config())
            if (sub.n > max_n) max_n = sub.n;
        T outer_bound = static_cast<T>(max_n * (max_n + 1)) + T{3};
        T nucleus_r = std::max(outer_bound * T{0.006}, T{0.04});
        model.nucleus = mesh::uv_sphere_mesh<T>(8, 4, nucleus_r);

        auto cfg = elem.electron_config();
        int total_subs = static_cast<int>(cfg.size());
        int sub_idx = 0;
        for (auto& sub : cfg) {
            Vec4f sub_color = orbital_palette(sub_idx, total_subs);
            for (int m = -static_cast<int>(sub.l); m <= static_cast<int>(sub.l); ++m) {
                auto cloud = sample_orbital_points<T>(sub.n, sub.l, m, num_points);

                if (cloud.positions.empty()) continue;

                model.orbitals.push_back({
                    sub.n, static_cast<int>(sub.l), m,
                    std::move(cloud),
                    sub_color
                });
            }
            ++sub_idx;
        }

        return model;
    }
};

} // namespace spatium::physics::atomic
