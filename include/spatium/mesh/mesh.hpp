#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/mesh/face_metrics.hpp>
#  include <array>
#  include <cmath>
#  include <cstdint>
#  include <format>
#  include <ostream>
#  include <ranges>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

template<Surface S>
struct Mesh {
    using PointT = typename S::PointType;
    using ScalarT = typename S::ScalarType;

    std::vector<PointT> vertices;
    std::vector<std::array<uint32_t, 3>> faces; // triangle indices

    [[nodiscard]] std::size_t vertex_count() const { return vertices.size(); }
    [[nodiscard]] std::size_t face_count() const { return faces.size(); }

    // Edge count via Euler: E = F * 3 / 2 (closed mesh)
    [[nodiscard]] std::size_t edge_count_approx() const { return faces.size() * 3 / 2; }

    // Get triangle at face index
    [[nodiscard]] std::array<PointT, 3> triangle(std::size_t i) const {
        auto [a, b, c] = faces[i];
        return {vertices[a], vertices[b], vertices[c]};
    }

    // Lazy range over face triangles: for (auto [a, b, c] : mesh.triangles()) { ... }
    [[nodiscard]] auto triangles() const {
        return faces | std::views::transform([this](const auto& f) {
            return std::array<PointT, 3>{vertices[f[0]], vertices[f[1]], vertices[f[2]]};
        });
    }

    // Total surface area
    ScalarT area(const S& /*space*/) const {
        ScalarT total{0};
        for (const auto& f : faces)
            total += face_area_gram(vertices[f[0]], vertices[f[1]], vertices[f[2]]);
        return total;
    }

    // Centroid (average of all vertices)
    [[nodiscard]] PointT centroid() const {
        PointT sum{};
        for (const auto& v : vertices) sum = PointT{sum + v};
        if (!vertices.empty())
            return PointT{sum / static_cast<ScalarT>(vertices.size())};
        return sum;
    }

    // Axis-aligned bounding box (returns {min_corner, max_corner})
    [[nodiscard]] std::pair<PointT, PointT> bounding_box() const {
        if (vertices.empty()) return {{}, {}};
        PointT lo = vertices[0], hi = vertices[0];
        for (std::size_t i = 1; i < vertices.size(); ++i) {
            for (std::size_t d = 0; d < PointT::size; ++d) {
                if (vertices[i][d] < lo[d]) lo[d] = vertices[i][d];
                if (vertices[i][d] > hi[d]) hi[d] = vertices[i][d];
            }
        }
        return {lo, hi};
    }

    // Stats summary
    friend std::ostream& operator<<(std::ostream& os, const Mesh& m) {
        return os << "Mesh{V=" << m.vertex_count()
                  << " F=" << m.face_count()
                  << " E\u2248" << m.edge_count_approx() << '}';
    }
};

} // namespace spatium::mesh

template<spatium::Surface S>
struct std::formatter<spatium::mesh::Mesh<S>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::mesh::Mesh<S>& m, auto& ctx) const {
        return std::format_to(ctx.out(), "Mesh{{V={} F={} E\u2248{}}}",
                              m.vertex_count(), m.face_count(), m.edge_count_approx());
    }
};
