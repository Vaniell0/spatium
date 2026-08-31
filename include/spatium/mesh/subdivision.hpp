#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/mesh.hpp>
#  include <unordered_map>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

// internal — do not use, no API stability. PairHash is the
// midpoint-edge cache key hasher used by subdivide_once(); it may
// be replaced with a different hashing strategy.
namespace detail {
struct PairHash {
    std::size_t operator()(std::pair<uint32_t, uint32_t> p) const {
        return std::hash<uint64_t>{}(
            (static_cast<uint64_t>(p.first) << 32) | p.second);
    }
};
} // namespace detail

template<Surface S>
Mesh<S> subdivide_once(const Mesh<S>& m, const S& space) {
    Mesh<S> result;
    result.vertices = m.vertices;

    std::unordered_map<std::pair<uint32_t, uint32_t>, uint32_t, detail::PairHash> midpoint_cache;

    auto get_midpoint = [&](uint32_t a, uint32_t b) -> uint32_t {
        auto key = std::minmax(a, b);
        if (auto it = midpoint_cache.find(key); it != midpoint_cache.end())
            return it->second;

        auto mid_ambient = (result.vertices[a] + result.vertices[b])
                           * typename S::ScalarType{0.5};
        auto mid_on_surface = space.project(mid_ambient);

        auto idx = static_cast<uint32_t>(result.vertices.size());
        result.vertices.push_back(mid_on_surface);
        midpoint_cache[key] = idx;
        return idx;
    };

    result.faces.reserve(m.faces.size() * 4);

    for (const auto& [i0, i1, i2] : m.faces) {
        uint32_t m01 = get_midpoint(i0, i1);
        uint32_t m12 = get_midpoint(i1, i2);
        uint32_t m20 = get_midpoint(i2, i0);
        result.faces.push_back({i0, m01, m20});
        result.faces.push_back({m01, i1, m12});
        result.faces.push_back({m20, m12, i2});
        result.faces.push_back({m01, m12, m20});
    }

    return result;
}

template<Surface S>
Mesh<S> subdivide(const Mesh<S>& m, const S& space, std::size_t levels) {
    auto result = m;
    for (std::size_t i = 0; i < levels; ++i)
        result = subdivide_once(result, space);
    return result;
}

// ── Adaptive subdivision by normal angle ──────────────────────
// Only subdivides faces where the angle between vertex normals exceeds
// the threshold. Handles T-junctions by splitting neighbor faces
// that share a subdivided edge.

template<Surface S>
Mesh<S> subdivide_adaptive(const Mesh<S>& m, const S& space,
                            typename S::ScalarType angle_threshold,
                            std::size_t max_levels = 4)
{
    using T = typename S::ScalarType;
    auto result = m;

    for (std::size_t level = 0; level < max_levels; ++level) {
        // Compute vertex normals
        std::vector<typename S::TangentVector> normals(result.vertices.size());
        for (std::size_t i = 0; i < result.vertices.size(); ++i)
            normals[i] = space.normal(result.vertices[i]);

        // Mark faces that need subdivision
        std::vector<bool> needs_refine(result.faces.size(), false);
        bool any = false;
        for (std::size_t fi = 0; fi < result.faces.size(); ++fi) {
            auto [a, b, c] = result.faces[fi];
            auto angle = [](const auto& n1, const auto& n2) {
                auto d = n1.dot(n2) / (n1.norm() * n2.norm());
                using std::acos;
                return acos(std::clamp(d, decltype(d){-1}, decltype(d){1}));
            };
            T max_a = std::max({angle(normals[a], normals[b]),
                                angle(normals[b], normals[c]),
                                angle(normals[a], normals[c])});
            if (max_a > angle_threshold) {
                needs_refine[fi] = true;
                any = true;
            }
        }

        if (!any) break;

        // Midpoint cache (shared across all faces for this level)
        std::unordered_map<std::pair<uint32_t, uint32_t>, uint32_t, detail::PairHash> midpoint_cache;
        Mesh<S> next;
        next.vertices = result.vertices;

        auto get_midpoint = [&](uint32_t a, uint32_t b) -> uint32_t {
            auto key = std::minmax(a, b);
            if (auto it = midpoint_cache.find(key); it != midpoint_cache.end())
                return it->second;
            auto mid = (next.vertices[a] + next.vertices[b]) * T{0.5};
            auto idx = static_cast<uint32_t>(next.vertices.size());
            next.vertices.push_back(space.project(mid));
            midpoint_cache[key] = idx;
            return idx;
        };

        // First pass: subdivide marked faces, track which edges are split
        for (std::size_t fi = 0; fi < result.faces.size(); ++fi) {
            if (needs_refine[fi]) {
                auto [i0, i1, i2] = result.faces[fi];
                uint32_t m01 = get_midpoint(i0, i1);
                uint32_t m12 = get_midpoint(i1, i2);
                uint32_t m20 = get_midpoint(i2, i0);
                next.faces.push_back({i0, m01, m20});
                next.faces.push_back({m01, i1, m12});
                next.faces.push_back({m20, m12, i2});
                next.faces.push_back({m01, m12, m20});
            }
        }

        // Second pass: handle unrefined faces (T-junction fix)
        for (std::size_t fi = 0; fi < result.faces.size(); ++fi) {
            if (needs_refine[fi]) continue;

            auto [a, b, c] = result.faces[fi];
            auto e01 = midpoint_cache.find(std::minmax(a, b));
            auto e12 = midpoint_cache.find(std::minmax(b, c));
            auto e20 = midpoint_cache.find(std::minmax(c, a));

            int splits = (e01 != midpoint_cache.end()) +
                         (e12 != midpoint_cache.end()) +
                         (e20 != midpoint_cache.end());

            if (splits == 0) {
                next.faces.push_back({a, b, c});
            } else if (splits == 3) {
                // All edges split → full subdivision
                uint32_t m01 = e01->second, m12 = e12->second, m20 = e20->second;
                next.faces.push_back({a, m01, m20});
                next.faces.push_back({m01, b, m12});
                next.faces.push_back({m20, m12, c});
                next.faces.push_back({m01, m12, m20});
            } else {
                // 1 or 2 splits → green refinement (split into 2 or 3 triangles)
                if (e01 != midpoint_cache.end() && splits == 1) {
                    uint32_t m = e01->second;
                    next.faces.push_back({a, m, c});
                    next.faces.push_back({m, b, c});
                } else if (e12 != midpoint_cache.end() && splits == 1) {
                    uint32_t m = e12->second;
                    next.faces.push_back({a, b, m});
                    next.faces.push_back({a, m, c});
                } else if (e20 != midpoint_cache.end() && splits == 1) {
                    uint32_t m = e20->second;
                    next.faces.push_back({a, b, m});
                    next.faces.push_back({b, c, m});
                } else {
                    // 2 splits → 3 triangles (find the unsplit edge)
                    if (e01 == midpoint_cache.end()) {
                        uint32_t m12 = e12->second, m20 = e20->second;
                        next.faces.push_back({a, b, m12});
                        next.faces.push_back({a, m12, m20});
                        next.faces.push_back({m20, m12, c});
                    } else if (e12 == midpoint_cache.end()) {
                        uint32_t m01 = e01->second, m20 = e20->second;
                        next.faces.push_back({a, m01, m20});
                        next.faces.push_back({m01, b, c});
                        next.faces.push_back({m01, c, m20});
                    } else {
                        uint32_t m01 = e01->second, m12 = e12->second;
                        next.faces.push_back({a, m01, c});
                        next.faces.push_back({m01, b, m12});
                        next.faces.push_back({m01, m12, c});
                    }
                }
            }
        }

        result = std::move(next);
    }

    return result;
}

} // namespace spatium::mesh
