#pragma once

// The mesh action space the substrate experiments share.
//
// Six operations over a mesh on Sphere<2>, chosen so the set is not
// degenerate: subdivision grows and cannot be undone, the transforms
// compose (which is why sequences collide at all), and `flip` is an
// involution. A set of only-growing operations would have no collisions by
// construction and would prove nothing.
//
// `fingerprint` is how two meshes are compared. Vertices are sorted first,
// so a mesh differing only in vertex order is the same mesh -- relabelling
// is not discovery.

#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/operations.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/spaces/sphere.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

namespace meshtoy {


using Space = spatium::Sphere<2, double>;
using MeshT = spatium::mesh::Mesh<Space>;

constexpr std::size_t kVertexCap = 5000;   // an op that exceeds this is skipped

struct Op {
    const char* name;
    MeshT (*apply)(const MeshT&, const Space&);
};

MeshT op_subdivide(const MeshT& m, const Space& s) {
    return spatium::mesh::subdivide_once(m, s);
}
MeshT op_scale_up(const MeshT& m, const Space&) {
    return spatium::mesh::transform(m, [](const auto& p) { return decltype(p){p * 1.5}; });
}
MeshT op_scale_down(const MeshT& m, const Space&) {
    return spatium::mesh::transform(m, [](const auto& p) { return decltype(p){p * (1.0 / 1.5)}; });
}
MeshT op_shift_x(const MeshT& m, const Space&) {
    return spatium::mesh::transform(m, [](const auto& p) {
        auto q = p; q[0] += 0.25; return q;
    });
}
MeshT op_shift_y(const MeshT& m, const Space&) {
    return spatium::mesh::transform(m, [](const auto& p) {
        auto q = p; q[1] += 0.25; return q;
    });
}
MeshT op_flip(const MeshT& m, const Space&) {
    return spatium::mesh::flip_normals(m);
}

const std::vector<Op>& ops() {
    static const std::vector<Op> v{
        {"subdivide", op_subdivide}, {"scale_up", op_scale_up},
        {"scale_down", op_scale_down}, {"shift_x", op_shift_x},
        {"shift_y", op_shift_y}, {"flip", op_flip},
    };
    return v;
}

// A mesh's identity, as a string so two of them can be compared exactly.
// `digits` < 0 means the full bit pattern; otherwise coordinates are
// rounded first, which is what decides whether near-misses count as the
// same mesh.
std::string fingerprint(const MeshT& m, int digits) {
    std::string out = std::to_string(m.vertex_count()) + ":" + std::to_string(m.face_count()) + "|";
    std::vector<std::string> coords;
    coords.reserve(m.vertex_count());
    for (const auto& v : m.vertices) {
        std::string s;
        for (int k = 0; k < 3; ++k) {
            double x = v[static_cast<std::size_t>(k)];
            if (digits >= 0) {
                const double scale = std::pow(10.0, digits);
                x = std::round(x * scale) / scale;
                if (x == 0.0) x = 0.0;   // fold -0.0 into 0.0
            }
            std::uint64_t bits = 0;
            std::memcpy(&bits, &x, sizeof(bits));
            s += std::to_string(bits) + ",";
        }
        coords.push_back(std::move(s));
    }
    // Sorted, so a mesh that differs only in vertex ordering is the same
    // mesh -- otherwise relabelling would be counted as discovery.
    std::sort(coords.begin(), coords.end());
    for (const auto& c : coords) out += c;
    return out;
}


}  // namespace meshtoy
