// Do different sequences of mesh operations reach the same mesh?
//
// This measures the assumption the whole "move the search onto the
// substrate" plan rests on, and which was until now only reasoning. Over
// numbers, a chain collapses hard: 9 765 625 paths at depth 10 reduce to
// 13 836 distinct states, a 706x collapse and an effective branching
// factor of 2.13 rather than 5, because `10+3+3` and `2*8` arrive at the
// same place. A search over a space that small needs no heuristic and
// learns nothing. The claim was that geometries do not collide that way.
//
// The claim as stated is too strong, and the structure says so before any
// measurement does: translations and scales compose, so `translate(a)`
// then `translate(b)` *is* `translate(a+b)` and two sequences must land on
// the same mesh. Subdivision is the opposite -- not invertible, changes
// the vertex count -- so it breaks collisions rather than creating them.
// The honest question is therefore not whether geometries collide but how
// much, against the 706x that numbers manage.
//
// Two hashes are reported rather than one. Exact bit patterns count two
// meshes as different if they differ in the last ulp, which flatters the
// hypothesis; rounding to a tolerance can merge meshes that are genuinely
// distinct, which flatters the opposite. Neither is the honest number on
// its own, so both are printed and the gap between them is itself the
// result.
//
// Build: part of the rsc tools target.

#include <spatium/mesh/mesh.hpp>
#include <cstring>
#include <spatium/mesh/operations.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/spaces/sphere.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <print>
#include <set>
#include <string>
#include <vector>

namespace {

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

}  // namespace

int main() {
    const Space space{};
    const MeshT start = spatium::mesh::icosahedron(space);

    std::println("Do different operation sequences reach the same mesh?");
    std::println("start: icosahedron on Sphere<2>, {} vertices; {} operations",
                 start.vertex_count(), ops().size());
    std::print("ops:");
    for (const auto& o : ops()) std::print(" {}", o.name);
    std::println("\n");

    std::println("  {:>5} | {:>12} | {:>10} | {:>12} | {:>10} | {:>12}",
                 "depth", "sequences", "exact", "collapse", "rounded 1e-9", "collapse");

    std::vector<MeshT> frontier{start};
    std::size_t sequences = 1;

    for (std::size_t depth = 1; depth <= 6; ++depth) {
        std::vector<MeshT> next;
        std::size_t attempted = 0;
        for (const auto& m : frontier) {
            for (const auto& o : ops()) {
                ++attempted;
                MeshT r = o.apply(m, space);
                if (r.vertex_count() > kVertexCap) continue;
                next.push_back(std::move(r));
            }
        }
        sequences = attempted;

        std::set<std::string> exact, rounded;
        for (const auto& m : next) {
            exact.insert(fingerprint(m, -1));
            rounded.insert(fingerprint(m, 9));
        }

        std::println("  {:>5} | {:>12} | {:>10} | {:>11.1f}x | {:>12} | {:>11.1f}x",
                     depth, sequences, exact.size(),
                     static_cast<double>(sequences) / static_cast<double>(exact.size()),
                     rounded.size(),
                     static_cast<double>(sequences) / static_cast<double>(rounded.size()));

        // Deduplicate the frontier by exact identity, the same way the
        // integer search deduplicated states -- otherwise the count grows
        // as the number of paths rather than of reachable meshes.
        std::set<std::string> seen;
        std::vector<MeshT> deduped;
        for (auto& m : next)
            if (seen.insert(fingerprint(m, -1)).second) deduped.push_back(std::move(m));
        frontier = std::move(deduped);
    }

    std::println("");
    std::println("Against numbers, for scale: five integer operations collapse 9 765 625");
    std::println("paths into 13 836 states by depth 10 -- 706x, effective branching 2.13");
    std::println("out of 5. The comparison that matters is the branching factor, not the");
    std::println("collapse, since the collapse compounds with depth.");
    return 0;
}
