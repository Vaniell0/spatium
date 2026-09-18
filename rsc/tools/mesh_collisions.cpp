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

#include "mesh_toy.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <print>
#include <set>
#include <string>
#include <vector>



int main() {
    const meshtoy::Space space{};
    const meshtoy::MeshT start = spatium::mesh::icosahedron(space);

    std::println("Do different operation sequences reach the same mesh?");
    std::println("start: icosahedron on Sphere<2>, {} vertices; {} operations",
                 start.vertex_count(), meshtoy::ops().size());
    std::print("ops:");
    for (const auto& o : meshtoy::ops()) std::print(" {}", o.name);
    std::println("\n");

    std::println("  {:>5} | {:>12} | {:>10} | {:>12} | {:>10} | {:>12}",
                 "depth", "sequences", "exact", "collapse", "rounded 1e-9", "collapse");

    std::vector<meshtoy::MeshT> frontier{start};
    std::size_t sequences = 1;

    for (std::size_t depth = 1; depth <= 6; ++depth) {
        std::vector<meshtoy::MeshT> next;
        std::size_t attempted = 0;
        for (const auto& m : frontier) {
            for (const auto& o : meshtoy::ops()) {
                ++attempted;
                meshtoy::MeshT r = o.apply(m, space);
                if (r.vertex_count() > meshtoy::kVertexCap) continue;
                next.push_back(std::move(r));
            }
        }
        sequences = attempted;

        std::set<std::string> exact, rounded;
        for (const auto& m : next) {
            exact.insert(meshtoy::fingerprint(m, -1));
            rounded.insert(meshtoy::fingerprint(m, 9));
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
        std::vector<meshtoy::MeshT> deduped;
        for (auto& m : next)
            if (seen.insert(meshtoy::fingerprint(m, -1)).second) deduped.push_back(std::move(m));
        frontier = std::move(deduped);
    }

    std::println("");
    std::println("Against numbers, for scale: five integer operations collapse 9 765 625");
    std::println("paths into 13 836 states by depth 10 -- 706x, effective branching 2.13");
    std::println("out of 5. The comparison that matters is the branching factor, not the");
    std::println("collapse, since the collapse compounds with depth.");
    return 0;
}
