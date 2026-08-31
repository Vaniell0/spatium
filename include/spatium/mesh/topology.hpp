#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/mesh/subdivision.hpp>  // detail::PairHash
#  include <algorithm>
#  include <cstdint>
#  include <memory>
#  include <span>
#  include <unordered_map>
#  include <utility>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

struct Edge {
    uint32_t v0, v1;  // ordered: v0 < v1
    std::array<uint32_t, 2> adj_faces;  // UINT32_MAX = no face
};

inline constexpr uint32_t no_vertex = UINT32_MAX;
inline constexpr uint32_t no_face = UINT32_MAX;

template<Surface S>
class MeshTopology {
public:
    // Zero-copy: caller manages lifetime via shared_ptr
    static MeshTopology build(std::shared_ptr<const Mesh<S>> mesh) {
        MeshTopology topo(std::move(mesh));
        topo.build_internal();
        return topo;
    }

    // Convenience: copies mesh internally (safe but uses more memory)
    static MeshTopology build(const Mesh<S>& mesh) {
        return build(std::make_shared<const Mesh<S>>(mesh));
    }

    std::span<const uint32_t> neighbors(uint32_t v) const {
        return vertex_neighbors_[v];
    }

    const Edge& edge(uint32_t e) const { return edges_[e]; }
    uint32_t edge_count() const { return static_cast<uint32_t>(edges_.size()); }
    uint32_t vertex_count() const { return static_cast<uint32_t>(vertex_neighbors_.size()); }

    const std::array<uint32_t, 3>& face_edges(uint32_t f) const { return face_edges_[f]; }

    bool has_boundary() const {
        for (auto& e : edges_)
            if (e.adj_faces[1] == no_face) return true;
        return false;
    }

    const Mesh<S>& mesh() const { return *mesh_; }

private:
    explicit MeshTopology(std::shared_ptr<const Mesh<S>> m) : mesh_(std::move(m)) {}

    void build_internal() {
        using Key = std::pair<uint32_t, uint32_t>;
        std::unordered_map<Key, uint32_t, detail::PairHash> edge_map;

        auto nv = static_cast<uint32_t>(mesh_->vertices.size());
        vertex_neighbors_.resize(nv);
        face_edges_.resize(mesh_->faces.size());

        for (uint32_t fi = 0; fi < mesh_->faces.size(); ++fi) {
            auto [a, b, c] = mesh_->faces[fi];
            std::array<uint32_t, 3> verts = {a, b, c};

            for (int k = 0; k < 3; ++k) {
                uint32_t u = verts[k], w = verts[(k + 1) % 3];
                auto key = std::minmax(u, w);

                auto [it, inserted] = edge_map.try_emplace(
                    key, static_cast<uint32_t>(edges_.size()));

                if (inserted) {
                    edges_.push_back({key.first, key.second, {fi, no_face}});
                    vertex_neighbors_[u].push_back(w);
                    vertex_neighbors_[w].push_back(u);
                } else {
                    edges_[it->second].adj_faces[1] = fi;
                }

                face_edges_[fi][k] = it->second;
            }
        }
    }

    std::shared_ptr<const Mesh<S>> mesh_;
    std::vector<Edge> edges_;
    std::vector<std::vector<uint32_t>> vertex_neighbors_;
    std::vector<std::array<uint32_t, 3>> face_edges_;
};

} // namespace spatium::mesh
