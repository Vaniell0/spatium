#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/spaces/euclidean.hpp>
#  include <array>
#  include <cstdint>
#  include <cstring>
#  include <filesystem>
#  include <fstream>
#  include <sstream>
#  include <string>
#  include <unordered_map>
#endif

SPATIUM_EXPORT namespace spatium::io {

using E3 = Euclidean<3>;

// ── Load STL (auto-detect binary vs ASCII) ────────────────────

inline Result<mesh::Mesh<E3>> load_stl(const std::filesystem::path& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file)
        return std::unexpected(Error(ErrorCode::InvalidArgument, "cannot open file"));

    // Peek: ASCII STL starts with "solid"
    char header[6]{};
    file.read(header, 5);
    file.seekg(0);

    bool is_ascii = (std::strncmp(header, "solid", 5) == 0);

    mesh::Mesh<E3> m;

    if (is_ascii) {
        std::string line;
        Vec3 normal{}, verts[3];
        int vi = 0;

        while (std::getline(file, line)) {
            std::istringstream ss(line);
            std::string word;
            ss >> word;

            if (word == "vertex") {
                double x, y, z;
                ss >> x >> y >> z;
                verts[vi++] = Vec3{x, y, z};
                if (vi == 3) {
                    auto base = static_cast<uint32_t>(m.vertices.size());
                    m.vertices.push_back(verts[0]);
                    m.vertices.push_back(verts[1]);
                    m.vertices.push_back(verts[2]);
                    m.faces.push_back({base, base + 1, base + 2});
                    vi = 0;
                }
            }
        }
    } else {
        // Binary STL: 80-byte header, uint32 triangle count, then triangles
        char hdr[80];
        file.read(hdr, 80);

        uint32_t num_tris = 0;
        file.read(reinterpret_cast<char*>(&num_tris), 4);

        for (uint32_t t = 0; t < num_tris; ++t) {
            float data[12];  // normal(3) + v0(3) + v1(3) + v2(3)
            file.read(reinterpret_cast<char*>(data), 48);
            uint16_t attr;
            file.read(reinterpret_cast<char*>(&attr), 2);

            auto base = static_cast<uint32_t>(m.vertices.size());
            for (int i = 0; i < 3; ++i)
                m.vertices.push_back(Vec3{
                    static_cast<double>(data[3 + i * 3]),
                    static_cast<double>(data[4 + i * 3]),
                    static_cast<double>(data[5 + i * 3])
                });
            m.faces.push_back({base, base + 1, base + 2});
        }
    }

    if (m.vertices.empty())
        return std::unexpected(Error(ErrorCode::InvalidArgument, "no triangles found"));

    return m;
}

// ── Save STL (binary) ─────────────────────────────────────────

inline Result<void> save_stl(const mesh::Mesh<E3>& m, const std::filesystem::path& path) {
    std::ofstream file(path, std::ios::binary);
    if (!file)
        return std::unexpected(Error(ErrorCode::InvalidArgument, "cannot open file for writing"));

    // 80-byte header
    char header[80]{};
    std::snprintf(header, 80, "Spatium STL export - %zu faces", m.faces.size());
    file.write(header, 80);

    auto num_tris = static_cast<uint32_t>(m.faces.size());
    file.write(reinterpret_cast<const char*>(&num_tris), 4);

    for (auto& [a, b, c] : m.faces) {
        auto e1 = m.vertices[b] - m.vertices[a];
        auto e2 = m.vertices[c] - m.vertices[a];
        auto n = Vec3{e1.cross(e2)}.normalized();

        float data[12] = {
            static_cast<float>(n[0]), static_cast<float>(n[1]), static_cast<float>(n[2]),
            static_cast<float>(m.vertices[a][0]), static_cast<float>(m.vertices[a][1]), static_cast<float>(m.vertices[a][2]),
            static_cast<float>(m.vertices[b][0]), static_cast<float>(m.vertices[b][1]), static_cast<float>(m.vertices[b][2]),
            static_cast<float>(m.vertices[c][0]), static_cast<float>(m.vertices[c][1]), static_cast<float>(m.vertices[c][2]),
        };
        file.write(reinterpret_cast<const char*>(data), 48);
        uint16_t attr = 0;
        file.write(reinterpret_cast<const char*>(&attr), 2);
    }

    return {};
}

} // namespace spatium::io
