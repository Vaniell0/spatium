#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/spaces/euclidean.hpp>
#  include <array>
#  include <charconv>
#  include <cstdint>
#  include <filesystem>
#  include <fstream>
#  include <sstream>
#  include <string>
#endif

SPATIUM_EXPORT namespace spatium::io {

using E3 = Euclidean<3>;

// ── Load OBJ ──────────────────────────────────────────────────

inline Result<mesh::Mesh<E3>> load_obj(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file)
        return std::unexpected(Error{ErrorCode::InvalidArgument, "cannot open file"});

    mesh::Mesh<E3> m;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream ss(line);
        std::string prefix;
        ss >> prefix;

        if (prefix == "v") {
            double x, y, z;
            if (!(ss >> x >> y >> z))
                return std::unexpected(Error{ErrorCode::InvalidArgument, "bad vertex line"});
            m.vertices.push_back(Vec3{x, y, z});
        } else if (prefix == "f") {
            // Support: "f 1 2 3", "f 1/2 3/4 5/6", "f 1/2/3 4/5/6 7/8/9", "f 1//3 4//6 7//9"
            std::vector<uint32_t> indices;
            std::string token;
            while (ss >> token) {
                auto slash = token.find('/');
                std::string_view idx_sv = token;
                if (slash != std::string::npos) idx_sv = std::string_view(token).substr(0, slash);
                int idx = 0;
                auto [ptr, ec] = std::from_chars(idx_sv.data(), idx_sv.data() + idx_sv.size(), idx);
                if (ec != std::errc{})
                    return std::unexpected(Error{ErrorCode::ParseError, "bad face index: " + token});
                // OBJ is 1-based; negative = relative
                if (idx > 0) indices.push_back(static_cast<uint32_t>(idx - 1));
                else if (idx < 0) indices.push_back(static_cast<uint32_t>(m.vertices.size() + idx));
            }
            // Triangulate face fan for polygons with >3 vertices
            for (std::size_t i = 1; i + 1 < indices.size(); ++i)
                m.faces.push_back({indices[0], indices[i], indices[i + 1]});
        }
    }

    if (m.vertices.empty())
        return std::unexpected(Error{ErrorCode::InvalidArgument, "no vertices found"});

    return m;
}

// ── Save OBJ ──────────────────────────────────────────────────

inline Result<void> save_obj(const mesh::Mesh<E3>& m, const std::filesystem::path& path) {
    std::ofstream file(path);
    if (!file)
        return std::unexpected(Error{ErrorCode::InvalidArgument, "cannot open file for writing"});

    file << "# Spatium OBJ export\n";
    file << "# Vertices: " << m.vertices.size() << " Faces: " << m.faces.size() << "\n\n";

    file << std::fixed;
    for (auto& v : m.vertices)
        file << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";

    file << "\n";
    for (auto& [a, b, c] : m.faces)
        file << "f " << (a + 1) << " " << (b + 1) << " " << (c + 1) << "\n";

    return {};
}

} // namespace spatium::io
