#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/mesh/mesh.hpp>
#  include <algorithm>
#  include <cmath>
#  include <cstdint>
#  include <format>
#  include <fstream>
#  include <functional>
#  include <span>
#  include <string>
#  include <utility>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io {

namespace svg_detail {
inline std::string xml_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            default: out += c;
        }
    }
    return out;
}
}  // namespace svg_detail

struct Svg {
    double width = 800;
    double height = 800;
    double scale = 100;
    double cx, cy;
    std::string content;

    inline Svg(double w = 800, double h = 800, double s = 100)
        : width(w), height(h), scale(s), cx(w / 2), cy(h / 2) {}

    inline double sx(double x) const { return cx + x * scale; }
    inline double sy(double y) const { return cy - y * scale; }

    inline Svg& line(double x1, double y1, double x2, double y2,
                     std::string stroke = "#cdd6f4", double w = 1.0) {
        content += std::format(
            R"(<line x1="{:.2f}" y1="{:.2f}" x2="{:.2f}" y2="{:.2f}" )"
            R"(stroke="{}" stroke-width="{:.1f}"/>)" "\n",
            sx(x1), sy(y1), sx(x2), sy(y2), stroke, w);
        return *this;
    }

    inline Svg& circle(double cx_w, double cy_w, double r,
                       std::string stroke = "#89b4fa", double w = 1.5) {
        content += std::format(
            R"(<circle cx="{:.2f}" cy="{:.2f}" r="{:.2f}" )"
            R"(fill="none" stroke="{}" stroke-width="{:.1f}"/>)" "\n",
            sx(cx_w), sy(cy_w), r * scale, stroke, w);
        return *this;
    }

    inline Svg& triangle(double x1, double y1, double x2, double y2,
                         double x3, double y3,
                         std::string fill = "#31324466",
                         std::string stroke = "#cdd6f4") {
        content += std::format(
            R"(<polygon points="{:.2f},{:.2f} {:.2f},{:.2f} {:.2f},{:.2f}" )"
            R"(fill="{}" stroke="{}" stroke-width="0.5"/>)" "\n",
            sx(x1), sy(y1), sx(x2), sy(y2), sx(x3), sy(y3), fill, stroke);
        return *this;
    }

    inline Svg& point(double x, double y, double r = 3,
                      std::string fill = "#f38ba8") {
        content += std::format(
            R"(<circle cx="{:.2f}" cy="{:.2f}" r="{:.1f}" fill="{}"/>)" "\n",
            sx(x), sy(y), r, fill);
        return *this;
    }

    inline Svg& text(double x, double y, const std::string& label,
                     std::string fill = "#cdd6f4", double size = 12) {
        content += std::format(
            R"(<text x="{:.2f}" y="{:.2f}" fill="{}" font-size="{:.0f}" )"
            R"(font-family="monospace">{}</text>)" "\n",
            sx(x), sy(y) - 5, fill, size, svg_detail::xml_escape(label));
        return *this;
    }

    inline Svg& polygon(const std::vector<std::pair<double, double>>& pts,
                        std::string fill = "#89b4fa44",
                        std::string stroke = "#89b4fa", double w = 1.5) {
        if (pts.empty()) return *this;
        std::string points_str;
        for (const auto& [x, y] : pts)
            points_str += std::format("{:.2f},{:.2f} ", sx(x), sy(y));
        content += std::format(
            R"(<polygon points="{}" fill="{}" stroke="{}" stroke-width="{:.1f}"/>)" "\n",
            points_str, fill, stroke, w);
        return *this;
    }

    // Template: stays in header
    template<typename MeshT>
    Svg& wireframe(const MeshT& m, std::string stroke = "#585b70") {
        for (const auto& [a, b, c] : m.faces) {
            const auto& va = m.vertices[a];
            const auto& vb = m.vertices[b];
            const auto& vc = m.vertices[c];
            double x0 = static_cast<double>(va[0]);
            double y0 = static_cast<double>(va[1]);
            double x1 = static_cast<double>(vb[0]);
            double y1 = static_cast<double>(vb[1]);
            double x2 = static_cast<double>(vc[0]);
            double y2 = static_cast<double>(vc[1]);
            triangle(x0, y0, x1, y1, x2, y2, "#18182466", stroke);
        }
        return *this;
    }

    inline bool save(const std::string& path) const {
        std::ofstream f(path);
        if (!f) return false;
        f << std::format(
            R"(<?xml version="1.0" encoding="UTF-8"?>)" "\n"
            R"(<svg xmlns="http://www.w3.org/2000/svg" )"
            R"(width="{}" height="{}" viewBox="0 0 {} {}">)" "\n"
            R"(<rect width="100%" height="100%" fill="#1e1e2e"/>)" "\n",
            width, height, width, height);
        f << content;
        f << "</svg>\n";
        return true;
    }

    // ── 3D Mesh Projection ────────────────────────────────────

    // Project 3D mesh to SVG wireframe
    template<Surface S, typename Proj>
    Svg& mesh_wireframe(const mesh::Mesh<S>& m, Proj&& proj,
                        std::string stroke = "#cdd6f4", double w = 0.5) {
        for (const auto& [a, b, c] : m.faces) {
            auto pa = proj(m.vertices[a]), pb = proj(m.vertices[b]), pc = proj(m.vertices[c]);
            double x0 = pa[0], y0 = pa[1];
            double x1 = pb[0], y1 = pb[1];
            double x2 = pc[0], y2 = pc[1];
            triangle(x0, y0, x1, y1, x2, y2, "none", stroke);
        }
        return *this;
    }

    // Project 3D mesh with filled faces (painter's algorithm)
    template<Surface S, typename Proj>
    Svg& mesh_filled(const mesh::Mesh<S>& m, Proj&& proj,
                     std::string fill = "#31324466",
                     std::string stroke = "#585b70") {
        // Sort faces by depth (back to front)
        auto faces_sorted = depth_sort(m, proj);
        for (auto fi : faces_sorted) {
            auto [a, b, c] = m.faces[fi];
            auto pa = proj(m.vertices[a]), pb = proj(m.vertices[b]), pc = proj(m.vertices[c]);
            triangle(pa[0], pa[1], pb[0], pb[1], pc[0], pc[1], fill, stroke);
        }
        return *this;
    }

    // Project 3D mesh with per-face scalar coloring
    template<Surface S, typename Proj>
    Svg& mesh_colored(const mesh::Mesh<S>& m, Proj&& proj,
                      std::span<const typename S::ScalarType> values,
                      std::function<std::string(double)> cmap = viridis_map(),
                      std::string stroke = "#585b7044") {
        using T = typename S::ScalarType;
        T vmin = *std::min_element(values.begin(), values.end());
        T vmax = *std::max_element(values.begin(), values.end());
        T range = vmax - vmin;
        if (range < T{1e-15}) range = T{1};

        auto faces_sorted = depth_sort(m, proj);
        for (auto fi : faces_sorted) {
            auto [a, b, c] = m.faces[fi];
            double t = static_cast<double>((values[a] + values[b] + values[c]) / T{3} - vmin) / static_cast<double>(range);
            auto pa = proj(m.vertices[a]), pb = proj(m.vertices[b]), pc = proj(m.vertices[c]);
            triangle(pa[0], pa[1], pb[0], pb[1], pc[0], pc[1], cmap(t), stroke);
        }
        return *this;
    }

    // ── Color maps ────────────────────────────────────────────

    static std::function<std::string(double)> viridis_map() {
        return [](double t) -> std::string {
            t = std::clamp(t, 0.0, 1.0);
            // Simplified viridis: purple → teal → yellow
            int r = static_cast<int>(std::clamp(68 + t * (253 - 68), 0.0, 255.0));
            int g = static_cast<int>(std::clamp(1 + t * (231 - 1), 0.0, 255.0));
            int b_val = static_cast<int>(std::clamp(84 + (0.5 - std::abs(t - 0.5)) * 2 * (150 - 84), 0.0, 255.0));
            return std::format("#{:02x}{:02x}{:02x}", r, g, b_val);
        };
    }

    static std::function<std::string(double)> heat_map() {
        return [](double t) -> std::string {
            t = std::clamp(t, 0.0, 1.0);
            int r = static_cast<int>(std::min(255.0, t * 2 * 255));
            int g = static_cast<int>(std::max(0.0, (t - 0.5) * 2 * 255));
            int b_val = static_cast<int>(std::max(0.0, (1.0 - t * 2) * 255));
            return std::format("#{:02x}{:02x}{:02x}", r, g, b_val);
        };
    }

private:
    // Depth-sort faces (back to front) using centroid Z after projection
    template<Surface S, typename Proj>
    std::vector<std::size_t> depth_sort(const mesh::Mesh<S>& m, Proj&& proj) const {
        std::vector<std::pair<double, std::size_t>> depths(m.faces.size());
        for (std::size_t i = 0; i < m.faces.size(); ++i) {
            auto [a, b, c] = m.faces[i];
            // Use original 3D centroid Z for depth (before projection flattens it)
            auto centroid = (m.vertices[a] + m.vertices[b] + m.vertices[c]) / 3.0;
            // Project and use the 3D z-coordinate of the centroid for sorting
            // Convention: higher z = closer to viewer (drawn last)
            double z = 0;
            if constexpr (S::PointType::size >= 3)
                z = static_cast<double>(centroid[2]);
            depths[i] = {z, i};
        }
        std::sort(depths.begin(), depths.end());  // ascending = far first
        std::vector<std::size_t> result(m.faces.size());
        for (std::size_t i = 0; i < m.faces.size(); ++i)
            result[i] = depths[i].second;
        return result;
    }
};

// ── Orthographic Projection ───────────────────────────────────

struct OrthoProjection {
    Vec3 forward;  // viewing direction (normalized)
    Vec3 up;       // up direction (normalized)
    Vec3 right_;

    OrthoProjection(Vec3 fwd = {0, 0, -1}, Vec3 u = {0, 1, 0})
        : forward(fwd.normalized()), up(u.normalized())
    {
        right_ = forward.cross(up).normalized();
        up = Vec3{right_.cross(forward)}.normalized();  // orthogonalize
    }

    template<typename P>
    Vec2 operator()(const P& p) const {
        return {
            static_cast<double>(p[0]) * right_[0] + static_cast<double>(p[1]) * right_[1] + static_cast<double>(p[2]) * right_[2],
            static_cast<double>(p[0]) * up[0] + static_cast<double>(p[1]) * up[1] + static_cast<double>(p[2]) * up[2]
        };
    }
};

// ── Convenience: one-shot mesh → SVG ──────────────────────────

template<Surface S>
Svg mesh_to_svg(const mesh::Mesh<S>& m,
                Vec3 view_dir = {0.5, -0.3, -1.0},
                Vec3 up = {0, 1, 0},
                double size = 800, double scale = 100) {
    Svg svg(size, size, scale);
    OrthoProjection proj(view_dir, up);
    svg.mesh_filled(m, proj);
    return svg;
}

template<Surface S>
Svg mesh_to_svg_colored(const mesh::Mesh<S>& m,
                        std::span<const typename S::ScalarType> values,
                        Vec3 view_dir = {0.5, -0.3, -1.0},
                        Vec3 up = {0, 1, 0},
                        double size = 800, double scale = 100) {
    Svg svg(size, size, scale);
    OrthoProjection proj(view_dir, up);
    svg.mesh_colored(m, proj, values);
    return svg;
}

} // namespace spatium::io
