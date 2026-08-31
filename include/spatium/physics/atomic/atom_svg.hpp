#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/physics/atomic/orbital.hpp>
#  include <spatium/physics/elements.hpp>
#  include <spatium/io/svg.hpp>
#  include <cmath>
#  include <format>
#  include <numbers>
#  include <string>
#endif

SPATIUM_EXPORT namespace spatium::physics::atomic {

// ── Orbital cross-section SVG ─────────────────────────────────
// Evaluates |ψ|² on a 2D grid in the XZ plane, renders as colored pixels.

inline io::Svg orbital_cross_section_svg(
    int n, int l, int m,
    double size = 600, int resolution = 100)
{
    double bound = n * n * 4.0 + 4.0;
    double step = 2.0 * bound / resolution;

    io::Svg svg(size, size, size / (2.0 * bound));

    // Evaluate density on grid
    double max_density = 0;
    std::vector<double> grid(resolution * resolution);

    for (int j = 0; j < resolution; ++j) {
        double z = -bound + (j + 0.5) * step;
        for (int i = 0; i < resolution; ++i) {
            double x = -bound + (i + 0.5) * step;
            // XZ plane (y=0)
            double d = orbital_density_xyz<double>(n, l, m, x, 0.0, z);
            grid[j * resolution + i] = d;
            max_density = std::max(max_density, d);
        }
    }

    if (max_density < 1e-30) return svg;

    // Render pixels
    auto cmap = io::Svg::viridis_map();
    double pixel_w = step;
    double pixel_h = step;

    for (int j = 0; j < resolution; ++j) {
        double z = -bound + (j + 0.5) * step;
        for (int i = 0; i < resolution; ++i) {
            double x = -bound + (i + 0.5) * step;
            double d = grid[j * resolution + i];
            if (d < max_density * 0.01) continue;  // skip near-zero

            double t = std::sqrt(d / max_density);  // sqrt for better contrast
            std::string color = cmap(t);

            // Draw filled rect
            double sx1 = svg.sx(x - pixel_w / 2);
            double sy1 = svg.sy(z + pixel_h / 2);
            double sw = pixel_w * svg.scale;
            double sh = pixel_h * svg.scale;

            svg.content += std::format(
                R"(<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}" fill="{}" stroke="none"/>)" "\n",
                sx1, sy1, sw, sh, color);
        }
    }

    // Label
    svg.text(-bound * 0.9, bound * 0.9,
             std::format("n={} l={} m={} ({})", n, l, m,
                         std::string(1, subshell_letter(l))));

    return svg;
}

// ── Periodic table SVG ────────────────────────────────────────

inline io::Svg periodic_table_svg(double width = 1400, double height = 900) {
    io::Svg svg(width, height, 1.0);
    svg.cx = 0; svg.cy = 0;  // raw pixel coords

    double cell_w = width / 19.0;
    double cell_h = height / 11.0;
    double margin = cell_w * 0.5;

    // Element category colors (Catppuccin Mocha palette)
    auto type_color = [](int Z) -> std::string {
        // Alkali metals
        if (Z == 3 || Z == 11 || Z == 19 || Z == 37 || Z == 55 || Z == 87)
            return "#89b4fa44";
        // Alkaline earth metals
        if (Z == 4 || Z == 12 || Z == 20 || Z == 38 || Z == 56 || Z == 88)
            return "#89b4fa33";
        // Transition metals
        if ((Z >= 21 && Z <= 30) || (Z >= 39 && Z <= 48) ||
            (Z >= 72 && Z <= 80) || (Z >= 104 && Z <= 112))
            return "#a6e3a133";
        // Lanthanides
        if (Z >= 57 && Z <= 71)
            return "#f9e2af33";
        // Actinides
        if (Z >= 89 && Z <= 103)
            return "#fab38733";
        // Noble gases
        if (Z == 2 || Z == 10 || Z == 18 || Z == 36 || Z == 54 || Z == 86 || Z == 118)
            return "#585b70";
        // Halogens
        if (Z == 9 || Z == 17 || Z == 35 || Z == 53 || Z == 85 || Z == 117)
            return "#cba6f744";
        // Metalloids
        if (Z == 5 || Z == 14 || Z == 32 || Z == 33 || Z == 51 || Z == 52 || Z == 84)
            return "#94e2d544";
        // Post-transition metals
        if (Z == 13 || Z == 31 || Z == 49 || Z == 50 || Z == 81 || Z == 82 || Z == 83 ||
            Z == 113 || Z == 114 || Z == 115 || Z == 116)
            return "#74c7ec33";
        // Nonmetals (H, C, N, O, P, S, Se)
        if (Z == 1 || Z == 6 || Z == 7 || Z == 8 || Z == 15 || Z == 16 || Z == 34)
            return "#45475a";
        return "#31324466";
    };

    // Full periodic table layout: {Z, row, col}
    struct Pos { int Z; int row; int col; };
    std::vector<Pos> layout;

    // Period 1
    layout.push_back({1, 0, 0});
    layout.push_back({2, 0, 17});
    // Period 2
    for (int Z = 3; Z <= 4; ++Z) layout.push_back({Z, 1, Z - 3});
    for (int Z = 5; Z <= 10; ++Z) layout.push_back({Z, 1, Z - 5 + 12});
    // Period 3
    for (int Z = 11; Z <= 12; ++Z) layout.push_back({Z, 2, Z - 11});
    for (int Z = 13; Z <= 18; ++Z) layout.push_back({Z, 2, Z - 13 + 12});
    // Period 4
    for (int Z = 19; Z <= 20; ++Z) layout.push_back({Z, 3, Z - 19});
    for (int Z = 21; Z <= 30; ++Z) layout.push_back({Z, 3, Z - 21 + 2});
    for (int Z = 31; Z <= 36; ++Z) layout.push_back({Z, 3, Z - 31 + 12});
    // Period 5
    for (int Z = 37; Z <= 38; ++Z) layout.push_back({Z, 4, Z - 37});
    for (int Z = 39; Z <= 48; ++Z) layout.push_back({Z, 4, Z - 39 + 2});
    for (int Z = 49; Z <= 54; ++Z) layout.push_back({Z, 4, Z - 49 + 12});
    // Period 6: s-block + d-block + p-block (skip lanthanides in main table)
    layout.push_back({55, 5, 0});
    layout.push_back({56, 5, 1});
    for (int Z = 72; Z <= 80; ++Z) layout.push_back({Z, 5, Z - 72 + 3});
    for (int Z = 81; Z <= 86; ++Z) layout.push_back({Z, 5, Z - 81 + 12});
    // Period 7: s-block + d-block + p-block (skip actinides in main table)
    layout.push_back({87, 6, 0});
    layout.push_back({88, 6, 1});
    for (int Z = 104; Z <= 112; ++Z) layout.push_back({Z, 6, Z - 104 + 3});
    for (int Z = 113; Z <= 118; ++Z) layout.push_back({Z, 6, Z - 113 + 12});
    // Lanthanides (row 8, cols 2-16)
    for (int Z = 57; Z <= 71; ++Z) layout.push_back({Z, 8, Z - 57 + 2});
    // Actinides (row 9, cols 2-16)
    for (int Z = 89; Z <= 103; ++Z) layout.push_back({Z, 9, Z - 89 + 2});

    for (auto& [Z, row, col] : layout) {
        double x = margin + col * cell_w;
        double y = margin + row * cell_h;

        // Cell background
        svg.content += std::format(
            R"(<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}" )"
            R"(fill="{}" stroke="#585b70" stroke-width="1" rx="3"/>)" "\n",
            x, y, cell_w - 2, cell_h - 2, type_color(Z));

        // Z number (top-left, small)
        svg.content += std::format(
            R"(<text x="{:.1f}" y="{:.1f}" fill="#a6adc8" font-size="10" font-family="monospace">{}</text>)" "\n",
            x + 4, y + 12, Z);

        // Symbol (center, large)
        auto& elem = element(Z);
        svg.content += std::format(
            R"(<text x="{:.1f}" y="{:.1f}" fill="#cdd6f4" font-size="18" font-family="monospace" )"
            R"(text-anchor="middle">{}</text>)" "\n",
            x + (cell_w - 2) / 2, y + cell_h / 2 + 4, elem.symbol);
    }

    // Lanthanide/actinide markers in main table
    auto draw_marker = [&](int row, const char* label) {
        double x = margin + 2 * cell_w;
        double y = margin + row * cell_h;
        svg.content += std::format(
            R"(<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}" )"
            R"(fill="none" stroke="#585b70" stroke-width="1" stroke-dasharray="4" rx="3"/>)" "\n",
            x, y, cell_w - 2, cell_h - 2);
        svg.content += std::format(
            R"(<text x="{:.1f}" y="{:.1f}" fill="#6c7086" font-size="9" font-family="monospace" )"
            R"(text-anchor="middle">{}</text>)" "\n",
            x + (cell_w - 2) / 2, y + cell_h / 2 + 3, label);
    };
    draw_marker(5, "57-71");
    draw_marker(6, "89-103");

    // Title
    svg.content += std::format(
        R"(<text x="{:.1f}" y="{:.1f}" fill="#cdd6f4" font-size="24" font-family="monospace" )"
        R"(text-anchor="middle">Periodic Table of Elements</text>)" "\n",
        width / 2, height - 10.0);

    return svg;
}

} // namespace spatium::physics::atomic
