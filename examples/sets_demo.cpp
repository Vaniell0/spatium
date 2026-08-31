// Finite set operations and a geometric Venn diagram of two overlapping
// circles, with intersection/symmetric-difference area readouts and an
// SVG render — the discrete/ domain's only demo.

#include "io_helpers.hpp"

#include <spatium/spatium.hpp>
#include <spatium/discrete/finite_set.hpp>
#include <spatium/io/svg.hpp>
#include <filesystem>
#include <format>
#include <iostream>
#include <numbers>
#include <print>
#include <string>
#include <string_view>

using namespace spatium;
using namespace spatium::discrete;

namespace {

void print_usage() {
    std::cout <<
        "Usage: sets_demo [--help] [--force] [--out-prefix PREFIX]\n"
        "  Finite set operations and a geometric Venn diagram of two\n"
        "  overlapping circles, with intersection / symmetric-difference\n"
        "  area readouts and an SVG render.\n"
        "  --help, -h         show this message\n"
        "  --force            overwrite existing output files\n"
        "  --out-prefix PFX   prefix prepended to every output filename\n"
        "  Outputs:           <prefix>venn_sets.svg\n";
}

}  // namespace

// Convert a Polygon<2> to SVG points
std::vector<std::pair<double, double>> to_svg_pts(const Polygon<2>& p) {
    std::vector<std::pair<double, double>> pts;
    for (const auto& v : p.vertices)
        pts.emplace_back(v[0], v[1]);
    return pts;
}

// Make a regular polygon (circle approximation) in 2D
Polygon<2> make_circle_poly(Vec2 center, double r, int n = 48) {
    std::vector<Vec2> verts;
    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * std::numbers::pi * i / n;
        verts.push_back(center + Vec2{r * std::cos(angle), r * std::sin(angle)});
    }
    return Polygon<2>{std::move(verts)};
}

int main(int argc, char** argv) {
    bool force = false;
    std::string out_prefix;
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--help" || a == "-h") { print_usage(); return 0; }
        if (a == "--force") { force = true; continue; }
        if (a == "--out-prefix" && i + 1 < argc) { out_prefix = argv[++i]; continue; }
        std::cerr << "unknown option: " << a << "\n";
        return 1;
    }

    // ── Part 1: Finite Set Operations ─────────────────────────
    std::println("=== Finite Set Operations ===\n");

    FiniteSet<int> A = {1, 2, 3, 4, 5};
    FiniteSet<int> B = {3, 4, 5, 6, 7};

    auto print_set = [](const std::string& name, const FiniteSet<int>& s) {
        std::print("  {} = {{", name);
        for (std::size_t i = 0; i < s.size(); ++i)
            std::print("{}{}", s.elements[i], i + 1 < s.size() ? ", " : "");
        std::println("}}  |{}| = {}", name, s.size());
    };

    print_set("A", A);
    print_set("B", B);
    std::println("");

    print_set("A ∪ B", A + B);
    print_set("A ∩ B", A & B);
    print_set("A ∖ B", A - B);
    print_set("B ∖ A", B - A);
    print_set("A △ B", A ^ B);

    std::println("\n  A ⊆ B? {}",  A <= B ? "yes" : "no");
    std::println("  (A ∩ B) ⊆ A? {}", (A & B) <= A ? "yes" : "no");

    std::println("\n  3 ∈ A? {}", A.contains(3) ? "yes" : "no");
    std::println("  6 ∈ A? {}", A.contains(6) ? "yes" : "no");

    // Power set of small set
    FiniteSet<int> C = {1, 2, 3};
    auto ps = C.power_set();
    std::println("\n  P({{1,2,3}}) has {} subsets", ps.size());

    // Cartesian product
    FiniteSet<char> X = {'a', 'b'};
    FiniteSet<int> Y = {1, 2, 3};
    auto prod = X.cartesian(Y);
    std::println("  {{a,b}} × {{1,2,3}} has {} pairs", prod.size());

    // ── Part 2: Geometric Venn Diagram ────────────────────────
    std::println("\n=== Geometric Venn Diagram ===\n");

    // Two overlapping circles as polygons
    auto circle_a = make_circle_poly(Vec2{-0.6, 0.0}, 1.2);
    auto circle_b = make_circle_poly(Vec2{0.6, 0.0}, 1.2);

    auto inter = intersection_region(circle_a, circle_b);

    std::println("  Circle A: center=(-0.6, 0), r=1.2, area={:.3f}", circle_a.area());
    std::println("  Circle B: center=(0.6, 0), r=1.2, area={:.3f}", circle_b.area());
    if (inter) {
        std::println("  A ∩ B: {} vertices, area={:.3f}", inter->size(), inter->area());
        auto sd = symmetric_difference_area<2, double>(circle_a, circle_b);
        if (sd) std::println("  A △ B: area={:.3f}", *sd);
    }

    // ── SVG Output ────────────────────────────────────────────
    spatium::io::Svg svg(600, 400, 120);

    // A ∖ B (blue, left crescent)
    svg.polygon(to_svg_pts(circle_a), "#89b4fa44", "#89b4fa", 2.0);

    // B ∖ A (green, right crescent)
    svg.polygon(to_svg_pts(circle_b), "#a6e3a144", "#a6e3a1", 2.0);

    // A ∩ B (purple, overlap)
    if (inter)
        svg.polygon(to_svg_pts(*inter), "#cba6f788", "#cba6f7", 2.0);

    // Labels
    svg.text(-1.2, 1.5, "A", "#89b4fa", 18);
    svg.text(0.9, 1.5, "B", "#a6e3a1", 18);
    if (inter)
        svg.text(-0.15, -0.1, std::format("{:.1f}", inter->area()), "#cba6f7", 14);

    // Element markers for A = {1,2,3,4,5}, B = {3,4,5,6,7}
    // A only: 1, 2
    svg.point(-1.2, 0.3, 4, "#89b4fa");
    svg.text(-1.15, 0.35, "1", "#89b4fa", 10);
    svg.point(-1.0, -0.3, 4, "#89b4fa");
    svg.text(-0.95, -0.25, "2", "#89b4fa", 10);

    // Intersection: 3, 4, 5
    svg.point(-0.15, 0.4, 4, "#cba6f7");
    svg.text(-0.1, 0.45, "3", "#cba6f7", 10);
    svg.point(0.0, 0.0, 4, "#cba6f7");
    svg.text(0.05, 0.05, "4", "#cba6f7", 10);
    svg.point(-0.1, -0.4, 4, "#cba6f7");
    svg.text(-0.05, -0.35, "5", "#cba6f7", 10);

    // B only: 6, 7
    svg.point(1.0, 0.3, 4, "#a6e3a1");
    svg.text(1.05, 0.35, "6", "#a6e3a1", 10);
    svg.point(1.2, -0.3, 4, "#a6e3a1");
    svg.text(1.25, -0.25, "7", "#a6e3a1", 10);

    std::string svg_path = out_prefix + "venn_sets.svg";
    if (spatium::examples::confirm_overwrite(svg_path, force)) {
        svg.save(svg_path);
        std::println("\n  Saved: {}", svg_path);
    }
}
