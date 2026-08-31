// Tour of Spatium: primitives, distances, transforms, morphisms,
// quaternions, sphere/hyperbolic geodesics, mesh subdivision,
// multiprecision, and SVG output — the broadest of the text demos.

#include "io_helpers.hpp"

#include <spatium/spatium.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/io/table.hpp>
#include <spatium/io/svg.hpp>
#include <filesystem>
#include <iostream>
#include <numbers>
#include <print>
#include <string>
#include <string_view>

using namespace spatium;
using namespace spatium::literals;

namespace {

void print_usage() {
    std::cout <<
        "Usage: showcase [--help] [--force] [--out-prefix PREFIX]\n"
        "  Tour of Spatium: primitives, distances, transforms, morphisms,\n"
        "  quaternions, sphere/hyperbolic geodesics, mesh subdivision,\n"
        "  multiprecision and SVG output.\n"
        "  --help, -h         show this message\n"
        "  --force            overwrite existing output files\n"
        "  --out-prefix PFX   prefix prepended to every output filename\n"
        "  Outputs:           <prefix>icosphere_l2.svg\n";
}

}  // namespace

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

    // ── Geometry Primitives ────────────────────────────────────
    section("Geometry Primitives");

    auto t = tri(Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0});
    Table("Property", "Value")
        .row("shape",     std::format("{}", t))
        .row("area",      std::format("{:.6f}", t.area()))
        .row("perimeter", std::format("{:.6f}", t.perimeter()))
        .row("normal",    std::format("{}", t.normal()))
        .row("centroid",  std::format("{}", t.centroid()))
        .print();

    // More primitives
    auto s = seg(Vec3{0, 0, 0}, Vec3{3, 4, 0});
    auto c = circle(Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1});
    auto d = disk(Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1});
    auto hex = poly<3>({
        Vec3{1, 0, 0}, Vec3{0.5, 0.866, 0}, Vec3{-0.5, 0.866, 0},
        Vec3{-1, 0, 0}, Vec3{-0.5, -0.866, 0}, Vec3{0.5, -0.866, 0}
    });

    Table("Shape", "Measure")
        .row("segment",     std::format("length = {:.4f}", s.length()))
        .row("circle",      std::format("circum = {:.4f}", c.circumference()))
        .row("disk",         std::format("area = {:.6f}  (pi = {:.6f})", d.area(), std::numbers::pi))
        .row("hexagon",     std::format("area = {:.6f}", hex.area()))
        .print();

    // ── Intersection Pipeline ──────────────────────────────────
    section("Intersection (pipe syntax)");

    auto r = *ray(Vec3{0.25, 0.25, 5.0}, Vec3{0, 0, -1});
    std::println("  ray:      {}", r);
    std::println("  triangle: {}", t);
    if (auto hit = r | t)
        std::println("  hit:      {}", *hit);

    auto p1 = *plane(Vec3{0, 0, 1}, Vec3{0, 0, 0});
    auto p2 = *plane(Vec3{0, 1, 0}, Vec3{0, 0, 0});
    if (auto l = p1 | p2)
        std::println("  {} | {} = {}", p1, p2, *l);

    // Ray-box intersection
    auto bx = box(Vec3{-1, -1, -1}, Vec3{1, 1, 1});
    auto rb = *ray(Vec3{-5, 0, 0}, Vec3{1, 0, 0});
    if (auto hit = rb | bx)
        std::println("  ray | box = {}", *hit);

    // ── Distance ───────────────────────────────────────────────
    section("Distance Functions");

    auto l1 = *line(Vec3{0, 0, 0}, Vec3{1, 0, 0});
    Vec3 point{0, 3, 4};
    Table("Query", "Distance")
        .row("point → line",     std::format("{:.4f}", distance(point, l1)))
        .row("point → triangle", std::format("{:.4f}", distance(point, t)))
        .row("point → box",      std::format("{:.4f}", distance(point, bx)))
        .print();

    // ── Clip (containment) ─────────────────────────────────────
    section("Clip (containment gate)");

    Vec3 inside{0.2, 0.2, 0};
    Vec3 outside{5, 5, 0};
    std::println("  clip({}, tri): {}", inside, clip(inside, t).has_value() ? "inside" : "outside");
    std::println("  clip({}, tri): {}", outside, clip(outside, t).has_value() ? "inside" : "outside");
    std::println("  clip({}, box): {}", Vec3{0,0,0}, clip(Vec3{0,0,0}, bx).has_value() ? "inside" : "outside");

    // ── Affine Transforms ──────────────────────────────────────
    section("Affine Transforms");

    auto chain = translate(10, 0, 0) * rotate_z(45_deg) * scale(2);
    auto p = pt<E3>(Vec3{1, 0, 0});
    std::println("  (1,0,0) | scale(2) * rot(45) * shift(10) = {}", chain(p).raw());
    std::println("  (1,0,0) | translate = {}", (p | translate(10, 0, 0)).raw());

    // ── Morphism Pipeline ──────────────────────────────────────
    section("Morphism Pipeline");

    auto scale   = morph<E3, E3>([](const Vec3& v) { return v * 2.0; });
    auto shift   = morph<E3, E3>([](const Vec3& v) { return v + Vec3{10, 0, 0}; });
    auto proj_xy = morph<E3, E2>([](const Vec3& v) -> Vec2 { return {v[0], v[1]}; });

    auto pipeline = scale | shift | proj_xy;
    auto origin = pt<E3>(Vec3{1, 2, 3});
    auto result = origin | pipeline;

    std::println("  {} | scale*2 | shift+10 | proj_xy = {}", origin, result);

    // identity
    auto id = identity<E3>();
    std::println("  {} | identity = {}", origin, (origin | id));

    // ── Quaternion ──────────────────────────────────────────────
    section("Quaternion Rotation");

    auto q = Quat::from_axis_angle(Vec3{0, 0, 1}, 90_deg);
    auto rotated = q.rotate(Vec3{1, 0, 0});
    std::println("  rotate (1,0,0) by 90 around Z: ({:.4f}, {:.4f}, {:.4f})",
                 rotated[0], rotated[1], rotated[2]);

    auto q2 = Quat::from_axis_angle(Vec3{0, 1, 0}, 180_deg);
    auto mid = Quat::slerp(q, q2, 0.5);
    std::println("  slerp(q1, q2, 0.5) norm = {:.4f}", mid.norm());

    // ── Sphere: Geodesics ──────────────────────────────────────
    section("Sphere S2 — Geodesics");

    S2 sphere;
    auto north = pt<S2>(Vec3{0, 0, 1});
    auto east  = pt<S2>(Vec3{1, 0, 0});

    Table("Metric", "Value")
        .row("d(N, E)", std::format("{:.10f}", north.distance_to(east, sphere)))
        .row("pi/2",    std::format("{:.10f}", std::numbers::pi / 2))
        .print();

    auto tangent = north.log(east, sphere);
    auto geo_mid = north.exp(tangent, 0.5, sphere);
    std::println("  geodesic midpoint: {}", geo_mid);

    // ── Hyperbolic Space ───────────────────────────────────────
    section("Hyperbolic H2");

    H2 hyp;
    auto o = pt<H2>(H2::origin());
    auto hq = pt<H2>(Vec3{std::cosh(2.0), std::sinh(2.0), 0.0});

    Table("Metric", "Value")
        .row("d(O, Q)", std::format("{:.10f}", o.distance_to(hq, hyp)))
        .row("expected", "2.0000000000")
        .print();

    // ── Verify Axioms ──────────────────────────────────────────
    section("Axiom Verification");

    auto m_ok = verify_metric(sphere,
        {Vec3{0,0,1}, Vec3{1,0,0}, Vec3{0,1,0}, Vec3{-1,0,0}});
    auto ip_ok = verify_inner_product(E3{},
        {Vec3{1,0,0}, Vec3{0,1,0}, Vec3{0,0,1}, Vec3{1,1,1}});
    auto el_ok = verify_exp_log(sphere,
        {Vec3{0,0,1}, Vec3{1,0,0}, Vec3{0,1,0}});

    Table("Space", "Axiom | Result")
        .row("S2",  std::format("metric        | {}", m_ok.passed ? "PASS" : "FAIL"))
        .row("E3",  std::format("inner product | {}", ip_ok.passed ? "PASS" : "FAIL"))
        .row("S2",  std::format("exp/log       | {}", el_ok.passed ? "PASS" : "FAIL"))
        .print();

    // ── Mesh Subdivision ───────────────────────────────────────
    section("Mesh Subdivision (Icosphere)");

    auto lod = LodChain<S2>::build(icosahedron(sphere), sphere, 4);
    auto true_area = 4.0 * std::numbers::pi;

    Table("LOD", "Faces | Vertices | Area | Error")
        .row("base", std::format("{:>5} | {:>8} | {:.6f} | {:.2e}",
             lod.at(0).face_count(), lod.at(0).vertex_count(),
             lod.at(0).area(sphere), std::abs(lod.at(0).area(sphere) - true_area)))
        .row("L1", std::format("{:>5} | {:>8} | {:.6f} | {:.2e}",
             lod.at(1).face_count(), lod.at(1).vertex_count(),
             lod.at(1).area(sphere), std::abs(lod.at(1).area(sphere) - true_area)))
        .row("L2", std::format("{:>5} | {:>8} | {:.6f} | {:.2e}",
             lod.at(2).face_count(), lod.at(2).vertex_count(),
             lod.at(2).area(sphere), std::abs(lod.at(2).area(sphere) - true_area)))
        .row("L3", std::format("{:>5} | {:>8} | {:.6f} | {:.2e}",
             lod.at(3).face_count(), lod.at(3).vertex_count(),
             lod.at(3).area(sphere), std::abs(lod.at(3).area(sphere) - true_area)))
        .row("L4", std::format("{:>5} | {:>8} | {:.6f} | {:.2e}",
             lod.at(4).face_count(), lod.at(4).vertex_count(),
             lod.at(4).area(sphere), std::abs(lod.at(4).area(sphere) - true_area)))
        .row("4pi", std::format("{:.10f}", true_area))
        .print();

    // ── Multiprecision ─────────────────────────────────────────
    section("Multiprecision (50-digit precision)");

    using boost::multiprecision::asin;
    Real50 pi50 = asin(Real50{1}) * Real50{2};

    int n = 100000;
    Real50 angle = Real50{2} * pi50 / Real50{n};
    Real50 poly_pi = Real50{n} / Real50{2} * sin(angle);
    Real50 error = abs(poly_pi - pi50);

    Table("", "")
        .row("pi (50-digit)", pi50.str(40))
        .row("100000-gon",    poly_pi.str(40))
        .row("error",         error.str(10))
        .print();

    // ── SVG Output ─────────────────────────────────────────────
    section("SVG Output");

    Svg svg(600, 600, 250);
    svg.circle(0, 0, 1, "#89b4fa", 1.5);
    svg.wireframe(lod.at(2), "#585b70");
    svg.point(0, 0, 4, "#f38ba8");
    svg.text(0.05, 0.05, "O");

    std::string svg_path = out_prefix + "icosphere_l2.svg";
    if (spatium::examples::confirm_overwrite(svg_path, force)) {
        if (svg.save(svg_path))
            std::println("  saved: {} ({} faces)", svg_path, lod.at(2).face_count());
    }

    std::println("\ndone.");
}
