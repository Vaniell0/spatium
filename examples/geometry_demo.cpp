// Console tour of geometry primitives, transforms, morphisms, quaternions,
// geodesics, hyperbolic distance, and axiom checks — no file output, no
// viewer. Smallest/plainest of the text demos; see showcase.cpp for the
// SVG-producing counterpart.

#include <spatium/spatium.hpp>
#include <format>
#include <iostream>
#include <print>
#include <string>
#include <string_view>

using namespace spatium;
using namespace spatium::literals;

namespace {

void print_usage() {
    std::cout <<
        "Usage: geometry_demo [--help] [--force] [--out-prefix PREFIX]\n"
        "  Console tour of geometry primitives, transforms, morphisms,\n"
        "  quaternions, geodesics, hyperbolic distance and axiom checks.\n"
        "  --help, -h         show this message\n"
        "  --force            accepted for CLI symmetry; demo writes no files\n"
        "  --out-prefix PFX   accepted for CLI symmetry; demo writes no files\n";
}

}  // namespace

int main(int argc, char** argv) {
    [[maybe_unused]] bool force = false;
    [[maybe_unused]] std::string out_prefix;
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--help" || a == "-h") { print_usage(); return 0; }
        if (a == "--force") { force = true; continue; }
        if (a == "--out-prefix" && i + 1 < argc) { out_prefix = argv[++i]; continue; }
        std::cerr << "unknown option: " << a << "\n";
        return 1;
    }

    // ── Factories: clean construction ─────────────────────────────
    auto t = tri(Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0});
    std::println("triangle area = {:.4f}, normal = {}", t.area(), t.normal());

    auto s = seg(Vec3{0, 0, 0}, Vec3{1, 1, 1});
    std::println("segment length = {:.4f}", s.length());

    auto c = circle(Vec3{0, 0, 0}, 2.0, Vec3{0, 0, 1});
    std::println("circle circumference = {:.4f}", c.circumference());

    auto d = disk(Vec3{0, 0, 0}, 2.0, Vec3{0, 0, 1});
    std::println("disk area = {:.4f}", d.area());

    auto pentagon = poly<3>({
        Vec3{1, 0, 0}, Vec3{0.309, 0.951, 0}, Vec3{-0.809, 0.588, 0},
        Vec3{-0.809, -0.588, 0}, Vec3{0.309, -0.951, 0}
    });
    std::println("pentagon area = {:.4f}", pentagon.area());

    auto subdivided = t.subdivide();
    std::println("after subdivision: {} triangles", subdivided.size());

    // ── Intersection via pipe operator ─────────────────────────────
    // ray(o, d) returns Result<Ray>; pipe works directly without *
    if (auto hit = ray(Vec3{0.25, 0.25, 5}, Vec3{0, 0, -1}) | t)
        std::println("ray | triangle = {}", *hit);

    auto a = box(Vec3{0, 0, 0}, Vec3{2, 2, 2});
    auto b = box(Vec3{1, 1, 1}, Vec3{3, 3, 3});
    std::println("boxes intersect: {}", a.intersects(b));
    if (auto overlap = a.intersection(b))
        std::println("overlap volume = {:.1f}", overlap->measure());

    auto p1 = *plane(Vec3{0, 0, 1}, Vec3{0, 0, 0});
    auto p2 = *plane(Vec3{0, 1, 0}, Vec3{0, 0, 0});
    if (auto l = p1 | p2)
        std::println("{} | {} = {}", p1, p2, *l);

    // plane from 3 points
    auto p3 = *plane(Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0});
    std::println("plane from 3 pts: normal = {}", p3.normal);

    // ── Distance functions ────────────────────────────────────────
    Vec3 point{2, 3, 0};
    auto l1 = *line(Vec3{0, 0, 0}, Vec3{1, 0, 0});
    std::println("point-to-line distance = {:.4f}", distance(point, l1));
    std::println("point-to-triangle distance = {:.4f}", distance(point, t));
    std::println("point-to-box distance = {:.4f}", distance(point, a));

    auto s2 = seg(Vec3{0, 0, 0}, Vec3{1, 0, 0});
    auto s3 = seg(Vec3{0, 1, 0}, Vec3{1, 1, 0});
    std::println("segment-segment distance = {:.4f}", distance(s2, s3));

    // ── Clip: containment gate ────────────────────────────────────
    Vec3 inside{0.2, 0.2, 0};
    Vec3 outside{2, 2, 0};
    std::println("clip inside triangle: {}", clip(inside, t).has_value());
    std::println("clip outside triangle: {}", clip(outside, t).has_value());
    std::println("clip inside box: {}", clip(Vec3{1, 1, 1}, a).has_value());

    // ── Result<T> error handling ──────────────────────────────────
    auto bad_ray = ray(Vec3{0, 0, 0}, Vec3{0, 0, 0}); // zero direction
    if (!bad_ray)
        std::println("ray error: {}", bad_ray.error().message);

    auto bad_plane = plane(Vec3{0, 0, 0}, Vec3{0, 0, 0}); // zero normal
    if (!bad_plane)
        std::println("plane error: {}", bad_plane.error().message);

    // ── Affine transforms ─────────────────────────────────────────
    auto combined = translate(5, 0, 0) * rotate_z(45_deg) * scale(2);
    auto p = pt<E3>(Vec3{1, 0, 0});
    auto transformed = combined(p);
    std::println("transform (1,0,0): {}", transformed.raw());

    // pipe syntax
    auto via_pipe = p | translate(5, 0, 0);
    std::println("pipe translate: {}", via_pipe.raw());

    // ── Morphism pipeline ──────────────────────────────────────────
    auto scale_m   = morph<E3, E3>([](const Vec3& v) { return v * 2.0; });
    auto shift_m   = morph<E3, E3>([](const Vec3& v) { return v + Vec3{10, 0, 0}; });
    auto proj_xy   = morph<E3, E2>([](const Vec3& v) -> Vec2 { return {v[0], v[1]}; });
    auto pipeline  = scale_m | shift_m | proj_xy;

    auto origin = pt<E3>(Vec3{1, 2, 3});
    auto result = origin | pipeline;
    std::println("(1,2,3) | scale*2 | shift+10 | proj_xy = {}", result.raw());

    // identity morphism
    auto id = identity<E3>();
    auto same = origin | id;
    std::println("identity: {} == {}", origin.raw(), same.raw());

    // ── Quaternion rotation ────────────────────────────────────────
    auto q = Quat::from_axis_angle(Vec3{0, 0, 1}, 90_deg);
    auto rotated = q.rotate(Vec3{1, 0, 0});
    std::println("quat rotate (1,0,0) by 90 around Z = ({:.4f}, {:.4f}, {:.4f})",
                 rotated[0], rotated[1], rotated[2]);
    std::println("quat: {}", q);

    // slerp (free function)
    auto q2 = Quat::from_axis_angle(Vec3{0, 0, 1}, 180_deg);
    auto mid_q = slerp(q, q2, 0.5);
    auto slerped = mid_q.rotate(Vec3{1, 0, 0});
    std::println("slerp(90, 180, 0.5) applied to (1,0,0) = ({:.4f}, {:.4f}, {:.4f})",
                 slerped[0], slerped[1], slerped[2]);

    // ── Sphere: geodesics ──────────────────────────────────────────
    S2 sphere;
    auto north = pt<S2>(Vec3{0, 0, 1});
    auto east  = pt<S2>(Vec3{1, 0, 0});

    std::println("sphere distance N->E = {:.4f} (pi/2 = {:.4f})",
                 north.distance_to(east, sphere), std::numbers::pi / 2);

    auto tangent = north.log(east, sphere);
    auto midpoint = north.exp(tangent, 0.5, sphere);
    std::println("geodesic midpoint = {}", midpoint.raw());

    // ── Hyperbolic space ───────────────────────────────────────────
    H2 hyp;
    auto o = pt<H2>(H2::origin());
    auto hq = pt<H2>(Vec3{std::cosh(2.0), std::sinh(2.0), 0.0});
    std::println("hyperbolic distance = {:.4f} (expected 2.0)", o.distance_to(hq, hyp));

    // ── Verify axioms ──────────────────────────────────────────────
    auto metric_ok = verify_metric(sphere,
        {Vec3{0,0,1}, Vec3{1,0,0}, Vec3{0,1,0}, Vec3{-1,0,0}});
    std::println("S2 metric axioms: {}{}", metric_ok.passed ? "PASS" : "FAIL",
                 metric_ok.passed ? "" : std::format(" ({})", metric_ok.message()));

    auto inner_ok = verify_inner_product(E3{},
        {Vec3{1,0,0}, Vec3{0,1,0}, Vec3{0,0,1}, Vec3{1,1,1}});
    std::println("E3 inner product axioms: {}", inner_ok.passed ? "PASS" : "FAIL");

    // ── Vec convenience ────────────────────────────────────────────
    auto zero = Vec3::zero();
    auto x_hat = Vec3::unit(0);
    auto y_hat = Vec3::unit(1);
    std::println("zero = {}, x = {}, y = {}", zero, x_hat, y_hat);

    auto a_v = Vec3{1, 2, 3};
    auto b_v = Vec3{4, 5, 6};
    std::println("lerp(a, b, 0.5) = {}", a_v.lerp(b_v, 0.5));
    std::println("a.distance_to(b) = {:.4f}", a_v.distance_to(b_v));
    std::println("a.norm() = {:.4f}", a_v.norm());

    auto reflected = Vec3{1, -1, 0}.reflect(Vec3{0, 1, 0});
    std::println("reflect (1,-1,0) across Y: {}", reflected);

    // ── Morphism inverse ───────────────────────────────────────────
    auto double_it = morph<E3, E3>(
        [](const Vec3& v) { return v * 2.0; },
        [](const Vec3& v) { return v * 0.5; }
    );
    std::println("morph has_inverse: {}", double_it.has_inverse());
    auto half_it = double_it.invert();
    auto test_pt = pt<E3>(Vec3{10, 20, 30});
    std::println("double({}) = {}", test_pt, (test_pt | double_it));
    std::println("half({}) = {}", test_pt, (test_pt | half_it));

    // ── Transform inverse ──────────────────────────────────────────
    auto tr = translate(5, 0, 0);
    if (auto inv = tr.inverse()) {
        auto fwd = tr(Vec3{1, 0, 0});
        auto back = (*inv)(fwd);
        std::println("translate(1,0,0) = {}, inverse = {}", fwd, back);
    }

    // ── Matrix trace ───────────────────────────────────────────────
    auto m = Matrix<double, 3, 3>::identity();
    std::println("I_3 trace = {:.0f}", m.trace());
}
