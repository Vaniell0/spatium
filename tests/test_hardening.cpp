#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/spatium.hpp>
#include <spatium/algebra/groups/so3.hpp>
#include <spatium/algebra/groups/se3.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/geometry/distance.hpp>
#include <numbers>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

// ── Vec degenerate ─────────────────────────────────────────────

TEST_CASE("Vec::normalized zero vector returns zero", "[hardening]") {
    Vec3 zero{0.0, 0.0, 0.0};
    auto n = zero.normalized();
    CHECK(n[0] == 0.0);
    CHECK(n[1] == 0.0);
    CHECK(n[2] == 0.0);
}

TEST_CASE("Vec compound operators", "[hardening]") {
    Vec3 v{1.0, 2.0, 3.0};
    v += Vec3{1.0, 1.0, 1.0};
    CHECK(v == Vec3{2.0, 3.0, 4.0});
    v -= Vec3{1.0, 1.0, 1.0};
    CHECK(v == Vec3{1.0, 2.0, 3.0});
    v *= 2.0;
    CHECK(v == Vec3{2.0, 4.0, 6.0});
    v /= 2.0;
    CHECK(v == Vec3{1.0, 2.0, 3.0});
}

// ── Sphere antipodal log_map ───────────────────────────────────

TEST_CASE("Sphere log_map antipodal returns valid tangent", "[hardening]") {
    S2 sphere;
    Vec3 north{0.0, 0.0, 1.0};
    Vec3 south{0.0, 0.0, -1.0};

    auto v = sphere.log_map(north, south);
    // Should have length = pi (half great circle)
    CHECK_THAT(v.norm(), WithinAbs(std::numbers::pi, 1e-6));
    // Should be perpendicular to north
    CHECK_THAT(v.dot(north), WithinAbs(0.0, 1e-8));
}

TEST_CASE("Sphere log_map same point returns zero", "[hardening]") {
    S2 sphere;
    Vec3 p{1.0, 0.0, 0.0};
    auto v = sphere.log_map(p, p);
    CHECK_THAT(v.norm(), WithinAbs(0.0, 1e-10));
}

// ── Euclidean is Surface ───────────────────────────────────────

TEST_CASE("Euclidean satisfies Surface", "[hardening]") {
    static_assert(Surface<E3>);
    E3 space;
    Vec3 p{1.0, 2.0, 3.0};
    CHECK(space.project(p) == p);
    auto n = space.normal(p);
    CHECK_THAT(n.norm(), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Mesh<Euclidean<3>> compiles", "[hardening]") {
    // Previously impossible: Euclidean didn't satisfy Surface
    E3 space;
    mesh::Mesh<E3> m;
    m.vertices = {Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}};
    m.faces = {{0, 1, 2}};
    CHECK(m.face_count() == 1);
    auto sub = mesh::subdivide_once(m, space);
    CHECK(sub.face_count() == 4);
}

// ── Polygon distance interior ──────────────────────────────────

TEST_CASE("Polygon distance returns 0 for interior point", "[hardening]") {
    Polygon2 sq{{{Vec2{0, 0}, Vec2{1, 0}, Vec2{1, 1}, Vec2{0, 1}}}};
    CHECK(sq.distance(Vec2{0.5, 0.5}) == 0.0);
}

// ── Morphism inverse composition ───────────────────────────────

TEST_CASE("Morphism composition preserves inverse", "[hardening]") {
    auto f = morph<E3, E3>(
        [](const Vec3& p) { return p * 2.0; },
        [](const Vec3& p) { return p * 0.5; }
    );
    auto g = morph<E3, E3>(
        [](const Vec3& p) { return p + Vec3{10, 0, 0}; },
        [](const Vec3& p) { return p - Vec3{10, 0, 0}; }
    );

    auto composed = g * f; // scale then shift
    REQUIRE(composed.inverse.has_value());

    // Forward: (1,0,0) → *2 → (2,0,0) → +10 → (12,0,0)
    auto fwd = composed(Vec3{1, 0, 0});
    CHECK(fwd == Vec3{12, 0, 0});

    // Inverse: (12,0,0) → -10 → (2,0,0) → /2 → (1,0,0)
    auto inv = (*composed.inverse)(Vec3{12, 0, 0});
    CHECK(inv == Vec3{1, 0, 0});
}

TEST_CASE("Morphism composition drops inverse if one is missing", "[hardening]") {
    auto f = morph<E3, E3>([](const Vec3& p) { return p * 2.0; }); // no inverse
    auto g = morph<E3, E3>(
        [](const Vec3& p) { return p + Vec3{1, 0, 0}; },
        [](const Vec3& p) { return p - Vec3{1, 0, 0}; }
    );
    auto composed = g * f;
    CHECK_FALSE(composed.inverse.has_value());
}

// ── Float instantiations ───────────────────────────────────────

TEST_CASE("Euclidean<3, float> works", "[hardening][float]") {
    Euclidean<3, float> space;
    Vec3f a{0.0f, 0.0f, 0.0f};
    Vec3f b{3.0f, 4.0f, 0.0f};
    auto d = space.distance(a, b);
    CHECK_THAT(static_cast<double>(d), WithinAbs(5.0, 1e-5));
    static_assert(EuclideanSpace<Euclidean<3, float>>);
}

TEST_CASE("Triangle<3, float> area", "[hardening][float]") {
    Triangle<3, float> t(
        Vec3f{0, 0, 0}, Vec3f{1, 0, 0}, Vec3f{0, 1, 0}
    );
    CHECK_THAT(static_cast<double>(t.area()), WithinAbs(0.5, 1e-5));
}

TEST_CASE("Box<3, float>", "[hardening][float]") {
    Box<3, float> b{Vec3f{0, 0, 0}, Vec3f{1, 1, 1}};
    CHECK(b.contains(Vec3f{0.5f, 0.5f, 0.5f}));
    CHECK_THAT(static_cast<double>(b.measure()), WithinAbs(1.0, 1e-5));
}

TEST_CASE("Sphere<2, float> distance", "[hardening][float]") {
    Sphere<2, float> sphere;
    Vec3f a{1.0f, 0.0f, 0.0f};
    Vec3f b{0.0f, 1.0f, 0.0f};
    auto d = sphere.distance(a, b);
    CHECK_THAT(static_cast<double>(d), WithinAbs(std::numbers::pi / 2, 1e-4));
}

// ── Degenerate geometry ────────────────────────────────────────

TEST_CASE("Triangle with zero area", "[hardening]") {
    // Collinear points
    Triangle3 t(Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{2, 0, 0});
    CHECK_THAT(t.area(), WithinAbs(0.0, 1e-10));
}

TEST_CASE("Triangle normal of degenerate triangle", "[hardening]") {
    Triangle3 t(Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{2, 0, 0});
    auto n = t.normal();
    // Should be zero vector (degenerate), not NaN
    CHECK(std::isfinite(n[0]));
    CHECK(std::isfinite(n[1]));
    CHECK(std::isfinite(n[2]));
}

TEST_CASE("Line from zero direction fails gracefully", "[hardening]") {
    auto l = Line3::from(Vec3{0, 0, 0}, Vec3{0, 0, 0});
    CHECK_FALSE(l.has_value());
    CHECK(l.error().code == ErrorCode::DegenerateInput);
}

TEST_CASE("Hyperplane from coincident points fails", "[hardening]") {
    auto p = Plane3::from_points(Vec3{1, 1, 1}, Vec3{1, 1, 1}, Vec3{1, 1, 1});
    CHECK_FALSE(p.has_value());
}

// ── Epsilon system ─────────────────────────────────────────────

TEST_CASE("epsilon<float> > epsilon<double>", "[hardening]") {
    CHECK(epsilon<float>() > epsilon<double>());
}

TEST_CASE("near_zero works", "[hardening]") {
    CHECK(near_zero(0.0));
    CHECK(near_zero(1e-15));
    CHECK_FALSE(near_zero(1.0));
}

// ── SE3 known-good test (not just roundtrip) ───────────────────

TEST_CASE("SE3 exp known-good: pure translation", "[hardening]") {
    algebra::SE3 se3;
    // Pure translation (zero rotation): exp should give identity + translation
    Vec<double, 6> xi{0, 0, 0, 1, 2, 3};
    auto T = se3.exp(xi);
    auto p = se3.act(T, Vec3{0, 0, 0});
    CHECK_THAT(p[0], WithinAbs(1.0, 1e-10));
    CHECK_THAT(p[1], WithinAbs(2.0, 1e-10));
    CHECK_THAT(p[2], WithinAbs(3.0, 1e-10));
}

TEST_CASE("SE3 exp known-good: 90deg Z rotation + translation", "[hardening]") {
    algebra::SE3 se3;
    // 90 degrees around Z, then check act on (1,0,0)
    auto R = se3.so3.rz(std::numbers::pi / 2);
    auto T = algebra::SE3::from_Rt(R, Vec3{5, 0, 0});
    auto p = se3.act(T, Vec3{1, 0, 0});
    // R*(1,0,0) = (0,1,0), then + (5,0,0) = (5,1,0)
    CHECK_THAT(p[0], WithinAbs(5.0, 1e-8));
    CHECK_THAT(p[1], WithinAbs(1.0, 1e-8));
    CHECK_THAT(p[2], WithinAbs(0.0, 1e-8));
}

// ── ProductSpace exp/log ───────────────────────────────────────

TEST_CASE("ProductSpace E1 x E1 exp/log", "[hardening]") {
    ProductSpace<Euclidean<1>, Euclidean<1>> r2;
    Vec<double, 2> p{1, 2};
    Vec<double, 2> q{4, 6};
    auto v = r2.log_map(p, q);
    auto recovered = r2.exp_map(p, v, 1.0);
    CHECK_THAT(r2.distance(q, recovered), WithinAbs(0.0, 1e-10));
}

// ── Matrix singular inverse ────────────────────────────────────

TEST_CASE("Matrix singular inverse returns error", "[hardening]") {
    Mat2 singular;
    singular(0, 0) = 1; singular(0, 1) = 2;
    singular(1, 0) = 2; singular(1, 1) = 4; // det = 0
    auto inv = singular.inverse();
    REQUIRE_FALSE(inv.has_value());
    CHECK(inv.error().code == ErrorCode::SingularMatrix);
}

// ── Segment-Segment 3D distance ────────────────────────────────

TEST_CASE("Segment-Segment 3D parallel", "[hardening]") {
    Segment3 s1{Vec3{0, 0, 0}, Vec3{1, 0, 0}};
    Segment3 s2{Vec3{0, 0, 3}, Vec3{1, 0, 3}};
    CHECK_THAT(geometry::distance(s1, s2), WithinAbs(3.0, 1e-10));
}

TEST_CASE("Segment-Segment 3D skew", "[hardening]") {
    Segment3 s1{Vec3{0, 0, 0}, Vec3{1, 0, 0}};
    Segment3 s2{Vec3{0.5, 0, 2}, Vec3{0.5, 1, 2}}; // above midpoint of s1, perpendicular
    CHECK_THAT(geometry::distance(s1, s2), WithinAbs(2.0, 1e-10));
}
