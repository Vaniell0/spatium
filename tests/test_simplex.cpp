#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/simplex.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

TEST_CASE("Simplex<3,0> (point)", "[simplex]") {
    Simplex<3, 0> point(Vec3{1.0, 2.0, 3.0});
    CHECK(point.centroid() == Vec3{1.0, 2.0, 3.0});
    CHECK(point.measure() == 1.0);
}

TEST_CASE("Simplex<3,1> (segment)", "[simplex]") {
    Simplex<3, 1> seg(Vec3{0.0, 0.0, 0.0}, Vec3{3.0, 4.0, 0.0});
    CHECK_THAT(seg.measure(), WithinAbs(5.0, 1e-12));
    CHECK(seg.centroid() == Vec3{1.5, 2.0, 0.0});
}

TEST_CASE("Simplex<3,2> (triangle) measure matches Triangle", "[simplex]") {
    Simplex<3, 2> tri(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0});
    CHECK_THAT(tri.measure(), WithinAbs(0.5, 1e-12));
}

TEST_CASE("Simplex<3,3> (tetrahedron) volume", "[simplex]") {
    // Unit tetrahedron: volume = 1/6
    Simplex<3, 3> tet(
        Vec3{0.0, 0.0, 0.0},
        Vec3{1.0, 0.0, 0.0},
        Vec3{0.0, 1.0, 0.0},
        Vec3{0.0, 0.0, 1.0}
    );
    CHECK_THAT(tet.measure(), WithinAbs(1.0 / 6.0, 1e-12));
}

TEST_CASE("Simplex face", "[simplex]") {
    Simplex<3, 2> tri(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0});

    // face(0) = remove vertex 0 → segment from v1 to v2
    auto f0 = tri.face(0);
    CHECK(f0[0] == Vec3{1.0, 0.0, 0.0});
    CHECK(f0[1] == Vec3{0.0, 1.0, 0.0});
}

TEST_CASE("Simplex barycentric at vertices", "[simplex]") {
    Simplex<3, 2> tri(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0});

    auto b0 = tri.barycentric(tri[0]);
    CHECK_THAT(b0[0], WithinAbs(1.0, 1e-10));
    CHECK_THAT(b0[1], WithinAbs(0.0, 1e-10));
    CHECK_THAT(b0[2], WithinAbs(0.0, 1e-10));
}

TEST_CASE("Simplex contains", "[simplex]") {
    Simplex<3, 2> tri(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0});
    CHECK(tri.contains(Vec3{0.1, 0.1, 0.0}));
    CHECK_FALSE(tri.contains(Vec3{1.0, 1.0, 0.0}));
}

TEST_CASE("Simplex bounding box", "[simplex]") {
    Simplex<3, 3> tet(
        Vec3{0.0, 0.0, 0.0},
        Vec3{1.0, 0.0, 0.0},
        Vec3{0.0, 1.0, 0.0},
        Vec3{0.0, 0.0, 1.0}
    );
    auto bb = tet.bounding_box();
    CHECK(bb.min_corner == Vec3{0.0, 0.0, 0.0});
    CHECK(bb.max_corner == Vec3{1.0, 1.0, 1.0});
}

TEST_CASE("Simplex concept satisfaction", "[simplex]") {
    static_assert(Shape<Simplex<3, 2>>);
    static_assert(ClosedShape<Simplex<3, 2>>);
    static_assert(Measurable<Simplex<3, 2>>);
    static_assert(Bounded<Simplex<3, 2>>);
    SUCCEED();
}
