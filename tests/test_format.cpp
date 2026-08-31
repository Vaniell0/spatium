#include <catch2/catch_test_macros.hpp>
#include <spatium/geometry/format.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <format>
#include <sstream>

using namespace spatium;
using namespace spatium::geometry;

TEST_CASE("Vec format", "[format]") {
    auto s = std::format("{}", Vec3{1.0, 2.0, 3.0});
    CHECK(s == "(1, 2, 3)");
}

TEST_CASE("Triangle format", "[format]") {
    auto t = tri(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0});
    auto s = std::format("{}", t);
    CHECK(s.starts_with("\u25b3["));
    CHECK(s.ends_with("]"));
}

TEST_CASE("Segment format", "[format]") {
    auto s = std::format("{}", Segment3{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0}});
    CHECK(s.find("\u2014") != std::string::npos);
}

TEST_CASE("Box format", "[format]") {
    auto s = std::format("{}", Box3{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0}});
    CHECK(s.starts_with("\u25a1["));
}

TEST_CASE("Hyperplane format", "[format]") {
    auto p = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 5.0});
    auto s = std::format("{}", p);
    CHECK(s.find("n=") != std::string::npos);
    CHECK(s.find("d=") != std::string::npos);
}

TEST_CASE("Line format", "[format]") {
    auto l = *Line3::from(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
    auto s = std::format("{}", l);
    CHECK(s.find("dir") != std::string::npos);
}

TEST_CASE("Ray format", "[format]") {
    auto r = *Ray3::from(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
    auto s = std::format("{}", r);
    CHECK(s.find("\u2192") != std::string::npos);
}

TEST_CASE("Circle format", "[format]") {
    auto s = std::format("{}", Circle3{Vec3{0.0, 0.0, 0.0}, 5.0, Vec3{0.0, 0.0, 1.0}});
    CHECK(s.find("r=") != std::string::npos);
}

TEST_CASE("Point format", "[format]") {
    auto p = pt<E3>(Vec3{1.0, 2.0, 3.0});
    auto s = std::format("{}", p);
    CHECK(s == "P(1, 2, 3)");
}

TEST_CASE("Stream output for Triangle", "[format]") {
    auto t = tri(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0});
    std::ostringstream oss;
    oss << t;
    CHECK(oss.str().starts_with("\u25b3["));
}

TEST_CASE("Polygon format", "[format]") {
    Polygon<2> p{{Vec2{0.0, 0.0}, Vec2{1.0, 0.0}, Vec2{0.5, 1.0}}};
    auto s = std::format("{}", p);
    CHECK(s.starts_with("Polygon{"));
    CHECK(s.ends_with("}"));
    CHECK(s.find("(0, 0)") != std::string::npos);
    CHECK(s.find("(1, 0)") != std::string::npos);
    CHECK(s.find("(0.5, 1)") != std::string::npos);
}

TEST_CASE("Disk format", "[format]") {
    Disk<3> d{Circle<3>{Vec3{1.0, 2.0, 3.0}, 4.0, Vec3{0.0, 0.0, 1.0}}};
    auto s = std::format("{}", d);
    CHECK(s.starts_with("Disk{"));
    CHECK(s.find("center=(1, 2, 3)") != std::string::npos);
    CHECK(s.find("r=4") != std::string::npos);
}

TEST_CASE("Simplex format", "[format]") {
    Simplex<3, 2> s{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    auto str = std::format("{}", s);
    CHECK(str.starts_with("Simplex2{"));
    CHECK(str.find("(0, 0, 0)") != std::string::npos);
    CHECK(str.find("(1, 0, 0)") != std::string::npos);
}

TEST_CASE("Quadric format", "[format]") {
    auto q = Quadric<double>::sphere(2.0);
    auto s = std::format("{}", q);
    CHECK(s.starts_with("Quadric{"));
    CHECK(s.find("Qxx=1") != std::string::npos);
    CHECK(s.find("Qww=-4") != std::string::npos);
}

TEST_CASE("Torus format", "[format]") {
    Torus<double> t{Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}, 1.5, 0.4};
    auto s = std::format("{}", t);
    CHECK(s.starts_with("Torus{"));
    CHECK(s.find("center=(0, 0, 0)") != std::string::npos);
    CHECK(s.find("axis=(0, 0, 1)") != std::string::npos);
    CHECK(s.find("R=1.5") != std::string::npos);
    CHECK(s.find("r=0.4") != std::string::npos);
}
