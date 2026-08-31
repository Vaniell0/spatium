// Phase 5 verification — spatial / discrete / io via `import`.
#include <catch2/catch_test_macros.hpp>
#include <sstream>

import spatium.core;
import spatium.algebra;
import spatium.spaces;
import spatium.geometry;
import spatium.spatial;
import spatium.discrete;
import spatium.io;

using spatium::Vec3;
using spatium::E3;
using spatium::geometry::tri;
using spatium::geometry::Triangle;
using spatium::spatial::BVH;
using spatium::discrete::FiniteSet;

TEST_CASE("module spatium.spatial: BVH ray hit", "[modules][phase5][spatial]") {
    std::vector<Triangle<3, double>> tris = {
        tri(Vec3{-1, -1, 0}, Vec3{1, -1, 0}, Vec3{0, 1, 0}),
    };
    auto bvh = BVH<Triangle<3, double>>::build(tris);
    spatium::geometry::Ray<3, double> r{Vec3{0, 0, -1}, Vec3{0, 0, 1}};
    auto hit = bvh.ray_cast(r);
    REQUIRE(hit.has_value());
}

TEST_CASE("module spatium.discrete: FiniteSet union/intersect", "[modules][phase5][discrete]") {
    FiniteSet<int> a{1, 2, 3, 4};
    FiniteSet<int> b{3, 4, 5, 6};
    auto u = a + b;
    REQUIRE(u.size() == 6);
    auto i = a & b;
    REQUIRE(i.size() == 2);
}

TEST_CASE("module spatium.io: Table renders headers and rows", "[modules][phase5][io]") {
    spatium::io::Table t("Key", "Value");
    t.row(std::string("alpha"), std::string("1"));
    t.row(std::string("beta"), std::string("2"));
    std::ostringstream os;
    t.print(os);
    auto s = os.str();
    REQUIRE(s.find("alpha") != std::string::npos);
    REQUIRE(s.find("beta") != std::string::npos);
    REQUIRE(s.find("Key") != std::string::npos);
}
