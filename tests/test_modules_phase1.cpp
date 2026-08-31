// Phase 1 verification — consume spatium.core and spatium.algebra purely
// through the modular API. Catch2 headers MUST precede any `import`
// (modules report §2.1 — std redefinition).
#include <catch2/catch_test_macros.hpp>
#include <cmath>     // std::abs in test bodies — Catch2 doesn't pull it
#include <expected>  // std::unexpected for Result construction

import spatium.core;
import spatium.algebra;

using spatium::Vec;
using spatium::Vec3;
using spatium::Quaternion;
using spatium::Complex;
using spatium::epsilon;
using spatium::approx_equal;
using spatium::Result;
using spatium::Error;
using spatium::ErrorCode;

TEST_CASE("module spatium.core: epsilon + approx_equal", "[modules][phase1][core]") {
    auto eps = epsilon<double>();
    REQUIRE(eps > 0.0);
    REQUIRE(approx_equal(1.0, 1.0 + eps / 2.0));
    REQUIRE_FALSE(approx_equal(1.0, 1.0 + 1e-3));
}

TEST_CASE("module spatium.core: Result + Error", "[modules][phase1][core]") {
    Result<int> ok = 42;
    REQUIRE(ok.has_value());
    REQUIRE(*ok == 42);

    Result<int> bad = std::unexpected(Error{ErrorCode::DegenerateInput, "bad"});
    REQUIRE_FALSE(bad.has_value());
    REQUIRE(bad.error().code == ErrorCode::DegenerateInput);
}

TEST_CASE("module spatium.algebra: Vec arithmetic", "[modules][phase1][algebra]") {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};
    Vec3 c = a + b;
    REQUIRE(c[0] == 5.0);
    REQUIRE(c[1] == 7.0);
    REQUIRE(c[2] == 9.0);

    REQUIRE(a.dot(b) == 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0);

    auto cross = Vec3{1.0, 0.0, 0.0}.cross(Vec3{0.0, 1.0, 0.0});
    REQUIRE(cross[0] == 0.0);
    REQUIRE(cross[1] == 0.0);
    REQUIRE(cross[2] == 1.0);
}

TEST_CASE("module spatium.algebra: Quaternion roundtrip", "[modules][phase1][algebra]") {
    auto q = Quaternion<double>::from_axis_angle({0.0, 0.0, 1.0}, 1.5707963267948966);
    auto v = q.rotate(Vec3{1.0, 0.0, 0.0});
    REQUIRE(std::abs(v[0]) < 1e-9);
    REQUIRE(std::abs(v[1] - 1.0) < 1e-9);
    REQUIRE(std::abs(v[2]) < 1e-9);
}

TEST_CASE("module spatium.algebra: Complex arithmetic", "[modules][phase1][algebra]") {
    Complex<double> z{3.0, 4.0};
    REQUIRE(z.magnitude() == 5.0);
    auto z2 = z * z;
    REQUIRE(z2.re == 3.0 * 3.0 - 4.0 * 4.0);
    REQUIRE(z2.im == 2.0 * 3.0 * 4.0);
}
