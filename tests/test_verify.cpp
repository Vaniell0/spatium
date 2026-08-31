#include <catch2/catch_test_macros.hpp>
#include <spatium/core/verify.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <array>
#include <cmath>

using namespace spatium;

TEST_CASE("Verify Euclidean metric axioms", "[verify]") {
    E3 space;
    std::array samples = {
        Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0},
        Vec3{1, 1, 1}, Vec3{-2, 3, 0.5},
    };
    auto result = verify_metric(space, std::span{samples});
    if (!result) FAIL(result.failure);
    CHECK(result);
}

TEST_CASE("Verify Euclidean inner product axioms", "[verify]") {
    E3 space;
    std::array samples = {
        Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1},
        Vec3{1, 2, 3}, Vec3{-1, 0.5, 2},
    };
    auto result = verify_inner_product(space, std::span{samples});
    if (!result) FAIL(result.failure);
    CHECK(result);
}

TEST_CASE("Verify Euclidean norm consistency", "[verify]") {
    E3 space;
    std::array samples = {
        Vec3{1, 0, 0}, Vec3{3, 4, 0}, Vec3{1, 1, 1},
    };
    auto result = verify_norm_consistency(space, std::span{samples});
    if (!result) FAIL(result.failure);
    CHECK(result);
}

TEST_CASE("Verify Euclidean exp/log roundtrip", "[verify]") {
    E3 space;
    std::array samples = {
        Vec3{0, 0, 0}, Vec3{1, 2, 3}, Vec3{-1, 0, 5},
    };
    auto result = verify_exp_log(space, std::span{samples});
    if (!result) FAIL(result.failure);
    CHECK(result);
}

TEST_CASE("Verify Sphere metric axioms", "[verify]") {
    S2 sphere;
    std::array samples = {
        Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1},
        Vec3{-1, 0, 0}, Vec3{0, -1, 0},
    };
    auto result = verify_metric(sphere, std::span{samples});
    if (!result) FAIL(result.failure);
    CHECK(result);
}

TEST_CASE("Verify Sphere exp/log roundtrip", "[verify]") {
    S2 sphere;
    // Points on unit sphere, not antipodal (log_map undefined for antipodes)
    std::array samples = {
        Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1},
    };
    auto result = verify_exp_log(sphere, std::span{samples});
    if (!result) FAIL(result.failure);
    CHECK(result);
}

TEST_CASE("Verify Hyperbolic metric axioms", "[verify]") {
    H2 hyp;
    std::array samples = {
        H2::origin(),
        Vec3{std::cosh(1.0), std::sinh(1.0), 0.0},
        Vec3{std::cosh(0.5), 0.0, std::sinh(0.5)},
        Vec3{std::cosh(2.0), std::sinh(2.0) * 0.6, std::sinh(2.0) * 0.8},
    };
    auto result = verify_metric(hyp, std::span{samples}, 1e-6);
    if (!result) FAIL(result.failure);
    CHECK(result);
}

TEST_CASE("Verify Hyperbolic exp/log roundtrip", "[verify]") {
    H2 hyp;
    std::array samples = {
        H2::origin(),
        Vec3{std::cosh(1.0), std::sinh(1.0), 0.0},
        Vec3{std::cosh(0.5), 0.0, std::sinh(0.5)},
    };
    auto result = verify_exp_log(hyp, std::span{samples}, 1e-5);
    if (!result) FAIL(result.failure);
    CHECK(result);
}
