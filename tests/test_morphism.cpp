#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/morphism.hpp>
#include <spatium/point.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/core/error.hpp>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("Morphism basic apply", "[morphism]") {
    auto scale = morph<E3, E3>([](const Vec3& p) { return p * 2.0; });
    auto result = scale(Vec3{1.0, 2.0, 3.0});
    CHECK(result == Vec3{2.0, 4.0, 6.0});
}

TEST_CASE("Morphism on typed Point", "[morphism]") {
    auto translate = morph<E3, E3>([](const Vec3& p) { return p + Vec3{10.0, 0.0, 0.0}; });
    auto p = pt<E3>(Vec3{1.0, 2.0, 3.0});
    auto q = translate(p);
    CHECK(q.raw() == Vec3{11.0, 2.0, 3.0});
}

TEST_CASE("Morphism pipe operator", "[morphism]") {
    auto scale = morph<E3, E3>([](const Vec3& p) { return p * 2.0; });
    auto p = pt<E3>(Vec3{1.0, 2.0, 3.0});

    auto q = p | scale;
    CHECK(q.raw() == Vec3{2.0, 4.0, 6.0});
}

TEST_CASE("Morphism composition with *", "[morphism]") {
    auto f = morph<E3, E3>([](const Vec3& p) { return p * 2.0; });
    auto g = morph<E3, E3>([](const Vec3& p) { return p + Vec3{1.0, 0.0, 0.0}; });

    // g * f = first scale, then translate
    auto gf = g * f;
    auto result = gf(Vec3{1.0, 0.0, 0.0});
    CHECK(result == Vec3{3.0, 0.0, 0.0});
}

TEST_CASE("Morphism pipe composition", "[morphism]") {
    auto scale = morph<E3, E3>([](const Vec3& p) { return p * 2.0; });
    auto translate = morph<E3, E3>([](const Vec3& p) { return p + Vec3{1.0, 0.0, 0.0}; });

    // f | g = first f, then g (pipe order)
    auto pipeline = scale | translate;
    auto p = pt<E3>(Vec3{1.0, 0.0, 0.0});
    auto q = p | pipeline;
    CHECK(q.raw() == Vec3{3.0, 0.0, 0.0});
}

TEST_CASE("Morphism chain: point | f | g", "[morphism]") {
    auto scale = morph<E3, E3>([](const Vec3& p) { return p * 3.0; });
    auto negate = morph<E3, E3>([](const Vec3& p) { return -p; });

    auto p = pt<E3>(Vec3{1.0, 2.0, 3.0});
    auto q = p | scale | negate;
    CHECK(q.raw() == Vec3{-3.0, -6.0, -9.0});
}

TEST_CASE("Morphism between spaces: E3 → E2 projection", "[morphism]") {
    auto proj_xy = morph<E3, E2>([](const Vec3& p) -> Vec2 {
        return Vec2{p[0], p[1]};
    });
    auto p = pt<E3>(Vec3{1.0, 2.0, 99.0});
    auto q = p | proj_xy;
    CHECK(q.raw() == Vec2{1.0, 2.0});
}

TEST_CASE("Identity morphism", "[morphism]") {
    auto id = identity<E3>();
    auto p = pt<E3>(Vec3{1.0, 2.0, 3.0});
    auto q = p | id;
    CHECK(q == p);
}

TEST_CASE("Morphism with inverse", "[morphism]") {
    auto scale = morph<E3, E3>(
        [](const Vec3& p) { return p * 2.0; },
        [](const Vec3& p) { return p * 0.5; }
    );
    CHECK(scale.inverse.has_value());
    auto result = (*scale.inverse)(Vec3{4.0, 6.0, 8.0});
    CHECK(result == Vec3{2.0, 3.0, 4.0});
}

// ── Result<Point> | Morphism pipe-unwrap ──────────────────────

TEST_CASE("Result<Point> pipe through morphism: success", "[morphism]") {
    auto scale = morph<E3, E3>([](const Vec3& p) { return p * 2.0; });
    Result<Point<E3>> ok = pt<E3>(Vec3{1.0, 2.0, 3.0});

    auto q = ok | scale;
    REQUIRE(q.has_value());
    CHECK(q->raw() == Vec3{2.0, 4.0, 6.0});
}

TEST_CASE("Result<Point> pipe through morphism: error propagates", "[morphism]") {
    auto scale = morph<E3, E3>([](const Vec3& p) { return p * 2.0; });
    Result<Point<E3>> err = std::unexpected(Error{ErrorCode::NoIntersection, "no hit"});

    auto q = err | scale;
    CHECK_FALSE(q.has_value());
}

// ── Generic Result<T> | Fn pipe ───────────────────────────────

TEST_CASE("Result<T> generic pipe: success propagates value", "[error][pipe]") {
    Result<int> r = 41;
    auto s = r | [](int x) { return x + 1; };
    REQUIRE(s.has_value());
    CHECK(*s == 42);
}

TEST_CASE("Result<T> generic pipe: failure short-circuits", "[error][pipe]") {
    Result<int> r = std::unexpected(Error{ErrorCode::InvalidArgument, "bad"});
    bool called = false;
    auto s = r | [&called](int x) { called = true; return x + 1; };
    CHECK_FALSE(s.has_value());
    CHECK_FALSE(called);
    CHECK(s.error().code == ErrorCode::InvalidArgument);
    CHECK(s.error().message == "bad");
}

TEST_CASE("Result<T> generic pipe: void closure yields Result<void>",
          "[error][pipe]") {
    Result<int> r = 7;
    int seen = 0;
    auto s = r | [&seen](int x) { seen = x; };
    static_assert(std::is_same_v<decltype(s), Result<void>>);
    CHECK(s.has_value());
    CHECK(seen == 7);
}

TEST_CASE("Result<T> generic pipe: closure returning Result flattens",
          "[error][pipe]") {
    auto checked_div = [](int x) -> Result<int> {
        if (x == 0) return std::unexpected(Error{ErrorCode::DegenerateInput, "zero"});
        return 100 / x;
    };
    Result<int> r = 4;
    auto s = r | checked_div;
    static_assert(std::is_same_v<decltype(s), Result<int>>);
    REQUIRE(s.has_value());
    CHECK(*s == 25);

    Result<int> z = 0;
    auto t = z | checked_div;
    CHECK_FALSE(t.has_value());
    CHECK(t.error().code == ErrorCode::DegenerateInput);
}

TEST_CASE("Result<T> generic pipe: chains compose left to right",
          "[error][pipe]") {
    Result<int> r = 3;
    auto s = r
           | [](int x) { return x * 2; }
           | [](int x) { return x + 1; };
    REQUIRE(s.has_value());
    CHECK(*s == 7);
}

TEST_CASE("Result<T> generic pipe: rvalue Result passes through", "[error][pipe]") {
    auto s = Result<int>{10} | [](int x) { return x * 3; };
    REQUIRE(s.has_value());
    CHECK(*s == 30);
}

// ── unwrap throws BadResult, not std::runtime_error ───────────

TEST_CASE("unwrap throws BadResult carrying the original Error",
          "[error][unwrap]") {
    Result<int> err =
        std::unexpected(Error{ErrorCode::SingularMatrix, "det=0"});
    try {
        (void)unwrap(err);
        FAIL("expected BadResult");
    } catch (const BadResult& e) {
        CHECK(e.code() == ErrorCode::SingularMatrix);
        CHECK(e.error().message == "det=0");
        CHECK(std::string{e.what()}.find("det=0") != std::string::npos);
    }
}

TEST_CASE("unwrap returns value on success", "[error][unwrap]") {
    Result<int> ok = 17;
    CHECK(unwrap(ok) == 17);
}

TEST_CASE("BadResult is also a std::exception", "[error][unwrap]") {
    Result<int> err = std::unexpected(Error{ErrorCode::ZeroNorm, "n=0"});
    try {
        (void)unwrap(std::move(err));
        FAIL("expected BadResult");
    } catch (const std::exception& e) {
        CHECK(std::string{e.what()}.find("n=0") != std::string::npos);
    }
}
