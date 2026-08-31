#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/vector.hpp>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("Vec construction", "[vector]") {
    Vec3 v{1.0, 2.0, 3.0};
    CHECK(v[0] == 1.0);
    CHECK(v[1] == 2.0);
    CHECK(v[2] == 3.0);

    Vec3 zero;
    CHECK(zero[0] == 0.0);
    CHECK(zero[1] == 0.0);
    CHECK(zero[2] == 0.0);
}

TEST_CASE("Vec arithmetic", "[vector]") {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};

    SECTION("addition") {
        auto c = a + b;
        CHECK(c[0] == 5.0);
        CHECK(c[1] == 7.0);
        CHECK(c[2] == 9.0);
    }

    SECTION("subtraction") {
        auto c = b - a;
        CHECK(c[0] == 3.0);
        CHECK(c[1] == 3.0);
        CHECK(c[2] == 3.0);
    }

    SECTION("negation") {
        auto c = -a;
        CHECK(c[0] == -1.0);
        CHECK(c[1] == -2.0);
        CHECK(c[2] == -3.0);
    }

    SECTION("scalar multiplication") {
        auto c = a * 2.0;
        CHECK(c[0] == 2.0);
        CHECK(c[1] == 4.0);
        CHECK(c[2] == 6.0);

        auto d = 3.0 * a;
        CHECK(d[0] == 3.0);
        CHECK(d[1] == 6.0);
        CHECK(d[2] == 9.0);
    }

    SECTION("scalar division") {
        auto c = a / 2.0;
        CHECK(c[0] == 0.5);
        CHECK(c[1] == 1.0);
        CHECK(c[2] == 1.5);
    }
}

TEST_CASE("Vec dot product", "[vector]") {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};
    CHECK(a.dot(b) == 32.0); // 4 + 10 + 18
}

TEST_CASE("Vec norm", "[vector]") {
    Vec3 v{3.0, 4.0, 0.0};
    CHECK(v.norm_squared() == 25.0);
    CHECK_THAT(v.norm(), WithinAbs(5.0, 1e-12));
}

TEST_CASE("Vec cross product", "[vector]") {
    Vec3 x{1.0, 0.0, 0.0};
    Vec3 y{0.0, 1.0, 0.0};
    auto z = x.cross(y);
    CHECK(z[0] == 0.0);
    CHECK(z[1] == 0.0);
    CHECK(z[2] == 1.0);
}

TEST_CASE("Vec normalized", "[vector]") {
    Vec3 v{3.0, 4.0, 0.0};
    auto u = v.normalized();
    CHECK_THAT(u.norm(), WithinAbs(1.0, 1e-12));
    CHECK_THAT(u[0], WithinAbs(0.6, 1e-12));
    CHECK_THAT(u[1], WithinAbs(0.8, 1e-12));
}

TEST_CASE("Vec constexpr arithmetic", "[vector]") {
    constexpr Vec3 a{1.0, 2.0, 3.0};
    constexpr Vec3 b{4.0, 5.0, 6.0};
    constexpr auto c = a + b;
    static_assert(c[0] == 5.0);
    static_assert(c[1] == 7.0);
    static_assert(c[2] == 9.0);

    constexpr auto d = a.dot(b);
    static_assert(d == 32.0);

    SUCCEED();
}

TEST_CASE("Vec equality", "[vector]") {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{1.0, 2.0, 3.0};
    Vec3 c{1.0, 2.0, 4.0};
    CHECK(a == b);
    CHECK(a != c);
}

TEST_CASE("Vec2 and Vec4", "[vector]") {
    Vec2 v2{1.0, 2.0};
    CHECK(v2.norm_squared() == 5.0);

    Vec4 v4{1.0, 1.0, 1.0, 1.0};
    CHECK(v4.norm_squared() == 4.0);
}

// ── Compile-time arithmetic guard ─────────────────────────────
// These static_asserts exercise the constexpr/`if consteval`
// branches in Vec (operator+/-/*///, dot, cross, norm_squared,
// lerp, reflect, abs, min, max, clamp, unit, equality). A
// regression that accidentally drops constexpr from any of these
// kernels — e.g. by introducing a non-constexpr helper or pulling
// in std::sqrt — collapses compilation here instead of silently
// shifting work to runtime.

namespace {

constexpr Vec3 ca{1.0, 2.0, 3.0};
constexpr Vec3 cb{4.0, 5.0, 6.0};

// dot / norm_squared
static_assert(ca.dot(cb) == 32.0);
static_assert(ca.dot(Vec3{1.0, 1.0, 1.0}) == 6.0);
static_assert(ca.norm_squared() == 14.0);

// cross — right-handed: x × y = z
static_assert(ca.cross(Vec3{0.0, 1.0, 0.0})[2] == 1.0);
static_assert(Vec3{1.0, 0.0, 0.0}.cross(Vec3{0.0, 1.0, 0.0})
                  == Vec3{0.0, 0.0, 1.0});
static_assert(Vec3{0.0, 1.0, 0.0}.cross(Vec3{1.0, 0.0, 0.0})
                  == Vec3{0.0, 0.0, -1.0});

// arithmetic — assignment from Vec to Vec must be constexpr
constexpr Vec3 cs = Vec3{ca + cb};
static_assert(cs == Vec3{5.0, 7.0, 9.0});
constexpr Vec3 csub = Vec3{cb - ca};
static_assert(csub == Vec3{3.0, 3.0, 3.0});
constexpr Vec3 cneg = Vec3{-ca};
static_assert(cneg == Vec3{-1.0, -2.0, -3.0});
constexpr Vec3 cmul = Vec3{ca * 2.0};
static_assert(cmul == Vec3{2.0, 4.0, 6.0});
constexpr Vec3 cdiv = Vec3{cb / 2.0};
static_assert(cdiv == Vec3{2.0, 2.5, 3.0});

// lerp — midpoint of (1,2,3) and (3,4,5) is (2,3,4)
static_assert(Vec3{1.0, 2.0, 3.0}.lerp(Vec3{3.0, 4.0, 5.0}, 0.5)
                  == Vec3{2.0, 3.0, 4.0});
// lerp endpoints
static_assert(ca.lerp(cb, 0.0) == ca);
static_assert(ca.lerp(cb, 1.0) == cb);

// reflect — bouncing (1,-1,0) off the +Y unit normal yields (1,1,0)
static_assert(Vec3{1.0, -1.0, 0.0}.reflect(Vec3{0.0, 1.0, 0.0})
                  == Vec3{1.0, 1.0, 0.0});

// component-wise: abs / min / max / clamp
static_assert(Vec3{-1.0, 2.0, -3.0}.abs() == Vec3{1.0, 2.0, 3.0});
static_assert(Vec3::min(ca, cb) == ca);
static_assert(Vec3::max(ca, cb) == cb);
static_assert(Vec3{0.0, 5.0, 10.0}.clamp(Vec3{1.0, 1.0, 1.0},
                                         Vec3{8.0, 8.0, 8.0})
                  == Vec3{1.0, 5.0, 8.0});

// factories
static_assert(Vec3::zero() == Vec3{0.0, 0.0, 0.0});
static_assert(Vec3::unit(0) == Vec3{1.0, 0.0, 0.0});
static_assert(Vec3::unit(1) == Vec3{0.0, 1.0, 0.0});
static_assert(Vec3::unit(2) == Vec3{0.0, 0.0, 1.0});

// compound assignment must keep its constexpr if-consteval branch
constexpr Vec3 caccum = []{
    Vec3 v{1.0, 2.0, 3.0};
    v += Vec3{4.0, 5.0, 6.0};
    v -= Vec3{2.0, 2.0, 2.0};
    v *= 2.0;
    v /= 2.0;
    return v;
}();
static_assert(caccum == Vec3{3.0, 5.0, 7.0});

// equality
static_assert(ca == Vec3{1.0, 2.0, 3.0});
static_assert(ca != cb);

} // namespace
