#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/matrix.hpp>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("Matrix identity", "[matrix]") {
    auto I = Mat3::identity();
    CHECK(I(0, 0) == 1.0);
    CHECK(I(1, 1) == 1.0);
    CHECK(I(2, 2) == 1.0);
    CHECK(I(0, 1) == 0.0);
}

TEST_CASE("Matrix arithmetic", "[matrix]") {
    Mat2 a;
    a(0, 0) = 1; a(0, 1) = 2;
    a(1, 0) = 3; a(1, 1) = 4;

    Mat2 b;
    b(0, 0) = 5; b(0, 1) = 6;
    b(1, 0) = 7; b(1, 1) = 8;

    auto c = a + b;
    CHECK(c(0, 0) == 6);
    CHECK(c(1, 1) == 12);

    auto s = a * 2.0;
    CHECK(s(0, 0) == 2);
    CHECK(s(1, 1) == 8);
}

TEST_CASE("Matrix multiply", "[matrix]") {
    auto I = Mat3::identity();
    Mat3 a;
    a(0, 0) = 1; a(0, 1) = 2; a(0, 2) = 3;
    a(1, 0) = 4; a(1, 1) = 5; a(1, 2) = 6;
    a(2, 0) = 7; a(2, 1) = 8; a(2, 2) = 9;

    auto b = I * a;
    CHECK(b(0, 0) == 1);
    CHECK(b(2, 2) == 9);
}

TEST_CASE("Matrix-vector multiply", "[matrix]") {
    auto I = Mat3::identity();
    Vec3 v{1, 2, 3};
    auto r = I * v;
    CHECK(r == v);

    Mat2 m;
    m(0, 0) = 2; m(0, 1) = 0;
    m(1, 0) = 0; m(1, 1) = 3;
    auto r2 = m * Vec2{1, 1};
    CHECK(r2 == Vec2{2, 3});
}

TEST_CASE("Matrix transpose", "[matrix]") {
    Mat2 m;
    m(0, 0) = 1; m(0, 1) = 2;
    m(1, 0) = 3; m(1, 1) = 4;
    auto t = m.transpose();
    CHECK(t(0, 0) == 1);
    CHECK(t(0, 1) == 3);
    CHECK(t(1, 0) == 2);
    CHECK(t(1, 1) == 4);
}

TEST_CASE("Matrix determinant 2x2", "[matrix]") {
    Mat2 m;
    m(0, 0) = 1; m(0, 1) = 2;
    m(1, 0) = 3; m(1, 1) = 4;
    CHECK_THAT(m.determinant(), WithinAbs(-2.0, 1e-12));
}

TEST_CASE("Matrix determinant 3x3", "[matrix]") {
    Mat3 m;
    m(0, 0) = 6; m(0, 1) = 1; m(0, 2) = 1;
    m(1, 0) = 4; m(1, 1) = -2; m(1, 2) = 5;
    m(2, 0) = 2; m(2, 1) = 8; m(2, 2) = 7;
    CHECK_THAT(m.determinant(), WithinAbs(-306.0, 1e-10));
}

TEST_CASE("Matrix determinant 4x4", "[matrix]") {
    Mat4 m;
    m(0,0)=1; m(0,1)=2; m(0,2)=3; m(0,3)=4;
    m(1,0)=5; m(1,1)=6; m(1,2)=7; m(1,3)=8;
    m(2,0)=2; m(2,1)=6; m(2,2)=4; m(2,3)=8;
    m(3,0)=3; m(3,1)=1; m(3,2)=1; m(3,3)=2;
    CHECK_THAT(m.determinant(), WithinAbs(72.0, 1e-8));
}

TEST_CASE("Matrix inverse 2x2", "[matrix]") {
    Mat2 m;
    m(0, 0) = 4; m(0, 1) = 7;
    m(1, 0) = 2; m(1, 1) = 6;
    auto inv = m.inverse();
    REQUIRE(inv.has_value());
    auto prod = m * *inv;
    auto I = Mat2::identity();
    for (int i = 0; i < 4; ++i)
        CHECK_THAT(prod.data[i], WithinAbs(I.data[i], 1e-12));
}

TEST_CASE("Matrix inverse 3x3", "[matrix]") {
    Mat3 m;
    m(0,0)=1; m(0,1)=2; m(0,2)=3;
    m(1,0)=0; m(1,1)=1; m(1,2)=4;
    m(2,0)=5; m(2,1)=6; m(2,2)=0;
    auto inv = m.inverse();
    REQUIRE(inv.has_value());
    auto prod = m * *inv;
    auto I = Mat3::identity();
    for (int i = 0; i < 9; ++i)
        CHECK_THAT(prod.data[i], WithinAbs(I.data[i], 1e-10));
}

TEST_CASE("Matrix col/row access", "[matrix]") {
    Mat2 m;
    m(0, 0) = 1; m(0, 1) = 2;
    m(1, 0) = 3; m(1, 1) = 4;
    CHECK(m.col(0) == Vec2{1, 3});
    CHECK(m.row(0) == Vec2{1, 2});
}

TEST_CASE("Matrix format", "[matrix]") {
    auto I = Mat2::identity();
    auto s = std::format("{}", I);
    CHECK(s.find("1") != std::string::npos);
}

// ── Compile-time arithmetic guard ─────────────────────────────
// Pins down Matrix's constexpr surface: identity, element access,
// addition, scalar multiplication, matrix*matrix, matrix*vector,
// transpose, trace, determinant (1x1, 2x2, 3x3), col/row access,
// equality. A regression that drops constexpr from any kernel —
// or that pulls a non-constexpr helper into them — fails to
// compile here instead of slipping through to runtime.

namespace {

// Build a 2x2 from scratch via a constexpr lambda. operator()(r,c)
// is mutating, so wrap construction in an immediately-invoked
// expression that stays inside constant evaluation.
constexpr Mat2 make_2x2(double a00, double a01, double a10, double a11) {
    Mat2 m;
    m(0, 0) = a00; m(0, 1) = a01;
    m(1, 0) = a10; m(1, 1) = a11;
    return m;
}

constexpr Mat3 make_3x3_diag(double d) {
    Mat3 m;
    m(0, 0) = d; m(1, 1) = d; m(2, 2) = d;
    return m;
}

constexpr auto Iden3 = Mat3::identity();
static_assert(Iden3(0, 0) == 1.0);
static_assert(Iden3(1, 1) == 1.0);
static_assert(Iden3(2, 2) == 1.0);
static_assert(Iden3(0, 1) == 0.0);
static_assert(Iden3(2, 0) == 0.0);

// Identity multiplied by itself is itself.
static_assert(Iden3 * Iden3 == Iden3);

// Identity acts as left/right neutral on a vector.
constexpr Vec3 mv = Iden3 * Vec3{1.0, 2.0, 3.0};
static_assert(mv == Vec3{1.0, 2.0, 3.0});

// Sum of two 2x2: [1,2;3,4] + [5,6;7,8] = [6,8;10,12]
constexpr auto m_a = make_2x2(1.0, 2.0, 3.0, 4.0);
constexpr auto m_b = make_2x2(5.0, 6.0, 7.0, 8.0);
constexpr auto m_sum = m_a + m_b;
static_assert(m_sum(0, 0) == 6.0 && m_sum(0, 1) == 8.0);
static_assert(m_sum(1, 0) == 10.0 && m_sum(1, 1) == 12.0);

// Difference, scalar multiplication.
constexpr auto m_diff = m_b - m_a;
static_assert(m_diff == make_2x2(4.0, 4.0, 4.0, 4.0));
constexpr auto m_scaled = m_a * 3.0;
static_assert(m_scaled == make_2x2(3.0, 6.0, 9.0, 12.0));
static_assert(2.0 * m_a == make_2x2(2.0, 4.0, 6.0, 8.0));

// 2x2 * 2x2: [1,2;3,4] * [5,6;7,8] = [19,22;43,50]
constexpr auto m_prod = m_a * m_b;
static_assert(m_prod == make_2x2(19.0, 22.0, 43.0, 50.0));

// 2x2 * Vec2: [1,2;3,4] * (1,1) = (3,7)
constexpr Vec2 mv2 = m_a * Vec2{1.0, 1.0};
static_assert(mv2 == Vec2{3.0, 7.0});

// Transpose
constexpr auto m_t = m_a.transpose();
static_assert(m_t == make_2x2(1.0, 3.0, 2.0, 4.0));

// Trace and determinants
static_assert(m_a.trace() == 5.0);
static_assert(m_a.determinant() == -2.0);
static_assert(make_3x3_diag(2.0).determinant() == 8.0);

// 1x1 determinant
constexpr Matrix<double, 1, 1> m1 = []{
    Matrix<double, 1, 1> m;
    m(0, 0) = 7.0;
    return m;
}();
static_assert(m1.determinant() == 7.0);

// col / row extraction
static_assert(m_a.col(0) == Vec2{1.0, 3.0});
static_assert(m_a.row(0) == Vec2{1.0, 2.0});

// Equality
static_assert(m_a == make_2x2(1.0, 2.0, 3.0, 4.0));
static_assert(m_a != m_b);

} // namespace
