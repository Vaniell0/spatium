#include <catch2/catch_test_macros.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/algebra/vector.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/geometry/triangle.hpp>
#include <limits>

using namespace spatium;
using namespace spatium::geometry;

TEST_CASE("Real50 satisfies Scalar", "[precision]") {
    static_assert(Scalar<Real50>);
    static_assert(Scalar<Real100>);
    SUCCEED();
}

TEST_CASE("Vec with Real50", "[precision]") {
    Vec<Real50, 3> a{Real50{1}, Real50{2}, Real50{3}};
    Vec<Real50, 3> b{Real50{4}, Real50{5}, Real50{6}};
    auto c = a + b;
    CHECK(c[0] == Real50{5});
    CHECK(c[1] == Real50{7});
    CHECK(c[2] == Real50{9});
}

TEST_CASE("Euclidean with Real50", "[precision]") {
    Euclidean<3, Real50> space;
    Vec<Real50, 3> a{Real50{0}, Real50{0}, Real50{0}};
    Vec<Real50, 3> b{Real50{3}, Real50{4}, Real50{0}};
    auto d = space.distance(a, b);
    CHECK(d == Real50{5});
}

TEST_CASE("Triangle with Real50 — high precision area", "[precision]") {
    // Right triangle with known area
    Triangle<3, Real50> tri(
        Vec<Real50, 3>{Real50{0}, Real50{0}, Real50{0}},
        Vec<Real50, 3>{Real50{1}, Real50{0}, Real50{0}},
        Vec<Real50, 3>{Real50{0}, Real50{1}, Real50{0}}
    );
    auto area = tri.area();
    // 0.5 with 50 digits of precision
    CHECK(area == Real50{"0.5"});
}

TEST_CASE("Real<N> generalizes Real50/Real100 to any digit count", "[precision]") {
    static_assert(std::is_same_v<Real50, Real<50>>);
    static_assert(std::is_same_v<Real100, Real<100>>);
    static_assert(Scalar<Real<17>>);

    using T = Real<17>;
    T a{1};
    T b{2};
    CHECK(a + b == T{3});
}

TEST_CASE("Real<N> reaches accuracy no 50- or 100-digit preset could represent", "[precision]") {
    // std::numeric_limits<>::digits10 is the exact, guaranteed contract --
    // proves digit count is a real compile-time parameter, not a choice
    // between two fixed presets.
    static_assert(std::numeric_limits<Real50>::digits10 == 50);
    static_assert(std::numeric_limits<Real100>::digits10 == 100);
    static_assert(std::numeric_limits<Real<300>>::digits10 == 300);

    // A runtime demonstration, not just the type-level contract: sqrt(2)
    // computed at Real<300> is accurate to a relative error near 1e-300 --
    // meaningless to even ask of a 50- or 100-digit type, since there
    // aren't that many digits to be accurate WITH. (cpp_dec_float's
    // internal guard digits mean naively probing "does adding 1e-60 to 1
    // get rounded away" is NOT a reliable way to test this -- the digit
    // count contract above and an actual high-precision computation are.)
    using T = Real<300>;
    using std::sqrt; using std::abs; // ADL: Real<N> provides its own
    T two{2};
    T root = sqrt(two);
    T residual = abs(root * root - two);
    CHECK(residual < T{"1e-290"});
}

TEST_CASE("Vec with an arbitrary (non-preset) Real<N> digit count", "[precision]") {
    using T = Real<37>;
    Vec<T, 3> a{T{1}, T{2}, T{3}};
    Vec<T, 3> b{T{4}, T{5}, T{6}};
    auto c = a + b;
    CHECK(c[0] == T{5});
    CHECK(c[1] == T{7});
    CHECK(c[2] == T{9});
}

TEST_CASE("High precision pi approximation via polygon", "[precision]") {
    // Inscribed regular N-gon in unit circle: area = N/2 * sin(2π/N)
    // As N→∞, area → π
    // With Real100, we can compute this very precisely

    using T = Real100;
    using boost::multiprecision::asin;

    // π via asin(1)*2
    T pi = asin(T{1}) * T{2};

    // 10000-gon inscribed in unit circle
    int n = 10000;
    T angle = T{2} * pi / T{n};
    T polygon_area = T{n} / T{2} * sin(angle);

    // Should be close to π
    T error = abs(polygon_area - pi);
    CHECK(error < T{"0.00001"});
}
