// Dimension as a value (core/dimension.hpp).

#include <catch2/catch_test_macros.hpp>
#include <spatium/core/concepts.hpp>
#include <spatium/core/dimension.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/l2.hpp>
#include <spatium/spaces/product.hpp>
#include <spatium/spaces/sphere.hpp>

using namespace spatium;

TEST_CASE("a number is a finite dimension, the old sentinels are the other two kinds", "[dimension]") {
    static_assert(Dimension(3).is_finite());
    static_assert(Dimension(kDynamic).is_dynamic());
    static_assert(Dimension(kInfinite).is_infinite());
    static_assert(Dimension::infinite().is_infinite() && !Dimension::infinite().is_finite());
    static_assert(Dimension::dynamic().is_dynamic());
    // read back as the number, so every comparison written against kDynamic still reads as it did
    static_assert(static_cast<std::size_t>(Dimension(3)) == 3);
    static_assert(Dimension::dynamic() == kDynamic || std::size_t(Dimension::dynamic()) == kDynamic);
    constexpr std::size_t as_old = Dimension(7);
    static_assert(as_old == 7);
    CHECK(Dimension(3).kind() == Dimension::Kind::Finite);
}

TEST_CASE("sums mean what they say: infinite absorbs, dynamic absorbs finite, finite adds", "[dimension]") {
    static_assert((Dimension(2) + Dimension(3)).is_finite() && std::size_t(Dimension(2) + Dimension(3)) == 5);
    static_assert((Dimension(2) + Dimension::infinite()).is_infinite());
    static_assert((Dimension::dynamic() + Dimension(2)).is_dynamic());
    static_assert((Dimension::dynamic() + Dimension::infinite()).is_infinite());
    CHECK(true);
}

TEST_CASE("the spaces written before it keep their dimension, and the new one has none to count", "[dimension]") {
    static_assert(Dimension(Sphere<2>::dimension).is_finite() && std::size_t(Sphere<2>::dimension) == 2);
    static_assert(EuclideanSpace<Euclidean<3>>);
    static_assert(!EuclideanSpace<L2Interval<double>>);
    static_assert(HilbertSpace<L2Interval<double>>);
    static_assert(L2Interval<double>::dimension.is_infinite());
    // a product with a function space is infinite-dimensional; with finite ones it adds
    static_assert(Dimension(ProductSpace<Sphere<2>, Euclidean<3>>::dimension).is_finite());
    static_assert(std::size_t(ProductSpace<Sphere<2>, Euclidean<3>>::dimension) == 5);
    CHECK(true);
}
