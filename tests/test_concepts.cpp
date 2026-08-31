#include <catch2/catch_test_macros.hpp>
#include <spatium/core/concepts.hpp>
#include <spatium/algebra/vector.hpp>

using namespace spatium;

// A minimal space that satisfies only Set
struct MinimalSet {
    using ScalarType = double;
    using PointType = Vec<double, 2>;
    static constexpr std::size_t dimension = 2;
};

// A space that satisfies MetricSpace but not VectorSpace
struct DiscreteMetric {
    using ScalarType = double;
    using PointType = Vec<double, 2>;
    static constexpr std::size_t dimension = 2;

    constexpr bool contains(const PointType&) const { return true; }

    constexpr ScalarType distance(const PointType& a, const PointType& b) const {
        return (a == b) ? 0.0 : 1.0;
    }
};

TEST_CASE("Scalar concept", "[concepts]") {
    static_assert(Scalar<double>);
    static_assert(Scalar<float>);
    static_assert(Scalar<int>);
    static_assert(!Scalar<std::string>);
    SUCCEED();
}

TEST_CASE("Set concept", "[concepts]") {
    static_assert(Set<MinimalSet>);
    SUCCEED();
}

TEST_CASE("MetricSpace without VectorSpace", "[concepts]") {
    static_assert(MetricSpace<DiscreteMetric>);
    static_assert(!VectorSpace<DiscreteMetric>);
    SUCCEED();
}

TEST_CASE("Concept hierarchy is strict refinement", "[concepts]") {
    // MetricSpace implies TopologicalSpace
    static_assert(TopologicalSpace<DiscreteMetric>);

    // InnerProductSpace implies NormedSpace implies MetricSpace
    // (tested via Euclidean in test_euclidean.cpp)
    SUCCEED();
}
