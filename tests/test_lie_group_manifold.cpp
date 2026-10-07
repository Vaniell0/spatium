// A Lie group as a space (spaces/lie_group.hpp): the algorithms written against exp,
// log and a metric -- the Frechet mean, geodesic, midpoint, the verifiers -- run on
// rotations and rigid motions with no copy of them for groups.
//
// The checks are paths that share nothing: the mean of two rotations about one axis
// against the half angle that is its closed form; the mean of rotated data against the
// rotated mean (equivariance); the Karcher condition the mean must satisfy, read off the
// logs; the geodesic midpoint against the one-parameter subgroup it lies on.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/core/access.hpp>
#include <spatium/core/verify.hpp>
#include <spatium/algebra/calculus.hpp>
#include <spatium/algebra/groups/se3.hpp>
#include <spatium/algebra/groups/so3.hpp>
#include <spatium/spaces/lie_group.hpp>
#include <cmath>
#include <numbers>
#include <vector>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {
using Rot = LieGroupManifold<algebra::SO3<>>;
using M3 = Matrix<double, 3, 3>;
const algebra::SO3<> so3;

M3 rot(double x, double y, double z) { return so3.exp(Vec<double, 3>{x, y, z}); }

double max_abs_diff(const M3& a, const M3& b) {
    double m = 0;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) m = std::max(m, std::abs(a(i, j) - b(i, j)));
    return m;
}
}  // namespace

TEST_CASE("a rotation group is a Riemannian manifold, and its distance is the angle", "[lie_group]") {
    static_assert(Manifold<Rot>);
    static_assert(RiemannianManifold<Rot>);
    static_assert(spaces::Riemannian<Rot>);
    const Rot space;
    const double angle = 0.9;
    // two rotations about one axis are `angle` apart
    CHECK_THAT(space.distance(rot(0, 0, 0.3), rot(0, 0, 0.3 + angle)), WithinAbs(angle, 1e-12));
    // the distance from the identity is the rotation angle
    CHECK_THAT(space.distance(so3.identity(), rot(0.2, -0.5, 0.4)),
               WithinAbs(std::sqrt(0.2 * 0.2 + 0.5 * 0.5 + 0.4 * 0.4), 1e-12));
}

TEST_CASE("exp and log invert each other on rotations", "[lie_group]") {
    const Rot space;
    const std::vector<M3> pts{rot(0, 0, 0), rot(0.4, 0.1, -0.3), rot(-0.2, 0.7, 0.5), rot(1.1, -0.4, 0.2)};
    const auto r = verify_exp_log(space, pts, 1e-9);
    INFO(r.message());
    CHECK(r.passed);
    const auto m = verify_metric(space, pts, 1e-9);
    INFO(m.message());
    CHECK(m.passed);
}

TEST_CASE("the Frechet mean of two rotations about one axis is the half angle", "[lie_group][symmetry]") {
    const Rot space;
    const M3 a = rot(0, 0.2, 0), b = rot(0, 0.9, 0);
    const auto mean = frechet_mean(space, std::vector<M3>{a, b}, a);
    CHECK(max_abs_diff(mean, rot(0, 0.55, 0)) < 1e-9);
    // the same point as the geodesic midpoint of the customization points
    CHECK(max_abs_diff(spaces::midpoint(space, a, b), rot(0, 0.55, 0)) < 1e-12);
}

TEST_CASE("the mean is equivariant: the mean of rotated data is the rotated mean", "[lie_group][symmetry]") {
    const Rot space;
    const std::vector<M3> data{rot(0.10, 0.20, -0.10), rot(-0.15, 0.05, 0.20), rot(0.20, -0.10, 0.05), rot(0.00, 0.30, -0.20)};
    const M3 g = rot(0.7, -0.2, 0.4);
    std::vector<M3> moved;
    for (const auto& r : data) moved.push_back(so3.compose(g, r));

    const auto mean = frechet_mean(space, data, data[0]);
    const auto mean_moved = frechet_mean(space, moved, moved[0]);
    CHECK(max_abs_diff(mean_moved, so3.compose(g, mean)) < 1e-8);
}

TEST_CASE("the mean satisfies the Karcher condition: the logs to the data sum to zero", "[lie_group][symmetry]") {
    const Rot space;
    const std::vector<M3> data{rot(0.3, 0.1, -0.2), rot(-0.2, 0.4, 0.1), rot(0.1, -0.3, 0.3), rot(0.0, 0.2, 0.2), rot(-0.1, -0.1, 0.0)};
    const auto mean = frechet_mean(space, data, data[0]);
    Vec<double, 3> sum{0, 0, 0};
    for (const auto& r : data) sum = sum + space.log_map(mean, r);
    CHECK_THAT(std::sqrt(sum.norm_squared()), WithinAbs(0.0, 1e-8));
}

TEST_CASE("rigid motions: exp and log invert each other", "[lie_group]") {
    using Pose = LieGroupManifold<algebra::SE3<>>;
    static_assert(Manifold<Pose>);
    const Pose space;
    const algebra::SE3<> se3;
    const auto p = se3.exp(Vec<double, 6>{0.1, -0.2, 0.3, 0.5, 0.1, -0.4});
    const auto q = se3.exp(Vec<double, 6>{-0.3, 0.1, 0.2, -0.2, 0.6, 0.3});
    const auto back = space.exp_map(p, space.log_map(p, q), 1.0);
    double m = 0;
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) m = std::max(m, std::abs(back(i, j) - q(i, j)));
    CHECK(m < 1e-9);
}
