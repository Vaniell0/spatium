// A net of one hidden tanh layer, the Adam step that trains it, and a pass of fitting.
//
// Each claim is held against a second route: the backward pass against Dual carried through the
// loss, weight by weight, and against a finite difference; the trained net against the function
// it was shown, on points it never saw; Adam's first step against its closed form; the pass of
// fitting against the losses its objective reported. Each check gets a case it must fail: a net
// trained at rate zero, and plain gradient descent, whose steps scale with the gradient.

#include <spatium/algebra/analysis/optimize.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/algebra/models/network.hpp>
#include <spatium/statistics/fit.hpp>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <random>
#include <vector>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

constexpr std::size_t inputs = 5, hidden = 8;
using Net = TanhNet<inputs, hidden>;

struct Point {
    Net::Input x;
    double target;
};

double truth(const Net::Input& x) { return std::tanh(2 * x[0] - x[1]); }

std::vector<Point> points(std::uint64_t seed, int count) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> u(-1, 1);
    std::vector<Point> out(static_cast<std::size_t>(count));
    for (Point& p : out) {
        for (double& v : p.x)
            v = u(rng);
        p.target = truth(p.x);
    }
    return out;
}

double mean_error(const Net& net, const std::vector<Point>& held_out) {
    double sum = 0;
    for (const Point& p : held_out)
        sum += net.loss(p.x, p.target);
    return sum / static_cast<double>(held_out.size());
}

} // namespace

TEST_CASE("the backward pass agrees with Dual and with a finite difference", "[network]") {
    std::mt19937_64 rng(2);
    std::normal_distribution<double> n(0, 0.7);
    using D = Dual<double>;
    for (int round = 0; round < 10; ++round) {
        Net net = Net::random(static_cast<std::uint64_t>(round + 1));
        for (double& v : net.weights())
            v += n(rng);
        Net::Input x{};
        for (double& v : x)
            v = n(rng);
        const double target = (round % 2) ? 0.6 : -0.3;
        const auto g = net.gradient(x, target);

        TanhNet<inputs, hidden, D>::Input dx{};
        for (std::size_t k = 0; k < inputs; ++k)
            dx[k] = D(x[k]);
        for (std::size_t i = 0; i < Net::weight_count; ++i) {
            TanhNet<inputs, hidden, D>::Weights dw;
            for (std::size_t k = 0; k < Net::weight_count; ++k)
                dw[k] = k == i ? D::variable(net.weights()[k]) : D(net.weights()[k]);
            const D loss = TanhNet<inputs, hidden, D>(dw).loss(dx, D(target));
            REQUIRE_THAT(g[i], WithinAbs(loss.deriv, 1e-12));

            Net up = net, down = net;
            up.weights()[i] += 1e-6, down.weights()[i] -= 1e-6;
            REQUIRE_THAT(g[i], WithinAbs((up.loss(x, target) - down.loss(x, target)) / 2e-6, 1e-6));
        }
    }
}

TEST_CASE("Adam's first step is the rate, signed, and its steps barely see the gradient's scale", "[network][adam]") {
    Adam<3> rule;
    std::array<double, 3> w{0, 0, 0};
    rule.step(w, {0.5, -2e-3, 7}, 0.1);
    CHECK_THAT(w[0], WithinAbs(-0.1, 1e-7));
    CHECK_THAT(w[1], WithinAbs(0.1, 1e-5)); // 2e-3 against epsilon 1e-8: off by 5e-6 of the rate
    CHECK_THAT(w[2], WithinAbs(-0.1, 1e-9));

    std::mt19937 rng(11);
    std::normal_distribution<double> n(0, 1);
    Adam<3> a, b;
    std::array<double, 3> wa{}, wb{}, plain_a{}, plain_b{};
    for (int i = 0; i < 50; ++i) {
        std::array<double, 3> g{n(rng), n(rng), n(rng)}, scaled = g;
        for (double& v : scaled)
            v *= 1000;
        a.step(wa, g, 0.01);
        b.step(wb, scaled, 0.01);
        for (std::size_t j = 0; j < 3; ++j) {
            plain_a[j] -= 0.01 * g[j];
            plain_b[j] -= 0.01 * scaled[j];
        }
    }
    for (std::size_t j = 0; j < 3; ++j) {
        CHECK_THAT(wa[j], WithinAbs(wb[j], 1e-6));
        CHECK(std::abs(plain_a[j] - plain_b[j]) > 0.1);
    }
}

TEST_CASE("fitted by Adam, the net learns the function it is shown", "[network][fit]") {
    const auto train = points(11, 6000);
    const auto held_out = points(99, 300);
    auto fitted = [&](double rate) {
        // random, not zero: with all weights zero every unit gets the same gradient and stays equal
        Net net = Net::random(1);
        Adam<Net::weight_count> rule;
        std::mt19937_64 rng(5);
        statistics::fit_pass(net.weights(), rule, std::span<const Point>(train), rate, rng,
                             [&](const Point& p) {
                                 return statistics::Fit<Net::Weights>{net.loss(p.x, p.target),
                                                                      net.gradient(p.x, p.target)};
                             });
        return mean_error(net, held_out);
    };
    REQUIRE(fitted(0.01) < 0.01);
    REQUIRE(fitted(0.0) > 0.1); // a net that never moves does not fit it
}

TEST_CASE("a pass visits every item once and returns the mean loss before its steps", "[fit]") {
    // a rule that moves nothing: the mean the pass returns is the plain mean of the losses
    struct Still {
        void step(std::array<double, 1>&, const std::array<double, 1>&, double) {}
    };
    std::array<double, 1> params{0};
    Still rule;
    std::vector<int> items{1, 2, 3, 4, 5, 6, 7};
    std::vector<int> seen;
    std::mt19937_64 rng(3);
    const double mean = statistics::fit_pass(params, rule, std::span<const int>(items), 1.0, rng,
                                             [&](const int& i) {
                                                 seen.push_back(i);
                                                 return statistics::Fit<std::array<double, 1>>{double(i), {0}};
                                             });
    CHECK(mean == 4.0);
    std::vector<int> sorted = seen;
    std::ranges::sort(sorted);
    CHECK(sorted == items);
    CHECK(seen != items); // drawn in an order, not walked in the given one
    // an empty sample is a pass with nothing to say
    CHECK(statistics::fit_pass(params, rule, std::span<const int>{}, 1.0, rng, [](const int&) {
              return statistics::Fit<std::array<double, 1>>{};
          }) == 0.0);
}
