// Entropy, cross-entropy and the Kullback-Leibler divergence, and distilling a softmax at a
// temperature.
//
// Each claim is held against a second route: the divergence against the cross-entropy less the
// entropy, computed apart; the divergence of two binomial distributions against n times that of
// one draw (it adds over independent draws), with the binomials built from
// discrete::binomial_coefficient; the distillation gradient written in softmax.hpp against Dual
// carried through statistics::kl_divergence. Each check gets a case it must fail: the divergence
// is not symmetric, and a q that gives zero to what p draws costs an infinite number of nats.

#include <spatium/algebra/analysis/optimize.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/algebra/models/softmax.hpp>
#include <spatium/discrete/combinatorics.hpp>
#include <spatium/statistics/information.hpp>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <random>
#include <vector>

using namespace spatium;
using namespace spatium::statistics;
using Catch::Matchers::WithinAbs;

namespace {

std::vector<double> random_distribution(std::mt19937& rng, std::size_t n) {
    std::uniform_real_distribution<double> u(0.01, 1);
    std::vector<double> p(n);
    double total = 0;
    for (double& v : p)
        total += v = u(rng);
    for (double& v : p)
        v /= total;
    return p;
}

std::vector<double> binomial(unsigned n, double a) {
    std::vector<double> p(n + 1);
    for (unsigned k = 0; k <= n; ++k)
        p[k] = double(*discrete::binomial_coefficient<double>(n, k)) * std::pow(a, k) * std::pow(1 - a, n - k);
    return p;
}

double span_kl(const std::vector<double>& p, const std::vector<double>& q) {
    return kl_divergence(std::span<const double>(p), std::span<const double>(q));
}

} // namespace

TEST_CASE("entropy: log n for n equal chances, zero for a certainty, one bit for a coin", "[statistics][information]") {
    for (std::size_t n : {1u, 2u, 7u, 64u}) {
        const std::vector<double> uniform(n, 1.0 / double(n));
        CHECK_THAT(entropy(std::span<const double>(uniform)), WithinAbs(std::log(double(n)), 1e-14));
    }
    const std::vector<double> sure{0, 0, 1, 0};
    CHECK(entropy(std::span<const double>(sure)) == 0.0);
    const std::vector<double> coin{0.5, 0.5};
    CHECK_THAT(in_bits(entropy(std::span<const double>(coin))), WithinAbs(1.0, 1e-15));
}

TEST_CASE("the divergence is the cross-entropy less the entropy, and not symmetric", "[statistics][information]") {
    std::mt19937 rng(9);
    int asymmetric = 0;
    for (int round = 0; round < 50; ++round) {
        const auto p = random_distribution(rng, 6), q = random_distribution(rng, 6);
        const std::span<const double> sp(p), sq(q);
        const double kl = kl_divergence(sp, sq);
        CHECK_THAT(kl, WithinAbs(cross_entropy(sp, sq) - entropy(sp), 1e-14));
        CHECK(kl >= 0);                                    // Gibbs
        CHECK(kl_divergence(sp, sp) == 0.0);
        asymmetric += std::abs(kl - kl_divergence(sq, sp)) > 1e-6;
    }
    CHECK(asymmetric > 40);

    // what p draws and q never would: infinite; what neither draws: nothing
    const std::vector<double> p{0.5, 0.5, 0}, q{1, 0, 0}, r{0.5, 0.25, 0.25};
    CHECK(std::isinf(span_kl(p, q)));
    CHECK(std::isinf(cross_entropy(std::span<const double>(p), std::span<const double>(q))));
    CHECK(std::isfinite(span_kl(p, r)));
    CHECK(std::isfinite(span_kl(q, p)));
}

TEST_CASE("the divergence adds over independent draws: two binomials are n times one draw", "[statistics][information]") {
    for (unsigned n : {1u, 5u, 20u, 60u})
        for (const auto& [a, b] : {std::pair{0.3, 0.6}, std::pair{0.5, 0.1}, std::pair{0.9, 0.85}}) {
            INFO("n = " << n << ", a = " << a << ", b = " << b);
            const double one = a * std::log(a / b) + (1 - a) * std::log((1 - a) / (1 - b));
            CHECK_THAT(span_kl(binomial(n, a), binomial(n, b)), WithinAbs(n * one, 1e-11 * (1 + n * one)));
        }
}

TEST_CASE("distilling a softmax: T^2 times the gradient at T is the derivative of T^2 KL", "[statistics][information][softmax]") {
    constexpr std::size_t K = 5;
    using Model = LinearSoftmax<K>;
    using D = Dual<double>;
    std::mt19937 rng(17);
    std::normal_distribution<double> n(0, 1);
    std::uniform_real_distribution<double> u(0, 1);
    for (double T : {1.0, 2.5, 8.0}) {
        INFO("T = " << T);
        Model teacher, student;
        for (double& w : teacher.weights())
            w = 2 * n(rng);
        for (double& w : student.weights())
            w = n(rng);
        std::vector<Model::Row> list(9);
        for (auto& r : list)
            for (double& v : r)
                v = u(rng);
        std::vector<double> p(list.size());
        teacher.distribution(list, p, T);

        const auto g = student.gradient(list, std::span<const double>(p), T);
        for (std::size_t j = 0; j < K; ++j) {
            LinearSoftmax<K, D>::Weights dw;
            for (std::size_t k = 0; k < K; ++k)
                dw[k] = k == j ? D::variable(student.weights()[k]) : D(student.weights()[k]);
            std::vector<LinearSoftmax<K, D>::Row> dlist(list.size());
            for (std::size_t i = 0; i < list.size(); ++i)
                for (std::size_t k = 0; k < K; ++k)
                    dlist[i][k] = D(list[i][k]);
            std::vector<D> q(list.size()), dp(list.size());
            LinearSoftmax<K, D>(dw).distribution(dlist, q, D(T));
            for (std::size_t i = 0; i < p.size(); ++i)
                dp[i] = D(p[i]);
            const D loss = D(T * T) * kl_divergence(std::span<const D>(dp), std::span<const D>(q));
            REQUIRE_THAT(T * T * g[j], WithinAbs(loss.deriv, 1e-11));
        }
        // the cross-entropy at T is the divergence plus the teacher's entropy
        std::vector<double> q(list.size());
        student.distribution(list, q, T);
        CHECK_THAT(student.cross_entropy(list, std::span<const double>(p), T),
                   WithinAbs(span_kl(p, q) + entropy(std::span<const double>(p)), 1e-13));
    }
}

TEST_CASE("a student distilled from a teacher's soft distribution chooses as the teacher does", "[statistics][information][softmax]") {
    constexpr std::size_t K = 6;
    using Model = LinearSoftmax<K>;
    const Model teacher(Model::Weights{3, -2, 1, 0, 4, -1});
    std::uniform_real_distribution<double> u(0, 1);
    auto lists = [&](std::uint64_t seed, int count) {
        std::mt19937 rng(static_cast<std::uint32_t>(seed));
        std::vector<std::vector<Model::Row>> out(static_cast<std::size_t>(count), std::vector<Model::Row>(10));
        for (auto& l : out)
            for (auto& r : l)
                for (double& v : r)
                    v = u(rng);
        return out;
    };
    auto distilled = [&](double T, double rate) {
        Model student;
        AdaGrad<K> rule;
        std::vector<double> p(10);
        for (const auto& l : lists(21, 1500)) {
            teacher.distribution(l, p, T);
            auto g = student.gradient(l, std::span<const double>(p), T);
            for (double& v : g)
                v *= T * T;
            rule.step(student.weights(), g, rate);
        }
        int same = 0;
        const auto held = lists(99, 400);
        for (const auto& l : held)
            same += student.choose(l) == teacher.choose(l);
        return double(same) / double(held.size());
    };
    CHECK(distilled(2.0, 0.5) > 0.95);
    CHECK(distilled(2.0, 0.0) < 0.3); // chance with ten candidates is 10%
}
