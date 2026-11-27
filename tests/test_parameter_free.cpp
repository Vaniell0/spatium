// Step rules with no rate (algebra/analysis/parameter_free.hpp).
//
// L-DoG solves a least-squares problem with no rate given, at the data's own scale and at a hundred
// times it; the case it must fail, a fixed-rate SGD good at one scale, diverges at the other. Muon-DoG
// fits a matrix the same way at both scales. The polynomial average of a constant is the constant, and
// its first value is the first point. The orthogonal factor at a shape given at run time is the one
// `Muon<N, R, C>` computes at a shape fixed at compile time, to the last bit, wide and tall.

#include <spatium/algebra/analysis/optimize.hpp>
#include <spatium/algebra/analysis/parameter_free.hpp>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <random>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

constexpr std::size_t rows = 30, cols = 5;

struct Problem {
    std::array<double, rows * cols> a{};
    std::array<double, rows> b{};
    double loss(const std::array<double, cols>& x) const {
        double s = 0;
        for (std::size_t i = 0; i < rows; ++i) {
            double r = -b[i];
            for (std::size_t j = 0; j < cols; ++j)
                r += a[i * cols + j] * x[j];
            s += r * r;
        }
        return s / rows;
    }
    std::array<double, cols> gradient(const std::array<double, cols>& x) const {
        std::array<double, cols> g{};
        for (std::size_t i = 0; i < rows; ++i) {
            double r = -b[i];
            for (std::size_t j = 0; j < cols; ++j)
                r += a[i * cols + j] * x[j];
            for (std::size_t j = 0; j < cols; ++j)
                g[j] += 2 * r * a[i * cols + j] / rows;
        }
        return g;
    }
};

Problem problem(double scale) {
    std::mt19937_64 rng(1);
    std::normal_distribution<double> n(0, 1);
    Problem p;
    for (double& v : p.a)
        v = scale * n(rng);
    for (double& v : p.b)
        v = scale * n(rng);
    return p;
}

// The least loss, by plain gradient steps at the rate that suits this problem (1 / the largest curvature
// bound), run long: the reference the rules are held to.
double best(const Problem& p) {
    double bound = 0;
    for (double v : p.a)
        bound += v * v;
    bound = 2 * bound / rows;
    std::array<double, cols> x{};
    for (int t = 0; t < 200000; ++t) {
        const auto g = p.gradient(x);
        for (std::size_t j = 0; j < cols; ++j)
            x[j] -= g[j] / bound;
    }
    return p.loss(x);
}

} // namespace

TEST_CASE("parameter-free: L-DoG reaches the least squares at any scale, a fixed rate does not", "[parameter_free]") {
    for (double scale : {1.0, 100.0}) {
        const Problem p = problem(scale);
        const double target = best(p);
        LayerDoG<> dog;
        std::array<double, cols> x{};
        for (int t = 0; t < 20000; ++t)
            dog.step(x, p.gradient(x));
        CHECK(p.loss(x) - target < 1e-3 * (1 + target));
    }
    // the case it must fail: a rate that works at scale 1 diverges at scale 100
    const auto sgd = [](const Problem& p, double rate) {
        std::array<double, cols> x{};
        for (int t = 0; t < 2000; ++t) {
            const auto g = p.gradient(x);
            for (std::size_t j = 0; j < cols; ++j)
                x[j] -= rate * g[j];
            if (!std::isfinite(p.loss(x)))
                return HUGE_VAL;
        }
        return p.loss(x);
    };
    const Problem small = problem(1.0), large = problem(100.0);
    CHECK(sgd(small, 0.05) - best(small) < 1e-3 * (1 + best(small)));
    CHECK(!(sgd(large, 0.05) < 1e6));
}

TEST_CASE("parameter-free: Muon-DoG fits a matrix at any scale with no rate", "[parameter_free]") {
    constexpr std::size_t r = 4, c = 6, samples = 40;
    for (double scale : {1.0, 100.0}) {
        std::mt19937_64 rng(2);
        std::normal_distribution<double> n(0, 1);
        std::array<double, r * c> truth{};
        for (double& v : truth)
            v = n(rng);
        std::array<std::array<double, c>, samples> xs{};
        for (auto& x : xs)
            for (double& v : x)
                v = scale * n(rng);
        const auto loss_grad = [&](const std::array<double, r * c>& w, std::array<double, r * c>* g) {
            double loss = 0;
            if (g)
                g->fill(0);
            for (const auto& x : xs)
                for (std::size_t i = 0; i < r; ++i) {
                    double e = 0;
                    for (std::size_t j = 0; j < c; ++j)
                        e += (w[i * c + j] - truth[i * c + j]) * x[j];
                    loss += e * e / samples;
                    if (g)
                        for (std::size_t j = 0; j < c; ++j)
                            (*g)[i * c + j] += 2 * e * x[j] / samples;
                }
            return loss;
        };
        MuonDoG<> rule;
        rule.blocks = {ParamBlock{.offset = 0, .rows = r, .cols = c}};
        std::array<double, r * c> w{}, g{};
        const double first = loss_grad(w, nullptr);
        for (int t = 0; t < 3000; ++t) {
            loss_grad(w, &g);
            rule.step(w, g);
        }
        CHECK(loss_grad(w, nullptr) < 1e-3 * first);
    }
}

TEST_CASE("parameter-free: the polynomial average", "[parameter_free]") {
    PolyAverage<> avg;
    const std::array<double, 3> point{1, 2, 3};
    avg.add(point);
    CHECK(avg.mean == std::vector<double>{1, 2, 3}); // the first value is the first point
    for (int t = 0; t < 50; ++t)
        avg.add(point);
    CHECK_THAT(avg.mean[2], WithinAbs(3, 1e-12)); // a constant stays itself
    // later points weigh more than earlier ones
    PolyAverage<> rising;
    for (int t = 1; t <= 100; ++t) {
        const std::array<double, 1> x{double(t)};
        rising.add(x);
    }
    CHECK(rising.mean[0] > 50.5);
}

TEST_CASE("parameter-free: the orthogonal factor at a run-time shape is Muon's", "[parameter_free][muon]") {
    std::mt19937_64 rng(3);
    std::normal_distribution<double> n(0, 1);
    std::array<double, 18> wide{}, tall{};
    for (double& v : wide)
        v = n(rng);
    for (double& v : tall)
        v = n(rng);
    const auto a = Muon<18, 3, 6>::orthogonalized(wide, 1e-7);
    const auto b = parameter_free_detail::orthogonalized<double>(wide, 3, 6, 1e-7);
    const auto c = Muon<18, 6, 3>::orthogonalized(tall, 1e-7);
    const auto d = parameter_free_detail::orthogonalized<double>(tall, 6, 3, 1e-7);
    for (std::size_t i = 0; i < 18; ++i) {
        CHECK(a[i] == b[i]);
        CHECK(c[i] == d[i]);
    }
    // the case it must fail: the same matrix read with the other shape is another factor
    const auto e = parameter_free_detail::orthogonalized<double>(wide, 6, 3, 1e-7);
    double apart = 0;
    for (std::size_t i = 0; i < 18; ++i)
        apart = std::max(apart, std::abs(a[i] - e[i]));
    CHECK(apart > 1e-3);
}

TEST_CASE("parameter-free: the rules take a model whose size is known only when it runs", "[parameter_free]") {
    // the same least squares as above, its parameters a vector: the steps equal those on the array
    const Problem p = problem(1.0);
    LayerDoG<> on_array, on_vector;
    std::array<double, cols> x{};
    std::vector<double> y(cols, 0.0);
    for (int t = 0; t < 500; ++t) {
        on_array.step(x, p.gradient(x));
        std::array<double, cols> at{};
        std::copy(y.begin(), y.end(), at.begin());
        const auto g = p.gradient(at);
        on_vector.step(y, g);
    }
    for (std::size_t j = 0; j < cols; ++j)
        CHECK(x[j] == y[j]);
}
