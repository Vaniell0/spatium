// Against Tight Inclusion (Wang et al., "A Large-scale Benchmark and an
// Inclusion-based Algorithm for Continuous Collision Detection", TOG 2021,
// github.com/Continuous-Collision-Detection/Tight-Inclusion, MIT) on the
// Sample-Queries set (github.com/Continuous-Collision-Detection/
// Sample-Queries): vertex-face and edge-edge queries with a ground truth,
// through Tight Inclusion's own vertexFaceCCD/edgeEdgeCCD with the
// benchmark's settings (tolerance 1e-6, no minimum separation, 10^6
// iterations) and through physics/mechanics/surface_ccd.hpp, which sees a
// vertex as a chart of no extent, an edge as a chart in one parameter and
// a triangle as a bilinear chart collapsed at one corner.
//
// Not part of the build. Tight Inclusion is built from its sources here
// with its spdlog logger replaced by a stub (it reports only errors):
//
//   g++ -std=c++23 -O3 -march=native -DNDEBUG -I<ti>/src -I<eigen3> -I<spatium>/include \
//       <spatium>/benchmarks/external/ti_compare.cpp <ti>/src/tight_inclusion/{ccd,interval,interval_root_finder,avx}.cpp -o ti_compare
//   ./ti_compare <Sample-Queries> [WIDTH] [PAIR_BUDGET]
//
// The numbers in the dataset are exact doubles -- at most 53 significant
// bits over a power of two -- so they are read without rationals.

#include <tight_inclusion/ccd.hpp>
#include <spatium/physics/mechanics/surface_ccd.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <print>
#include <sstream>
#include <string>
#include <vector>

namespace mech = spatium::physics::mechanics;
using SV = spatium::Vec<double, 3>;
namespace fs = std::filesystem;

struct Query {
    std::array<SV, 8> x;
    bool truth;
};

static std::vector<Query> read(const fs::path& file) {
    std::ifstream in(file);
    std::vector<Query> out;
    std::string line;
    Query q{};
    int row = 0;
    while (std::getline(in, line)) {
        // Numerators and denominators can pass 2^64; each value has at most
        // 53 significant bits over a power of two, so a long double reads
        // both exactly and their quotient is the double itself.
        std::array<long double, 7> v{};
        std::stringstream ss(line);
        std::string cell;
        for (int k = 0; k < 7 && std::getline(ss, cell, ','); ++k) v[k] = std::strtold(cell.c_str(), nullptr);
        q.x[row] = SV{double(v[0] / v[1]), double(v[2] / v[3]), double(v[4] / v[5])};
        if (row == 0) q.truth = v[6] != 0;
        if (++row == 8) { out.push_back(q); row = 0; }
    }
    return out;
}

// A point that moves from a to b.
static mech::MovingChart<double> vertex(const SV& a, const SV& b) {
    mech::MovingChart<double> c;
    const SV d{b - a};
    c.p = [a](double, double) { return a; };
    c.w = [d](double, double) { return d; };
    c.u1 = c.v1 = 0;
    return c;
}
// An edge from e0 to e1 at the start, f0 to f1 at the end.
static mech::MovingChart<double> edge(const SV& e0, const SV& e1, const SV& f0, const SV& f1) {
    mech::MovingChart<double> c;
    const SV d0{f0 - e0}, d1{f1 - e1};
    c.p = [e0, e1](double u, double) { return SV{e0 + (e1 - e0) * u}; };
    c.w = [d0, d1](double u, double) { return SV{d0 + (d1 - d0) * u}; };
    c.bp = {SV{e1 - e0}.norm(), 0, 0, 0};
    c.bw = {SV{d1 - d0}.norm(), 0, 0, 0};
    c.v1 = 0;
    return c;
}
// A triangle a, b, c as x(u, v) = a + u ((1 - v)(b - a) + v (c - a)).
static mech::MovingChart<double> triangle(const SV& a, const SV& b, const SV& c, const SV& a1, const SV& b1,
                                          const SV& c1) {
    mech::MovingChart<double> m;
    const SV da{a1 - a}, db{b1 - b}, dc{c1 - c};
    const auto tri = [](const SV& a, const SV& b, const SV& c, double u, double v) {
        return SV{a + (SV{b - a} * (1 - v) + SV{c - a} * v) * u};
    };
    m.p = [=](double u, double v) { return tri(a, b, c, u, v); };
    m.w = [=](double u, double v) { return tri(da, db, dc, u, v); };
    m.bp = {std::max(SV{b - a}.norm(), SV{c - a}.norm()), SV{c - b}.norm(), 0, 0};
    m.bw = {std::max(SV{db - da}.norm(), SV{dc - da}.norm()), SV{dc - db}.norm(), 0, 0};
    return m;
}

int main(int argc, char** argv) {
    const fs::path root = argc > 1 ? argv[1] : "Sample-Queries";
    const double width = argc > 2 ? std::atof(argv[2]) : 1e-6;
    // Tight Inclusion's cap is 10^6 iterations; ours counts pairs bounded.
    const std::size_t budget = argc > 3 ? std::size_t(std::atoll(argv[3])) : 1000000;
    using clk = std::chrono::steady_clock;
    for (const char* kind : {"vertex-face", "edge-edge"}) {
        const bool vf = std::string(kind) == "vertex-face";
        std::size_t n = 0, positives = 0, fn_ti = 0, fp_ti = 0, fn_us = 0, fp_us = 0;
        double s_ti = 0, s_us = 0;
        std::vector<fs::path> files;
        for (const auto& e : fs::recursive_directory_iterator(root))
            if (e.is_regular_file() && e.path().extension() == ".csv" && e.path().parent_path().filename() == kind)
                files.push_back(e.path());
        std::sort(files.begin(), files.end());
        for (const auto& f : files) {
            const double f_ti = s_ti, f_us = s_us;
            const std::size_t f_n = n;
            std::size_t f_pairs = 0, f_worst = 0;
            for (const auto& q : read(f)) {
                const auto& x = q.x;
                const auto E = [](const SV& v) { return ticcd::Vector3(v[0], v[1], v[2]); };
                double toi = 0, out_tol = 0;
                auto t0 = clk::now();
                const bool ti = vf ? ticcd::vertexFaceCCD(E(x[0]), E(x[1]), E(x[2]), E(x[3]), E(x[4]), E(x[5]), E(x[6]),
                                                          E(x[7]), {-1, -1, -1}, 0.0, toi, 1e-6, 1.0, 1000000, out_tol)
                                   : ticcd::edgeEdgeCCD(E(x[0]), E(x[1]), E(x[2]), E(x[3]), E(x[4]), E(x[5]), E(x[6]),
                                                        E(x[7]), {-1, -1, -1}, 0.0, toi, 1e-6, 1.0, 1000000, out_tol);
                auto t1 = clk::now();
                const auto a = vf ? vertex(x[0], x[4]) : edge(x[0], x[1], x[4], x[5]);
                const auto b = vf ? triangle(x[1], x[2], x[3], x[5], x[6], x[7]) : edge(x[2], x[3], x[6], x[7]);
                const auto r = mech::first_contact(a, b, width, budget);
                auto t2 = clk::now();
                s_ti += std::chrono::duration<double>(t1 - t0).count();
                s_us += std::chrono::duration<double>(t2 - t1).count();
                ++n;
                positives += q.truth;
                fn_ti += q.truth && !ti;
                fp_ti += !q.truth && ti;
                fn_us += q.truth && !r.hit;
                fp_us += !q.truth && r.hit;
                f_pairs += r.pairs;
                f_worst = std::max(f_worst, r.pairs);
            }
            const std::size_t k = n - f_n;
            if (k) std::println(stderr, "{}: {} queries, TI {:.2f} us, ours {:.2f} us, pairs {:.0f} mean {} worst",
                                f.parent_path().parent_path().filename().string() + "/" + f.filename().string(), k,
                                1e6 * (s_ti - f_ti) / k, 1e6 * (s_us - f_us) / k, double(f_pairs) / k, f_worst);
        }
        std::println("{}: {} queries, {} positive", kind, n, positives);
        std::println("  Tight Inclusion: {:.3f} us a query, {} false negatives, {} false positives", 1e6 * s_ti / n, fn_ti,
                     fp_ti);
        std::println("  ours:            {:.3f} us a query, {} false negatives, {} false positives", 1e6 * s_us / n, fn_us,
                     fp_us);
    }
}
