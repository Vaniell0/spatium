// The moments of zeta on the critical line against the conjecture of Conrey, Farmer, Keating,
// Rubinstein and Snaith: for blocks [C, D] of t, the integral of |zeta(1/2 + it)|^{2k} there beside
// the integral of P_k(log(t / 2 pi)), the polynomial their recipe gives, and the ratio. P_2 is
// Heath-Brown's, proved; P_3 is conjectured. Coefficients from "Integral moments of L-functions"
// (Proc. London Math. Soc. 2005), (5.1.5) and (5.1.6); their Table 5.1.2 has the sixth moment on
// blocks of 50000 up to 1.45e6, and [0, 2350000].
//
//   zeta_moments <k> <from> <to> <block> [threads]
//
// Above t = 1000 each block is integrated by the fixed rule (`zeta_moment_fixed`), and one unit of
// it, the middle one, again with twice its points: the rule's own convergence. Z itself is held
// apart, against the proved Z and against tables of zeros. The adaptive rule is no check up there:
// Z carries a rounding of about 1e-15 t, and at t near 1e8 the rule either stopped early on an
// estimate the oscillation fooled (9e-5 off at a tolerance loosened to that rounding) or never met
// one below it. Blocks go from the highest down, so the dearest are not left for last. One line per
// block as it is done: C D conjecture reality ratio evaluations seconds, the relative difference of
// the two rules on the checked unit, and the bound on what the cut of the Riemann-Siegel formula costs
// the block (`error` of `zeta_moment_fixed`), relative.
#include <spatium/algebra/zeta.hpp>

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <numbers>
#include <thread>
#include <vector>

namespace {

// the coefficients of P_k(x), highest power first
std::vector<long double> polynomial(int k) {
    switch (k) {
    case 1: return {1, 2 * std::numbers::egamma_v<long double>};
    case 2: return {0.0506605918211688857219397316048638L, 0.69886988487897996984709628427658502L, 2.425962198846682004756575310160663L,
                    3.227907964901254764380689851274668L, 1.312424385961669226168440066229978L};
    case 3: return {0.000005708527034652788398376841445252313L, 0.00040502133088411440331215332025984L, 0.011072455215246998350410400826667L,
                    0.14840073080150272680851401518774L, 1.0459251779054883439385323798059L, 3.984385094823534724747964073429L,
                    8.60731914578120675614834763629L, 10.274330830703446134183009522L, 6.59391302064975810465713392L,
                    0.9165155076378930590178543L};
    default: return {};
    }
}

// the integral of P(log(t / 2 pi)) from 0 to T: with t = 2 pi e^x it is 2 pi e^x Q(x), Q + Q' = P,
// so Q is P minus its derivative plus its second, and so on
long double conjecture_to(const std::vector<long double>& p, long double T) {
    if (T <= 0) return 0;
    const long double x = std::log(T / (2 * std::numbers::pi_v<long double>));
    std::vector<long double> d(p.rbegin(), p.rend());   // lowest power first
    long double q = 0;
    for (int sign = 1; !d.empty(); sign = -sign) {
        long double v = 0;
        for (std::size_t i = d.size(); i-- > 0;) v = v * x + d[i];
        q += sign * v;
        for (std::size_t i = 1; i < d.size(); ++i) d[i - 1] = d[i] * static_cast<long double>(i);
        d.pop_back();
    }
    return 2 * std::numbers::pi_v<long double> * std::exp(x) * q;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr, "zeta_moments <k> <from> <to> <block> [threads]\n");
        return 2;
    }
    const int k = std::atoi(argv[1]);
    const double from = std::atof(argv[2]), to = std::atof(argv[3]), block = std::atof(argv[4]);
    const int threads = argc > 5 ? std::atoi(argv[5]) : 1;
    const auto p = polynomial(k);
    if (p.empty() || !(from < to) || !(block > 0) || threads < 1) {
        std::fprintf(stderr, "k is 1, 2 or 3; from < to; block > 0\n");
        return 2;
    }
    const long blocks = static_cast<long>(std::ceil((to - from) / block));
    std::atomic<long> next{0};
    std::mutex out;
    std::vector<std::thread> pool;
    for (int i = 0; i < threads; ++i)
        pool.emplace_back([&] {
            for (long b; (b = next++) < blocks;) {
                const double c = from + (blocks - 1 - b) * block, d = std::min(to, c + block);
                const auto start = std::chrono::steady_clock::now();
                const bool fixed = c >= 1000;
                const auto m = fixed ? spatium::zeta_moment_fixed(k, c, d) : spatium::zeta_moment(k, c, d);
                double routes = 0;
                if (fixed) {
                    const double u = std::floor((c + d) / 2);
                    const auto one = spatium::zeta_moment_fixed(k, u, u + 1);
                    const double other = spatium::zeta_moment_fixed(k, u, u + 1, 2 * static_cast<int>(one.evaluations)).value;
                    routes = (one.value - other) / other;
                }
                const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
                const long double conj = conjecture_to(p, d) - conjecture_to(p, c);
                const std::lock_guard lock(out);
                std::printf("%.0f %.0f %.4Lf %.4f %.9Lf %ld %.1f %.2e %.2e\n", c, d, conj, m.value, m.value / conj, m.evaluations, seconds, routes, m.error / m.value);
                std::fflush(stdout);
            }
        });
    for (auto& t : pool) t.join();
}
