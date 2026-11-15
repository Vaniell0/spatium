// Z against Odlyzko's tables of zeros (www-users.cse.umn.edu/~odlyzko/zeta_tables): between two
// neighbouring zeros of a table Z must change sign once, so its sign at the midpoints alternates, or a
// zero is missing or extra; each zero, refined by the Illinois method to 1e-12, is set beside the
// table's. Out of the build and of CI: it needs GCC's __float128 for the phases past 10^12.
//
//   g++ -std=gnu++23 -O3 -march=native -Iinclude scripts/external/zeta/odlyzko_check.cpp -lquadmath -o odlyzko_check
//   ./odlyzko_check first zeros1                        the first 100000 zeros, riemann_siegel_z
//   ./odlyzko_check high zeros3 267653395647            zeros 10^12 + 1 ..., phases in double-double
//
// Measured 2026-11-13: the first 100000 all found, none extra, at most 2.99e-9 from the table (its
// stated accuracy 3e-9), 21 s; the 9998 between the first and last of zeros3 likewise, at most 2.53e-9
// (stated 1e-8), 18 minutes on one core.
#include <spatium/algebra/zeta.hpp>
#include <quadmath.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
using I = spatium::Interval<double>;

// t log n and theta(t) to about 32 digits: t and log n as pairs of doubles (log n from a quad), the
// product by fma, the phase reduced by 2 pi as a pair, the cosine of what is left in double
struct DD { double hi, lo; };
static DD dd(__float128 x) { const double h = (double)x; return {h, (double)(x - h)}; }
static std::vector<DD> logs;
static double z_dd(__float128 t) {
    const __float128 pi = M_PIq;
    const long N = (long)sqrtq(t / (2 * pi));
    while ((long)logs.size() <= N) logs.push_back(dd(logq((__float128)logs.size())));
    const __float128 theta = t / 2 * logq(t / (2 * pi)) - t / 2 - pi / 8 + 1 / (48 * t) + 7 / (5760 * t * t * t);
    const DD T = dd(t), th = dd(theta), tp = dd(2 * pi);
    double sum = 0;
    for (long n = 1; n <= N; ++n) {
        const DD l = logs[n];
        const double p = T.hi * l.hi, e = std::fma(T.hi, l.hi, -p) + (T.hi * l.lo + T.lo * l.hi);
        // theta - p - e, as hi + lo
        const double s = th.hi - p, bb = s - th.hi, err = (th.hi - (s - bb)) - (p + bb);
        const double lo = err + th.lo - e;
        const double k = std::nearbyint(s / tp.hi);
        const double r = std::fma(-k, tp.hi, s) + (lo - k * tp.lo);
        sum += std::cos(r) / std::sqrt((double)n);
    }
    return 2 * sum + spatium::zeta_detail::riemann_siegel_rest((double)t);
}

static double z_first(double t) { return t >= 1000 ? spatium::riemann_siegel_z(t) : spatium::hardy_z(I(t)).mid(); }

template<class F, class T> T refine(F f, T a, T b) {   // Illinois on a sign change
    double fa = f(a), fb = f(b);
    int side = 0;
    for (int i = 0; i < 100 && b - a > 1e-12; ++i) {
        T c = (a * fb - b * fa) / (fb - fa);
        if (!(c > a && c < b)) c = (a + b) / 2;
        const double fc = f(c);
        if ((fc > 0) == (fa > 0)) { a = c, fa = fc; if (side == -1) fb /= 2; side = -1; }
        else { b = c, fb = fc; if (side == 1) fa /= 2; side = 1; }
    }
    return (a * fb - b * fa) / (fb - fa);
}

int main(int argc, char** argv) {
    const std::string mode = argv[1];
    std::ifstream in(argv[2]);
    std::vector<std::string> rows;
    for (std::string s; in >> s;) if (s.find_first_not_of("0123456789.") == std::string::npos && s.find('.') != std::string::npos) rows.push_back(s);
    const __float128 base = argc > 3 ? strtoflt128(argv[3], nullptr) : 0;
    std::vector<__float128> g;
    for (const auto& s : rows) g.push_back(base + strtoflt128(s.c_str(), nullptr));
    const auto Z = [&](__float128 t) { return mode == "high" ? z_dd(t) : z_first((double)t); };
    long missing = 0, beyond = 0;
    double worst = 0;
    // past 10^12 the table starts mid-stream: its first and last zeros only bound the others
    const std::size_t first = mode == "high" ? 1 : 0, last = mode == "high" ? g.size() - 1 : g.size();
    __float128 prev_mid = first ? (g[0] + g[1]) / 2 : g[0] - (g[1] - g[0]) / 2;
    double prev_sign = Z(prev_mid);
    for (std::size_t i = first; i < last; ++i) {
        const __float128 mid = i + 1 < g.size() ? (g[i] + g[i + 1]) / 2 : g[i] + (g[i] - g[i - 1]) / 2;
        const double here = Z(mid);
        if ((here > 0) == (prev_sign > 0)) { ++missing; std::printf("no sign change around zero %zu at %.10f\n", i + 1, (double)(g[i] - base)); }
        else {
            const __float128 ours = refine(Z, prev_mid, mid);
            const double d = std::abs((double)(ours - g[i]));
            worst = std::max(worst, d);
            if (d > (mode == "high" ? 1e-8 : 3e-9)) ++beyond;
        }
        prev_mid = mid, prev_sign = here;
    }
    std::printf("%s: %zu zeros of the table, %ld without a sign change, %ld past its stated accuracy, largest difference %.2e\n",
                mode.c_str(), last - first, missing, beyond, worst);
}
