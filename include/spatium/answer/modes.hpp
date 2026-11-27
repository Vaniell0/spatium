#pragma once

// Damped sinusoids in data that carries noise: how many the data holds, by description length, and their
// frequencies and decay times with standard errors.
//
// `synthesize_law` (law.hpp) fits a sum of terms whose only free numbers are the coefficients in front of them, so
// a resonance, exp(-x / tau) cos(2 pi f x), whose f and tau sit inside the operations, is out of its reach. Here
// the model is a sum of N modes
//
//   y(x) = sum_j exp(-x / tau_j) (a_j cos(2 pi f_j x) + b_j sin(2 pi f_j x)),
//
// the amplitudes linear and the f_j, tau_j nonlinear. For N modes the amplitudes are eliminated at each trial of
// the frequencies and decay times (variable projection), so the search is over 2N numbers: mode N + 1 starts from
// the best point of a grid of f and tau with the first N held, then all 2(N + 1) numbers are refined by the
// simplex method of Nelder and Mead. The count is the one of least description length
//
//   n ln(RSS / n)  +  4 N ln n,
//
// four fitted numbers to a mode. The noise is the caller's, as in law.hpp: a sigma at each point, or a covariance
// whitened by its Cholesky factor. The standard errors are those of Gauss-Newton at the optimum (the Jacobian of
// the whitened model in all 4N numbers, scaled by RSS / (n - 4N)); they are NaN where that Jacobian is singular.
//
// The grid and the simplex find a local optimum of a surface with several. The grid is what keeps the start near
// the global one, so its range must hold the modes and its step be finer than the line's width, f_step below
// 1 / (2 pi tau) or so; a mode narrower than a grid step can be missed.

#include <spatium/answer/law.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <numbers>
#include <optional>
#include <vector>

namespace spatium::answer {

struct ModeSpec {
    std::size_t max_modes = 2;
    // the grid the next mode starts from: frequencies from lo to hi in steps, decay times geometric from lo to hi
    double frequency_lo = 0, frequency_hi = 0, frequency_step = 1;
    double decay_lo = 0, decay_hi = 0;
    std::size_t decay_count = 25;
    std::size_t refine_iterations = 800;   // of the simplex, for each count of modes
    std::vector<double> covariance;        // as in LawSpec: when given, sigma is not read
};

struct Mode {
    double frequency = 0, decay = 0;         // f and tau, in the units of x
    double cosine = 0, sine = 0;             // the amplitudes a and b
    double frequency_error = std::numeric_limits<double>::quiet_NaN();
    double decay_error = std::numeric_limits<double>::quiet_NaN();
};

struct Modes {
    std::vector<Mode> modes;                 // the count of least length
    std::vector<double> lengths, rss;        // by count 0 .. max_modes (empty when the input was refused)
    double rms = 0;                          // of the chosen residual in units of sigma
};

namespace law_detail {

// the simplex method of Nelder and Mead on f over n numbers, the start's neighbours 5% away in each (0.00025 where
// one is zero), until both the simplex and the values on it are within xatol and fatol or the iterations are out
inline std::vector<double> nelder_mead(const std::function<double(const std::vector<double>&)>& f, std::vector<double> x0,
                                       std::size_t iterations, double xatol, double fatol, double* value) {
    const std::size_t n = x0.size();
    std::vector<std::vector<double>> sim(n + 1, x0);
    for (std::size_t i = 0; i < n; ++i) sim[i + 1][i] = x0[i] != 0 ? 1.05 * x0[i] : 0.00025;
    std::vector<double> fv(n + 1);
    for (std::size_t i = 0; i <= n; ++i) fv[i] = f(sim[i]);
    const auto order = [&] {
        std::vector<std::size_t> idx(n + 1);
        for (std::size_t i = 0; i <= n; ++i) idx[i] = i;
        std::ranges::sort(idx, [&](std::size_t a, std::size_t b) { return fv[a] < fv[b]; });
        auto s2 = sim;
        auto f2 = fv;
        for (std::size_t i = 0; i <= n; ++i) sim[i] = s2[idx[i]], fv[i] = f2[idx[i]];
    };
    order();
    for (std::size_t it = 0; it < iterations; ++it) {
        double dx = 0, df = 0;
        for (std::size_t i = 1; i <= n; ++i) {
            df = std::max(df, std::abs(fv[i] - fv[0]));
            for (std::size_t j = 0; j < n; ++j) dx = std::max(dx, std::abs(sim[i][j] - sim[0][j]));
        }
        if (dx <= xatol && df <= fatol) break;
        std::vector<double> centre(n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) centre[j] += sim[i][j] / double(n);
        const auto along = [&](double t) {
            std::vector<double> p(n);
            for (std::size_t j = 0; j < n; ++j) p[j] = centre[j] + t * (sim[n][j] - centre[j]);
            return p;
        };
        const auto xr = along(-1.0);
        const double fr = f(xr);
        bool shrink = false;
        if (fr < fv[0]) {
            const auto xe = along(-2.0);
            const double fe = f(xe);
            if (fe < fr) sim[n] = xe, fv[n] = fe;
            else sim[n] = xr, fv[n] = fr;
        } else if (fr < fv[n - 1]) {
            sim[n] = xr, fv[n] = fr;
        } else {
            const bool outside = fr < fv[n];
            const auto xc = along(outside ? -0.5 : 0.5);
            const double fc = f(xc);
            if (fc <= (outside ? fr : fv[n])) sim[n] = xc, fv[n] = fc;
            else shrink = true;
        }
        if (shrink)
            for (std::size_t i = 1; i <= n; ++i) {
                for (std::size_t j = 0; j < n; ++j) sim[i][j] = sim[0][j] + 0.5 * (sim[i][j] - sim[0][j]);
                fv[i] = f(sim[i]);
            }
        order();
    }
    if (value) *value = fv[0];
    return sim[0];
}

// the inverse of a symmetric positive definite matrix by its Cholesky factor; empty when it has none
inline std::vector<double> inverse_spd(const std::vector<double>& a, std::size_t k) {
    const auto l = cholesky_lower(a, k);
    if (l.empty()) return {};
    std::vector<double> inv(k * k, 0.0);
    for (std::size_t c = 0; c < k; ++c) {   // solve A x = e_c: L z = e_c, then L^T x = z
        std::vector<double> z(k, 0.0);
        for (std::size_t i = 0; i < k; ++i) {
            double s = i == c ? 1.0 : 0.0;
            for (std::size_t m = 0; m < i; ++m) s -= l[i * k + m] * z[m];
            z[i] = s / l[i * k + i];
        }
        for (std::size_t i = k; i-- > 0;) {
            double s = z[i];
            for (std::size_t m = i + 1; m < k; ++m) s -= l[m * k + i] * inv[m * k + c];
            inv[i * k + c] = s / l[i * k + i];
        }
    }
    return inv;
}

// the two columns of one mode at x
inline void mode_columns(const std::vector<double>& x, double f, double tau, std::vector<double>& c, std::vector<double>& s) {
    c.resize(x.size()), s.resize(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double e = std::exp(-x[i] / tau), w = 2 * std::numbers::pi * f * x[i];
        c[i] = e * std::cos(w), s[i] = e * std::sin(w);
    }
}

} // namespace law_detail

inline Modes synthesize_modes(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& sigma,
                              const ModeSpec& spec) {
    using namespace law_detail;
    const std::size_t n = x.size();
    Modes out;
    if (n == 0 || y.size() != n || !(spec.frequency_step > 0) || !(spec.frequency_hi >= spec.frequency_lo) ||
        !(spec.decay_lo > 0) || !(spec.decay_hi >= spec.decay_lo) || spec.decay_count == 0)
        return out;
    for (std::size_t i = 0; i < n; ++i)
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return out;
    const auto w = weigh(n, sigma, spec.covariance);
    if (!w) return out;
    std::vector<double> b = y;
    w->apply(b);
    const double ln_n = std::log(double(n));
    const auto length = [&](double rss, std::size_t modes) {
        return double(n) * std::log(std::max(rss, 1e-300) / double(n)) + 4.0 * double(modes) * ln_n;
    };

    // the whitened columns of the modes at p = (f1, tau1, f2, tau2, ...), and the fit of the amplitudes to them
    const auto columns_of = [&](const std::vector<double>& p) {
        std::vector<std::vector<double>> cols;
        std::vector<double> c, s;
        for (std::size_t j = 0; j + 1 < p.size(); j += 2) {
            mode_columns(x, p[j], p[j + 1], c, s);
            w->apply(c), w->apply(s);
            cols.push_back(c), cols.push_back(s);
        }
        return cols;
    };
    const auto rss_of = [&](const std::vector<double>& p) {
        for (std::size_t j = 1; j < p.size(); j += 2)
            if (!(p[j] > 0)) return std::numeric_limits<double>::infinity();   // a decay time is positive
        return fit(columns_of(p), b).rss;
    };

    double zero = 0;
    for (double v : b) zero += v * v;
    out.rss.push_back(zero), out.lengths.push_back(length(zero, 0));
    std::vector<double> p;
    std::vector<Fit> fits{Fit{}};
    std::vector<std::vector<double>> params{{}};
    for (std::size_t m = 1; m <= spec.max_modes && n > 4 * m; ++m) {
        // the best start for mode m on the grid, the first m - 1 held
        std::optional<std::pair<double, std::pair<double, double>>> best;
        const double ratio = spec.decay_count > 1 ? std::pow(spec.decay_hi / spec.decay_lo, 1.0 / double(spec.decay_count - 1)) : 1.0;
        for (double f = spec.frequency_lo; f <= spec.frequency_hi + 1e-9 * spec.frequency_step; f += spec.frequency_step) {
            double tau = spec.decay_lo;
            for (std::size_t k = 0; k < spec.decay_count; ++k, tau *= ratio) {
                auto q = p;
                q.push_back(f), q.push_back(tau);
                const double r = rss_of(q);
                if (!best || r < best->first) best = {r, {f, tau}};
            }
        }
        if (!best || !std::isfinite(best->first)) break;
        p.push_back(best->second.first), p.push_back(best->second.second);
        double value = 0;
        p = nelder_mead(rss_of, p, spec.refine_iterations, 0.05, 1e-4, &value);
        const Fit f = fit(columns_of(p), b);
        if (!std::isfinite(f.rss)) break;
        out.rss.push_back(f.rss), out.lengths.push_back(length(f.rss, m));
        fits.push_back(f), params.push_back(p);
    }
    std::size_t chosen = 0;
    for (std::size_t m = 1; m < out.lengths.size(); ++m)
        if (out.lengths[m] < out.lengths[chosen]) chosen = m;
    out.rms = std::sqrt(out.rss[chosen] / double(n));
    if (chosen == 0) return out;

    // the modes with their amplitudes, and the errors from the Jacobian of the whitened model in all 4 N numbers
    const auto& q = params[chosen];
    const auto& f = fits[chosen];
    const std::size_t k = 4 * chosen;
    std::vector<std::vector<double>> jac;
    for (std::size_t j = 0; j < chosen; ++j) {
        const double fj = q[2 * j], tj = q[2 * j + 1];
        std::vector<double> c0, s0, c1, s1, c2, s2;
        mode_columns(x, fj, tj, c0, s0);
        mode_columns(x, fj * (1 + 1e-3), tj, c1, s1);
        mode_columns(x, fj, tj * (1 + 1e-3), c2, s2);
        std::vector<double> df(n), dt(n);
        for (std::size_t i = 0; i < n; ++i) {
            const double m0 = f.c[2 * j] * c0[i] + f.c[2 * j + 1] * s0[i];
            df[i] = (f.c[2 * j] * c1[i] + f.c[2 * j + 1] * s1[i] - m0) / (1e-3 * fj);
            dt[i] = (f.c[2 * j] * c2[i] + f.c[2 * j + 1] * s2[i] - m0) / (1e-3 * tj);
        }
        w->apply(df), w->apply(dt);
        jac.push_back(df), jac.push_back(dt);
    }
    for (const auto& c : columns_of(q)) jac.push_back(c);
    std::vector<double> normal(k * k);
    for (std::size_t a = 0; a < k; ++a)
        for (std::size_t c = 0; c < k; ++c) {
            double s = 0;
            for (std::size_t i = 0; i < n; ++i) s += jac[a][i] * jac[c][i];
            normal[a * k + c] = s;
        }
    const auto inv = inverse_spd(normal, k);
    const double s2 = f.rss / double(n - k);
    for (std::size_t j = 0; j < chosen; ++j) {
        Mode m{.frequency = q[2 * j], .decay = q[2 * j + 1], .cosine = f.c[2 * j], .sine = f.c[2 * j + 1]};
        if (!inv.empty()) {
            m.frequency_error = std::sqrt(s2 * inv[(2 * j) * k + 2 * j]);
            m.decay_error = std::sqrt(s2 * inv[(2 * j + 1) * k + 2 * j + 1]);
        }
        out.modes.push_back(m);
    }
    return out;
}

} // namespace spatium::answer
