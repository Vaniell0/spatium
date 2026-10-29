#pragma once

// Synthesis of a particular analytic solution on a region, and of a piecewise one where a single
// expression does not reach the tolerance.
//
// The notation is the library's own: a candidate is a field of the scene IR (io/field.hpp) in
// one variable, so what is found is evaluated by the IR's interpreter and lowered by its emitters
// with nothing written for it. The conditions of the job -- the reference on an interval, the
// tolerance, the operations allowed and how large an expression may grow -- are the set the
// search runs over.
//
// On one interval: every expression of the allowed operations up to `max_size` nodes, over the
// interval's own variable t = (x - centre) * 2 / width, is built from the smaller ones (bottom up), and an expression whose values at the sample points equal an
// earlier one's is dropped -- the paths collapse to distinct states, which is what keeps the
// enumeration finite. The answer is a sum of at most `max_terms` of them with coefficients, the
// terms picked greedily (the one that best explains what is left, the coefficients refitted by
// least squares each time, and each chosen term swapped for any other while that lowers the
// error at the samples, since a greedy pick cannot take a term back), accepted when its error on
// a dense grid of the interval is within the tolerance. That is a certificate of the strength "sample": the error is measured, not bounded.
//
// Piecewise: an interval on which nothing is accepted is split at its middle and each half is
// synthesised alone, down to `min_width`; what remains unreached is reported as such, with the
// interval. The baseline every learned policy for where to split has to beat.

#include <spatium/io/field.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <numbers>
#include <optional>
#include <vector>

namespace rsc {

using Expr = spatium::io::build::Field<double>;

enum class Unary { Sqrt, Sin, Cos };

struct SynthesisSpec {
    double lo = 0, hi = 1;
    double tolerance = 1e-6;  // on |value - reference|, scaled by |reference| when `relative`
    bool relative = true;
    std::size_t max_size = 5; // nodes of one term
    std::size_t max_terms = 4;
    std::size_t samples = 48;   // where terms are compared and fitted
    std::size_t checks = 2001;  // where the answer is certified
    std::vector<Unary> unary{Unary::Sqrt};
    std::size_t max_candidates = 4000;
    std::size_t exhaustive_pairs = 200000; // all sets of one and two terms are tried when there are at most this many
};

struct Term {
    Expr expr;
    std::size_t size = 1;
    std::vector<double> at; // values at the sample points
};

struct Solution {
    double lo = 0, hi = 0;
    Expr expr;
    std::vector<double> coefficients;
    std::size_t terms = 0, nodes = 0;
    double error = 0;     // the largest error on the check grid, in the spec's sense
    bool accepted = false;
};

namespace synthesis_detail {

inline std::vector<double> grid(double lo, double hi, std::size_t n) {
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; ++i) // Chebyshev points: dense at the ends, where errors grow
        x[i] = 0.5 * (lo + hi) - 0.5 * (hi - lo) * std::cos(std::numbers::pi * (double(i) + 0.5) / double(n));
    return x;
}

inline std::vector<double> values(const Expr& e, const std::vector<double>& x) {
    std::vector<double> out(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
        out[i] = e(x[i], 0.0);
    return out;
}

inline bool usable(const std::vector<double>& v) {
    return std::ranges::all_of(v, [](double a) { return std::isfinite(a) && std::abs(a) < 1e12; });
}

// A key that two expressions with the same values share: the values rounded to 1e-10 of their
// scale.
inline std::vector<long long> key(const std::vector<double>& v) {
    double scale = 0;
    for (double a : v)
        scale = std::max(scale, std::abs(a));
    std::vector<long long> k(v.size());
    for (std::size_t i = 0; i < v.size(); ++i)
        k[i] = scale > 0 ? std::llround(v[i] / scale * 1e10) : 0;
    k.push_back(std::llround(std::log2(scale + 1e-300)));
    return k;
}

// Least squares by modified Gram-Schmidt on the columns: no normal equations, whose condition
// number is the square of the basis's, and no size fixed when compiling.
inline std::optional<std::vector<double>> least_squares(const std::vector<std::vector<double>>& basis,
                                                        const std::vector<double>& y) {
    const std::size_t k = basis.size(), n = y.size();
    std::vector<std::vector<double>> q = basis;
    std::vector<std::vector<double>> r(k, std::vector<double>(k, 0.0));
    for (std::size_t j = 0; j < k; ++j) {
        for (std::size_t i = 0; i < j; ++i) {
            double d = 0;
            for (std::size_t t = 0; t < n; ++t)
                d += q[i][t] * q[j][t];
            r[i][j] = d;
            for (std::size_t t = 0; t < n; ++t)
                q[j][t] -= d * q[i][t];
        }
        double norm = 0;
        for (double v : q[j])
            norm += v * v;
        norm = std::sqrt(norm);
        if (!(norm > 1e-13 * std::sqrt(double(n))))
            return std::nullopt; // the new column adds nothing the others do not span
        r[j][j] = norm;
        for (double& v : q[j])
            v /= norm;
    }
    std::vector<double> c(k);
    for (std::size_t j = k; j-- > 0;) {
        double d = 0;
        for (std::size_t t = 0; t < n; ++t)
            d += q[j][t] * y[t];
        for (std::size_t i = j + 1; i < k; ++i)
            d -= r[j][i] * c[i];
        c[j] = d / r[j][j];
    }
    return c;
}

} // namespace synthesis_detail

// Every distinct expression of up to spec.max_size nodes over the interval's variable t and x, the
// constant 1, + - * / and the spec's unary operations, by their values at `x`.
inline std::vector<Term> enumerate(const SynthesisSpec& spec, const std::vector<double>& x) {
    using namespace synthesis_detail;
    std::vector<std::vector<Term>> by_size(spec.max_size + 1);
    std::map<std::vector<long long>, bool> seen;
    std::size_t count = 0;
    auto offer = [&](Expr e, std::size_t size) {
        if (count >= spec.max_candidates)
            return;
        auto v = values(e, x);
        if (!usable(v))
            return;
        if (!seen.emplace(key(v), true).second)
            return;
        by_size[size].push_back({std::move(e), size, std::move(v)});
        ++count;
    };
    // the variable of the interval: t = (x - centre) * 2 / width on [-1, 1], a field of x like any
    // other, so a term in t is a term in x; on a narrow interval powers of x are all but constant
    // and only t keeps them apart
    const double centre = 0.5 * (spec.lo + spec.hi), scale = 2 / (spec.hi - spec.lo);
    offer((Expr::u() - Expr(centre)) * Expr(scale), 1);
    offer(Expr::u(), 1); // and x itself, for the forms that are global (sin x, 4 / b)
    offer(Expr(1.0), 1);
    for (std::size_t s = 2; s <= spec.max_size; ++s) {
        for (const Term& t : std::vector<Term>(by_size[s - 1]))
            for (Unary op : spec.unary)
                offer(op == Unary::Sqrt ? sqrt(t.expr) : op == Unary::Sin ? sin(t.expr) : cos(t.expr), s);
        for (std::size_t a = 1; a + 1 < s; ++a) {
            const std::size_t b = s - 1 - a;
            const auto left = by_size[a], right = by_size[b];
            for (const Term& l : left)
                for (const Term& r : right) {
                    if (a <= b) { // + and * commute: one order is enough
                        offer(l.expr + r.expr, s);
                        offer(l.expr * r.expr, s);
                    }
                    offer(l.expr - r.expr, s);
                    offer(l.expr / r.expr, s);
                }
        }
    }
    std::vector<Term> out;
    for (auto& level : by_size)
        for (Term& t : level)
            out.push_back(std::move(t));
    return out;
}

// One interval: the sum of at most max_terms terms that best explains the reference, accepted
// or not by its error on the check grid.
inline Solution synthesize(const std::function<double(double)>& reference, const SynthesisSpec& spec) {
    using namespace synthesis_detail;
    const auto x = grid(spec.lo, spec.hi, spec.samples);
    std::vector<double> y(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
        y[i] = reference(x[i]);
    const auto terms = enumerate(spec, x);
    const auto checks = grid(spec.lo, spec.hi, spec.checks);
    std::vector<double> truth(checks.size());
    for (std::size_t i = 0; i < checks.size(); ++i)
        truth[i] = reference(checks[i]);

    auto error_of = [&](const Expr& e) {
        double worst = 0;
        for (std::size_t i = 0; i < checks.size(); ++i) {
            const double d = std::abs(e(checks[i], 0.0) - truth[i]);
            worst = std::max(worst, spec.relative ? d / std::max(std::abs(truth[i]), 1e-300) : d);
        }
        return worst;
    };

    Solution best{spec.lo, spec.hi, Expr(0.0), {}, 0, 1, 1e300, false};
    // the squared error at the sample points of the fit on `set`, or infinity
    auto sample_error = [&](const std::vector<std::size_t>& set) {
        std::vector<std::vector<double>> cols;
        for (std::size_t t : set)
            cols.push_back(terms[t].at);
        const auto c = least_squares(cols, y);
        if (!c)
            return 1e300;
        double e = 0;
        for (std::size_t i = 0; i < y.size(); ++i) {
            double fit = 0;
            for (std::size_t j = 0; j < set.size(); ++j)
                fit += (*c)[j] * cols[j][i];
            e += (y[i] - fit) * (y[i] - fit);
        }
        return e;
    };

    // Sets of one and two terms are tried in full (they are cheap and a greedy pick misses the
    // pair whose members each explain little alone); larger sets grow greedily from the best pair.
    std::vector<std::size_t> seed;
    if (terms.size() * terms.size() <= spec.exhaustive_pairs * 2) {
        double top = 1e300;
        for (std::size_t a = 0; a < terms.size(); ++a)
            for (std::size_t b = a; b < terms.size(); ++b) {
                std::vector<std::size_t> set = a == b ? std::vector<std::size_t>{a} : std::vector<std::size_t>{a, b};
                if (set.size() > spec.max_terms)
                    continue;
                const double e = sample_error(set);
                if (e < top * (1 - 1e-12) || (e <= top && set.size() < seed.size()))
                    top = e, seed = std::move(set);
            }
    }
    std::vector<std::size_t> chosen;
    std::vector<double> residual = y;
    for (std::size_t k = seed.empty() ? 1 : seed.size(); k <= spec.max_terms; ++k) {
        if (k == seed.size() && chosen.empty()) {
            chosen = seed;
        } else {
        // the term that best explains what is left, by the cosine of the angle to the residual
        std::size_t pick = terms.size();
        double top = -1;
        for (std::size_t t = 0; t < terms.size(); ++t) {
            if (std::ranges::find(chosen, t) != chosen.end())
                continue;
            double dot = 0, nt = 0;
            for (std::size_t i = 0; i < y.size(); ++i)
                dot += terms[t].at[i] * residual[i], nt += terms[t].at[i] * terms[t].at[i];
            const double score = nt > 0 ? std::abs(dot) / std::sqrt(nt) : -1;
            if (score > top)
                top = score, pick = t;
        }
        if (pick == terms.size())
            break;
        chosen.push_back(pick);
        }
        // the greedy pick cannot take a term back: swap each chosen one for any candidate while
        // that lowers the error at the samples
        double current = sample_error(chosen);
        for (bool improved = true; improved;) {
            improved = false;
            for (std::size_t j = 0; j < chosen.size(); ++j)
                for (std::size_t t = 0; t < terms.size(); ++t) {
                    if (std::ranges::find(chosen, t) != chosen.end())
                        continue;
                    auto trial = chosen;
                    trial[j] = t;
                    const double e = sample_error(trial);
                    if (e < current * (1 - 1e-12)) {
                        chosen = std::move(trial), current = e, improved = true;
                    }
                }
        }
        std::vector<std::vector<double>> basis;
        for (std::size_t t : chosen)
            basis.push_back(terms[t].at);
        const auto c = least_squares(basis, y);
        if (!c)
            break;
        Expr e(0.0);
        std::size_t nodes = 0;
        for (std::size_t j = 0; j < chosen.size(); ++j) {
            e = j == 0 ? Expr((*c)[j]) * terms[chosen[j]].expr : e + Expr((*c)[j]) * terms[chosen[j]].expr;
            nodes += terms[chosen[j]].size + 2; // the coefficient and its product, and a sum
        }
        for (std::size_t i = 0; i < y.size(); ++i) {
            double fit = 0;
            for (std::size_t j = 0; j < chosen.size(); ++j)
                fit += (*c)[j] * basis[j][i];
            residual[i] = y[i] - fit;
        }
        const double err = error_of(e);
        if (err < best.error)
            best = {spec.lo, spec.hi, e, *c, chosen.size(), nodes, err, err <= spec.tolerance};
        if (best.accepted)
            break;
    }
    return best;
}

// Piecewise: each interval that nothing accepts is halved, down to min_width; returns the pieces
// in order, unaccepted ones marked so.
inline std::vector<Solution> synthesize_piecewise(const std::function<double(double)>& reference,
                                                  SynthesisSpec spec, double min_width) {
    std::vector<Solution> out;
    std::function<void(double, double)> go = [&](double lo, double hi) {
        spec.lo = lo, spec.hi = hi;
        Solution s = synthesize(reference, spec);
        if (s.accepted || hi - lo <= min_width) {
            out.push_back(std::move(s));
            return;
        }
        const double mid = 0.5 * (lo + hi);
        go(lo, mid);
        go(mid, hi);
    };
    go(spec.lo, spec.hi);
    return out;
}

} // namespace rsc
