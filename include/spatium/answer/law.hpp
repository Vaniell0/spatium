#pragma once

// A law from data that carries noise: the shortest description of y(x), not the first sum that
// comes within a tolerance.
//
// `synthesize` (synthesis.hpp) answers a question whose values are exact, and stops at the first
// sum of terms within its tolerance. Data measured or computed with a remainder is another job:
// any sum of enough terms fits it, and the law is the one that says the most in the fewest words.
// So the terms are the same -- every expression of the allowed operations over x up to a size,
// those with equal values merged (`enumerate`) -- but a sum of them is fitted by least squares
// weighted by the noise each point carries (`sigma`), and chosen by its description length:
//
//   n ln(RSS / n)  +  ln n for each fitted coefficient  +  its nodes times ln(symbols) for each term,
//
// a term written out symbol by symbol from x, 1, + - * / and the unary operations, so naming a short
// term costs less than naming a long one. Sets of one and two terms are tried in full among the first
// ones, and of three among the shortest.
//
// Each coefficient comes with its standard error, and is then read as a named constant where it is
// one: a combination, with small whole or fractional coefficients, of 1, Euler's gamma, log 2pi,
// pi, 1/pi, log(2pi)/pi, gamma/pi and log 2, within three standard errors. A combination is taken
// only when one that cheap is unlikely to land that close by chance -- the count of combinations
// no costlier, times the width of the window, over the spread of their values, below 1/100 --
// since among a hundred thousand combinations one lands near any number. A coefficient with no
// such reading keeps its digits.

#include <spatium/answer/synthesis.hpp>
#include <spatium/io/field_text.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

namespace spatium::answer {

struct LawSpec {
    std::size_t max_size = 5;      // nodes of one term
    std::size_t max_terms = 3;
    std::vector<Unary> unary{Unary::Log, Unary::Sqrt};
    std::size_t max_candidates = 4000;
    std::size_t pair_pool = 600;   // the first terms among which every pair is tried
    std::size_t triple_pool = 60;  // and every triple
};

// A coefficient read as a named constant: its value and how it is written.
struct NamedConstant {
    double value = 0;
    std::string text;
    double cost = 0;       // bits to write it
    double by_chance = 0;  // the expected count of combinations that cheap landing as close
};

struct Law {
    Expr expr;                                // the law, a recognised coefficient as its exact value
    std::vector<Expr> terms;
    std::vector<double> coefficients, errors; // fitted, and their standard errors
    std::vector<std::optional<NamedConstant>> named;
    double length = std::numeric_limits<double>::infinity();   // its description length
    double rms = 0;                           // of the residual in units of sigma
    std::string text;                         // the law as a person writes it, constants by name
};

namespace law_detail {

struct Fit {
    std::vector<double> c, se;
    double rss = std::numeric_limits<double>::infinity();
};

// Least squares by modified Gram-Schmidt, with the standard errors from R: (A^T A)^-1 = R^-1 R^-T.
inline Fit fit(const std::vector<std::vector<double>>& cols, const std::vector<double>& y) {
    const std::size_t k = cols.size(), n = y.size();
    Fit out;
    if (k == 0 || n <= k) return out;
    auto q = cols;
    std::vector<std::vector<double>> r(k, std::vector<double>(k, 0.0));
    for (std::size_t j = 0; j < k; ++j) {
        for (std::size_t i = 0; i < j; ++i) {
            double d = 0;
            for (std::size_t t = 0; t < n; ++t) d += q[i][t] * q[j][t];
            r[i][j] = d;
            for (std::size_t t = 0; t < n; ++t) q[j][t] -= d * q[i][t];
        }
        double norm = 0;
        for (double v : q[j]) norm += v * v;
        norm = std::sqrt(norm);
        double scale = 0;
        for (double v : cols[j]) scale += v * v;
        if (!(norm > 1e-10 * std::sqrt(scale))) return out;   // a column the others nearly span
        r[j][j] = norm;
        for (double& v : q[j]) v /= norm;
    }
    out.c.assign(k, 0.0);
    for (std::size_t j = k; j-- > 0;) {
        double d = 0;
        for (std::size_t t = 0; t < n; ++t) d += q[j][t] * y[t];
        for (std::size_t i = j + 1; i < k; ++i) d -= r[j][i] * out.c[i];
        out.c[j] = d / r[j][j];
    }
    out.rss = 0;
    for (std::size_t t = 0; t < n; ++t) {
        double f = 0;
        for (std::size_t j = 0; j < k; ++j) f += out.c[j] * cols[j][t];
        out.rss += (y[t] - f) * (y[t] - f);
    }
    // R^-1 by back substitution; the variance of c_j is s^2 times the squared norm of its row
    std::vector<std::vector<double>> inv(k, std::vector<double>(k, 0.0));
    for (std::size_t j = 0; j < k; ++j) {
        inv[j][j] = 1 / r[j][j];
        for (std::size_t i = j; i-- > 0;) {
            double d = 0;
            for (std::size_t m = i + 1; m <= j; ++m) d += r[i][m] * inv[m][j];
            inv[i][j] = -d / r[i][i];
        }
    }
    const double s2 = out.rss / static_cast<double>(n - k);
    out.se.assign(k, 0.0);
    for (std::size_t i = 0; i < k; ++i) {
        double v = 0;
        for (std::size_t j = i; j < k; ++j) v += inv[i][j] * inv[i][j];
        out.se[i] = std::sqrt(s2 * v);
    }
    return out;
}

// a named number: what it is written above and below a fraction's bar
struct Basis {
    const char* over;
    const char* under;
    double value;
};

inline const std::vector<Basis>& basis() {
    static const std::vector<Basis> b{
        {"", "", 1.0},
        {"gamma", "", std::numbers::egamma},
        {"log(2 pi)", "", std::log(2 * std::numbers::pi)},
        {"pi", "", std::numbers::pi},
        {"", "pi", 1 / std::numbers::pi},
        {"log(2 pi)", "pi", std::log(2 * std::numbers::pi) / std::numbers::pi},
        {"gamma", "pi", std::numbers::egamma / std::numbers::pi},
        {"log(2)", "", std::numbers::ln2},
    };
    return b;
}

struct Combination {
    double value, cost;
    std::string text;
};

// p / q times a basis element, as a person writes it: 3 gamma, log(2 pi)/(2 pi), 1/2
inline std::string scaled(long long p, long long q, const Basis& what, bool first) {
    const std::string sign = p < 0 ? (first ? "-" : " - ") : (first ? "" : " + ");
    const long long a = std::abs(p);
    const std::string over = std::string(what.over), under = std::string(what.under);
    std::string top = over.empty() ? std::to_string(a) : a == 1 ? over : std::format("{} {}", a, over);
    std::string bottom = under.empty() ? (q == 1 ? "" : std::to_string(q)) : q == 1 ? under : std::format("({} {})", q, under);
    return sign + (bottom.empty() ? top : top + "/" + bottom);
}

// Every combination of one to three basis elements: one with a fraction p / q, |p| <= 12 and q in
// 1, 2, 3, 4, 6, 8; two with |p| <= 6, q in 1, 2, 4; three with whole |p| <= 3. Its cost is the
// bits to name its elements and write their coefficients.
inline const std::vector<Combination>& combinations() {
    static const std::vector<Combination> all = [] {
        std::vector<Combination> out;
        const auto& b = basis();
        const double name = std::log2(static_cast<double>(b.size()));
        const auto bits = [](long long p, long long q) { return std::log2(2.0 * double(std::abs(p)) + 1) + std::log2(double(q)); };
        const auto reduced = [](long long p, long long q) { return std::gcd(std::abs(p), q) == 1; };
        for (std::size_t i = 0; i < b.size(); ++i)
            for (long long q : {1, 2, 3, 4, 6, 8})
                for (long long p = -12; p <= 12; ++p)
                    if (p != 0 && reduced(p, q))
                        out.push_back({double(p) / double(q) * b[i].value, name + bits(p, q), scaled(p, q, b[i], true)});
        for (std::size_t i = 0; i < b.size(); ++i)
            for (std::size_t j = i + 1; j < b.size(); ++j)
                for (long long q1 : {1, 2, 4})
                    for (long long p1 = -6; p1 <= 6; ++p1)
                        for (long long q2 : {1, 2, 4})
                            for (long long p2 = -6; p2 <= 6; ++p2)
                                if (p1 && p2 && reduced(p1, q1) && reduced(p2, q2))
                                    out.push_back({double(p1) / double(q1) * b[i].value + double(p2) / double(q2) * b[j].value,
                                                   2 * name + bits(p1, q1) + bits(p2, q2),
                                                   scaled(p1, q1, b[i], true) + scaled(p2, q2, b[j], false)});
        for (std::size_t i = 0; i < b.size(); ++i)
            for (std::size_t j = i + 1; j < b.size(); ++j)
                for (std::size_t m = j + 1; m < b.size(); ++m)
                    for (long long p1 = -3; p1 <= 3; ++p1)
                        for (long long p2 = -3; p2 <= 3; ++p2)
                            for (long long p3 = -3; p3 <= 3; ++p3)
                                if (p1 && p2 && p3)
                                    out.push_back({double(p1) * b[i].value + double(p2) * b[j].value + double(p3) * b[m].value,
                                                   3 * name + bits(p1, 1) + bits(p2, 1) + bits(p3, 1),
                                                   scaled(p1, 1, b[i], true) + scaled(p2, 1, b[j], false) +
                                                       scaled(p3, 1, b[m], false)});
        std::ranges::sort(out, {}, &Combination::cost);
        return out;
    }();
    return all;
}

// The cheapest combination within three standard errors of c, if one that cheap is unlikely to land
// there by chance.
inline std::optional<NamedConstant> named(double c, double se) {
    if (!(se > 0) || !std::isfinite(c)) return std::nullopt;
    const auto& all = combinations();
    double lo = all.front().value, hi = lo;
    for (const auto& k : all) lo = std::min(lo, k.value), hi = std::max(hi, k.value);
    const double window = 6 * se, spread = hi - lo;
    for (std::size_t i = 0; i < all.size(); ++i) {
        if (std::abs(all[i].value - c) > 3 * se) continue;
        // the combinations no costlier than this one
        std::size_t cheaper = i + 1;
        while (cheaper < all.size() && all[cheaper].cost <= all[i].cost) ++cheaper;
        const double chance = double(cheaper) * window / spread;
        if (chance >= 0.01) return std::nullopt;
        return NamedConstant{all[i].value, all[i].text, all[i].cost, chance};
    }
    return std::nullopt;
}

} // namespace law_detail

inline Law synthesize_law(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& sigma,
                          const LawSpec& spec = {}) {
    using namespace law_detail;
    namespace bd = spatium::io::build;
    const std::size_t n = x.size();
    SynthesisSpec s{.lo = *std::ranges::min_element(x), .hi = *std::ranges::max_element(x), .max_size = spec.max_size,
                    .unary = spec.unary, .max_candidates = spec.max_candidates, .scaled_variable = false};
    const std::vector<Term> terms = enumerate(s, x);
    // each term's column and y, divided by the noise at each point
    std::vector<std::vector<double>> col(terms.size(), std::vector<double>(n));
    std::vector<double> b(n);
    for (std::size_t i = 0; i < n; ++i) b[i] = y[i] / sigma[i];
    for (std::size_t t = 0; t < terms.size(); ++t)
        for (std::size_t i = 0; i < n; ++i) col[t][i] = terms[t].at[i] / sigma[i];
    const double ln_n = std::log(double(n));
    // a term written out, a symbol a node, each one of x, 1, + - * / and the unary operations
    const double symbol = std::log(6.0 + double(spec.unary.size()));
    const auto length = [&](const std::vector<std::size_t>& set, const Fit& f) {
        double l = double(n) * std::log(std::max(f.rss, 1e-300) / double(n)) + double(set.size()) * ln_n;
        for (std::size_t t : set) l += double(terms[t].size) * symbol;
        return l;
    };
    Law best;
    std::vector<std::size_t> chosen;
    Fit chosen_fit;
    std::vector<std::pair<double, std::vector<std::size_t>>> pairs_by_length;
    const auto offer = [&](std::vector<std::size_t> set) {
        if (set.size() > spec.max_terms) return;
        std::vector<std::vector<double>> cols;
        for (std::size_t t : set) cols.push_back(col[t]);
        Fit f = fit(cols, b);
        if (!std::isfinite(f.rss)) return;
        const double l = length(set, f);
        if (set.size() == 2) pairs_by_length.emplace_back(l, set);
        if (l < best.length) best.length = l, chosen = std::move(set), chosen_fit = std::move(f);
    };
    const std::size_t pairs = std::min(terms.size(), spec.pair_pool), triples = std::min(terms.size(), spec.triple_pool);
    for (std::size_t a = 0; a < pairs; ++a) {
        offer({a});
        for (std::size_t c = a + 1; c < pairs; ++c) offer({a, c});
    }
    if (spec.max_terms >= 3) {
        for (std::size_t a = 0; a < triples; ++a)
            for (std::size_t c = a + 1; c < triples; ++c)
                for (std::size_t d = c + 1; d < triples; ++d) offer({a, c, d});
        // the best pairs, each with a third term from the whole pool: a law of three whose terms are
        // not all among the shortest -- x log x, x and 1 for the count of zeros
        std::ranges::sort(pairs_by_length, {}, &std::pair<double, std::vector<std::size_t>>::first);
        const auto top = std::vector(pairs_by_length.begin(), pairs_by_length.begin() + std::min<std::ptrdiff_t>(20, std::ssize(pairs_by_length)));
        for (const auto& [l, pair] : top)
            for (std::size_t d = 0; d < pairs; ++d)
                if (d != pair[0] && d != pair[1]) offer({pair[0], pair[1], d});
    }
    if (chosen.empty()) return best;
    best.coefficients = chosen_fit.c, best.errors = chosen_fit.se;
    best.rms = std::sqrt(chosen_fit.rss / double(n));
    best.named.assign(chosen.size(), std::nullopt);
    for (std::size_t j = 0; j < chosen.size(); ++j) best.terms.push_back(terms[chosen[j]].expr);
    // One coefficient at a time, the surest first: once it is fixed at its named value the others
    // are fitted again without it, and their errors shrink where its term and theirs were nearly
    // the same column -- x and x log x over a few decades.
    for (bool more = true; more;) {
        more = false;
        std::vector<std::size_t> free;
        std::vector<double> rest = b;
        for (std::size_t j = 0; j < chosen.size(); ++j) {
            if (!best.named[j]) { free.push_back(j); continue; }
            for (std::size_t i = 0; i < n; ++i) rest[i] -= best.named[j]->value * col[chosen[j]][i];
        }
        if (free.empty()) break;
        std::vector<std::vector<double>> cols;
        for (std::size_t j : free) cols.push_back(col[chosen[j]]);
        const Fit f = fit(cols, rest);
        if (!std::isfinite(f.rss)) break;
        std::optional<NamedConstant> surest;
        std::size_t at = 0;
        for (std::size_t m = 0; m < free.size(); ++m) {
            best.coefficients[free[m]] = f.c[m], best.errors[free[m]] = f.se[m];
            if (auto k = named(f.c[m], f.se[m]); k && (!surest || k->by_chance < surest->by_chance)) surest = k, at = free[m];
        }
        if (surest) best.named[at] = surest, more = true;
    }
    // the law: each coefficient its named value where it has one, and the text with its names
    std::string text;
    for (std::size_t j = 0; j < chosen.size(); ++j) {
        const double c = best.named[j] ? best.named[j]->value : best.coefficients[j];
        const Expr term = Expr(c) * best.terms[j];
        best.expr = j == 0 ? term : best.expr + term;
        const std::string coefficient = best.named[j] ? "(" + best.named[j]->text + ")" : std::format("{:.6g}", c);
        text += (j ? " + " : "") + coefficient + " * " + bd::to_text(best.terms[j], "x").value_or("?");
    }
    best.text = text;
    return best;
}

} // namespace spatium::answer
