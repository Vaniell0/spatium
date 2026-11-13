#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/interval.hpp>
#  include <spatium/core/rational.hpp>
#  include <spatium/io/field.hpp>
#  include <spatium/io/field_parse.hpp>
#  include <spatium/io/field_pod.hpp>
#  include <spatium/io/field_text.hpp>
#  include <algorithm>
#  include <cctype>
#  include <cmath>
#  include <cstdint>
#  include <map>
#  include <numbers>
#  include <optional>
#  include <string>
#  include <string_view>
#  include <utility>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io::build {

// ── A field written in its shortest form ─────────────────────────
//
// `shortest(f)` is f written as briefly as this file knows how, or f itself: a rewrite is kept
// only when what a page shows of it is no longer (`shown_length`, the marks of its LaTeX) -- of two
// as long, the one written here, a sum with its terms in order and a root out of a denominator,
// when its text is no longer either, so sqrt(50) stays and x^0.5 is not x^(1 / 2) -- and its value
// holds (`same_values`). The rewrite reads f as a sum
// of terms, each a fraction times a product of atoms raised to whole powers -- a name, pi, the
// square root of a squarefree whole number, a function of the language or a power with its
// arguments written shortest in turn, a denominator that is no single term -- with the fractions
// exact (core/rational.hpp). So like terms are collected, (x + x) * x + x^2 is 3 * x^2 and
// sin(x) * cos(x) + sin(x) * cos(x) is 2 * sin(x) * cos(x); terms over one denominator share it,
// (1 / x - log(x) / x) / x is (1 - log(x)) / x^2; a root is taken out of its square factors and its
// square is a number, sqrt(24) + sqrt(150) is 7 * sqrt(6); and a sine or cosine of a multiple of
// pi with a denominator of 1, 2, 3, 4 or 6 is its exact value, sin(pi / 6) is 1/2.
//
// What makes a rewrite safe is the check, not the algebra: collecting terms cancels x / x to 1,
// which has a value at 0 where x / x has none, so the values are compared where either has none
// too, and such a rewrite is refused whole. A fraction that does not fit in 64 bits, a power past
// 64 or a sum past 256 terms leaves f as it is.

namespace detail {

namespace shortest_impl {

using F = Field<double>;

// a product of atoms, each with its whole power, sorted by atom
using Monomial = std::vector<std::pair<int, int>>;
using Poly = std::map<Monomial, Rational>;

// What an atom is, and so where it stands in a product: a constant before a name, a name before a
// function.
enum class Kind { Pi, Radical, Name, Other };

struct Atom {
    F field;
    std::string key;
    Kind kind = Kind::Other;
    long long square = 0;   // a radical's: sqrt(square)^2 is square
};

class Reader {
public:
    std::vector<Atom> atoms;
    bool failed = false;

    std::optional<Poly> read(const F& whole) {
        const F f = whole.reachable_from_root();
        std::vector<Poly> at(f.size());
        for (std::size_t i = 0; i < f.size() && !failed; ++i) at[i] = op(f, i, at);
        if (failed) return std::nullopt;
        return at.back();
    }

    F write(const Poly& p) const {
        if (p.empty()) return F(0.0);
        // one denominator for every term: each atom to the largest power any term divides by
        std::map<int, int> below;
        for (const auto& [m, c] : p)
            for (const auto& [a, e] : m)
                if (e < 0) below[a] = std::max(below[a], -e);
        std::vector<std::pair<Monomial, Rational>> terms;
        for (const auto& [m, c] : p) {
            std::map<int, int> e;
            for (const auto& [a, k] : m) e[a] += k;
            for (const auto& [a, k] : below) e[a] += k;
            Monomial up;
            for (const auto& [a, k] : e)
                if (k != 0) up.emplace_back(a, k);
            terms.emplace_back(std::move(up), c);
        }
        // the highest degree in the names first, then by text
        const auto degree = [&](const Monomial& m) {
            int d = 0;
            for (const auto& [a, k] : m) d += atoms[static_cast<std::size_t>(a)].kind >= Kind::Name ? k : 0;
            return d;
        };
        // a number after a root or pi, as sqrt(2) + 2 is written
        std::ranges::stable_sort(terms, [&](const auto& x, const auto& y) {
            const int dx = degree(x.first), dy = degree(y.first);
            return dx != dy ? dx > dy : !x.first.empty() && y.first.empty();
        });
        // a term led by a minus does not lead the sum, where one is not
        if (terms.size() > 1 && terms.front().second < Rational(0)) {
            const auto first = std::ranges::find_if(terms, [](const auto& t) { return t.second > Rational(0); });
            if (first != terms.end()) std::rotate(terms.begin(), first, first + 1);
        }
        F sum;
        bool started = false;
        for (const auto& [m, c] : terms) {
            const bool minus = started && c < Rational(0);
            const F t = term(m, minus ? -c : c);
            sum = !started ? t : minus ? std::move(sum) - t : std::move(sum) + t;
            started = true;
        }
        if (below.empty()) return sum;
        Monomial down;
        for (const auto& [a, k] : below) down.emplace_back(a, k);
        return std::move(sum) / product(down);
    }

private:
    static constexpr std::size_t most_terms = 256;

    void fail() { failed = true; }

    int atom(F field, Kind kind = Kind::Other, long long square = 0) {
        std::string key = to_text(field).value_or("");
        if (key.empty()) key = "#" + std::to_string(atoms.size());   // no notation: an atom of its own
        for (std::size_t k = 0; k < atoms.size(); ++k)
            if (atoms[k].key == key) return static_cast<int>(k);
        atoms.push_back({std::move(field), std::move(key), kind, square});
        return static_cast<int>(atoms.size() - 1);
    }

    static Poly constant(Rational c) {
        Poly p;
        if (c != Rational(0)) p[{}] = c;
        return p;
    }

    Poly of_atom(int a, int power = 1) {
        Poly p;
        p[Monomial{{a, power}}] = Rational(1);
        return normal(p);
    }

    // a radical's power brought below 2 -- sqrt(6)^2 is 6, 1 / sqrt(2) is sqrt(2) / 2 -- and zero
    // terms dropped
    Poly normal(const Poly& p) {
        Poly out;
        for (auto [m, c] : p) {
            Monomial kept;
            for (const auto& [a, e] : m) {
                const auto& at = atoms[static_cast<std::size_t>(a)];
                int k = e;
                if (at.kind == Kind::Radical) {
                    const int r = ((k % 2) + 2) % 2, half = (k - r) / 2;
                    for (int j = 0; j < std::abs(half); ++j) c = half > 0 ? c * Rational(at.square) : c / Rational(at.square);
                    k = r;
                }
                if (k != 0) kept.emplace_back(a, k);
            }
            if (!c.ok()) return fail(), Poly{};
            Rational& slot = out[kept];
            slot += c;
            if (!slot.ok()) return fail(), Poly{};
            if (slot == Rational(0)) out.erase(kept);
        }
        if (out.size() > most_terms) fail();
        return out;
    }

    Poly add(const Poly& x, const Poly& y, Rational sign = Rational(1)) {
        Poly s = x;
        for (const auto& [m, c] : y) {
            Rational& slot = s[m];
            slot += sign * c;
        }
        return normal(s);
    }

    Poly mul(const Poly& x, const Poly& y) {
        if (x.size() * y.size() > most_terms) return fail(), Poly{};
        Poly s;
        for (const auto& [m1, c1] : x)
            for (const auto& [m2, c2] : y) {
                std::map<int, int> e;
                for (const auto& [a, k] : m1) e[a] += k;
                for (const auto& [a, k] : m2) e[a] += k;
                Monomial m;
                for (const auto& [a, k] : e)
                    if (k != 0) m.emplace_back(a, k);
                Rational& slot = s[m];
                slot += c1 * c2;
            }
        return normal(s);
    }

    // a single term to a whole power; none where p is a sum
    std::optional<Poly> term_power(const Poly& p, int n) {
        if (p.size() != 1) return std::nullopt;
        const auto& [m, c] = *p.begin();
        Rational r(1);
        for (int k = 0; k < std::abs(n); ++k) r = n > 0 ? r * c : r / c;
        Monomial out;
        for (const auto& [a, e] : m) out.emplace_back(a, e * n);
        return normal(Poly{{out, r}});
    }

    static std::optional<Rational> as_constant(const Poly& p) {
        if (p.empty()) return Rational(0);
        if (p.size() == 1 && p.begin()->first.empty()) return p.begin()->second;
        return std::nullopt;
    }

    // c * pi, its c
    std::optional<Rational> times_pi(const Poly& p) const {
        if (p.empty()) return Rational(0);
        if (p.size() != 1) return std::nullopt;
        const auto& [m, c] = *p.begin();
        if (m.size() == 1 && m[0].second == 1 && atoms[static_cast<std::size_t>(m[0].first)].kind == Kind::Pi) return c;
        return std::nullopt;
    }

    // sin(c pi), or cos(c pi) as sin((c + 1/2) pi), for c of denominator 1, 2, 3, 4 or 6
    std::optional<Poly> exact_sine(Rational c) {
        const long long d = c.den();
        if (d != 1 && d != 2 && d != 3 && d != 4 && d != 6) return std::nullopt;
        long long k = (c.num() * (12 / d)) % 24;   // c in twelfths of pi, on [0, 24)
        if (k < 0) k += 24;
        const int sign = k >= 12 ? -1 : 1;
        k %= 12;
        if (k > 6) k = 12 - k;   // sin(pi - x) = sin(x)
        Poly p;
        switch (k) {
            case 0: break;
            case 1: case 5: return std::nullopt;   // pi / 12: no value this file writes
            case 2: p = constant(Rational(1, 2)); break;
            case 3: p = mul(constant(Rational(1, 2)), of_atom(radical(2))); break;
            case 4: p = mul(constant(Rational(1, 2)), of_atom(radical(3))); break;
            case 6: p = constant(Rational(1)); break;
        }
        return sign > 0 ? p : mul(constant(Rational(-1)), p);
    }

    int radical(long long m) { return atom(sqrt(F(static_cast<double>(m))), Kind::Radical, m); }

    // the q-th root of a fraction where it is one: 8 to 2, 27/8 to 3/2
    static std::optional<Rational> exact_root(Rational c, long long q) {
        if (q < 2 || q > 64 || c < Rational(0)) return std::nullopt;
        const auto root = [q](long long n) -> std::optional<long long> {
            const auto r = static_cast<long long>(std::llround(std::pow(static_cast<double>(n), 1.0 / static_cast<double>(q))));
            for (long long t = std::max(0LL, r - 1); t <= r + 1; ++t) {
                __int128 x = 1;
                for (long long k = 0; k < q && x <= n; ++k) x *= t;
                if (x == n) return t;
            }
            return std::nullopt;
        };
        const auto n = root(c.num()), d = root(c.den());
        if (!n || !d) return std::nullopt;
        return Rational(*n, *d);
    }

    // sqrt of a fraction p / q: sqrt(p q) / q, its square factors taken out
    std::optional<Poly> root_of(Rational c) {
        if (c < Rational(0)) return std::nullopt;
        const __int128 n = static_cast<__int128>(c.num()) * c.den();
        if (n > (__int128(1) << 62)) return std::nullopt;
        long long rest = static_cast<long long>(n), out = 1;
        for (long long k = 2; k * k <= rest; ++k)
            while (rest % (k * k) == 0) rest /= k * k, out *= k;
        const Poly whole = constant(Rational(out, c.den()));
        return rest == 1 ? whole : mul(whole, of_atom(radical(rest)));
    }

    Poly op(const F& f, std::size_t i, const std::vector<Poly>& at) {
        const auto& n = f.op(i);
        const auto as_field = [&](std::uint32_t j) { return write(at[j]); };
        switch (n.op) {
            case Op::Const: {
                if (std::abs(n.value) == std::numbers::pi)
                    return mul(constant(Rational(n.value > 0 ? 1 : -1)), of_atom(atom(F(std::numbers::pi), Kind::Pi)));
                const Rational r(n.value);
                return r.ok() ? constant(r) : of_atom(atom(f.subfield(i)));
            }
            case Op::U: case Op::V: case Op::Symbol: return of_atom(atom(f.subfield(i), Kind::Name));
            case Op::Add: return add(at[n.a], at[n.b]);
            case Op::Sub: return add(at[n.a], at[n.b], Rational(-1));
            case Op::Mul: return mul(at[n.a], at[n.b]);
            case Op::Div: {
                if (const auto inverse = term_power(at[n.b], -1)) return mul(at[n.a], *inverse);
                return mul(at[n.a], of_atom(atom(as_field(n.b)), -1));
            }
            case Op::Pow: {
                const auto k = as_constant(at[n.b]);
                if (k && k->den() == 1 && std::abs(k->num()) <= 64) {
                    const int e = static_cast<int>(k->num());
                    if (const auto t = term_power(at[n.a], e)) return *t;
                    if (e > 0) {
                        Poly p = constant(Rational(1));
                        for (int j = 0; j < e && !failed; ++j) p = mul(p, at[n.a]);
                        return p;
                    }
                    return of_atom(atom(as_field(n.a)), e);
                }
                // a number to a fraction p / q: its q-th root where it is a fraction, a square root
                // out of its square factors where q is 2 -- 8^(1/3) is 2, 2^(1/2) is sqrt(2)
                if (const auto c = as_constant(at[n.a]); k && c && *c > Rational(0) && std::abs(k->num()) <= 64) {
                    const int e = static_cast<int>(k->num());
                    if (const auto r = exact_root(*c, k->den()))
                        if (const auto t = term_power(constant(*r), e)) return *t;
                    if (k->den() == 2)
                        if (const auto r = root_of(*c))
                            if (const auto t = term_power(*r, e)) return *t;
                }
                return of_atom(atom(pow(as_field(n.a), as_field(n.b))));
            }
            case Op::Sqrt:
                if (const auto c = as_constant(at[n.a]))
                    if (const auto r = root_of(*c)) return *r;
                return of_atom(atom(sqrt(as_field(n.a))));
            case Op::Sin: case Op::Cos:
                if (const auto c = times_pi(at[n.a]))
                    if (const auto s = exact_sine(n.op == Op::Sin ? *c : *c + Rational(1, 2))) return *s;
                break;
            default: break;
        }
        // a function of the language, its arguments written shortest; anything else as it was
        if (const Function* fn = function_of_op(n.op)) {
            std::vector<F> args{as_field(n.a)};
            if (fn->arity > 1) args.push_back(as_field(n.b));
            return of_atom(atom(apply(*fn, std::move(args))));
        }
        return of_atom(atom(f.subfield(i)));
    }

    // the atoms of m multiplied left to right after `first`, or alone when there is none
    F product(const Monomial& m, std::optional<F> first = std::nullopt) const {
        std::vector<std::pair<int, int>> order(m.begin(), m.end());
        std::ranges::stable_sort(order, [&](const auto& x, const auto& y) {
            const auto &a = atoms[static_cast<std::size_t>(x.first)], &b = atoms[static_cast<std::size_t>(y.first)];
            return a.kind != b.kind ? a.kind < b.kind : a.key < b.key;
        });
        F out = first ? *first : F(1.0);
        bool started = first.has_value();
        for (const auto& [a, e] : order) {
            const F& base = atoms[static_cast<std::size_t>(a)].field;
            F power = base;
            for (int k = 1; k < e; ++k) power = std::move(power) * base;   // one op read again: written base^e
            out = started ? std::move(out) * power : power;
            started = true;
        }
        return out;
    }

    // c times the product, c's denominator written below it
    F term(const Monomial& m, Rational c) const {
        const long long p = c.num(), q = c.den();
        F top;
        if (m.empty()) top = F(static_cast<double>(p));
        else if (p == 1) top = product(m);
        else top = product(m, F(static_cast<double>(p)));
        return q == 1 ? top : std::move(top) / F(static_cast<double>(q));
    }
};

// what a value is compared at: as many points as names, spread over both signs and through 0
inline constexpr double sample_values[] = {0.37, -1.21, 2.05, -0.58, 1.3, 0.0, 3.7, -2.4, 0.91, 1.0, -1.0, 0.5, 6.2, -4.75};

// a command of LaTeX as the marks it shows: none for \left and \right, one for a fraction's bar, a
// root, a dot or a Greek letter, a function's name by its letters
inline std::size_t command_length(std::string_view name) {
    if (name == "left" || name == "right" || name == ",") return 0;
    if (name == "frac" || name == "sqrt" || name == "cdot" || name.size() <= 2 || name == "infty") return name == "ln" ? 2 : 1;
    static constexpr std::string_view greek[] = {"alpha", "beta", "gamma", "delta", "epsilon", "theta", "lambda", "sigma", "omega", "phi", "psi"};
    return std::ranges::find(greek, name) != std::end(greek) ? 1 : name.size();
}

inline bool same_value(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
    if (std::isinf(a) || std::isinf(b)) return a == b;
    return std::abs(a - b) <= 1e-9 * std::max({1.0, std::abs(a), std::abs(b)});
}

} // namespace shortest_impl

} // namespace detail

// What a page shows of a formula, counted: the marks of its LaTeX, a command as the marks it draws,
// and braces, spaces, ^ and _ as none. (x + 1)^2 is 6 and x^2 + 2x + 1 is 7.
inline std::size_t shown_length(std::string_view latex) {
    std::size_t n = 0;
    for (std::size_t i = 0; i < latex.size();) {
        const char c = latex[i];
        if (c == '\\') {
            std::size_t end = i + 1;
            while (end < latex.size() && std::isalpha(static_cast<unsigned char>(latex[end]))) ++end;
            if (end == i + 1) end = std::min(latex.size(), i + 2);   // \, \{ : one mark or none
            n += detail::shortest_impl::command_length(latex.substr(i + 1, end - i - 1));
            i = end;
            continue;
        }
        if (c != '{' && c != '}' && c != ' ' && c != '^' && c != '_') ++n;
        ++i;
    }
    return n;
}

// Whether g has f's values: at fourteen points, every name of either given a value from a fixed
// list, the two agree where both have a value (to 1e-9 of their size) and have none at the same
// points; then on a small box about each point, computed in intervals, the two overlap.
inline bool same_values(const Field<double>& f, const Field<double>& g) {
    namespace si = detail::shortest_impl;
    std::vector<std::string> names = f.free_symbols();
    for (auto& n : g.free_symbols())
        if (std::ranges::find(names, n) == names.end()) names.push_back(std::move(n));
    const auto pf = lower(f), pg = lower(g);
    if (!pf || !pg) return false;
    constexpr std::size_t count = std::size(si::sample_values);
    for (std::size_t k = 0; k < count; ++k) {
        const auto value = [&](std::size_t j) { return si::sample_values[(k + 5 * j) % count]; };
        const auto slots = [&](const auto& pod) {
            std::vector<double> s(pod->symbols.size(), 0.0);
            for (std::size_t j = 0; j < names.size(); ++j)
                for (std::size_t m = 0; m < s.size(); ++m)
                    if (pod->symbols[m] == names[j]) s[m] = value(j);
            return s;
        };
        const double u = value(names.size()), v = value(names.size() + 1);
        if (!si::same_value(interpret(*pf, u, v, slots(pf)), interpret(*pg, u, v, slots(pg)))) return false;
    }
    using I = spatium::Interval<double>;
    if (beyond<I>(*pf) || beyond<I>(*pg)) return true;   // no box to compute in: the points are the check
    const auto in_f = converted<I>(*pf), in_g = converted<I>(*pg);
    for (std::size_t k = 0; k < count; ++k) {
        const auto box = [&](std::size_t j) { const double x = si::sample_values[(k + 5 * j) % count]; return I{x - 1e-6, x + 1e-6}; };
        const auto slots = [&](const auto& pod) {
            std::vector<I> s(pod->symbols.size(), I{0});
            for (std::size_t j = 0; j < names.size(); ++j)
                for (std::size_t m = 0; m < s.size(); ++m)
                    if (pod->symbols[m] == names[j]) s[m] = box(j);
            return s;
        };
        const spatium::Undecided undecided;
        const I a = interpret(in_f, box(names.size()), box(names.size() + 1), slots(pf));
        const I b = interpret(in_g, box(names.size()), box(names.size() + 1), slots(pg));
        if (undecided.raised() || a.is_empty() || b.is_empty()) continue;
        if (a.hi < b.lo || b.hi < a.lo) return false;
    }
    return true;
}

inline Field<double> shortest(const Field<double>& f) {
    detail::shortest_impl::Reader r;
    const auto p = r.read(f);
    if (!p) return f;
    const Field<double> g = r.write(*p);
    const auto before = to_latex(f), after = to_latex(g);
    if (!before || !after || *after == *before) return f;
    const std::size_t b = shown_length(*before), a = shown_length(*after);
    if (a > b || (a == b && to_text(g)->size() > to_text(f)->size())) return f;
    return same_values(f, g) ? g : f;
}

} // namespace spatium::io::build
