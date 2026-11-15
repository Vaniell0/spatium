#pragma once

// A definition of a document given its value: the one door the calculator, a scene and the
// gallery go through.
//
// `evaluate(doc, "x", chosen)` puts the parameters in `chosen` in as numbers, and what the
// definition still reads then decides what its value is. Reading no parameter, it is a number: a
// value evaluated, a root solved in closed form (every root, or those on its bracket) or found on
// its bracket, an integral computed. Reading one, it is a function of it on the parameter's
// interval: a value is its own expression, and a root or an integral is synthesised and
// certified (query.hpp). Reading more, it is refused, since an answer is a function of one
// parameter. So
//
//   p in [0.5, 4]
//   x: x^2 - p = 0 on [0, 10]
//
// gives x as sqrt(p) with nothing chosen, as the calculator shows it, and as the number sqrt(2)
// with p = 2 chosen, as a scene draws it.
//
// Every route that can give the definition its value is tried and kept as a candidate, with its
// cost in evaluations of the definition's field or the reason it gives none, and one is taken.
// Which one is decided in one place, `take`, by a fixed order of preference -- a closed form
// before a root on a bracket, since its steps are exact and it finds every root. The candidates
// not taken stay in the record, so what was chosen among, and at what cost, can be read back
// from it, not only the choice.
//
// How it computes is the engine's (`Engine`), apart from what is asked: the arithmetic -- a
// double; a float, as a device computes; a fraction, in which 1/3 + 1/6 is 1/2; a decimal, in
// which ten times 0.1 is 1; Real50 -- the tolerance of a synthesis, and whether an arithmetic
// that cannot compute something refuses (`strict`) or hands it to a double. Each candidate says
// which arithmetic it computed in, so a record holds the engine with the answer. The arithmetic
// is the interpreter's own parameter (io/field_pod.hpp), one program run in each.

#include <spatium/answer/closed_form.hpp>
#include <spatium/answer/query.hpp>
#include <spatium/algebra/quadrature.hpp>
#include <spatium/algebra/roots.hpp>
#include <spatium/algebra/series.hpp>
#include <spatium/core/decimal.hpp>
#include <spatium/core/interval.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/core/rational.hpp>
#include <spatium/io/document.hpp>
#include <spatium/io/field_pod.hpp>
#include <spatium/io/field_shortest.hpp>
#include <spatium/io/steps_text.hpp>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace spatium::answer {

enum class Route { AtInfinity, Evaluate, Expression, Derivative, ClosedForm, Linear, Intervals, Enclosed, Bracket, Integral, Synthesis };

enum class Arithmetic { Double, Float, Rational, Decimal, Real50 };

inline const char* arithmetic_name(Arithmetic a) {
    switch (a) {
        case Arithmetic::Double: return "double";
        case Arithmetic::Float: return "float";
        case Arithmetic::Rational: return "rational";
        case Arithmetic::Decimal: return "decimal";
        case Arithmetic::Real50: return "real50";
    }
    return "?";
}

inline constexpr Arithmetic arithmetics[] = {Arithmetic::Double, Arithmetic::Float, Arithmetic::Rational,
                                             Arithmetic::Decimal, Arithmetic::Real50};

namespace detail {

// Why this build cannot compute in `a`; empty when it can.
inline std::string unavailable(Arithmetic a) {
    if (a == Arithmetic::Decimal && !SPATIUM_HAS_DECIMAL)
        return "this build has no decimal: its compiler has no decimal64 (GCC's <decimal/decimal>)";
#if !(defined(SPATIUM_HAS_BOOST_MULTIPRECISION) && SPATIUM_HAS_BOOST_MULTIPRECISION)
    if (a == Arithmetic::Real50) return "this build has no real50: it was built without Boost.Multiprecision";
#endif
    return {};
}

} // namespace detail

// How the engine computes, apart from what is asked. The defaults are a double and a synthesis
// tolerance of 1e-8.
struct Engine {
    Arithmetic arithmetic = Arithmetic::Double;
    bool strict = false;       // an arithmetic that cannot compute a value refuses, instead of a double
    double tolerance = 1e-8;   // of a synthesis
};

// A chain of steps that checks a route's answer by another way over the same solution: Vieta's
// formulas over a closed form's coefficients and roots.
struct Check {
    std::string name;
    std::vector<spatium::io::build::Step> steps;
};

struct Candidate {
    Route route = Route::Evaluate;
    std::size_t cost = 0;        // evaluations of the definition's field
    std::string refused;         // why the route gives no value; empty when it gives one
    std::vector<double> values;  // the numbers it gives, when it gives numbers
    Arithmetic arithmetic = Arithmetic::Double;   // what it computed in
    std::string text;            // its number as that arithmetic holds it, when that is not a double (`21/4`)
    std::vector<spatium::io::build::Step> steps;   // what the route did, as it ran
    std::vector<Check> checks;
};

// An inequality's solution on its box, by the interval method: the axis in order, each part a
// piece where the inequality is proved to hold or to fail, or a point between two -- a root of
// its left side less its right, proved alone in its interval, or a break where it has no value --
// or a piece nothing was proved of.
struct Axis {
    enum class Kind { Holds, Fails, Root, Break, Undecided };
    struct Part {
        Kind kind = Kind::Undecided;
        double lo = 0, hi = 0;
        std::string text;   // a root's value: a fraction where it is one, else its proved digits
    };
    std::string relation;
    std::vector<Part> parts;   // adjacent, from the box's low end to its high end
    bool line = false;         // the box is every double: its ends are infinite
};

// The set where an inequality holds, as LaTeX: its intervals joined, an end in a bracket where it
// belongs to the set -- a root of a relation that holds at equality, or a bracket's end -- and in a
// parenthesis where it does not; the empty set, or every real number.
inline std::string set_latex(const Axis& a, const std::string& name) {
    namespace bd = spatium::io::build;
    const auto number = [](double v) { return bd::detail::number_text(v, bd::detail::TextStyle::Latex); };
    const auto value = [](const std::string& t) {
        const auto slash = t.find('/');
        if (slash == std::string::npos) return t;
        const bool minus = t.starts_with('-');
        return std::string(minus ? "-" : "") + "\\frac{" + t.substr(minus, slash - minus) + "}{" + t.substr(slash + 1) + "}";
    };
    const bool at_equality = a.relation.size() == 2;
    using K = Axis::Kind;
    std::vector<std::string> spans;
    std::string lo;
    bool open = false, whole = true;
    for (std::size_t i = 0; i < a.parts.size(); ++i) {
        const auto& p = a.parts[i];
        const K prev = i ? a.parts[i - 1].kind : K::Fails, next = i + 1 < a.parts.size() ? a.parts[i + 1].kind : K::Fails;
        if (p.kind == K::Holds) {
            if (open) continue;
            open = true;
            if (i == 0) lo = a.line ? "(-\\infty" : "[" + number(p.lo);
            else if (prev == K::Root) lo = (at_equality ? "[" : "(") + value(a.parts[i - 1].text);
            else if (prev == K::Break && !a.parts[i - 1].text.empty()) lo = "(" + value(a.parts[i - 1].text);
            else lo = "(" + number(p.lo);
            continue;
        }
        whole = false;
        if (p.kind == K::Root && at_equality) {
            if (open && next == K::Holds) continue;   // a root inside the set, as x^2 >= 0 has at 0
            if (open) spans.push_back(lo + ", " + value(p.text) + "]"), open = false;
            else if (next != K::Holds) spans.push_back("\\{" + value(p.text) + "\\}");
            continue;
        }
        if (open) spans.push_back(lo + ", " + (!p.text.empty() ? value(p.text) : number(p.lo)) + ")"), open = false;
    }
    if (open) spans.push_back(lo + ", " + (a.line ? std::string("\\infty)") : number(a.parts.back().hi) + "]"));
    const std::string x = bd::detail::name_latex(name) + " \\in ";
    if (spans.empty()) return x + "\\varnothing";
    if (whole && a.line) return x + "\\mathbb{R}";
    std::string out;
    for (const auto& sp : spans) out += (out.empty() ? "" : " \\cup ") + sp;
    return x + out;
}

// A system's roots on every double of each unknown, each proved alone in its box (roots_on in N
// unknowns), and the boxes where one may still be.
struct SystemRoots {
    std::vector<std::string> unknowns;
    std::vector<std::vector<spatium::Interval<double>>> roots, undecided;
    std::vector<std::vector<spatium::Interval<double>>> at_roots;   // each equation on each root's box: holds 0
    std::size_t pieces = 0, excluded = 0;                            // the space was cut into; proved to hold none
};

struct Evaluation {
    std::string name;
    std::vector<std::string> free;       // the parameter the value is a function of, when there is one
    std::vector<Candidate> candidates;   // in the order of preference
    std::optional<std::size_t> taken;    // the candidate whose value this is; none when every one refused
    std::string refused;                 // why there is no value, when there is none
    // the value, as the route taken gives it
    std::vector<double> values;          // the taken candidate's: a number, or the real roots a closed form keeps
    Expr written;                        // a value's expression as written, and as the function it is when free
    // the expression a page shows, where its shortest form (io/field_shortest.hpp) shows less than
    // `written`: a function's, a derivative's, or the exact form of a number, 7 sqrt(6)
    std::optional<Expr> shortest;
    std::optional<ClosedForm> closed;
    std::optional<LinearSystem> linear;   // a system's solution, every unknown's value
    std::optional<Axis> axis;             // an inequality's solution
    std::optional<SystemRoots> system;    // a system's roots, where it is not linear
    std::optional<double> error;         // of an integral, as estimated
    std::optional<Answer> answer;        // of a synthesis, in pieces with their certificate
    // a number's digits proved over the reals (core/interval.hpp), or why they are not
    std::optional<spatium::Proved> proved;
    std::string unproved;
    std::optional<spatium::Rational> exact;   // the value as a fraction, where every operation is one on fractions
    // a number's computation, a step an operation, as it ran (io/steps_text.hpp)
    std::vector<spatium::io::build::Step> steps;
    // an equation's roots on its bracket, or on the doubles' whole range: each proved alone in its
    // interval, and where one may still be
    std::optional<spatium::Roots<double>> enclosed;
    // a value that reads an infinity, as its limit with every infinity one x growing without bound
    std::optional<spatium::LimitResult<double>> limit;
    std::string parameter;               // the name the answer's pieces read as U
    // a derivative: the variable it is in, and the point it is taken at, where it is a number there
    std::string variable;
    std::optional<double> at;

    // the route of the value, or of the first refusal when there is none
    std::optional<Route> route() const {
        if (taken) return candidates[*taken].route;
        if (!candidates.empty()) return candidates.front().route;
        return std::nullopt;
    }
};

namespace detail {

template<class F>
std::string shortest(F x) {
    char buf[32];
    const auto [end, ec] = std::to_chars(buf, buf + sizeof buf, x);
    return ec == std::errc{} ? std::string(buf, end) : "nan";
}

// A number of arithmetic T as a double, as text, and whether it holds a number at all.
template<class T> double as_double(const T& x) {
    if constexpr (requires { x.to_double(); }) return x.to_double();
    else return static_cast<double>(x);
}
template<class T> std::string as_text(const T& x) {
    if constexpr (requires { x.text(); }) return x.text();
    else if constexpr (std::is_floating_point_v<T>) return shortest(x);   // a float's own shortest digits
    else return x.str();
}
template<class T> bool holds(const T& x) {
    if constexpr (requires { x.ok(); }) return x.ok();
    else if constexpr (requires { x.finite(); }) return x.finite();
    else return std::isfinite(as_double(x));
}

// A field with nothing left to read, computed in arithmetic T.
template<spatium::Scalar T>
Candidate value_in(const spatium::io::build::PodField<double>& pod, Arithmetic a) {
    namespace bd = spatium::io::build;
    Candidate c{Route::Evaluate, 1, {}, {}, a, {}};
    if (const auto op = bd::beyond<T>(pod)) {
        c.refused = std::format("{} has no {} in it", arithmetic_name(a), bd::op_name(*op));
        return c;
    }
    const T v = bd::interpret(bd::converted<T>(pod), T{0}, T{0});
    if (!holds(v)) {
        c.refused = a == Arithmetic::Rational ? "the fraction does not fit in 64 bits, or divides by zero"
                                              : std::format("the value is not a number in {}", arithmetic_name(a));
        return c;
    }
    c.values = {as_double(v)};
    c.text = as_text(v);
    return c;
}

inline Candidate value_in(const Expr& f, Arithmetic a) {
    namespace bd = spatium::io::build;
    Candidate c{Route::Evaluate, 1, unavailable(a), {}, a, {}};
    if (!c.refused.empty()) return c;
    const auto pod = bd::lower(f);
    if (!pod) {
        c.refused = pod.error().message;
        return c;
    }
    switch (a) {
        case Arithmetic::Double: break;
        case Arithmetic::Float: return value_in<float>(*pod, a);
        case Arithmetic::Rational: return value_in<spatium::Rational>(*pod, a);
        case Arithmetic::Decimal:
#if SPATIUM_HAS_DECIMAL
            return value_in<spatium::Decimal>(*pod, a);
#else
            break;
#endif
        case Arithmetic::Real50:
#if defined(SPATIUM_HAS_BOOST_MULTIPRECISION) && SPATIUM_HAS_BOOST_MULTIPRECISION
            return value_in<spatium::Real50>(*pod, a);
#else
            break;
#endif
    }
    return value_in<double>(*pod, a);
}

// A number with nothing left to read, enclosed: computed in an interval of B, each constant held
// as it was written, so the digits its ends share are digits of what was written, not of a double
// near it. None are proved when a comparison on the way was undecided, since the program then
// followed a guessed branch. False when the interval is empty: the value is no real number, as
// 1/0 is not; an infinity read in is one, exactly.
template<spatium::Scalar B>
bool prove(const spatium::io::build::PodField<double>& pod, Evaluation& out) {
    namespace bd = spatium::io::build;
    using I = spatium::Interval<B>;
    if (const auto op = bd::beyond<I>(pod)) {
        out.unproved = std::format("an interval has no {} in it", bd::op_name(*op));
        return true;
    }
    const spatium::Undecided undecided;
    const I v = bd::interpret(bd::converted<I>(pod), I{0}, I{0});
    const double lo = static_cast<double>(v.lo), hi = static_cast<double>(v.hi);
    if (undecided.raised()) out.unproved = "a comparison on the way was undecided, or a root or logarithm reached below zero";
    else if (v.is_empty()) return out.unproved = "the value is not a number", false;
    else if (lo == hi && std::isinf(lo)) out.proved = spatium::Proved{0, lo > 0 ? "inf" : "-inf", std::abs(lo), true};
    else out.proved = spatium::proved(v);
    return true;
}

inline bool prove(const Expr& f, Arithmetic a, Evaluation& out) {
    namespace bd = spatium::io::build;
    const auto pod = bd::lower(f);
    if (!pod) {
        out.unproved = pod.error().message;
        return true;
    }
    // where every operation is one on fractions, the value is the fraction it is, exactly: the
    // twelve a page of exercises asks for, not 12.000000000000002 and its proved digits. Not where
    // it reads pi, which is no fraction: the double nearest it read as one made 2 pi "exactly"
    // 3141592653589793/500000000000000
    const bool reads_pi = std::ranges::any_of(pod->ops, [](const auto& n) {
        return static_cast<bd::Op>(n.code) == bd::Op::Const && std::abs(n.value) == std::numbers::pi;
    });
    if (!reads_pi && !bd::beyond<spatium::Rational>(*pod)) {
        const auto r = bd::interpret(bd::converted<spatium::Rational>(*pod), spatium::Rational(0), spatium::Rational(0));
        if (r.ok()) {
            const auto decimal = r.decimal();
            out.exact = r;
            out.proved = spatium::Proved{decimal ? static_cast<int>(std::ranges::count_if(*decimal, [](char c) { return c >= '0' && c <= '9'; })) : 0,
                                         decimal ? *decimal : r.text(), std::abs(r.to_double()), true};
            return true;
        }
    }
#if defined(SPATIUM_HAS_BOOST_MULTIPRECISION) && SPATIUM_HAS_BOOST_MULTIPRECISION
    if (a == Arithmetic::Real50) return prove<spatium::Real50>(*pod, out);
#endif
    (void)a;
    return prove<double>(*pod, out);
}

// what a page shows of an expression: its shortest form, where it is not the expression itself
inline std::optional<Expr> shortest_of(const Expr& e) {
    Expr s = spatium::io::build::shortest(e);
    if (spatium::io::build::to_text(s).value_or("") == spatium::io::build::to_text(e).value_or("")) return std::nullopt;
    return s;
}

inline std::string plain(double v) {
    namespace bd = spatium::io::build;
    return bd::detail::number_text(v, bd::detail::TextStyle::Plain);
}

// The steps of a search for every root: where it searched, and how the space was cut.
inline void searched(Candidate& c, std::string latex, std::size_t pieces, std::size_t excluded, std::size_t undecided) {
    c.steps.emplace_back("searched", std::move(latex));
    c.steps.emplace_back("pieces cut where undecided: with no root + with one + undecided",
                         std::format("{} = {} + {} + {}", pieces, excluded, pieces - excluded - undecided, undecided));
}

// The roots of a field of one unknown `x` on a box, each proved: computed in an interval and in a
// dual number of intervals, which the lowered program runs in unchanged.
inline Candidate enclosed_roots(const Expr& f, spatium::Interval<double> box, const std::string& x, Evaluation& out) {
    namespace bd = spatium::io::build;
    using I = spatium::Interval<double>;
    using D = spatium::Dual<I>;
    Candidate c{Route::Enclosed, 0, {}, {}, Arithmetic::Double, {}};
    const auto pod = bd::lower(f);
    if (!pod) {
        c.refused = pod.error().message;
        return c;
    }
    if (const auto op = bd::beyond<I>(*pod)) {
        c.refused = std::format("an interval has no {} in it", bd::op_name(*op));
        return c;
    }
    const auto in_interval = bd::converted<I>(*pod);
    const auto in_dual = bd::converted<D>(*pod);
    const auto field = [&]<class X>(const X& x) {
        if constexpr (std::is_same_v<X, I>) return bd::interpret(in_interval, x, I{0});
        else return bd::interpret(in_dual, x, D{I{0}});
    };
    auto found = spatium::roots_on(field, box);
    c.cost = found.evaluations;
    for (const I& r : found.roots) c.values.push_back(r.mid());
    const auto range = [](const I& x) { return std::format("[{}, {}]", shortest(x.lo), shortest(x.hi)); };
    if (found.roots.empty()) {
        if (found.complete()) c.refused = std::format("no real root on {}: every piece of it is proved to hold none", range(box));
        else {
            c.refused = std::format("no root proved on {}; one may lie in {}", range(box), range(found.undecided.front()));
            if (found.undecided.size() > 1) c.refused += std::format(" or {} more pieces", found.undecided.size() - 1);
        }
    }
    const std::string name = spatium::io::build::detail::name_latex(x);
    const bool line = box.hi == std::numeric_limits<double>::max() && box.lo == -std::numeric_limits<double>::max();
    searched(c, line ? name + " \\in \\mathbb{R}" : std::format("{} \\in \\left[{}, {}\\right]", name, plain(box.lo), plain(box.hi)),
             found.pieces, found.excluded, found.undecided.size());
    for (const I& r : found.roots) {
        const auto p = spatium::proved(r);
        c.steps.emplace_back("Newton's step lands inside: one root", std::format("{} \\approx {}", name, p.digits > 0 ? p.text : std::string("0")));
    }
    out.enclosed = std::move(found);
    return c;
}

// The roots of a system of N equations in N unknowns on every double of each, each proved (N from
// 2 to 4): each equation lowered, run in intervals and in dual numbers, its unknowns put in its slots.
template<std::size_t N = 2>
Candidate enclosed_system(const std::vector<Expr>& equations, const std::vector<std::string>& unknowns, Evaluation& out) {
    namespace bd = spatium::io::build;
    using I = spatium::Interval<double>;
    using D = spatium::Dual<I>;
    Candidate c{Route::Enclosed, 0, {}, {}, Arithmetic::Double, {}};
    if constexpr (N > 4) {
        c.refused = "a system of more than four unknowns";
        return c;
    } else {
        if (unknowns.size() != N) return enclosed_system<N + 1>(equations, unknowns, out);
        std::vector<bd::PodField<I>> in_interval;
        std::vector<bd::PodField<D>> in_dual;
        std::vector<std::array<int, N>> slot;   // where each unknown is in each equation's slots, or -1
        for (const auto& e : equations) {
            const auto pod = bd::lower(e);
            if (!pod) return c.refused = pod.error().message, c;
            if (const auto op = bd::beyond<I>(*pod)) return c.refused = std::format("an interval has no {} in it", bd::op_name(*op)), c;
            std::array<int, N> at;
            for (std::size_t j = 0; j < N; ++j) {
                const auto s = std::ranges::find(pod->symbols, unknowns[j]);
                at[j] = s == pod->symbols.end() ? -1 : static_cast<int>(s - pod->symbols.begin());
            }
            in_interval.push_back(bd::converted<I>(*pod)), in_dual.push_back(bd::converted<D>(*pod)), slot.push_back(at);
        }
        const auto f = [&]<class X>(const std::array<X, N>& v) {
            std::array<X, N> y;
            for (std::size_t k = 0; k < N; ++k) {
                const auto& program = [&]() -> const auto& { if constexpr (std::is_same_v<X, I>) return in_interval[k]; else return in_dual[k]; }();
                std::vector<X> slots(program.symbols.size(), X(0));
                for (std::size_t j = 0; j < N; ++j)
                    if (slot[k][j] >= 0) slots[static_cast<std::size_t>(slot[k][j])] = v[j];
                y[k] = bd::interpret(program, X(0), X(0), slots);
            }
            return y;
        };
        const double big = std::numeric_limits<double>::max();
        spatium::Box<double, N> box;
        box.fill(I{-big, big});
        const spatium::RootsN<double, N> found = spatium::roots_on(f, box, spatium::Budget{.evaluations = 65536});
        c.cost = found.evaluations;
        SystemRoots s{unknowns, {}, {}, {}, found.pieces, found.excluded};
        for (const auto& r : found.roots) {
            s.roots.emplace_back(r.begin(), r.end());
            const auto y = f(r);
            s.at_roots.emplace_back(y.begin(), y.end());
        }
        for (const auto& r : found.undecided) s.undecided.emplace_back(r.begin(), r.end());
        if (s.roots.empty())
            c.refused = found.complete() ? "no real root: every piece of the space is proved to hold none"
                                         : std::format("no root proved; one may lie in {} pieces", found.undecided.size());
        for (const auto& r : s.roots)
            for (const I& v : r) c.values.push_back(v.mid());
        const auto name = [](const std::string& u) { return bd::detail::name_latex(u); };
        std::string names;
        for (const auto& u : unknowns) names += (names.empty() ? "" : ", ") + name(u);
        searched(c, names + " \\in \\mathbb{R}", s.pieces, s.excluded, s.undecided.size());
        for (const auto& r : s.roots) {
            std::string at;
            for (std::size_t j = 0; j < r.size(); ++j) at += (j ? ",\\ " : "") + name(unknowns[j]) + " \\approx " + spatium::proved(r[j]).text;
            c.steps.emplace_back("Krawczyk's step lands inside: one root", std::move(at));
        }
        out.system = std::move(s);
        return c;
    }
}

// The fraction with the smallest denominator in [a, b], by continued fractions: -2 in
// [-2.0000000005, -1.9999999999999998]; none past 64 bits or 40 steps.
inline std::optional<spatium::Rational> simplest_in(double a, double b, int depth = 0) {
    using spatium::Rational;
    if (!(a <= b) || depth > 40 || !std::isfinite(a) || !std::isfinite(b)) return std::nullopt;
    if (a <= 0 && 0 <= b) return Rational(0);
    if (b < 0) {
        const auto r = simplest_in(-b, -a, depth);
        return r ? std::optional(Rational(0) - *r) : std::nullopt;
    }
    if (const double c = std::ceil(a); c <= b) return c < 0x1p62 ? std::optional(Rational(static_cast<long long>(c), 1)) : std::nullopt;
    const double f = std::floor(a);
    const auto r = simplest_in(1 / (b - f), 1 / (a - f), depth + 1);
    if (!r || *r == Rational(0)) return std::nullopt;
    return Rational(static_cast<long long>(f), 1) + Rational(1) / *r;
}

// An inequality f rel 0 in `x` on a box, by the interval method (Axis). `written` is f as written,
// its constants not folded, where a root that is a short fraction is checked to be one exactly.
inline Candidate interval_method(const Expr& f, const Expr& written, const std::string& rel, spatium::Interval<double> box,
                                 const std::string& x, Evaluation& out) {
    namespace bd = spatium::io::build;
    using I = spatium::Interval<double>;
    using D = spatium::Dual<I>;
    using spatium::Truth;
    Candidate c{Route::Intervals, 0, {}, {}, Arithmetic::Double, {}};
    const auto pod = bd::lower(f);
    if (!pod) return c.refused = pod.error().message, c;
    if (const auto op = bd::beyond<I>(*pod)) return c.refused = std::format("an interval has no {} in it", bd::op_name(*op)), c;
    const auto in_interval = bd::converted<I>(*pod);
    const auto in_dual = bd::converted<D>(*pod);
    const auto field = [&]<class X>(const X& x) {
        if constexpr (std::is_same_v<X, I>) return bd::interpret(in_interval, x, I{0});
        else return bd::interpret(in_dual, x, D{I{0}});
    };
    // f rel 0 on a piece. Where f has no value on it the piece is cut, as roots_on cuts it, and
    // left undecided once it is one magnitude with no value at either end: a break, found below
    const auto one_scale = [](const I& x) { return (0 < x.lo && x.hi <= 16 * x.lo) || (x.hi < 0 && 16 * x.hi <= x.lo); };
    const auto holds = [&](const spatium::Box<double, 1>& b) {
        const I y = field(b[0]), zero{0};
        if (y.is_empty())
            return spatium::Decided{Truth::Undecided, one_scale(b[0]) && field(I{b[0].lo}).is_empty() && field(I{b[0].hi}).is_empty()};
        if (rel == ">") return spatium::Decided{less(zero, y)};
        if (rel == ">=") return spatium::Decided{less_equal(zero, y)};
        if (rel == "<") return spatium::Decided{less(y, zero)};
        return spatium::Decided{less_equal(y, zero)};
    };
    auto proof = spatium::prove_on(spatium::Box<double, 1>{box}, holds);
    c.cost = proof.evaluations;
    std::ranges::sort(proof.pieces, {}, [](const auto& p) { return p.first[0].lo; });
    // a root or a break at a short fraction, checked in fractions: f is zero there, or has no
    // value -- 2, not 2.000000000000000; -2, not -1.9999999999999998
    const auto pod_written = bd::lower(written);
    const auto exact = [&](const I& r, bool root) -> std::string {
        if (!pod_written || bd::beyond<spatium::Rational>(*pod_written)) return {};
        const auto found = simplest_in(r.lo, r.hi);
        if (!found || !found->ok() || found->den() > 1000000) return {};
        const spatium::Rational q = *found;
        const auto v = bd::interpret(bd::converted<spatium::Rational>(*pod_written), q, spatium::Rational(0));
        return (root ? v.ok() && v == spatium::Rational(0) : !v.ok()) ? q.text() : std::string();
    };
    Axis axis{rel, {}, box.lo == -std::numeric_limits<double>::max() && box.hi == std::numeric_limits<double>::max()};
    const auto add = [&](Axis::Kind k, double lo, double hi, std::string text = {}) {
        if (!axis.parts.empty() && axis.parts.back().kind == k && k != Axis::Kind::Root && k != Axis::Kind::Break)
            axis.parts.back().hi = hi;
        else axis.parts.push_back({k, lo, hi, std::move(text)});
    };
    // the pieces nothing was proved of, run together, are where it changes: a root, a break, or undecided
    std::vector<I> runs;
    for (const auto& [b, t] : proof.pieces) {
        if (t != Truth::Undecided) {
            add(t == Truth::Yes ? Axis::Kind::Holds : Axis::Kind::Fails, b[0].lo, b[0].hi);
            continue;
        }
        if (!axis.parts.empty() && axis.parts.back().kind == Axis::Kind::Undecided && axis.parts.back().hi == b[0].lo) {
            axis.parts.back().hi = b[0].hi;
            continue;
        }
        axis.parts.push_back({Axis::Kind::Undecided, b[0].lo, b[0].hi, {}});
    }
    // past the doubles at the ends of the line the value overflows (inf / inf): an end piece nothing
    // was proved of takes the verdict of the piece inside it, until the tails are read at infinity
    if (axis.line && axis.parts.size() > 1) {
        const double past = 0x1p1020;
        if (axis.parts.front().kind == Axis::Kind::Undecided && axis.parts.front().hi <= -past)
            axis.parts[1].lo = axis.parts.front().lo, axis.parts.erase(axis.parts.begin());
        if (axis.parts.size() > 1 && axis.parts.back().kind == Axis::Kind::Undecided && axis.parts.back().lo >= past)
            axis.parts[axis.parts.size() - 2].hi = axis.parts.back().hi, axis.parts.pop_back();
    }
    for (auto& p : axis.parts) {
        if (p.kind != Axis::Kind::Undecided) continue;
        // a break at a short fraction, where f in fractions has no value: 1/(x + 2) at -2
        if (auto at = exact(I{p.lo, p.hi}, false); !at.empty()) {
            p.kind = Axis::Kind::Break, p.text = std::move(at);
            continue;
        }
        const auto found = spatium::roots_on(field, I{p.lo, p.hi}, spatium::Budget{.evaluations = 256});
        c.cost += found.evaluations;
        if (found.complete() && found.roots.size() == 1) {
            const I& r = found.roots.front();
            const auto proved = spatium::proved(r);
            const std::string fraction = exact(r, true);
            p = {Axis::Kind::Root, r.lo, r.hi, !fraction.empty() ? fraction : proved.digits > 0 ? proved.text : std::format("{}", r.mid())};
            c.values.push_back(r.mid());
        } else if (found.complete() && found.roots.empty()) {
            p.kind = Axis::Kind::Break, p.text = exact(I{p.lo, p.hi}, false);
        }
    }
    // its steps: each end of the set, a root or a break, then the set
    const std::string name = bd::detail::name_latex(x);
    for (const auto& part : axis.parts)
        if (part.kind == Axis::Kind::Root)
            c.steps.emplace_back("an end: a root, proved alone by Newton's step", std::format("{} \\approx {}", name, part.text));
        else if (part.kind == Axis::Kind::Break)
            c.steps.emplace_back("a break: no value", std::format("{} \\in \\left[{}, {}\\right]", name, plain(part.lo), plain(part.hi)));
    c.steps.emplace_back("where it holds", set_latex(axis, x));
    out.axis = std::move(axis);
    return c;
}

// A derivative its rules write (io/field_derivative.hpp) against the one a dual number computes
// from the function by its own rules (algebra/dual.hpp), both in intervals, at each point: where
// both have a value they must overlap. A point where they do not is a fault in one of the two.
struct DualAgreement {
    std::size_t compared = 0;          // points where both have a value
    double difference = 0;             // the largest gap between the two midpoints there
    std::optional<double> disagrees;   // the first point where they do not overlap
    std::string refused;               // why nothing could be compared
};

inline DualAgreement agree_with_dual(const Expr& f, const Expr& df, const std::string& x, const std::vector<double>& points) {
    namespace bd = spatium::io::build;
    using I = spatium::Interval<double>;
    using D = spatium::Dual<I>;
    DualAgreement out;
    const auto pf = bd::lower(f), pd = bd::lower(df);
    if (!pf || !pd) return out.refused = (!pf ? pf : pd).error().message, out;
    if (const auto op = bd::beyond<D>(*pf)) return out.refused = std::format("a dual number has no {} in it", bd::op_name(*op)), out;
    if (const auto op = bd::beyond<I>(*pd)) return out.refused = std::format("an interval has no {} in it", bd::op_name(*op)), out;
    const auto in_dual = bd::converted<D>(*pf);
    const auto in_interval = bd::converted<I>(*pd);
    for (const double p : points) {
        const spatium::Undecided undecided;
        std::vector<D> fs(pf->symbols.size(), D{I{0}});
        for (std::size_t k = 0; k < fs.size(); ++k)
            if (pf->symbols[k] == x) fs[k] = D{I{p}, I{1}};
        std::vector<I> ds(pd->symbols.size(), I{0});
        for (std::size_t k = 0; k < ds.size(); ++k)
            if (pd->symbols[k] == x) ds[k] = I{p};
        const I by_dual = bd::interpret(in_dual, D{I{0}}, D{I{0}}, fs).deriv;
        const I by_rules = bd::interpret(in_interval, I{0}, I{0}, ds);
        const auto valued = [](const I& v) { return !v.is_empty() && std::isfinite(v.lo) && std::isfinite(v.hi); };
        if (undecided.raised() || !valued(by_dual) || !valued(by_rules)) continue;
        ++out.compared;
        out.difference = std::max(out.difference, std::abs(by_dual.mid() - by_rules.mid()));
        if ((by_rules.hi < by_dual.lo || by_dual.hi < by_rules.lo) && !out.disagrees) out.disagrees = p;
    }
    if (!out.compared) out.refused = "neither has a value at any point compared";
    return out;
}

} // namespace detail

// ── The routes, as a table ────────────────────────────────────────
//
// A route is a row: its name, how it finds its value in the words a page shows, what it applies
// to, and how it computes -- the candidates it adds to an evaluation, each with the steps it ran.
// The rows are in the order of preference `take` reads: a closed form before every root proved,
// and that before a root on its bracket.

namespace detail {

// What a route is given: the document, the definition asked for, unfolded, and how to compute.
struct RouteInput {
    const spatium::io::build::Document& doc;
    const spatium::io::build::Statement& d;
    const detail::Unfolded& u;
    const Chosen& chosen;
    const Engine& engine;
};

struct RouteRow {
    Route route;
    const char* name;
    const char* caption;
    bool (*applies)(const RouteInput&, const Evaluation&);   // the evaluation so far: a route may follow another's refusal
    void (*run)(const RouteInput&, Evaluation&);
};

using Kind = spatium::io::build::StatementKind;

inline bool is_value(const RouteInput& in) { return in.d.kind == Kind::Value || in.d.kind == Kind::Parameter; }
// a root or an integral reading no parameter: a number, found by a solver
inline bool is_number(const RouteInput& in) { return !is_value(in) && in.u.free.empty(); }
inline bool is_equation(const RouteInput& in) { return is_number(in) && in.d.kind == Kind::Root && in.d.relation == "="; }
inline bool is_root(const RouteInput& in) { return is_equation(in) && in.d.also.empty(); }
inline double number(const Expr& e) { return e(0.0, 0.0); }
inline spatium::Interval<double> box_of(const RouteInput& in) {
    const double big = std::numeric_limits<double>::max();
    return spatium::io::build::has_bracket(in.d) ? spatium::Interval<double>{number(in.u.lo), number(in.u.hi)} : spatium::Interval<double>{-big, big};
}

// Whether a field reads an infinity: inf or -inf among the constants of its program.
inline bool reads_infinity(const Expr& f) {
    namespace bd = spatium::io::build;
    const auto pod = bd::lower(f);
    return pod && std::ranges::any_of(pod->ops, [](const auto& n) { return static_cast<bd::Op>(n.code) == bd::Op::Const && std::isinf(n.value); });
}

// The leading term of a limit at infinity as LaTeX, f ~ c x^k: k = -exponent / ramification.
inline std::string leading_latex(const spatium::LimitResult<double>& l) {
    namespace bd = spatium::io::build;
    if (l.coefficient == 0) return "0";
    const std::string c = bd::detail::number_text(l.coefficient, bd::detail::TextStyle::Latex);
    if (l.exponent == 0) return c;
    const int k = -l.exponent, r = l.ramification;
    const std::string power = r == 1 || k % r == 0 ? std::to_string(k / r) : std::format("{}/{}", k, r);
    const std::string x = power == "1" ? "x" : std::format("x^{{{}}}", power);
    return (l.coefficient == 1 ? "" : l.coefficient == -1 ? "-" : c + "\\,") + x;
}

// a value that reads an infinity: its limit as x grows without bound, every inf being that one
// x, computed as a Laurent series in 1/x (algebra/series.hpp) -- arithmetic on the chart at
// infinity, where no infinity enters an operation
inline void at_infinity(const RouteInput& in, Evaluation& out) {
    namespace bd = spatium::io::build;
    using S = spatium::Series<double, 8>;
    out.written = in.d.value.field;
    Candidate c{Route::AtInfinity, 0, {}, {}, Arithmetic::Double, {}};
    const auto pod = bd::lower(in.u.f);
    if (!pod) return c.refused = pod.error().message, out.candidates.push_back(std::move(c));
    if (const auto op = bd::beyond<S>(*pod))
        return c.refused = std::format("a series has no {} in it", bd::op_name(*op)), out.candidates.push_back(std::move(c));
    const auto f = [&](const S& x) {
        auto program = bd::converted<S>(*pod);
        for (std::size_t i = 0; i < program.ops.size(); ++i)
            if (const double v = pod->ops[i].value; static_cast<bd::Op>(pod->ops[i].code) == bd::Op::Const && std::isinf(v))
                program.ops[i].value = v > 0 ? x : -x;
        ++c.cost;
        return bd::interpret(program, S{0}, S{0});
    };
    const auto limit = spatium::limit(f, spatium::AtInfinity<double>{});
    if (!limit.determined()) {
        c.refused = "with every inf one x growing without bound, the value is no Laurent series in 1/x: an exponential, a "
                    "logarithm or a non-integer power of the infinity, or a division by zero";
        return out.candidates.push_back(std::move(c));
    }
    c.values = {limit.value};
    const std::string written = *bd::to_latex(out.written);
    c.steps.emplace_back("every \\infty one x, as x grows without bound: the leading term", written + " \\sim " + leading_latex(limit));
    c.steps.emplace_back("its limit", std::format("{} \\to {}", written, bd::detail::number_text(limit.value, bd::detail::TextStyle::Latex)));
    out.limit = limit;
    out.candidates.push_back(std::move(c));
}

// a value with nothing left to read: in the engine's arithmetic, and where it cannot, in a double
// unless it is strict; proved, and its computation a step an operation
inline void evaluate_value(const RouteInput& in, Evaluation& out) {
    out.written = in.d.kind == Kind::Value ? in.d.value.field : Expr::symbol(in.d.name);
    if (in.engine.arithmetic == Arithmetic::Double) {
        out.candidates.push_back({Route::Evaluate, 1, {}, {number(in.u.f)}, Arithmetic::Double, {}});
    } else {
        out.candidates.push_back(value_in(in.u.f, in.engine.arithmetic));
        if (!out.candidates.back().refused.empty() && !in.engine.strict)
            out.candidates.push_back({Route::Evaluate, 1, {}, {number(in.u.f)}, Arithmetic::Double, {}});
    }
    const bool real = prove(in.u.f, in.engine.arithmetic, out);
    // the exact form of a number that is no fraction: sqrt(24) + sqrt(150) is 7 sqrt(6)
    if (real && !out.exact) out.shortest = shortest_of(in.u.f);
    if (auto steps = spatium::io::build::evaluation_steps(in.u.f)) out.steps = std::move(*steps);
    // no real number, though a double gives inf or NaN: refused, naming the operation that has none
    if (!real) {
        std::string why = "it has no value";
        if (!out.steps.empty() && out.steps.back().lead == "No value")
            why += std::format(": the {} {} has none", out.steps.back().label, out.steps.back().note);
        for (auto& c : out.candidates)
            if (c.route == Route::Evaluate && c.refused.empty()) c.refused = why, c.values.clear();
    }
    // one operation is the value itself, not a chain
    if (out.steps.size() > 1)
        for (auto& c : out.candidates)
            if (c.route == Route::Evaluate) c.steps = out.steps;
}

inline void expression(const RouteInput& in, Evaluation& out) {
    out.written = in.u.f;   // the function, exact, in the parameter's name
    out.shortest = shortest_of(out.written);
    out.candidates.push_back({Route::Expression, 0, {}, {}, Arithmetic::Double, {}});
}

// a derivative: the function its rules write, checked against a dual number at sixteen points of
// its variable's interval (of [-2, 2] for a variable of its own); at a point, its value there,
// proved in intervals as a value is, and checked at that point
inline void derivative(const RouteInput& in, Evaluation& out) {
    namespace bd = spatium::io::build;
    Candidate c{Route::Derivative, 0, {}, {}, Arithmetic::Double, {}};
    const std::string& x = in.d.variable;
    out.variable = x;
    const bd::Statement* parameter = in.doc.find(x);
    if (parameter && parameter->kind != Kind::Parameter) parameter = nullptr;
    if (bd::has_bracket(in.d)) {
        if (const auto open = in.u.lo.free_symbols(); !open.empty())
            return c.refused = std::format("the point reads {}; it is a number", open.front()), out.candidates.push_back(std::move(c));
        out.at = number(in.u.lo);
    } else if (const auto chosen = std::ranges::find(in.chosen, x, &Chosen::value_type::first); chosen != in.chosen.end()) {
        out.at = chosen->second;
    }
    const Expr& g = in.u.f;
    for (const auto& n : g.free_symbols())
        if (n != x) return c.refused = std::format("the derivative reads {} as well as {}; give it a value", n, x), out.candidates.push_back(std::move(c));
    std::vector<double> points;
    if (out.at) {
        points = {*out.at};
    } else {
        const double lo = parameter ? number(parameter->lo.field) : -2.0, hi = parameter ? number(parameter->hi.field) : 2.0;
        for (int k = 0; k < 16; ++k) points.push_back(lo + (hi - lo) * (k + 0.5) / 16);
    }
    const DualAgreement agreed = agree_with_dual(in.u.written, g, x, points);
    c.cost = 2 * points.size();
    const std::string name = bd::detail::name_latex(x);
    const std::string written = std::format("\\frac{{d}}{{d{}}}\\left({}\\right)", name, *bd::to_latex(in.d.value.field));
    c.steps.emplace_back("each operation by its rule, the chain rule between them", written + " = " + *bd::to_latex(g));
    if (agreed.disagrees) {
        c.refused = std::format("the rules and a dual number disagree on the derivative at {} = {}", x, *agreed.disagrees);
        return out.candidates.push_back(std::move(c));
    }
    if (agreed.compared) {
        const std::string where = out.at ? std::format("at {} = {}", x, plain(*out.at))
                                         : std::format("at {} points of [{}, {}]", agreed.compared, plain(points.front()), plain(points.back()));
        c.checks.push_back({"a dual number", {spatium::io::build::Step(
            std::format("the derivative a dual number computes from the function, in intervals, {}: they overlap", where),
            std::format("\\max \\left|\\Delta\\right| = {}", bd::detail::number_text(agreed.difference, bd::detail::TextStyle::Latex)))}});
    }
    out.written = g;
    if ((out.shortest = shortest_of(g)))
        c.steps.emplace_back("like terms collected, over one denominator", "= " + *bd::to_latex(*out.shortest));
    if (!out.at) return out.candidates.push_back(std::move(c));
    const Expr there = g.bind(x, Expr(*out.at));
    c.values = {number(there)};
    if (!prove(there, in.engine.arithmetic, out)) {
        c.refused = std::format("the derivative has no value at {} = {}", x, plain(*out.at));
        c.values.clear();
    } else {
        c.steps.emplace_back(std::format("at {} = {}", x, plain(*out.at)),
                             std::format("{} = {}", *bd::to_latex(g), bd::detail::number_text(c.values.front(), bd::detail::TextStyle::Latex)));
    }
    out.candidates.push_back(std::move(c));
}

// a root, as fractions when the engine asks for them: the roots as they are, or why they are not
inline void closed_form(const RouteInput& in, Evaluation& out) {
    const bool exact = in.engine.arithmetic == Arithmetic::Rational;
    const bool bracket = spatium::io::build::has_bracket(in.d);
    const auto box = box_of(in);
    Candidate c{Route::ClosedForm, 0, {}, {}, Arithmetic::Double, {}};
    // the roots as they are are always sought: a page shows them beside their decimals
    if (auto solved = solve_closed_form(in.u.f, in.d.name, &in.u.written); !solved) {
        c.refused = solved.error().message;
    } else {
        c.cost = solved->evaluations;
        for (const auto& z : solved->solution.roots)
            if (z.im == 0 && (!bracket || (box.lo <= z.re && z.re <= box.hi))) c.values.push_back(z.re);
        if (bracket && c.values.empty()) c.refused = std::format("no root of the closed form is on [{}, {}]", in.d.lo.text, in.d.hi.text);
        c.steps = solved->steps;
        // Vieta's formulas: the equation and its coefficients, then the sum and the product of the roots
        if (solved->check && c.steps.size() >= 2)
            c.checks.push_back({"Vieta's formulas", {c.steps[0], c.steps[1], spatium::io::build::to_step(*solved->check, solved->unknown)}});
        if (exact && solved->exact) {
            c.arithmetic = Arithmetic::Rational, c.text = solved->exact->text;
            // a fraction is the root itself, not the double's rounding of it
            if (!solved->exact->fractions.empty()) {
                c.values.clear();
                for (const auto& r : solved->exact->fractions)
                    if (!bracket || (box.lo <= r.to_double() && r.to_double() <= box.hi)) c.values.push_back(r.to_double());
            }
        } else if (exact) {
            Candidate refused{Route::ClosedForm, c.cost, solved->inexact, {}, Arithmetic::Rational, {}};
            if (in.engine.strict) c = std::move(refused);
            else out.candidates.push_back(std::move(refused));
        }
        out.closed = std::move(*solved);
    }
    out.candidates.push_back(std::move(c));
}

// a system: every unknown together, by elimination where its equations are linear
inline std::vector<std::string> unknowns_of(const RouteInput& in) {
    std::vector<std::string> unknowns{in.d.name};
    for (const auto& e : in.d.also) unknowns.push_back(e.name);
    return unknowns;
}

inline void linear(const RouteInput& in, Evaluation& out) {
    Candidate c{Route::Linear, 0, {}, {}, Arithmetic::Double, {}};
    if (auto solved = solve_linear_system(in.u.system, unknowns_of(in), &in.u.system_written); !solved) {
        c.refused = solved.error().message;
    } else {
        c.cost = solved->evaluations, c.values = solved->values, c.steps = solved->steps;
        if (!solved->exact.empty()) c.arithmetic = Arithmetic::Rational;
        out.linear = std::move(*solved);
    }
    out.candidates.push_back(std::move(c));
}

inline bool not_linear(const Evaluation& out) {
    return !out.candidates.empty() && out.candidates.back().route == Route::Linear &&
           out.candidates.back().refused.find("not linear") != std::string::npos;
}

// every root on the bracket, or on the doubles' whole range without one, each proved; in intervals
// of doubles, so not where a strict engine asks for another arithmetic alone. A system, where it
// is not linear, on every double of each unknown.
inline void enclosed(const RouteInput& in, Evaluation& out) {
    if (!in.d.also.empty()) return out.candidates.push_back(enclosed_system(in.u.system, unknowns_of(in), out));
    out.candidates.push_back(enclosed_roots(in.u.f.bind(in.u.bound, Expr::u()), box_of(in), in.d.name, out));
}

inline void bracket(const RouteInput& in, Evaluation& out) {
    std::size_t evaluations = 0;
    const Expr of_u = in.u.f.bind(in.u.bound, Expr::u());
    const auto box = box_of(in);
    const auto root = spatium::find_root([&](double x) { return ++evaluations, of_u(x, 0.0); }, box.lo, box.hi);
    Candidate c{Route::Bracket, evaluations, root ? std::string() : root.error().message, {}, Arithmetic::Double, {}};
    if (root) {
        c.values = {*root};
        c.steps.emplace_back("bisection to the last bit", std::format("{} = {}", spatium::io::build::detail::name_latex(in.d.name), plain(*root)));
    }
    out.candidates.push_back(std::move(c));
}

// an inequality: the set where it holds, by the interval method
inline void intervals(const RouteInput& in, Evaluation& out) {
    out.candidates.push_back(interval_method(in.u.f.bind(in.d.name, Expr::u()), in.u.written.bind(in.d.name, Expr::u()), in.d.relation,
                                             box_of(in), in.d.name, out));
}

inline void integral(const RouteInput& in, Evaluation& out) {
    std::size_t evaluations = 0;
    const Expr of_u = in.u.f.bind(in.u.bound, Expr::u());
    const auto r = spatium::integrate([&](double x) { return ++evaluations, of_u(x, 0.0); },
                                      spatium::Finite<double>{number(in.u.lo), number(in.u.hi)});
    Candidate c{Route::Integral, evaluations, {}, {}, Arithmetic::Double, {}};
    if (r.trusted()) c.values = {r.value}, out.error = r.error_estimate;
    else c.refused = std::format("the integral of {} did not converge", in.d.name);
    out.candidates.push_back(std::move(c));
}

// a root or an integral reading one parameter: a function of it, synthesised and certified
inline void synthesis(const RouteInput& in, Evaluation& out) {
    Candidate c{Route::Synthesis, 0, {}, {}, Arithmetic::Double, {}};
    if (auto q = query_of(in.doc, in.d.name, in.chosen); !q) {
        c.refused = q.error().message;
    } else {
        q->tolerance = in.engine.tolerance;
        Answer a = answer(*q);
        c.cost = a.evaluations, c.refused = a.refused;
        out.parameter = q->parameter;
        out.answer = std::move(a);
    }
    out.candidates.push_back(std::move(c));
}

} // namespace detail

inline constexpr detail::RouteRow route_table[] = {
    {Route::AtInfinity, "at infinity",
     "As a limit. Every inf in the value is one x growing without bound, and the value is computed as a Laurent "
     "series in 1/x -- arithmetic on the chart at infinity, where no infinity enters an operation -- whose leading "
     "term is the limit: a number, or an infinity with its order. So inf times 0 is 0 and inf - inf is 0, the "
     "limit of x - x. An exponential or a logarithm of the infinity is no such series, and is refused here.",
     // in a series of doubles, so not where a strict engine asks for another arithmetic alone
     [](const detail::RouteInput& in, const Evaluation&) {
         return in.d.kind == detail::Kind::Value && in.u.free.empty() && detail::reads_infinity(in.u.f) &&
                (!in.engine.strict || in.engine.arithmetic == Arithmetic::Double);
     },
     detail::at_infinity},
    {Route::Evaluate, "evaluate",
     "Evaluated as written, with nothing rewritten: in fractions where every operation is one on "
     "fractions, else in intervals, each number held as it was typed and every operation rounding "
     "outward, so the true value lies between the ends; the digits both ends share are proved.",
     [](const detail::RouteInput& in, const Evaluation&) { return detail::is_value(in) && in.u.free.empty(); }, detail::evaluate_value},
    {Route::Expression, "expression", "The expression as written, a function of the parameter it reads.",
     [](const detail::RouteInput& in, const Evaluation&) { return detail::is_value(in) && !in.u.free.empty(); }, detail::expression},
    {Route::Derivative, "derivative",
     "Differentiated by rules: each operation by its own, the chain rule between them, a function of the "
     "language by the rule its row of the table writes. The result is checked against the derivative a dual "
     "number computes from the function by its own rules, both in intervals, at sixteen points of the "
     "variable's interval, or at the point it is taken at, where its value is proved as a value's is.",
     [](const detail::RouteInput& in, const Evaluation&) { return in.d.kind == detail::Kind::Derivative; }, detail::derivative},
    {Route::ClosedForm, "closed form",
     "Solved in closed form. The coefficients are read from the equation's own values and the degree "
     "is checked at more points; the roots come from the quadratic formula in the form that keeps the "
     "small root, and Vieta's formulas check them by a second route.",
     [](const detail::RouteInput& in, const Evaluation&) { return detail::is_root(in); }, detail::closed_form},
    {Route::Linear, "linear system",
     "Solved by elimination. Each equation's coefficients are read from its own values and that it is "
     "linear is checked at more points; in fractions where every coefficient is one, so the solution "
     "is exact, else in doubles.",
     [](const detail::RouteInput& in, const Evaluation&) { return detail::is_equation(in) && !in.d.also.empty(); }, detail::linear},
    {Route::Intervals, "interval method",
     "The interval method. The axis is cut where the inequality is undecided on a piece, computed in "
     "intervals, until every piece is proved to hold it or not; between a piece where it holds and one "
     "where it does not is a root of its left side less its right, proved alone there by Newton's "
     "interval step and an end of the set, or a break where it has no value.",
     [](const detail::RouteInput& in, const Evaluation&) { return detail::is_number(in) && in.d.kind == detail::Kind::Root && in.d.relation != "="; },
     detail::intervals},
    {Route::Enclosed, "every root, proved",
     "Every root on the box, each proved. The box is cut in two where it is undecided; a piece is "
     "excluded where the equation's interval leaves out zero, and holds exactly one root where Newton's "
     "step in intervals, with the derivative enclosed by a dual number, lands inside it. Where the "
     "derivative may vanish, at a double root, nothing is proved. With no bracket the box is every double.",
     [](const detail::RouteInput& in, const Evaluation& out) {
         if (!in.d.also.empty()) return detail::is_equation(in) && detail::not_linear(out);
         return detail::is_root(in) && (!in.engine.strict || in.engine.arithmetic == Arithmetic::Double);
     },
     detail::enclosed},
    {Route::Bracket, "root on a bracket", "A root where the equation changes sign on its bracket, by bisection to the last bit.",
     [](const detail::RouteInput& in, const Evaluation&) { return detail::is_root(in) && spatium::io::build::has_bracket(in.d); }, detail::bracket},
    {Route::Integral, "integral", "The integral by the library's quadrature, its error estimated, not proved.",
     [](const detail::RouteInput& in, const Evaluation&) { return detail::is_number(in) && in.d.kind == detail::Kind::Integral; }, detail::integral},
    {Route::Synthesis, "synthesis",
     "Synthesised. The library computes the question's value at each parameter by its own solvers; the "
     "answer is the shortest expression in the parameter that matches it, searched bottom up, or a "
     "Chebyshev polynomial where none of a few terms does, checked against those values on a dense grid.",
     [](const detail::RouteInput& in, const Evaluation&) {
         return !detail::is_value(in) && in.d.kind != detail::Kind::Derivative && !in.u.free.empty();
     },
     detail::synthesis},
};

static_assert([] {
    for (std::size_t i = 0; i < std::size(route_table); ++i)
        if (route_table[i].route != static_cast<Route>(i)) return false;
    return true;
}(), "route_table is in the order of Route");

inline const char* route_name(Route r) { return route_table[static_cast<std::size_t>(r)].name; }
// How a route finds its value, in words a page shows under "How it was found".
inline const char* route_caption(Route r) { return route_table[static_cast<std::size_t>(r)].caption; }

// The choice among the candidates, made in this one place: by the rule of the table's order, the
// first that gives a value.
inline std::optional<std::size_t> take(const std::vector<Candidate>& candidates) {
    for (std::size_t i = 0; i < candidates.size(); ++i)
        if (candidates[i].refused.empty()) return i;
    return std::nullopt;
}

inline Evaluation evaluate(const spatium::io::build::Document& doc, std::string_view name, const Chosen& chosen = {},
                           const Engine& engine = {}) {
    namespace bd = spatium::io::build;
    Evaluation out;
    out.name = std::string(name);
    const bd::Statement* d = doc.find(name);
    if (!d) {
        out.refused = std::format("the document defines no '{}'", name);
        return out;
    }
    if (d->kind == bd::StatementKind::Shape) {
        out.refused = std::format("{} is a shape, which a scene builds (to_trace) and nothing evaluates", d->name);
        return out;
    }
    auto u = detail::unfolded(doc, *d, chosen);
    if (!u) {
        out.refused = u.error().message;
        return out;
    }
    out.free = u->free;
    if (out.free.size() > 1) {
        std::string names;
        for (const auto& n : out.free) names += (names.empty() ? "" : ", ") + n;
        out.refused = std::format("{} reads the parameters {}; give all of them but one a value", d->name, names);
        return out;
    }
    const detail::RouteInput in{doc, *d, *u, chosen, engine};
    for (const detail::RouteRow& row : route_table)
        if (row.applies(in, out)) row.run(in, out);
    out.taken = take(out.candidates);
    if (out.taken) {
        out.values = out.candidates[*out.taken].values;
        return out;
    }
    // each refusal under its route's name, unless it already starts with it
    for (const auto& c : out.candidates) {
        const std::string route = std::string(route_name(c.route)) + ": ";
        out.refused += std::string(out.refused.empty() ? "" : "; ") +
                       (out.candidates.size() > 1 && !c.refused.starts_with(route) ? route : "") + c.refused;
    }
    return out;
}

} // namespace spatium::answer
