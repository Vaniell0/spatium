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

#include <rsc/spatium/closed_form.hpp>
#include <rsc/spatium/query.hpp>
#include <spatium/algebra/quadrature.hpp>
#include <spatium/algebra/roots.hpp>
#include <spatium/core/decimal.hpp>
#include <spatium/core/interval.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/core/rational.hpp>
#include <spatium/io/document.hpp>
#include <spatium/io/field_pod.hpp>
#include <spatium/io/steps_text.hpp>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace rsc {

enum class Route { Evaluate, Expression, ClosedForm, Linear, Intervals, Enclosed, Bracket, Integral, Synthesis };

inline constexpr Route routes[] = {Route::Evaluate,  Route::Expression, Route::ClosedForm,
                                   Route::Linear,    Route::Intervals,  Route::Enclosed,
                                   Route::Bracket,   Route::Integral,   Route::Synthesis};

inline const char* route_name(Route r) {
    switch (r) {
        case Route::Evaluate: return "evaluate";
        case Route::Expression: return "expression";
        case Route::ClosedForm: return "closed form";
        case Route::Linear: return "linear system";
        case Route::Intervals: return "interval method";
        case Route::Enclosed: return "every root, proved";
        case Route::Bracket: return "root on a bracket";
        case Route::Integral: return "integral";
        case Route::Synthesis: return "synthesis";
    }
    return "?";
}

// How a route finds its value, in words a page shows under "How it was found".
inline const char* route_caption(Route r) {
    switch (r) {
        case Route::Evaluate:
            return "Evaluated as written, with nothing rewritten: in fractions where every operation is one on "
                   "fractions, else in intervals, each number held as it was typed and every operation rounding "
                   "outward, so the true value lies between the ends; the digits both ends share are proved.";
        case Route::Expression: return "The expression as written, a function of the parameter it reads.";
        case Route::ClosedForm:
            return "Solved in closed form. The coefficients are read from the equation's own values and the degree "
                   "is checked at more points; the roots come from the quadratic formula in the form that keeps the "
                   "small root, and Vieta's formulas check them by a second route.";
        case Route::Linear:
            return "Solved by elimination. Each equation's coefficients are read from its own values and that it is "
                   "linear is checked at more points; in fractions where every coefficient is one, so the solution "
                   "is exact, else in doubles.";
        case Route::Intervals:
            return "The interval method. The axis is cut where the inequality is undecided on a piece, computed in "
                   "intervals, until every piece is proved to hold it or not; between a piece where it holds and one "
                   "where it does not is a root of its left side less its right, proved alone there by Newton's "
                   "interval step and an end of the set, or a break where it has no value.";
        case Route::Enclosed:
            return "Every root on the box, each proved. The box is cut in two where it is undecided; a piece is "
                   "excluded where the equation's interval leaves out zero, and holds exactly one root where Newton's "
                   "step in intervals, with the derivative enclosed by a dual number, lands inside it. Where the "
                   "derivative may vanish, at a double root, nothing is proved. With no bracket the box is every double.";
        case Route::Bracket: return "A root where the equation changes sign on its bracket, by bisection to the last bit.";
        case Route::Integral: return "The integral by the library's quadrature, its error estimated, not proved.";
        case Route::Synthesis:
            return "Synthesised. The library computes the question's value at each parameter by its own solvers; the "
                   "answer is the shortest expression in the parameter that matches it, searched bottom up, or a "
                   "Chebyshev polynomial where none of a few terms does, checked against those values on a dense grid.";
    }
    return "";
}

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

// Why this build cannot compute in `a`; empty when it can.
inline std::string unavailable(Arithmetic a) {
    if (a == Arithmetic::Decimal && !SPATIUM_HAS_DECIMAL)
        return "this build has no decimal: its compiler has no decimal64 (GCC's <decimal/decimal>)";
#if !(defined(SPATIUM_HAS_BOOST_MULTIPRECISION) && SPATIUM_HAS_BOOST_MULTIPRECISION)
    if (a == Arithmetic::Real50) return "this build has no real50: it was built without Boost.Multiprecision";
#endif
    return {};
}

// How the engine computes, apart from what is asked. The defaults are a double and a synthesis
// tolerance of 1e-8.
struct Engine {
    Arithmetic arithmetic = Arithmetic::Double;
    bool strict = false;       // an arithmetic that cannot compute a value refuses, instead of a double
    double tolerance = 1e-8;   // of a synthesis
};

struct Candidate {
    Route route = Route::Evaluate;
    std::size_t cost = 0;        // evaluations of the definition's field
    std::string refused;         // why the route gives no value; empty when it gives one
    std::vector<double> values;  // the numbers it gives, when it gives numbers
    Arithmetic arithmetic = Arithmetic::Double;   // what it computed in
    std::string text;            // its number as that arithmetic holds it, when that is not a double (`21/4`)
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
    std::string parameter;               // the name the answer's pieces read as U

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
// followed a guessed branch.
template<spatium::Scalar B>
void prove(const spatium::io::build::PodField<double>& pod, Evaluation& out) {
    namespace bd = spatium::io::build;
    using I = spatium::Interval<B>;
    if (const auto op = bd::beyond<I>(pod)) {
        out.unproved = std::format("an interval has no {} in it", bd::op_name(*op));
        return;
    }
    const spatium::Undecided undecided;
    const I v = bd::interpret(bd::converted<I>(pod), I{0}, I{0});
    if (undecided.raised()) out.unproved = "a comparison on the way was undecided, or a root or logarithm reached below zero";
    else if (v.is_empty()) out.unproved = "the value is not a number";
    else out.proved = spatium::proved(v);
}

inline void prove(const Expr& f, Arithmetic a, Evaluation& out) {
    namespace bd = spatium::io::build;
    const auto pod = bd::lower(f);
    if (!pod) {
        out.unproved = pod.error().message;
        return;
    }
    // where every operation is one on fractions, the value is the fraction it is, exactly: the
    // twelve a page of exercises asks for, not 12.000000000000002 and its proved digits
    if (!bd::beyond<spatium::Rational>(*pod)) {
        const auto r = bd::interpret(bd::converted<spatium::Rational>(*pod), spatium::Rational(0), spatium::Rational(0));
        if (r.ok()) {
            const auto decimal = r.decimal();
            out.exact = r;
            out.proved = spatium::Proved{decimal ? static_cast<int>(std::ranges::count_if(*decimal, [](char c) { return c >= '0' && c <= '9'; })) : 0,
                                         decimal ? *decimal : r.text(), std::abs(r.to_double()), true};
            return;
        }
    }
#if defined(SPATIUM_HAS_BOOST_MULTIPRECISION) && SPATIUM_HAS_BOOST_MULTIPRECISION
    if (a == Arithmetic::Real50) return prove<spatium::Real50>(*pod, out);
#endif
    (void)a;
    prove<double>(*pod, out);
}

// The roots of a field of one unknown on a box, each proved: computed in an interval and in a dual
// number of intervals, which the lowered program runs in unchanged.
inline Candidate enclosed_roots(const Expr& f, spatium::Interval<double> box, Evaluation& out) {
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

// An inequality f rel 0 on a box, by the interval method (Axis). `written` is f as written, its
// constants not folded, where a root that is a short fraction is checked to be one exactly.
inline Candidate interval_method(const Expr& f, const Expr& written, const std::string& rel, spatium::Interval<double> box, Evaluation& out) {
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
    out.axis = std::move(axis);
    return c;
}

} // namespace detail

// The route whose value is taken: the first, in the order of preference, that gives one.
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
    const auto number = [](const Expr& e) { return e(0.0, 0.0); };
    const auto finish = [&] {
        out.taken = take(out.candidates);
        if (out.taken) {
            out.values = out.candidates[*out.taken].values;
            return;
        }
        // each refusal under its route's name, unless it already starts with it
        for (const auto& c : out.candidates) {
            const std::string name = std::string(route_name(c.route)) + ": ";
            out.refused += std::string(out.refused.empty() ? "" : "; ") +
                           (out.candidates.size() > 1 && !c.refused.starts_with(name) ? name : "") + c.refused;
        }
    };

    if (d->kind == bd::StatementKind::Value || d->kind == bd::StatementKind::Parameter) {
        out.written = d->kind == bd::StatementKind::Value ? d->value.field : Expr::symbol(d->name);
        if (out.free.empty() && engine.arithmetic == Arithmetic::Double) {
            out.candidates.push_back({Route::Evaluate, 1, {}, {number(u->f)}, Arithmetic::Double, {}});
        } else if (out.free.empty()) {
            // in the engine's arithmetic, and where it cannot, in a double unless it is strict
            out.candidates.push_back(detail::value_in(u->f, engine.arithmetic));
            if (!out.candidates.back().refused.empty() && !engine.strict)
                out.candidates.push_back({Route::Evaluate, 1, {}, {number(u->f)}, Arithmetic::Double, {}});
        } else {
            out.written = u->f;   // the function, exact, in the parameter's name
            out.candidates.push_back({Route::Expression, 0, {}, {}, Arithmetic::Double, {}});
        }
        if (out.free.empty()) {
            detail::prove(u->f, engine.arithmetic, out);
            if (auto steps = spatium::io::build::evaluation_steps(u->f)) out.steps = std::move(*steps);
        }
        finish();
        return out;
    }

    if (!out.free.empty()) {
        Candidate c{Route::Synthesis, 0, {}, {}, Arithmetic::Double, {}};
        if (auto q = query_of(doc, name, chosen); !q) {
            c.refused = q.error().message;
        } else {
            q->tolerance = engine.tolerance;
            Answer a = answer(*q);
            c.cost = a.evaluations, c.refused = a.refused;
            out.parameter = q->parameter;
            out.answer = std::move(a);
        }
        out.candidates.push_back(std::move(c));
        finish();
        return out;
    }

    std::size_t evaluations = 0;
    const Expr of_u = u->f.bind(u->bound, Expr::u());
    const auto field = [&](double x) {
        ++evaluations;
        return of_u(x, 0.0);
    };
    if (d->kind == bd::StatementKind::Integral) {
        const auto r = spatium::integrate(field, spatium::Finite<double>{number(u->lo), number(u->hi)});
        Candidate c{Route::Integral, evaluations, {}, {}, Arithmetic::Double, {}};
        if (r.trusted()) c.values = {r.value}, out.error = r.error_estimate;
        else c.refused = std::format("the integral of {} did not converge", d->name);
        out.candidates.push_back(std::move(c));
        finish();
        return out;
    }

    // an inequality: the set where it holds, by the interval method
    if (d->relation != "=") {
        const bool bracket = bd::has_bracket(*d);
        const double big = std::numeric_limits<double>::max();
        out.candidates.push_back(detail::interval_method(u->f.bind(d->name, Expr::u()), u->written.bind(d->name, Expr::u()), d->relation,
            bracket ? spatium::Interval<double>{number(u->lo), number(u->hi)} : spatium::Interval<double>{-big, big}, out));
        finish();
        return out;
    }

    // a system: every unknown together, by elimination where its equations are linear
    if (!d->also.empty()) {
        Candidate c{Route::Linear, 0, {}, {}, Arithmetic::Double, {}};
        std::vector<std::string> unknowns{d->name};
        for (const auto& e : d->also) unknowns.push_back(e.name);
        if (auto solved = solve_linear_system(u->system, unknowns, &u->system_written); !solved) {
            c.refused = solved.error().message;
        } else {
            c.cost = solved->evaluations, c.values = solved->values;
            if (!solved->exact.empty()) c.arithmetic = Arithmetic::Rational;
            out.linear = std::move(*solved);
        }
        const bool not_linear = c.refused.find("not linear") != std::string::npos;
        out.candidates.push_back(std::move(c));
        // where it is not linear, every root on every double, each proved
        if (not_linear) out.candidates.push_back(detail::enclosed_system(u->system, unknowns, out));
        finish();
        return out;
    }

    // a root: in closed form where its equation is one, and on its bracket where it has one
    const bool bracket = bd::has_bracket(*d);
    const double lo = bracket ? number(u->lo) : 0, hi = bracket ? number(u->hi) : 0;
    {
        // as fractions when the engine asks for them: the roots as they are, or why they are not
        const bool exact = engine.arithmetic == Arithmetic::Rational;
        Candidate c{Route::ClosedForm, 0, {}, {}, Arithmetic::Double, {}};
        // the roots as they are are always sought: a page shows them beside their decimals
        if (auto solved = solve_closed_form(u->f, d->name, &u->written); !solved) {
            c.refused = solved.error().message;
        } else {
            c.cost = solved->evaluations;
            for (const auto& z : solved->solution.roots)
                if (z.im == 0 && (!bracket || (lo <= z.re && z.re <= hi))) c.values.push_back(z.re);
            if (bracket && c.values.empty()) c.refused = std::format("no root of the closed form is on [{}, {}]", d->lo.text, d->hi.text);
            if (exact && solved->exact) {
                c.arithmetic = Arithmetic::Rational, c.text = solved->exact->text;
                // a fraction is the root itself, not the double's rounding of it
                if (!solved->exact->fractions.empty()) {
                    c.values.clear();
                    for (const auto& r : solved->exact->fractions)
                        if (!bracket || (lo <= r.to_double() && r.to_double() <= hi)) c.values.push_back(r.to_double());
                }
            } else if (exact) {
                Candidate refused{Route::ClosedForm, c.cost, solved->inexact, {}, Arithmetic::Rational, {}};
                if (engine.strict) c = std::move(refused);
                else out.candidates.push_back(std::move(refused));
            }
            out.closed = std::move(*solved);
        }
        out.candidates.push_back(std::move(c));
    }
    // every root on the bracket, or on the doubles' whole range without one, each proved; in
    // intervals of doubles, so not where a strict engine asks for another arithmetic alone
    if (!engine.strict || engine.arithmetic == Arithmetic::Double)
        out.candidates.push_back(detail::enclosed_roots(
            of_u, bracket ? spatium::Interval<double>{lo, hi}
                          : spatium::Interval<double>{-std::numeric_limits<double>::max(), std::numeric_limits<double>::max()},
            out));
    if (bracket) {
        evaluations = 0;
        const auto root = spatium::find_root(field, lo, hi);
        Candidate c{Route::Bracket, evaluations, root ? std::string() : root.error().message, {}, Arithmetic::Double, {}};
        if (root) c.values = {*root};
        out.candidates.push_back(std::move(c));
    }
    finish();
    return out;
}

} // namespace rsc
