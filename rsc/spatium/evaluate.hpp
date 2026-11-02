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
#include <spatium/core/precision.hpp>
#include <spatium/core/rational.hpp>
#include <spatium/io/document.hpp>
#include <spatium/io/field_pod.hpp>

#include <charconv>
#include <cmath>
#include <cstddef>
#include <format>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace rsc {

enum class Route { Evaluate, Expression, ClosedForm, Bracket, Integral, Synthesis };

inline const char* route_name(Route r) {
    switch (r) {
        case Route::Evaluate: return "evaluate";
        case Route::Expression: return "expression";
        case Route::ClosedForm: return "closed form";
        case Route::Bracket: return "root on a bracket";
        case Route::Integral: return "integral";
        case Route::Synthesis: return "synthesis";
    }
    return "?";
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
    std::optional<double> error;         // of an integral, as estimated
    std::optional<Answer> answer;        // of a synthesis, in pieces with their certificate
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
        for (const auto& c : out.candidates)
            out.refused += std::string(out.refused.empty() ? "" : "; ") +
                           (out.candidates.size() > 1 ? std::string(route_name(c.route)) + ": " : "") + c.refused;
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

    // a root: in closed form where its equation is one, and on its bracket where it has one
    const bool bracket = bd::has_bracket(*d);
    const double lo = bracket ? number(u->lo) : 0, hi = bracket ? number(u->hi) : 0;
    {
        Candidate c{Route::ClosedForm, 0, {}, {}, Arithmetic::Double, {}};
        if (auto solved = solve_closed_form(u->f, d->name); !solved) {
            c.refused = solved.error().message;
            if (!bracket) c.refused += std::format(", so {} needs a bracket, as in {}: ... on [a, b]", d->name, d->name);
        } else {
            c.cost = solved->evaluations;
            for (const auto& z : solved->solution.roots)
                if (z.im == 0 && (!bracket || (lo <= z.re && z.re <= hi))) c.values.push_back(z.re);
            if (bracket && c.values.empty()) c.refused = std::format("no root of the closed form is on [{}, {}]", d->lo.text, d->hi.text);
            out.closed = std::move(*solved);
        }
        out.candidates.push_back(std::move(c));
    }
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
