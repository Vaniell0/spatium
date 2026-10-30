#pragma once

// A question with a free symbol, answered by synthesis instead of by rewriting.
//
// The question is written as a field of the scene IR in two symbols: U is the unknown (or the
// variable of integration) and V the parameter the answer is a function of. What the question
// means is computed by the library -- a root on a bracket for `solve`, the integration door for
// `integrate` -- and that number, for each value of the parameter, is the reference. The answer
// is a field of the parameter synthesised by synthesis.hpp on the parameter's interval, in pieces
// where one expression does not reach the tolerance, with the error of every piece measured
// against the reference. No expression is ever transformed: where the operations hold a closed
// form it is found, where they do not the pieces say so.
//
// The question can be written in named symbols instead (x * x - p, with `unknown = "x"` and
// `parameter = "p"`): they are bound to U and V before anything is evaluated, and a symbol the
// question reads but does not name is a question not yet asked, so it is refused by name.
//
// The answer's field reads the parameter as U, its own variable.
//
// A question can also be written as text, the way a person types it: `ask` reads
// `solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]` (or `integrate exp(p * log(x)) for x in
// [0, 1], p in [0, 4]`) into the Query with those names, so what it answers is the question
// built in code from the same symbols, to the last bit.

#include <rsc/spatium/synthesis.hpp>
#include <spatium/io/field_parse.hpp>
#include <spatium/algebra/quadrature.hpp>
#include <spatium/algebra/roots.hpp>

#include <cmath>
#include <format>
#include <limits>
#include <optional>
#include <string>
#include <string_view>

namespace rsc {

enum class Ask { Solve, Integrate };

struct Query {
    Ask ask = Ask::Solve;
    Expr f;                // f(U = x, V = p), or f of the symbols named below
    std::string unknown{};   // the symbol bound to U; empty when f already reads U
    std::string parameter{}; // the symbol bound to V; empty when f already reads V
    double x_lo = 0, x_hi = 1; // the bracket of the root, or the limits of the integral
    // a limit that moves with the parameter, as a field of V (the parameter): Si(x) is the integral
    // of sin(t) / t over [0, x]; when set it is read instead of the number beside it
    std::optional<Expr> x_lo_at{}, x_hi_at{};
    double p_lo = 0, p_hi = 1; // where the answer is wanted
    double tolerance = 1e-8;
    bool relative = true;
    std::vector<Unary> unary{Unary::Sqrt};
    std::size_t max_terms = 4;
    double min_width = 1e-6;   // of a piece, in the parameter
    std::size_t polynomial_degree = 32;   // a piece no short expression answers is a Chebyshev polynomial up to this
};

struct Answer {
    std::vector<Solution> pieces;
    bool answered = false;   // every piece accepted
    std::string refused;     // why there is no answer at all, when there is none
};

// The meaning of the question at one value of the parameter, or NaN where it has none (a
// bracket with no change of sign; an integral the door does not trust).
inline double reference(const Query& q, double p) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double lo = q.x_lo_at ? (*q.x_lo_at)(0.0, p) : q.x_lo, hi = q.x_hi_at ? (*q.x_hi_at)(0.0, p) : q.x_hi;
    if (!(lo < hi)) return nan;
    if (q.ask == Ask::Solve) {
        const auto root = spatium::find_root([&](double x) { return q.f(x, p); }, lo, hi);
        return root ? *root : nan;
    }
    const auto r = spatium::integrate([&](double x) { return q.f(x, p); }, spatium::Finite<double>{lo, hi});
    return r.trusted() ? r.value : nan;
}

// The question with its named symbols bound to U and V.
inline Query bound(Query q) {
    if (!q.unknown.empty()) q.f = q.f.bind(q.unknown, Expr::u());
    if (!q.parameter.empty()) {
        q.f = q.f.bind(q.parameter, Expr::v());
        if (q.x_lo_at) q.x_lo_at = q.x_lo_at->bind(q.parameter, Expr::v());
        if (q.x_hi_at) q.x_hi_at = q.x_hi_at->bind(q.parameter, Expr::v());
    }
    q.unknown.clear();
    q.parameter.clear();
    return q;
}

inline Answer answer(const Query& asked) {
    Answer out;
    const Query q = bound(asked);
    if (const auto open = q.f.free_symbols(); !open.empty()) {
        out.refused = "the question reads the symbol '" + open.front() + "', which it does not name";
        return out;
    }
    // a question that has no meaning somewhere on the parameter's interval is refused, not fitted
    for (int i = 0; i <= 16; ++i) {
        const double p = q.p_lo + (q.p_hi - q.p_lo) * i / 16.0;
        if (!std::isfinite(reference(q, p))) {
            out.refused = "the question has no value at parameter " + std::to_string(p);
            return out;
        }
    }
    SynthesisSpec spec{.lo = q.p_lo, .hi = q.p_hi, .tolerance = q.tolerance, .relative = q.relative,
                       .max_terms = q.max_terms, .unary = q.unary, .polynomial_degree = q.polynomial_degree};
    out.pieces = synthesize_piecewise([&](double p) { return reference(q, p); }, spec, q.min_width);
    out.answered = std::ranges::all_of(out.pieces, [](const Solution& s) { return s.accepted; });
    return out;
}

// `solve f [= g] for x in [a, b], p in [c, d]` or `integrate f for x in [a, b], p in [c, d]`:
// the first name is the unknown (or the variable of integration) on its bracket, the second the
// parameter on the interval where the answer is wanted. An equation f = g is f - g, rewritten
// only where no value changes (so f = 0 is f). A bound is a number, or an expression in no name
// (sqrt(2), 1 / 3, pi); a bound of the unknown may read the parameter, so `integrate sin(t) / t
// for t in [0, x], x in [0.5, 10]` is the sine integral Si(x). The rest of the Query keeps its
// defaults.
inline spatium::Result<Query> ask(std::string_view text) {
    using spatium::Error, spatium::ErrorCode;
    spatium::io::build::FieldReader r(text, "", "", {"solve", "integrate", "for", "in"});   // every other name is a symbol
    Query q;
    if (r.take_word("solve")) q.ask = Ask::Solve;
    else if (r.take_word("integrate")) q.ask = Ask::Integrate;
    else return r.fail("'solve' or 'integrate'");
    auto f = r.expression();
    if (!f) return std::unexpected(f.error());
    if (q.ask == Ask::Solve && r.take('=')) {
        auto g = r.expression();
        if (!g) return std::unexpected(g.error());
        f = (std::move(*f) - *g).simplified();
    }
    q.f = std::move(*f);
    if (!r.take_word("for")) return r.fail(q.ask == Ask::Solve ? "'=' or 'for'" : "'for'");

    // a bound is an expression; one of the unknown's bounds may read the parameter, named after it
    const auto bound_of = [&]() -> spatium::Result<Expr> { return r.expression(); };
    const auto number = [](const Expr& e) { return e(0.0, 0.0); };
    Expr ulo, uhi;
    if (const std::string name(r.word()); name.empty()) return r.fail("a name");
    else {
        q.unknown = name;
        r.take_word(name);
    }
    if (!r.take_word("in")) return r.fail("'in'");
    if (!r.take('[')) return r.fail("'['");
    if (auto a = bound_of(); !a) return std::unexpected(a.error()); else ulo = *a;
    if (!r.take(',')) return r.fail("','");
    if (auto b = bound_of(); !b) return std::unexpected(b.error()); else uhi = *b;
    if (!r.take(']')) return r.fail("']'");
    if (!r.take(',')) return r.fail("',' and the parameter the answer is a function of");
    if (const std::string name(r.word()); name.empty()) return r.fail("a name");
    else {
        q.parameter = name;
        r.take_word(name);
    }
    if (!r.take_word("in")) return r.fail("'in'");
    if (!r.take('[')) return r.fail("'['");
    Expr plo, phi;
    if (auto a = bound_of(); !a) return std::unexpected(a.error()); else plo = *a;
    if (!r.take(',')) return r.fail("','");
    if (auto b = bound_of(); !b) return std::unexpected(b.error()); else phi = *b;
    if (!r.take(']')) return r.fail("']'");
    if (q.parameter == q.unknown)
        return std::unexpected(Error(ErrorCode::ParseError,
            "ask: '" + q.unknown + "' is named twice, as the unknown and as the parameter"));
    for (const Expr* e : {&plo, &phi})
        if (const auto open = e->free_symbols(); !open.empty())
            return std::unexpected(Error(ErrorCode::ParseError,
                "ask: a bound of " + q.parameter + " reads the name '" + open.front() + "'; it is a number"));
    for (const Expr* e : {&ulo, &uhi})
        for (const auto& name : e->free_symbols())
            if (name != q.parameter)
                return std::unexpected(Error(ErrorCode::ParseError,
                    "ask: a bound of " + q.unknown + " reads the name '" + name + "'; it may read only the parameter " +
                    q.parameter));
    q.p_lo = number(plo), q.p_hi = number(phi);
    if (!(q.p_lo < q.p_hi))
        return std::unexpected(Error(ErrorCode::ParseError,
            std::format("ask: {} in [{}, {}] is empty", q.parameter, q.p_lo, q.p_hi)));
    const auto at = [&](const Expr& e, double p) { return e.bind(q.parameter, Expr(p))(0.0, 0.0); };
    if (ulo.free_symbols().empty()) q.x_lo = number(ulo); else q.x_lo_at = ulo;
    if (uhi.free_symbols().empty()) q.x_hi = number(uhi); else q.x_hi_at = uhi;
    // the unknown's interval must hold something wherever the answer is wanted
    for (int i = 0; i <= 16; ++i) {
        const double p = q.p_lo + (q.p_hi - q.p_lo) * i / 16.0;
        const double lo = at(ulo, p), hi = at(uhi, p);
        if (!(lo < hi))
            return std::unexpected(Error(ErrorCode::ParseError,
                q.x_lo_at || q.x_hi_at
                    ? std::format("ask: {} in [{}, {}] is empty at {} = {}", q.unknown, lo, hi, q.parameter, p)
                    : std::format("ask: {} in [{}, {}] is empty", q.unknown, lo, hi)));
    }
    if (!r.done()) return r.fail("the end");
    return q;
}

} // namespace rsc
