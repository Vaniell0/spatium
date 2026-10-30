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

#include <rsc/spatium/synthesis.hpp>
#include <spatium/algebra/quadrature.hpp>
#include <spatium/algebra/roots.hpp>

#include <cmath>
#include <limits>
#include <string>

namespace rsc {

enum class Ask { Solve, Integrate };

struct Query {
    Ask ask = Ask::Solve;
    Expr f;                // f(U = x, V = p), or f of the symbols named below
    std::string unknown{};   // the symbol bound to U; empty when f already reads U
    std::string parameter{}; // the symbol bound to V; empty when f already reads V
    double x_lo = 0, x_hi = 1; // the bracket of the root, or the limits of the integral
    double p_lo = 0, p_hi = 1; // where the answer is wanted
    double tolerance = 1e-8;
    bool relative = true;
    std::vector<Unary> unary{Unary::Sqrt};
    std::size_t max_terms = 4;
    double min_width = 1e-6;   // of a piece, in the parameter
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
    if (q.ask == Ask::Solve) {
        const auto root = spatium::find_root([&](double x) { return q.f(x, p); }, q.x_lo, q.x_hi);
        return root ? *root : nan;
    }
    const auto r = spatium::integrate([&](double x) { return q.f(x, p); },
                                      spatium::Finite<double>{q.x_lo, q.x_hi});
    return r.trusted() ? r.value : nan;
}

// The question with its named symbols bound to U and V.
inline Query bound(Query q) {
    if (!q.unknown.empty()) q.f = q.f.bind(q.unknown, Expr::u());
    if (!q.parameter.empty()) q.f = q.f.bind(q.parameter, Expr::v());
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
                       .max_terms = q.max_terms, .unary = q.unary};
    out.pieces = synthesize_piecewise([&](double p) { return reference(q, p); }, spec, q.min_width);
    out.answered = std::ranges::all_of(out.pieces, [](const Solution& s) { return s.accepted; });
    return out;
}

} // namespace rsc
