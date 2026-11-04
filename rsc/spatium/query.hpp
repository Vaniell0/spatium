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
// A question is also a document (`query_of`): the parameter `p in [0.5, 4]` and the root
// `x: x^2 - p = 0 on [0, 10]` are the Query with those names, so what it answers is the question
// built in code from the same symbols, to the last bit.

#include <rsc/spatium/synthesis.hpp>
#include <spatium/io/document.hpp>
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
    std::size_t evaluations = 0;   // of the question's field, over every reference it took: the cost
};

// The meaning of the question at one value of the parameter, or NaN where it has none (a
// bracket with no change of sign; an integral the door does not trust).
// `evaluations`, when given, counts the evaluations of the field it took.
inline double reference(const Query& q, double p, std::size_t* evaluations = nullptr) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double lo = q.x_lo_at ? (*q.x_lo_at)(0.0, p) : q.x_lo, hi = q.x_hi_at ? (*q.x_hi_at)(0.0, p) : q.x_hi;
    if (!(lo < hi)) return nan;
    const auto f = [&](double x) {
        if (evaluations) ++*evaluations;
        return q.f(x, p);
    };
    if (q.ask == Ask::Solve) {
        const auto root = spatium::find_root(f, lo, hi);
        return root ? *root : nan;
    }
    const auto r = spatium::integrate(f, spatium::Finite<double>{lo, hi});
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
        if (!std::isfinite(reference(q, p, &out.evaluations))) {
            out.refused = asked.parameter.empty() ? std::format("the question has no value at parameter {}", p)
                                               : std::format("the question has no value at {} = {}", asked.parameter, p);
            return out;
        }
    }
    SynthesisSpec spec{.lo = q.p_lo, .hi = q.p_hi, .tolerance = q.tolerance, .relative = q.relative,
                       .max_terms = q.max_terms, .unary = q.unary, .polynomial_degree = q.polynomial_degree};
    out.pieces = synthesize_piecewise([&](double p) { return reference(q, p, &out.evaluations); }, spec, q.min_width);
    out.answered = std::ranges::all_of(out.pieces, [](const Solution& s) { return s.accepted; });
    return out;
}

// A document and the name asked of it: a question as a page writes it, its rows the document's
// lines (field_mathjson.hpp), so a question and a document ask the same thing.
struct AskedDocument {
    spatium::io::build::Document document;
    std::string wanted;
};

// Parameters given a value: a scene's, or a slider's. Each names a parameter of the document and
// lies on its interval.
using Chosen = std::vector<std::pair<std::string, double>>;

namespace detail {

// A definition with the values it reads put in place, latest first so a value reading an earlier
// one unfolds, and then the chosen parameters: what is left are the parameters it is a function of.
struct Unfolded {
    Expr f;                          // a value's expression, a root's left side less its right, an integrand
    // a root's left side less its right as written: `f` before its constants are folded, which
    // rounds 1/3 and sqrt(2) to doubles; what an exact route reads
    Expr written;
    Expr lo, hi;                     // the bracket or the limits, where written
    std::string bound;               // the name it binds itself: a root's unknown, an integral's variable
    std::vector<std::string> free;   // the parameters it reads and is not given, in the order it reads them
};

inline spatium::Result<Unfolded> unfolded(const spatium::io::build::Document& doc, const spatium::io::build::Statement& d,
                                          const Chosen& chosen) {
    namespace bd = spatium::io::build;
    using spatium::Error, spatium::ErrorCode;
    const auto refuse = [](std::string m) { return std::unexpected(Error(ErrorCode::InvalidArgument, std::move(m))); };
    const auto number = [](const Expr& e) { return e(0.0, 0.0); };
    for (const auto& [name, value] : chosen) {
        const bd::Statement* p = doc.find(name);
        if (!p || p->kind != bd::StatementKind::Parameter)
            return refuse(std::format("{} is given a value, and the document has no parameter {}", name, name));
        if (!(number(p->lo.field) <= value && value <= number(p->hi.field)))
            return refuse(std::format("{} = {} is not on its interval [{}, {}]", name, value, p->lo.text, p->hi.text));
    }
    const auto unfold = [&](Expr e) {
        for (auto it = doc.statements.rbegin(); it != doc.statements.rend(); ++it)
            if (it->kind == bd::StatementKind::Value) e = e.bind(it->name, it->value.field);
        for (const auto& [name, value] : chosen) e = e.bind(name, Expr(value));
        return e;
    };
    Unfolded out;
    switch (d.kind) {
        case bd::StatementKind::Root:
            out.written = unfold(d.value.field - d.rhs.field);
            out.f = unfold((d.value.field - d.rhs.field).simplified()), out.bound = d.name;
            break;
        case bd::StatementKind::Integral: out.f = unfold(d.value.field), out.bound = d.variable; break;
        case bd::StatementKind::Parameter: out.f = unfold(Expr::symbol(d.name)); break;
        default: out.f = unfold(d.value.field); break;
    }
    out.lo = unfold(d.lo.field), out.hi = unfold(d.hi.field);
    for (const Expr* e : {&out.f, &out.lo, &out.hi})
        for (const auto& n : e->free_symbols()) {
            if (n == out.bound || std::ranges::find(out.free, n) != out.free.end()) continue;
            const bd::Statement* p = doc.find(n);
            if (!p || p->kind != bd::StatementKind::Parameter)
                return refuse(std::format("{} reads {}, which is not a parameter or a value", d.name, n));
            out.free.push_back(n);
        }
    return out;
}

} // namespace detail

// The question a document asks of one of its names: a root on a bracket or an integral that
// reads one parameter once the `chosen` ones are put in, with the values it reads put in place.
// The Query is the one `answer` synthesises.
inline spatium::Result<Query> query_of(const spatium::io::build::Document& doc, std::string_view wanted, const Chosen& chosen = {}) {
    namespace bd = spatium::io::build;
    using spatium::Error, spatium::ErrorCode;
    const auto refuse = [](std::string m) { return std::unexpected(Error(ErrorCode::InvalidArgument, m)); };
    const bd::Statement* d = doc.find(wanted);
    if (!d) return refuse(std::format("the document defines no '{}'", wanted));
    if (d->kind != bd::StatementKind::Root && d->kind != bd::StatementKind::Integral)
        return refuse(std::format("{} is a {}; a question asks for a root or an integral", d->name, bd::statement_kind_name(d->kind)));
    if (d->kind == bd::StatementKind::Root && !bd::has_bracket(*d))
        return refuse(std::format("{} has no bracket; a question is answered on one, as in on [a, b]", d->name));
    auto u = detail::unfolded(doc, *d, chosen);
    if (!u) return refuse(u.error().message);
    Query q;
    q.ask = d->kind == bd::StatementKind::Root ? Ask::Solve : Ask::Integrate;
    q.f = u->f;
    q.unknown = u->bound;
    const Expr &ulo = u->lo, &uhi = u->hi;
    if (u->free.size() != 1)
        return refuse(std::format("{} reads {} parameters; a question is answered as a function of one", d->name, u->free.size()));
    const bd::Statement& p = *doc.find(u->free.front());
    q.parameter = p.name;
    const auto number = [](const Expr& e) { return e(0.0, 0.0); };
    q.p_lo = number(p.lo.field), q.p_hi = number(p.hi.field);
    if (!(q.p_lo < q.p_hi))
        return std::unexpected(Error(ErrorCode::ParseError, std::format("{} in [{}, {}] is empty", q.parameter, q.p_lo, q.p_hi)));
    for (const Expr* e : {&ulo, &uhi})
        for (const auto& name : e->free_symbols())
            if (name != q.parameter)
                return std::unexpected(Error(ErrorCode::ParseError,
                    "a bound of " + q.unknown + " reads the name '" + name + "'; it may read only the parameter " + q.parameter));
    if (ulo.free_symbols().empty()) q.x_lo = number(ulo); else q.x_lo_at = ulo;
    if (uhi.free_symbols().empty()) q.x_hi = number(uhi); else q.x_hi_at = uhi;
    // the unknown's interval must hold something wherever the answer is wanted
    const auto at = [&](const Expr& e, double v) { return e.bind(q.parameter, Expr(v))(0.0, 0.0); };
    for (int i = 0; i <= 16; ++i) {
        const double v = q.p_lo + (q.p_hi - q.p_lo) * i / 16.0;
        const double lo = at(ulo, v), hi = at(uhi, v);
        if (!(lo < hi))
            return std::unexpected(Error(ErrorCode::ParseError,
                q.x_lo_at || q.x_hi_at
                    ? std::format("{} in [{}, {}] is empty at {} = {}", q.unknown, lo, hi, q.parameter, v)
                    : std::format("{} in [{}, {}] is empty", q.unknown, lo, hi)));
    }
    return q;
}

} // namespace rsc
