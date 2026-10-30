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

#include <rsc/spatium/closed_form.hpp>
#include <rsc/spatium/query.hpp>
#include <spatium/algebra/quadrature.hpp>
#include <spatium/algebra/roots.hpp>
#include <spatium/io/document.hpp>

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

struct Candidate {
    Route route = Route::Evaluate;
    std::size_t cost = 0;        // evaluations of the definition's field
    std::string refused;         // why the route gives no value; empty when it gives one
    std::vector<double> values;  // the numbers it gives, when it gives numbers
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

// The route whose value is taken: the first, in the order of preference, that gives one.
inline std::optional<std::size_t> take(const std::vector<Candidate>& candidates) {
    for (std::size_t i = 0; i < candidates.size(); ++i)
        if (candidates[i].refused.empty()) return i;
    return std::nullopt;
}

inline Evaluation evaluate(const spatium::io::build::Document& doc, std::string_view name, const Chosen& chosen = {},
                           double tolerance = 1e-8) {
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
        if (out.free.empty()) {
            out.candidates.push_back({Route::Evaluate, 1, {}, {number(u->f)}});
        } else {
            out.written = u->f;   // the function, exact, in the parameter's name
            out.candidates.push_back({Route::Expression, 0, {}, {}});
        }
        finish();
        return out;
    }

    if (!out.free.empty()) {
        Candidate c{Route::Synthesis, 0, {}, {}};
        if (auto q = query_of(doc, name, chosen); !q) {
            c.refused = q.error().message;
        } else {
            q->tolerance = tolerance;
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
        Candidate c{Route::Integral, evaluations, {}, {}};
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
        Candidate c{Route::ClosedForm, 0, {}, {}};
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
        Candidate c{Route::Bracket, evaluations, root ? std::string() : root.error().message, {}};
        if (root) c.values = {*root};
        out.candidates.push_back(std::move(c));
    }
    finish();
    return out;
}

} // namespace rsc
