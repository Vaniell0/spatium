#pragma once

// An equation the library solves in closed form, typed as text, with the steps it was solved by.
//
// `2x^2 - 3x - 5 = 0` is read into a field in one unknown; its coefficients are the field's
// values at x = 0, 1 and -1 (c = f(0), b = (f(1) - f(-1)) / 2, a = (f(1) + f(-1)) / 2 - c), and
// that it is a polynomial of degree at most two is checked, not assumed: the quadratic with those
// coefficients must give the field's values at five more points, or the equation is refused by
// name. No expression is rewritten. The roots come from `solve_quadratic_steps`, the solver
// itself, so the steps a page shows are the ones the answer was computed through; the last step
// checks the two roots by Vieta's formulas, a second route that says when the first lost digits.
//
// This is the closed-form side of the calculator; a question with a parameter whose answer has
// to be found goes to `ask` and synthesis (query.hpp).

#include <spatium/algebra/polynomial.hpp>
#include <spatium/io/field_parse.hpp>
#include <spatium/io/steps_text.hpp>

#include <array>
#include <cmath>
#include <format>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace rsc {

struct ClosedForm {
    std::string unknown;
    spatium::QuadraticSteps<double> solution;
    std::optional<spatium::VietaCheck<double>> check;   // of the roots, when there are two
    std::vector<spatium::io::build::Step> steps;        // the solution, and the check last
    std::size_t evaluations = 0;                        // of the equation's field, the route's cost
};

// The equation f = 0 in its one unknown, solved in closed form.
inline spatium::Result<ClosedForm> solve_closed_form(const spatium::io::build::Field<double>& f, const std::string& x) {
    using spatium::Error, spatium::ErrorCode;
    namespace bd = spatium::io::build;
    if (const auto names = f.free_symbols(); names.size() != 1 || names.front() != x)
        return std::unexpected(Error(ErrorCode::InvalidArgument,
            std::format("closed form: the equation must read exactly one unknown, it reads {}", names.size())));
    std::size_t evaluations = 0;
    const auto at = [&](double v) { return ++evaluations, f.bind(x, bd::Field<double>(v))(0.0, 0.0); };
    const double c = at(0), up = at(1), down = at(-1);
    const double a = (up + down) / 2 - c, b = (up - down) / 2;
    for (double v : {2.0, -2.0, 0.5, 3.0, -7.0}) {
        const double want = at(v), got = (a * v + b) * v + c;
        const double scale = std::max({std::abs(want), std::abs(a * v * v), std::abs(b * v), std::abs(c), 1.0});
        if (!(std::abs(want - got) <= 1e-12 * scale))
            return std::unexpected(Error(ErrorCode::InvalidArgument,
                std::format("closed form: the equation is not a polynomial of degree at most two in {} "
                            "(at {} = {} it is {}, a quadratic through x = 0, 1, -1 gives {})", x, x, v, want, got)));
    }
    ClosedForm out{x, spatium::solve_quadratic_steps(a, b, c), std::nullopt, {}, evaluations};
    out.check = spatium::check_vieta(out.solution);
    out.steps = bd::to_steps(out.solution, x);
    if (out.check) out.steps.push_back(bd::to_step(*out.check, x));
    return out;
}

inline spatium::Result<ClosedForm> solve_closed_form(std::string_view text) {
    using spatium::Error, spatium::ErrorCode;
    namespace bd = spatium::io::build;
    bd::FieldReader r(text, "", "");
    auto lhs = r.expression();
    if (!lhs) return std::unexpected(lhs.error());
    if (!r.take('=')) return r.fail("'='");
    auto rhs = r.expression();
    if (!rhs) return std::unexpected(rhs.error());
    if (!r.done()) return r.fail("the end");
    const auto f = std::move(*lhs) - *rhs;
    const auto names = f.free_symbols();
    return solve_closed_form(f, names.size() == 1 ? names.front() : std::string());
}

} // namespace rsc
