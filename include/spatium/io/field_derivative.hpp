#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/io/field.hpp>
#  include <spatium/io/field_parse.hpp>
#  include <array>
#  include <cassert>
#  include <cstdint>
#  include <format>
#  include <optional>
#  include <string>
#  include <string_view>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io::build {

// ── The derivative of a field, as a field ────────────────────────
//
// `derivative(f, "x")` is f differentiated in the symbol x, written as another field: an op at a
// time, in order, each op's derivative from its operands' values and their derivatives by the
// chain rule. The rules of + - * / and of a power are here; a function of the language carries
// its own in its row of the table (`Function::derivative`, written in the language over its
// arguments a and b: the sine's is cos(a)), so a function added to the table says its derivative
// where it says its name, and the page's reference shows it from the same row.
//
// An op that reads nothing that reads x has derivative zero and adds no term, and a derivative of
// exactly 1 is no factor: d/dx (2 * x) is 2, not 0 * x + 2 * 1, since 0 * x is not 0 where x is
// infinite and nothing would fold it. Constants are left as written, not folded: 1/3 stays a
// quotient, so an exact arithmetic computes the derivative of x^3 / 3 at 1 as 1. A comparison (`less`)
// and a gather are constant between their steps, so their derivative is 0 wherever there is one;
// noise and a net that read x are refused by name, since no rule for them is written.
//
// Where the derivative has no value -- abs at 0, sqrt at 0 -- the field has none either: a / abs(a)
// is NaN at 0, 1 / (2 sqrt(a)) infinite. A dual number computes the same derivative by other rules
// (algebra/dual.hpp), which is what an answer checks this one against (answer/evaluate.hpp).

namespace detail {

// The rule of a function in its argument k, read once from the table: its arguments the symbols
// "#a" and "#b", which no text of the language can write, so binding them to operands that read a
// symbol named a or b binds nothing twice.
inline const Field<double>& table_rule(const Function& fn, int k) {
    static const auto rules = [] {
        std::array<std::array<Field<double>, 2>, functions.size()> out{};
        for (std::size_t i = 0; i < functions.size(); ++i)
            for (int j = 0; j < functions[i].arity; ++j) {
                auto r = parse(functions[i].derivative[static_cast<std::size_t>(j)], "", "");
                assert(r && "a function's derivative in the table does not read");   // held by the tests
                out[i][static_cast<std::size_t>(j)] = r->bind("a", Field<double>::symbol("#a")).bind("b", Field<double>::symbol("#b"));
            }
        return out;
    }();
    return rules[static_cast<std::size_t>(&fn - functions.data())][static_cast<std::size_t>(k)];
}

inline bool is_constant(const Field<double>& f, double c) {
    return f.size() == 1 && f.op(0).op == Op::Const && f.op(0).value == c;
}

// a product, with a factor of exactly 1 left out
inline Field<double> times(const Field<double>& a, const Field<double>& b) {
    if (is_constant(a, 1)) return b;
    if (is_constant(b, 1)) return a;
    return a * b;
}

} // namespace detail

inline Result<Field<double>> derivative(const Field<double>& f, std::string_view x) {
    using F = Field<double>;
    const auto& ops = f.ops();
    // each op's derivative, none where it is zero because the op reads nothing that reads x
    std::vector<std::optional<F>> d(ops.size());
    std::vector<std::optional<F>> values(ops.size());
    const auto value = [&](std::uint32_t j) -> const F& {
        if (!values[j]) values[j] = f.subfield(j);
        return *values[j];
    };
    const auto sum = [](std::optional<F> a, std::optional<F> b) -> std::optional<F> {
        if (a && b) return std::move(*a) + *b;
        return a ? std::move(a) : std::move(b);
    };
    for (std::size_t i = 0; i < ops.size(); ++i) {
        const auto& n = ops[i];
        const auto self = static_cast<std::uint32_t>(i);
        const std::optional<F>& da = arity(n.op) > 0 ? d[n.a] : d[self];
        const std::optional<F>& db = arity(n.op) > 1 ? d[n.b] : d[self];
        switch (n.op) {
            case Op::Symbol:
                if (n.name == x) d[i] = F(1.0);
                break;
            case Op::Add: d[i] = sum(da, db); break;
            case Op::Sub:
                d[i] = sum(da, db ? std::optional<F>(F(-1.0) * *db) : std::nullopt);
                break;
            case Op::Mul:
                d[i] = sum(da ? std::optional<F>(detail::times(*da, value(n.b))) : std::nullopt,
                           db ? std::optional<F>(detail::times(value(n.a), *db)) : std::nullopt);
                break;
            case Op::Div:
                // (a / b)' = (a' - (a / b) b') / b
                if (da || db) d[i] = sum(da, db ? std::optional<F>(F(-1.0) * detail::times(value(self), *db)) : std::nullopt).value() / value(n.b);
                break;
            case Op::Pow:
                // (a^b)' = b a^(b - 1) a' + a^b log(a) b'
                d[i] = sum(da ? std::optional<F>(detail::times(value(n.b) * pow(value(n.a), value(n.b) - F(1.0)), *da)) : std::nullopt,
                           db ? std::optional<F>(detail::times(value(self) * log(value(n.a)), *db)) : std::nullopt);
                break;
            case Op::Gather: break;   // constant between the steps of its index
            case Op::Noise:
            case Op::Net:
                if (d[n.a] || d[n.b] || d[n.c])
                    return std::unexpected(Error(ErrorCode::NotImplemented,
                        std::format("derivative in {}: {} reads {}, and no rule for its derivative is written", x, op_name(n.op), x)));
                break;
            default:
                if (const Function* fn = function_of_op(n.op))
                    for (int k = 0; k < fn->arity; ++k) {
                        const std::optional<F>& dk = k == 0 ? d[n.a] : d[n.b];
                        const F& rule = detail::table_rule(*fn, k);
                        if (!dk || detail::is_constant(rule, 0)) continue;
                        F factor = rule.bind("#a", value(n.a));
                        if (fn->arity > 1) factor = factor.bind("#b", value(n.b));
                        d[i] = sum(std::move(d[i]), detail::times(factor, *dk));
                    }
                break;   // a constant, an input other than x, an opaque leaf: they read no symbol
        }
    }
    return d.back() ? *d.back() : F(0.0);
}

} // namespace spatium::io::build
