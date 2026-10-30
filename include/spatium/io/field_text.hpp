#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/io/field.hpp>
#  include <algorithm>
#  include <array>
#  include <cmath>
#  include <cstdint>
#  include <format>
#  include <limits>
#  include <string>
#  include <string_view>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io::build {

// ── A field as text, and as LaTeX ─────────────────────────────────
//
// What a person reads of a field: its ops written out in the usual notation, the inputs under
// names the caller gives (u and v by default), a symbol under its own name, a constant in the
// shortest digits that read back to the same double. A synthesised answer is a field, and this
// is how it is read in code (`to_text`) and on a page (`to_latex`, which a page renders as it
// is). Both are one walk over the ops the root reads; they differ only in how each op is
// spelled. Parentheses go where precedence needs them and nowhere else, so u / (v * 2) keeps
// its parentheses and u + v * 2 has none; the text is read back by `parse` (field_parse.hpp)
// left to right, so in it u + (v + 1) keeps them too, while LaTeX, read as mathematics, drops them. The text is the field as built, not simplified, with
// one exception that changes no value: a factor of exactly 1 is left out, since x * 1 is x to the
// last bit for every x. A product of an op with itself is written as a power (x^2, x^3: the same op read
// again, which `parse` reads back into the same products). A field with an opaque leaf has no text and says why.

namespace detail {

enum class TextStyle { Plain, Latex };

// How tightly an expression binds, for deciding parentheses: a sum, a product, a power, then an
// atom (a name, a number, a call). A negative number binds like a sum, so x * (-2) keeps its
// parentheses.
enum class Binding { Sum = 1, Product = 2, Power = 3, Atom = 4 };

struct Spelled {
    std::string text;
    Binding binds = Binding::Atom;
    bool one = false;        // the constant 1 exactly, which a product leaves out
    bool negative = false;   // a negative constant, which a sum writes as a subtraction
    bool fraction = false;   // a LaTeX \frac, which a power's base puts in parentheses
    // a product of one op with itself, k times, multiplied left to right: written base^k
    std::uint32_t base = std::numeric_limits<std::uint32_t>::max();
    int power = 1;
};

inline std::string grouped(const Spelled& e, Binding need, TextStyle style) {
    if (e.binds >= need) return e.text;
    return style == TextStyle::Latex ? "\\left(" + e.text + "\\right)" : "(" + e.text + ")";
}

inline std::string number_text(double x, TextStyle style) {
    if (style == TextStyle::Plain) return std::format("{}", x);
    if (std::isnan(x)) return "\\mathrm{NaN}";
    if (std::isinf(x)) return x > 0 ? "\\infty" : "-\\infty";
    std::string s = std::format("{}", x);
    const auto e = s.find('e');
    if (e == std::string::npos) return s;
    // 1.8e-05 is 1.8 * 10^{-5}: no '+', no leading zeros, and no mantissa of exactly 1
    std::string exponent = s.substr(e + 1);
    std::string sign;
    if (!exponent.empty() && (exponent[0] == '+' || exponent[0] == '-')) {
        if (exponent[0] == '-') sign = "-";
        exponent.erase(0, 1);
    }
    exponent.erase(0, std::min(exponent.find_first_not_of('0'), exponent.size() - 1));
    const std::string mantissa = s.substr(0, e);
    const std::string power = std::format("10^{{{}{}}}", sign, exponent);
    if (mantissa == "1") return power;
    if (mantissa == "-1") return "-" + power;
    return std::format("{} \\cdot {}", mantissa, power);
}

// A name in LaTeX: a Greek letter by its command, one letter as it is, a longer name upright,
// and what follows the first '_' as a subscript (r_0 is r with subscript 0).
inline std::string name_latex(std::string_view name) {
    static constexpr std::array<std::string_view, 30> greek{
        "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa",
        "lambda", "mu", "nu", "xi", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi",
        "omega", "Gamma", "Delta", "Theta", "Lambda", "Sigma", "Phi", "Omega"};
    const auto bar = name.find('_');
    const std::string_view stem = name.substr(0, bar);
    std::string out;
    if (std::ranges::find(greek, stem) != greek.end()) out = std::format("\\{}", stem);
    else if (stem.size() == 1) out = std::string(stem);
    else out = std::format("\\mathrm{{{}}}", stem);
    if (bar != std::string_view::npos) out += std::format("_{{{}}}", name.substr(bar + 1));
    return out;
}

inline Spelled spell_op(const FieldOp<double>& n, const std::vector<Spelled>& at, const std::string& u,
                        const std::string& v, TextStyle style) {
    const bool tex = style == TextStyle::Latex;
    // The text is read back by `parse`, left to right, so a + (b + c) keeps its parentheses there:
    // floating-point addition does not reassociate. A page reads LaTeX as mathematics, where it does.
    const Binding right_of_sum = tex ? Binding::Sum : Binding::Product;
    const Binding right_of_product = tex ? Binding::Product : Binding::Power;
    const auto arg = [&](std::uint32_t i) -> const Spelled& { return at[i]; };
    const auto call = [&](std::string_view plain, std::string_view latex, const std::string& inside) {
        return tex ? std::format("{}\\left({}\\right)", latex, inside) : std::format("{}({})", plain, inside);
    };
    switch (n.op) {
        case Op::Const: {
            const Binding b = std::signbit(n.value) ? Binding::Sum : Binding::Atom;
            return {number_text(n.value, style), b, n.value == 1.0, std::signbit(n.value) && !std::isnan(n.value)};
        }
        case Op::U:      return {tex ? name_latex(u) : u};
        case Op::V:      return {tex ? name_latex(v) : v};
        case Op::Symbol: return {tex ? name_latex(n.name) : n.name};
        case Op::Add: {
            // a + (-c) is a - c to the last bit, and reads as one
            const auto& b = arg(n.b);
            if (b.negative) return {arg(n.a).text + " - " + b.text.substr(1), Binding::Sum};
            return {arg(n.a).text + " + " + grouped(b, right_of_sum, style), Binding::Sum};
        }
        case Op::Sub:
            return {arg(n.a).text + " - " + grouped(arg(n.b), Binding::Product, style), Binding::Sum};
        case Op::Mul: {
            if (arg(n.a).one) return arg(n.b);
            if (arg(n.b).one) return arg(n.a);
            // x * x is one op read twice in the pool, and x * x * x the same op once more on the
            // right: written x^2 and x^3, which parse reads back into the same products. A square
            // of a square stays (x^2)^2, since x^4 would be read as x * x * x * x, another order.
            const std::uint32_t base = n.a == n.b ? n.a : (arg(n.a).base == n.b ? n.b : arg(n.a).base);
            if (n.a == n.b || arg(n.a).base == n.b) {
                const int k = n.a == n.b ? 2 : arg(n.a).power + 1;
                const Spelled& b = arg(base);
                // a fraction, or a base that has a superscript of its own (x^{1}, a coordinate), is put in
                // parentheses in LaTeX, where x^{1}^{2} is not a formula
                const bool raised = tex && (b.fraction || b.text.find('^') != std::string::npos);
                const std::string under = raised ? "\\left(" + b.text + "\\right)" : grouped(b, Binding::Atom, style);
                Spelled out{tex ? std::format("{}^{{{}}}", under, k) : std::format("{}^{}", under, k), Binding::Power};
                out.base = base;
                out.power = k;
                return out;
            }
            return {grouped(arg(n.a), Binding::Product, style) + (tex ? " \\cdot " : " * ") +
                        grouped(arg(n.b), right_of_product, style),
                    Binding::Product};
        }
        case Op::Div:
            if (tex) {
                Spelled out{std::format("\\frac{{{}}}{{{}}}", arg(n.a).text, arg(n.b).text)};
                out.fraction = true;
                return out;
            }
            return {grouped(arg(n.a), Binding::Product, style) + " / " + grouped(arg(n.b), Binding::Power, style),
                    Binding::Product};
        case Op::Min:  return {call("min", "\\min", arg(n.a).text + ", " + arg(n.b).text)};
        case Op::Max:  return {call("max", "\\max", arg(n.a).text + ", " + arg(n.b).text)};
        case Op::Sin:  return {call("sin", "\\sin", arg(n.a).text)};
        case Op::Cos:  return {call("cos", "\\cos", arg(n.a).text)};
        case Op::Log:  return {call("log", "\\ln", arg(n.a).text)};
        case Op::Exp:  return {call("exp", "\\exp", arg(n.a).text)};
        case Op::Sqrt: return {tex ? std::format("\\sqrt{{{}}}", arg(n.a).text) : call("sqrt", "", arg(n.a).text)};
        case Op::Less:
            // 1 where a < b and 0 elsewhere: an Iverson bracket on a page
            return {tex ? std::format("\\left[{} < {}\\right]", arg(n.a).text, arg(n.b).text)
                        : call("less", "", arg(n.a).text + ", " + arg(n.b).text)};
        case Op::Noise:
            return {call("noise", "\\operatorname{noise}",
                         arg(n.a).text + ", " + arg(n.b).text + ", " + arg(n.c).text)};
        case Op::Gather:
            return {tex ? std::format("\\operatorname{{gather}}_{{{}}}\\left({}\\right)", n.k, arg(n.a).text)
                        : std::format("gather[{}]({})", n.k, arg(n.a).text)};
        case Op::Coord:  return {tex ? std::format("x^{{{}}}", n.k) : std::format("x[{}]", n.k)};
        case Op::Origin: return {tex ? std::format("o_{{{}}}", n.k) : std::format("origin[{}]", n.k)};
        case Op::Point:  return {tex ? std::format("p_{{{}}}", n.k) : std::format("point[{}]", n.k)};
        case Op::Hash:   return {tex ? std::format("h_{{{}}}", n.k) : std::format("hash[{}]", n.k)};
        case Op::Id:     return {tex ? std::string("\\mathrm{id}") : std::string("id")};
        case Op::Opaque: break;   // refused before the walk
    }
    return {"?"};
}

inline Result<std::string> spell(const Field<double>& whole, const std::string& u, const std::string& v,
                                 TextStyle style) {
    const Field<double> f = whole.reachable_from_root();
    std::vector<Spelled> at(f.size());
    for (std::size_t i = 0; i < f.size(); ++i) {
        const auto& n = f.op(i);
        if (n.op == Op::Opaque)
            return std::unexpected(Error(ErrorCode::InvalidArgument,
                std::format("to_text: op {} is Opaque -- a closure has no notation", i)));
        at[i] = spell_op(n, at, u, v, style);
    }
    return at.back().text;
}

} // namespace detail

inline Result<std::string> to_text(const Field<double>& f, std::string u = "u", std::string v = "v") {
    return detail::spell(f, u, v, detail::TextStyle::Plain);
}

inline Result<std::string> to_latex(const Field<double>& f, std::string u = "u", std::string v = "v") {
    return detail::spell(f, u, v, detail::TextStyle::Latex);
}

} // namespace spatium::io::build
