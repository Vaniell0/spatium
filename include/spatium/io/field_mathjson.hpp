#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/io/field_parse.hpp>
#  include <spatium/io/json.hpp>
#  include <algorithm>
#  include <array>
#  include <cctype>
#  include <cmath>
#  include <format>
#  include <optional>
#  include <string>
#  include <string_view>
#  include <utility>
#  include <vector>
#endif

// ── An expression written as MathJSON ─────────────────────────────
//
// MathJSON is the JSON a LaTeX editor writes for a formula (MathLive's, and the Compute Engine
// behind it): a number, a symbol, or a function as an array whose first item names it, such as
// ["Add", ["Square", "x"], ["Negate", 2]]. A page reads what a person types with that editor and
// sends the MathJSON; this reads it as the text of the language (field_parse.hpp) and hands that
// text to the language's own reader, so a formula means one thing whichever way it came in and is
// refused for the same reasons. The text is kept as written, so a document read from MathJSON
// prints back as the language.
//
// What MathJSON can say and the language cannot -- a complex number, a derivative, a set -- is
// refused with the name it was given and where in the formula it stood ("at [2][1]"), which is
// what a page needs to point at it.
//
// A symbol is written {"sym": "x"}; a bare string is left to the caller, which in a document is
// the text of the language.

SPATIUM_EXPORT namespace spatium::io::build {

// A MathJSON function the language has, and how it is written there: `fn` a function of the
// language by name, `op` an operator joining the arguments, or neither, a form written out below.
struct MathJsonHead {
    std::string_view name;   // as MathJSON names it
    int min_args, max_args;  // -1: any number
    std::string_view caption;
};

inline constexpr std::array<MathJsonHead, 23> mathjson_heads{{
    {"Add", 2, -1, "a + b + ..."},
    {"Subtract", 2, 2, "a - b"},
    {"Negate", 1, 1, "-a"},
    {"Multiply", 2, -1, "a * b * ..."},
    {"Divide", 2, 2, "a / b"},
    {"Rational", 2, 2, "a / b"},
    {"Half", 1, 1, "a / 2"},
    {"Power", 2, 2, "a^b"},
    {"Square", 1, 1, "a^2"},
    {"Sqrt", 1, 1, "sqrt(a)"},
    {"Root", 2, 2, "a^(1/n)"},
    {"Exp", 1, 1, "exp(a)"},
    {"Ln", 1, 1, "log(a)"},
    {"Log", 1, 2, "log(a) / log(10), or log(a) / log(b)"},
    {"Sin", 1, 1, "sin(a)"},
    {"Cos", 1, 1, "cos(a)"},
    {"Tan", 1, 1, "sin(a) / cos(a)"},
    {"Abs", 1, 1, "abs(a)"},
    {"Max", 2, 2, "max(a, b)"},
    {"Min", 2, 2, "min(a, b)"},
    {"Delimiter", 1, 1, "(a)"},
    {"Equal", 2, 2, "a = b, the equation of a root"},
    {"Error", 0, -1, "what the editor could not read"},
}};

// MathJSON's named constants that the language has, as the language writes them
struct MathJsonConstant {
    std::string_view name;
    std::string_view text;
};

inline constexpr std::array<MathJsonConstant, 4> mathjson_constants{{
    {"Pi", "pi"},
    {"ExponentialE", "exp(1)"},
    {"PositiveInfinity", "inf"},
    {"NaN", "nan"},
}};

namespace detail {

// How tightly a text holds together, so a part is put in parentheses only where it must be:
// a sum, a product, a minus in front, a power, and a number, name or call.
enum class Tight { Sum, Product, Minus, Power, Atom };

struct MathText {
    std::string text;
    Tight tight;
};

inline std::unexpected<Error> mathjson_error(const std::string& path, std::string message) {
    return std::unexpected(Error(ErrorCode::ParseError, path.empty() ? message : std::format("at {}: {}", path, message)));
}

inline std::string wrapped(const MathText& t, Tight at_least) {
    return t.tight < at_least ? "(" + t.text + ")" : t.text;
}

inline MathText mathjson_number(double x) {
    if (std::isnan(x)) return MathText{"nan", Tight::Atom};
    if (std::isinf(x)) return MathText{x > 0 ? "inf" : "-inf", x > 0 ? Tight::Atom : Tight::Minus};
    std::string s = std::format("{}", std::abs(x));
    return MathText{x < 0 ? "-" + s : s, x < 0 ? Tight::Minus : Tight::Atom};
}

// A number MathJSON writes as a string, {"num": "..."}: the digits of a decimal, kept as written
// so 0.1 stays the shortest decimal it is (the engine's exact arithmetics read it that way).
inline Result<MathText> mathjson_digits(std::string_view s, const std::string& path) {
    if (s == "NaN") return MathText{"nan", Tight::Atom};
    if (s == "+Infinity" || s == "Infinity") return MathText{"inf", Tight::Atom};
    if (s == "-Infinity") return MathText{"-inf", Tight::Minus};
    std::string_view body = s;
    const bool negative = !body.empty() && body.front() == '-';
    if (negative || (!body.empty() && body.front() == '+')) body.remove_prefix(1);
    std::size_t k = 0, digits = 0;
    while (k < body.size() && std::isdigit(static_cast<unsigned char>(body[k]))) ++k, ++digits;
    if (k < body.size() && body[k] == '.') {
        ++k;
        while (k < body.size() && std::isdigit(static_cast<unsigned char>(body[k]))) ++k, ++digits;
    }
    if (digits > 0 && k < body.size() && (body[k] == 'e' || body[k] == 'E')) {
        std::size_t e = k + 1;
        if (e < body.size() && (body[e] == '+' || body[e] == '-')) ++e;
        const std::size_t first = e;
        while (e < body.size() && std::isdigit(static_cast<unsigned char>(body[e]))) ++e;
        if (e > first) k = e;
    }
    if (digits == 0 || k != body.size())
        return mathjson_error(path, std::format("'{}' is not a number the language writes (a repeating decimal or another base is not)", s));
    return MathText{(negative ? "-" : "") + std::string(body), negative ? Tight::Minus : Tight::Atom};
}

inline Result<MathText> mathjson_symbol(std::string_view name, const std::string& path) {
    for (const auto& c : mathjson_constants)
        if (c.name == name) return MathText{std::string(c.text), Tight::Atom};
    const bool word = !name.empty() && (std::isalpha(static_cast<unsigned char>(name.front())) || name.front() == '_') &&
                      std::ranges::all_of(name, [](char c) { return std::isalnum(static_cast<unsigned char>(c)) || c == '_'; });
    if (!word)
        return mathjson_error(path, std::format("the symbol '{}' is not a name the language reads (letters, digits and '_', not first a digit)", name));
    for (const auto& c : constants)
        if (c.name == name) return MathText{std::string(name), Tight::Atom};
    if (name == "ImaginaryUnit" || name == "ComplexInfinity")
        return mathjson_error(path, std::format("'{}' is a complex number; the language computes with real numbers", name));
    return MathText{std::string(name), Tight::Atom};
}

inline Result<MathText> mathjson_text(const JsonValue& v, const std::string& path);

// The same text with its sign turned over, when that is simply written: a sum shows `a - b`
// for a part that is a minus, a negative number, or a product led by one.
inline std::optional<MathText> mathjson_flipped(const JsonValue& v, const std::string& path) {
    if (v.is_number() && v.as_number() < 0) return mathjson_number(-v.as_number());
    if (const JsonValue* num = v.find("num"); v.is_object() && num && num->is_string() && num->as_string().starts_with('-')) {
        auto t = mathjson_digits(std::string_view(num->as_string()).substr(1), path);
        if (t) return *t;
    }
    if (!v.is_array() || v.as_array().empty() || !v.as_array()[0].is_string()) return std::nullopt;
    const auto& a = v.as_array();
    if (a[0].as_string() == "Negate" && a.size() == 2) {
        auto t = mathjson_text(a[1], path + "[1]");
        if (t) return *t;
    }
    if (a[0].as_string() == "Multiply" && a.size() >= 3 && a[1].is_number() && a[1].as_number() < 0) {
        JsonValue::ArrayType rest(a.begin(), a.end());
        if (a[1].as_number() == -1) rest.erase(rest.begin() + 1);
        else rest[1] = JsonValue(-a[1].as_number());
        const JsonValue flipped = rest.size() == 2 ? rest[1] : JsonValue(std::move(rest));
        auto t = mathjson_text(flipped, path);
        if (t) return *t;
    }
    return std::nullopt;
}

inline Result<MathText> mathjson_text(const JsonValue& v, const std::string& path) {
    if (v.is_number()) return mathjson_number(v.as_number());
    if (v.is_object()) {
        if (const JsonValue* n = v.find("num"); n && n->is_string()) return mathjson_digits(n->as_string(), path);
        if (const JsonValue* s = v.find("sym"); s && s->is_string()) return mathjson_symbol(s->as_string(), path);
        if (const JsonValue* f = v.find("fn"); f && f->is_array()) return mathjson_text(*f, path);
        return mathjson_error(path, "an object that is not {\"num\"}, {\"sym\"} or {\"fn\"}");
    }
    if (v.is_string()) return mathjson_symbol(v.as_string(), path);
    if (!v.is_array() || v.as_array().empty() || !v.as_array()[0].is_string())
        return mathjson_error(path, "not a number, a symbol or a function");

    const auto& a = v.as_array();
    const std::string head = a[0].as_string();
    const int n = static_cast<int>(a.size()) - 1;
    const auto known = std::ranges::find_if(mathjson_heads, [&](const MathJsonHead& h) { return h.name == head; });
    if (known == mathjson_heads.end())
        return mathjson_error(path, std::format("'{}' is not in the language", head));
    if (head == "Error") {
        std::string what;
        for (int k = 1; k <= n; ++k)
            if (a[k].is_string()) what += (what.empty() ? "" : ", ") + a[k].as_string();
        return mathjson_error(path, what.empty() ? std::string("the editor could not read the formula")
                                                 : std::format("the editor could not read the formula ({})", what));
    }
    if (head == "Equal") return mathjson_error(path, "an equation stands only as a whole formula, the equation of a root");
    if (n < known->min_args || (known->max_args >= 0 && n > known->max_args))
        return mathjson_error(path, known->min_args == known->max_args
                                        ? std::format("'{}' takes {} argument{}, not {}", head, known->min_args, known->min_args == 1 ? "" : "s", n)
                                        : std::format("'{}' takes {} or more arguments, not {}", head, known->min_args, n));

    std::vector<MathText> arg;
    for (int k = 1; k <= n; ++k) {
        auto t = mathjson_text(a[k], std::format("{}[{}]", path, k));
        if (!t) return std::unexpected(t.error());
        arg.push_back(std::move(*t));
    }
    const auto call = [&](std::string_view fn) {
        std::string s(fn);
        s += '(';
        for (std::size_t k = 0; k < arg.size(); ++k) s += (k ? ", " : "") + arg[k].text;
        return MathText{s + ")", Tight::Atom};
    };
    const auto quotient = [](const MathText& p, const MathText& q) {
        return MathText{wrapped(p, Tight::Product) + " / " + wrapped(q, Tight::Power), Tight::Product};
    };
    const auto power = [](const MathText& b, const MathText& e) {
        return MathText{wrapped(b, Tight::Atom) + "^" + wrapped(e, Tight::Atom), Tight::Power};
    };

    if (head == "Add") {
        std::string s = wrapped(arg[0], Tight::Sum);
        for (int k = 1; k < n; ++k) {
            const auto flipped = mathjson_flipped(a[k + 1], std::format("{}[{}]", path, k + 1));
            s += flipped ? " - " + wrapped(*flipped, Tight::Product) : " + " + wrapped(arg[k], Tight::Product);
        }
        return MathText{s, Tight::Sum};
    }
    if (head == "Subtract") return MathText{wrapped(arg[0], Tight::Sum) + " - " + wrapped(arg[1], Tight::Product), Tight::Sum};
    if (head == "Negate") return MathText{"-" + wrapped(arg[0], Tight::Power), Tight::Minus};
    if (head == "Multiply") {
        std::string s = wrapped(arg[0], Tight::Product);
        for (int k = 1; k < n; ++k) s += " * " + wrapped(arg[k], Tight::Minus);
        return MathText{s, Tight::Product};
    }
    if (head == "Divide" || head == "Rational") return quotient(arg[0], arg[1]);
    if (head == "Half") return quotient(arg[0], MathText{"2", Tight::Atom});
    if (head == "Power") return power(arg[0], arg[1]);
    if (head == "Square") return power(arg[0], MathText{"2", Tight::Atom});
    if (head == "Root") return power(arg[0], quotient(MathText{"1", Tight::Atom}, arg[1]));
    if (head == "Sqrt") return call("sqrt");
    if (head == "Exp") return call("exp");
    if (head == "Ln") return call("log");
    if (head == "Log") {
        const MathText base = n == 2 ? arg[1] : MathText{"10", Tight::Atom};
        return quotient(MathText{"log(" + arg[0].text + ")", Tight::Atom}, MathText{"log(" + base.text + ")", Tight::Atom});
    }
    if (head == "Sin") return call("sin");
    if (head == "Cos") return call("cos");
    if (head == "Tan") return quotient(MathText{"sin(" + arg[0].text + ")", Tight::Atom}, MathText{"cos(" + arg[0].text + ")", Tight::Atom});
    if (head == "Abs") return call("abs");
    if (head == "Max") return call("max");
    if (head == "Min") return call("min");
    return MathText{"(" + arg[0].text + ")", Tight::Atom};   // Delimiter
}

} // namespace detail

// A MathJSON expression as the text of the language, which the language's reader then reads.
inline Result<std::string> mathjson_text(const JsonValue& expression) {
    auto t = detail::mathjson_text(expression, "");
    if (!t) return std::unexpected(t.error());
    return std::move(t->text);
}

// A MathJSON equation, ["Equal", a, b], as the texts of its two sides.
inline Result<std::pair<std::string, std::string>> mathjson_equation(const JsonValue& equation) {
    const auto& a = equation.as_array();   // empty when it is not an array
    if (a.size() != 3 || !a[0].is_string() || a[0].as_string() != "Equal")
        return detail::mathjson_error("", "an equation is [\"Equal\", a, b]");
    auto l = detail::mathjson_text(a[1], "[1]");
    if (!l) return std::unexpected(l.error());
    auto r = detail::mathjson_text(a[2], "[2]");
    if (!r) return std::unexpected(r.error());
    return std::pair{std::move(l->text), std::move(r->text)};
}

} // namespace spatium::io::build
