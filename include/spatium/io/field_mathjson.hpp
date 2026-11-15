#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/io/field_parse.hpp>
#  include <spatium/io/field_text.hpp>
#  include <spatium/io/json.hpp>
#  include <algorithm>
#  include <array>
#  include <cctype>
#  include <charconv>
#  include <cmath>
#  include <format>
#  include <limits>
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
// sends the MathJSON; this reads the tree straight into the field the language's reader builds for
// the same formula typed (field_parse.hpp), op for op, so a formula means one thing whichever way it
// came in. Where a document needs its text, the language's printer writes it (field_text.hpp):
// there is one spelling of the language, the library's.
//
// What MathJSON can say and the language cannot -- a complex number, a derivative, a set -- is
// refused with the name it was given and where in the formula it stood ("at [2][1]"), which is
// what a page needs to point at it.
//
// A symbol is written {"sym": "x"}, or as a bare string inside a formula; a field of a document
// that is a bare string is the text of the language. What the editor's own parser writes before
// any rewriting is read too: a product written without a sign (`2x`, `InvisibleOperator`) and a
// mixed number (`1\frac{13}{20}`, the same head before a fraction of whole numbers), the
// logarithms in base 2 and 10 (`Lb`, `Lg`), and `e` as the number e.
//
// A row of the editor is one line of a document, read by the forms of a statement
// (`mathjson_row_text`, document.hpp).

SPATIUM_EXPORT namespace spatium::io::build {

// A MathJSON head that is not a function of the language: a form (a sum, a power, a fraction) or a
// composition the language writes out (`Tan` is sin / cos). A function the language has is read by
// the head its row of `functions` names (field_parse.hpp), so the two tables cannot disagree.
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
    {"InvisibleOperator", 2, -1, "a b ..., a product written without a sign"},
    {"Divide", 2, 2, "a / b"},
    {"Rational", 2, 2, "a / b"},
    {"Half", 1, 1, "a / 2"},
    {"Power", 2, 2, "a^b"},
    {"Square", 1, 1, "a^2"},
    {"Root", 2, 2, "a^(1/n)"},
    {"Log", 1, 2, "log(a) / log(10), or log(a) / log(b)"},
    {"Lb", 1, 1, "log(a) / log(2)"},
    {"Lg", 1, 1, "log(a) / log(10)"},
    {"Tan", 1, 1, "sin(a) / cos(a)"},
    {"Boole", 1, 1, "[a < b], an Iverson bracket: less(a, b)"},
    {"Delimiter", 1, 1, "(a)"},
    {"Equal", 2, 2, "a = b, the equation of a root"},
    {"Less", 2, 2, "a < b, an inequality"},
    {"LessEqual", 2, 2, "a <= b, an inequality"},
    {"Greater", 2, 2, "a > b, an inequality"},
    {"GreaterEqual", 2, 2, "a >= b, an inequality"},
    {"Error", 0, -1, "what the editor could not read"},
}};

// MathJSON's named constants that the language has, as the language writes them
struct MathJsonConstant {
    std::string_view name;
    std::string_view text;
};

inline constexpr std::array<MathJsonConstant, 5> mathjson_constants{{
    {"Pi", "pi"},
    {"ExponentialE", "exp(1)"},
    {"e", "exp(1)"},   // the editor's parser writes e as typed; the language has no constant e
    {"PositiveInfinity", "inf"},
    {"NaN", "nan"},
}};

namespace detail {

inline std::unexpected<Error> mathjson_error(const std::string& path, std::string message) {
    return std::unexpected(Error(ErrorCode::ParseError, path.empty() ? message : std::format("at {}: {}", path, message)));
}

// A number MathJSON writes as a string, {"num": "..."}: the digits of a decimal, read as the
// language's reader reads them, so 0.1 is the double the text 0.1 is.
inline Result<Field<double>> mathjson_digits(std::string_view s, const std::string& path) {
    if (s == "NaN") return Field<double>(std::numeric_limits<double>::quiet_NaN());
    if (s == "+Infinity" || s == "Infinity") return Field<double>(std::numeric_limits<double>::infinity());
    if (s == "-Infinity") return Field<double>(-std::numeric_limits<double>::infinity());
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
    double x = 0;
    const auto [end, ec] = std::from_chars(body.data(), body.data() + body.size(), x);
    if (digits == 0 || k != body.size() || ec != std::errc{} || end != body.data() + body.size())
        return mathjson_error(path, std::format("'{}' is not a number the language writes (a repeating decimal or another base is not)", s));
    return Field<double>(negative ? -x : x);
}

inline Result<Field<double>> mathjson_symbol(std::string_view name, const std::string& path) {
    for (const auto& c : mathjson_constants)
        if (c.name == name) return parse(c.text, "", "");
    const bool word = !name.empty() && (std::isalpha(static_cast<unsigned char>(name.front())) || name.front() == '_') &&
                      std::ranges::all_of(name, [](char c) { return std::isalnum(static_cast<unsigned char>(c)) || c == '_'; });
    if (!word)
        return mathjson_error(path, std::format("the symbol '{}' is not a name the language reads (letters, digits and '_', not first a digit)", name));
    for (const auto& c : constants)
        if (c.name == name) return Field<double>(c.value);
    if (name == "ImaginaryUnit" || name == "ComplexInfinity")
        return mathjson_error(path, std::format("'{}' is a complex number; the language computes with real numbers", name));
    return Field<double>::symbol(std::string(name));
}

inline Result<Field<double>> mathjson_field(const JsonValue& v, const std::string& path);

// The same item with its sign turned over, where a person writes it so: a sum reads `a - b` for a
// part that is a minus, a negative number, or a product led by one -- the ops the language's reader
// builds for that text, so a formula typed and the same formula from the editor are one field.
inline std::optional<Result<Field<double>>> mathjson_flipped(const JsonValue& v, const std::string& path) {
    if (v.is_number() && v.as_number() < 0) return Result<Field<double>>(Field<double>(-v.as_number()));
    if (const JsonValue* num = v.is_object() ? v.find("num") : nullptr; num && num->is_string() && num->as_string().starts_with('-'))
        return mathjson_digits(std::string_view(num->as_string()).substr(1), path);
    if (!v.is_array() || v.as_array().empty() || !v.as_array()[0].is_string()) return std::nullopt;
    const auto& a = v.as_array();
    if (a[0].as_string() == "Negate" && a.size() == 2) return mathjson_field(a[1], path + "[1]");
    if ((a[0].as_string() == "Multiply" || a[0].as_string() == "InvisibleOperator") && a.size() >= 3 && a[1].is_number() &&
        a[1].as_number() < 0) {
        JsonValue::ArrayType rest(a.begin(), a.end());
        if (a[1].as_number() == -1) rest.erase(rest.begin() + 1);
        else rest[1] = JsonValue(-a[1].as_number());
        return mathjson_field(rest.size() == 2 ? rest[1] : JsonValue(std::move(rest)), path);
    }
    return std::nullopt;
}

// A whole power from 1 to 64, written as a number: x^n, which the reader takes as n factors
inline std::optional<int> mathjson_power(const JsonValue& e) {
    const JsonValue* num = e.is_object() ? e.find("num") : &e;
    double k = std::numeric_limits<double>::quiet_NaN();
    if (num && num->is_number()) k = num->as_number();
    else if (num && num->is_string() && std::ranges::all_of(num->as_string(), [](char c) { return c >= '0' && c <= '9'; }) &&
             !num->as_string().empty() && num->as_string().size() <= 2)
        k = std::stoi(num->as_string());
    if (k >= 1 && k <= 64 && std::floor(k) == k) return static_cast<int>(k);
    return std::nullopt;
}

// whether a MathJSON item is a whole number written as one: the parts of a mixed number
inline bool mathjson_whole(const JsonValue& x) {
    const JsonValue* num = x.is_object() ? x.find("num") : &x;
    if (num && num->is_string()) {
        const auto& t = num->as_string();
        return !t.empty() && std::ranges::all_of(t, [](char c) { return c >= '0' && c <= '9'; });
    }
    return num && num->is_number() && num->as_number() >= 0 && std::floor(num->as_number()) == num->as_number();
}

// A MathJSON expression read into a field, the ops the language's reader builds for the same text:
// a head of the table is an op of the field or a composition of them (`Tan` is sin / cos), and
// what the language cannot say is refused with its name and where in the formula it stood.
inline Result<Field<double>> mathjson_field(const JsonValue& v, const std::string& path) {
    if (v.is_number()) return Field<double>(v.as_number());
    if (v.is_object()) {
        if (const JsonValue* n = v.find("num"); n && n->is_string()) return mathjson_digits(n->as_string(), path);
        if (const JsonValue* s = v.find("sym"); s && s->is_string()) return mathjson_symbol(s->as_string(), path);
        if (const JsonValue* f = v.find("fn"); f && f->is_array()) return mathjson_field(*f, path);
        return mathjson_error(path, "an object that is not {\"num\"}, {\"sym\"} or {\"fn\"}");
    }
    if (v.is_string()) return mathjson_symbol(v.as_string(), path);
    if (!v.is_array() || v.as_array().empty() || !v.as_array()[0].is_string())
        return mathjson_error(path, "not a number, a symbol or a function");

    const auto& a = v.as_array();
    const std::string head = a[0].as_string();
    const int n = static_cast<int>(a.size()) - 1;
    // a function of the language, by the head its row names
    if (const auto f = std::ranges::find(functions, std::string_view(head), &Function::mathjson); !head.empty() && f != functions.end()) {
        if (n != f->arity)
            return mathjson_error(path, std::format("'{}' takes {} argument{}, not {}", head, f->arity, f->arity == 1 ? "" : "s", n));
        std::vector<Field<double>> arg;
        for (int k = 1; k <= n; ++k) {
            auto x = mathjson_field(a[k], std::format("{}[{}]", path, k));
            if (!x) return x;
            arg.push_back(std::move(*x));
        }
        return apply(*f, std::move(arg));
    }
    if (head == "Open") return mathjson_error(path, "an open end; an interval here holds its ends, [a, b]");
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
    // the bracket to_latex writes for less(a, b), [a < b], read back: 1 where a < b holds, else 0
    if (head == "Boole" && n == 1) {
        const auto& rel = a[1].as_array();   // empty when it is not an array
        const std::string inner = rel.size() == 3 && rel[0].is_string() ? rel[0].as_string() : "";
        if (inner != "Less" && inner != "Greater")
            return mathjson_error(path + "[1]", "an Iverson bracket here holds a strict inequality, [a < b] or [a > b]");
        auto x = mathjson_field(rel[1], path + "[1][1]");
        if (!x) return x;
        auto y = mathjson_field(rel[2], path + "[1][2]");
        if (!y) return y;
        return inner == "Less" ? less(std::move(*x), *y) : less(std::move(*y), *x);
    }
    if (head == "Equal" || head.starts_with("Less") || head.starts_with("Greater"))
        return mathjson_error(path, "an equation or an inequality stands only as a whole formula, of a root or of a set");
    if (n < known->min_args || (known->max_args >= 0 && n > known->max_args))
        return mathjson_error(path, known->min_args == known->max_args
                                        ? std::format("'{}' takes {} argument{}, not {}", head, known->min_args, known->min_args == 1 ? "" : "s", n)
                                        : std::format("'{}' takes {} or more arguments, not {}", head, known->min_args, n));
    // a whole number written right before a fraction of whole numbers is a mixed number, 1 13/20,
    // as a page of exercises writes it: the sum, not the product the bare juxtaposition would be
    const bool mixed = head == "InvisibleOperator" && n == 2 && mathjson_whole(a[1]) && a[2].is_array() &&
                       a[2].as_array().size() == 3 && a[2].as_array()[0].is_string() &&
                       (a[2].as_array()[0].as_string() == "Divide" || a[2].as_array()[0].as_string() == "Rational") &&
                       mathjson_whole(a[2].as_array()[1]) && mathjson_whole(a[2].as_array()[2]);

    std::vector<Field<double>> arg;
    for (int k = 1; k <= n; ++k) {
        auto f = mathjson_field(a[k], std::format("{}[{}]", path, k));
        if (!f) return f;
        arg.push_back(std::move(*f));
    }
    const auto fold = [&](auto join) {
        Field<double> out = arg[0];
        for (std::size_t k = 1; k < arg.size(); ++k) out = join(std::move(out), arg[k]);
        return out;
    };
    const auto log_of = [](Field<double> x, double base) { return log(std::move(x)) / log(Field<double>(base)); };
    if (mixed) return std::move(arg[0]) + arg[1];
    if (head == "Add") {
        Field<double> out = arg[0];
        for (int k = 1; k < n; ++k) {
            auto flipped = mathjson_flipped(a[k + 1], std::format("{}[{}]", path, k + 1));
            if (flipped && !*flipped) return std::move(*flipped);
            out = flipped ? std::move(out) - **flipped : std::move(out) + arg[k];
        }
        return out;
    }
    if (head == "Subtract") return std::move(arg[0]) - arg[1];
    // -2 is a number, -x a product with -1, as the reader takes a typed minus
    if (head == "Negate" && (a[1].is_number() || (a[1].is_object() && a[1].find("num")))) return Field<double>(-arg[0](0.0, 0.0));
    if (head == "Negate") return Field<double>(-1.0) * arg[0];
    if (head == "Multiply" || head == "InvisibleOperator")
        return fold([](Field<double> x, const Field<double>& y) { return std::move(x) * y; });
    if (head == "Divide" || head == "Rational") return std::move(arg[0]) / arg[1];
    if (head == "Half") return std::move(arg[0]) / Field<double>(2.0);
    if (head == "Square") return arg[0] * arg[0];
    // e to a power is exp of it, as \exp is: the two ways of writing it read the same
    if (head == "Power" && same_value(arg[0], exp(Field<double>(1.0)))) return exp(std::move(arg[1]));
    if (head == "Power") {
        // a whole power from 1 to 64 is a product of that many factors, as the reader takes x^n
        if (const auto k = mathjson_power(a[2])) {
            Field<double> out = arg[0];
            for (int i = 1; i < *k; ++i) out = std::move(out) * arg[0];
            return out;
        }
        return pow(std::move(arg[0]), arg[1]);
    }
    if (head == "Root") return pow(std::move(arg[0]), Field<double>(1.0) / arg[1]);
    if (head == "Log") return n == 2 ? log(std::move(arg[0])) / log(std::move(arg[1])) : log_of(std::move(arg[0]), 10);
    if (head == "Lb") return log_of(std::move(arg[0]), 2);
    if (head == "Lg") return log_of(std::move(arg[0]), 10);
    if (head == "Tan") return sin(arg[0]) / cos(arg[0]);
    return std::move(arg[0]);   // Delimiter
}

// a side of a statement, as the text the language's printer writes for it
inline Result<std::string> mathjson_side(const JsonValue& v, const std::string& path) {
    auto f = mathjson_field(v, path);
    if (!f) return std::unexpected(f.error());
    return to_text(*f, "", "");
}

} // namespace detail

// A MathJSON expression as the text of the language: the field it reads as, printed.
inline Result<std::string> mathjson_text(const JsonValue& expression) { return detail::mathjson_side(expression, ""); }

// A MathJSON equation, ["Equal", a, b], as the texts of its two sides.
inline Result<std::pair<std::string, std::string>> mathjson_equation(const JsonValue& equation) {
    const auto& a = equation.as_array();   // empty when it is not an array
    if (a.size() != 3 || !a[0].is_string() || a[0].as_string() != "Equal")
        return detail::mathjson_error("", "an equation is [\"Equal\", a, b]");
    auto l = detail::mathjson_side(a[1], "[1]");
    if (!l) return std::unexpected(l.error());
    auto r = detail::mathjson_side(a[2], "[2]");
    if (!r) return std::unexpected(r.error());
    return std::pair{std::move(*l), std::move(*r)};
}

namespace detail {

// a name a row defines: a symbol the language reads as a name, not one of its constants
inline Result<std::string> mathjson_name(const JsonValue& v, const std::string& path) {
    const JsonValue* sym = v.is_object() ? v.find("sym") : nullptr;
    const JsonValue& s = sym ? *sym : v;
    if (!s.is_string()) return mathjson_error(path, "a name is a symbol, such as x");
    auto t = mathjson_symbol(s.as_string(), path);
    if (!t) return std::unexpected(t.error());
    const auto names = t->free_symbols();
    if (names.size() != 1 || names.front() != s.as_string())
        return mathjson_error(path, std::format("'{}' is a constant, not a name to define", s.as_string()));
    return std::string(s.as_string());
}

} // namespace detail

} // namespace spatium::io::build
