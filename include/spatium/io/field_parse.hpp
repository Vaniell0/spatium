#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/io/field.hpp>
#  include <algorithm>
#  include <array>
#  include <cctype>
#  include <charconv>
#  include <format>
#  include <limits>
#  include <numbers>
#  include <optional>
#  include <string>
#  include <string_view>
#  include <utility>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io::build {

// ── Text read into a field ────────────────────────────────────────
//
// The inverse of `to_text`: a field written in the usual notation is read into the same ops, so
// a question typed by a person and an answer printed by the library are one language. Sums and
// products associate to the left, as `to_text` assumes when it leaves their parentheses out; a
// unary minus binds tighter than a product and looser than a power, so -x^2 is -(x * x); x^n with
// a whole n from 1 to 64 is the product of n factors of x, multiplied left to right, which is
// the same field as writing them out. The names given for the inputs (u and v by default) read
// as U and V, any other name as a symbol, `inf`, `nan` and `pi` as numbers; a minus in front of a
// number is part of the number (-0 is negative zero), in front of anything else it is a product
// with -1, which negates every double exactly. The functions are the field's own: sin, cos, log,
// exp and sqrt of one argument, min, max and less of two. A product after a number or a closing
// parenthesis may leave out its '*' (2x^2, (x - 1)(x + 2)), which `to_text` never does. A text
// that is not a field is refused at the column where it stops being one, with what was expected
// there.
//
// What `to_text` writes for a field reads back to the same values bit for bit, and to the same
// text again; when the field has no factor of exactly 1 and adds no negative constant (the two
// things the printer writes more simply), to the same ops.

// The reader itself, for a text in which a field is one part: it reads an expression from where
// it stands up to the first thing that cannot continue it, and a caller reads the words around
// it (a question's `solve ... for x in [0, 1]` is read this way).
class FieldReader {
public:
    // `reserved` are words of the surrounding text (a question's `for` and `in`), which are not names
    FieldReader(std::string_view text, std::string u, std::string v, std::vector<std::string> reserved = {})
        : s_(text), u_(std::move(u)), v_(std::move(v)), reserved_(std::move(reserved)) {}

    // A sum, from the current position up to the first thing that cannot continue it.
    Result<Field<double>> expression() {
        auto x = product();
        if (!x) return x;
        while (true) {
            const char c = peek();
            if (c != '+' && c != '-') return x;
            ++at_;
            auto y = product();
            if (!y) return y;
            *x = c == '+' ? std::move(*x) + *y : std::move(*x) - *y;
        }
    }

    // The next word, without taking it: a name of letters, digits and '_', or empty.
    std::string_view word() {
        skip();
        std::size_t end = at_;
        while (end < s_.size() && is_name_char(s_[end], end == at_)) ++end;
        return s_.substr(at_, end - at_);
    }

    bool take_word(std::string_view w) {
        if (word() != w) return false;
        at_ += w.size();
        return true;
    }

    bool take(char c) {
        if (peek() != c) return false;
        ++at_;
        return true;
    }

    bool done() { return peek() == '\0'; }

    std::unexpected<Error> fail(std::string_view expected) {
        skip();
        const std::string found = at_ < s_.size() ? std::format("'{}'", found_here()) : "the end";
        return std::unexpected(Error(ErrorCode::ParseError,
            std::format("parse: at column {}, expected {}, found {}", at_ + 1, expected, found)));
    }

private:
    static bool is_name_char(char c, bool first) {
        const auto b = static_cast<unsigned char>(c);
        return std::isalpha(b) || c == '_' || (!first && std::isdigit(b));
    }

    void skip() {
        while (at_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[at_]))) ++at_;
    }

    char peek() {
        skip();
        return at_ < s_.size() ? s_[at_] : '\0';
    }

    std::string_view found_here() {
        const std::string_view w = word();
        return w.empty() ? s_.substr(at_, 1) : w;
    }

    Result<Field<double>> product() {
        auto x = signed_power();
        if (!x) return x;
        while (true) {
            const char c = peek();
            const bool implied = c != '*' && c != '/' && follows_number_or_group() &&
                                 (c == '(' || (!word().empty() && std::ranges::find(reserved_, word()) == reserved_.end()));
            if (c != '*' && c != '/' && !implied) return x;
            if (!implied) ++at_;
            auto y = signed_power();
            if (!y) return y;
            *x = c == '/' ? std::move(*x) / *y : std::move(*x) * *y;
        }
    }

    // A product may be written without its '*' after a number or a closing parenthesis, the way
    // a person types one: 2x^2 is 2 * x^2 and (x - 1)(x + 2) a product of two. After a name it may
    // not -- u v is refused -- since a name followed by a name is more often a typo than a product.
    bool follows_number_or_group() const {
        std::size_t k = at_;
        while (k > 0 && std::isspace(static_cast<unsigned char>(s_[k - 1]))) --k;
        if (k == 0) return false;
        const char last = s_[k - 1];
        return last == ')' || std::isdigit(static_cast<unsigned char>(last)) || last == '.';
    }

    Result<Field<double>> signed_power() {
        if (!take('-')) return power();
        // -2 is a number, -2^2 is -(2 * 2)
        const std::size_t from = at_;
        if (const auto n = number(); n && peek() != '^') return Field<double>(-*n);
        at_ = from;
        auto x = signed_power();
        if (!x) return x;
        return Field<double>(-1.0) * *x;
    }

    Result<Field<double>> power() {
        auto x = atom();
        if (!x || !take('^')) return x;
        skip();
        const std::size_t from = at_;
        int n = 0;
        const auto [end, ec] = std::from_chars(s_.data() + at_, s_.data() + s_.size(), n);
        const std::size_t to = static_cast<std::size_t>(end - s_.data());
        if (ec != std::errc{} || n < 1 || n > 64 || (to < s_.size() && (s_[to] == '.' || s_[to] == 'e'))) {
            at_ = from;
            return fail("a whole power from 1 to 64 (x^a for another a is exp(a * log(x)))");
        }
        at_ = to;
        Field<double> out = *x;
        for (int i = 1; i < n; ++i) out = std::move(out) * *x;
        return out;
    }

    // A literal: digits in the shortest form `to_text` writes them, or inf, or nan.
    std::optional<double> number() {
        skip();
        if (take_word("inf")) return std::numeric_limits<double>::infinity();
        if (take_word("nan")) return std::numeric_limits<double>::quiet_NaN();
        if (take_word("pi")) return std::numbers::pi;
        if (at_ >= s_.size() || !(std::isdigit(static_cast<unsigned char>(s_[at_])) || s_[at_] == '.'))
            return std::nullopt;
        double x = 0;
        const auto [end, ec] = std::from_chars(s_.data() + at_, s_.data() + s_.size(), x);
        if (ec != std::errc{}) return std::nullopt;
        at_ = static_cast<std::size_t>(end - s_.data());
        return x;
    }

    Result<Field<double>> atom() {
        if (const auto n = number()) return Field<double>(*n);
        if (take('(')) {
            auto x = expression();
            if (!x) return x;
            if (!take(')')) return fail("')'");
            return x;
        }
        const std::string_view name = word();
        if (name.empty() || std::ranges::find(reserved_, name) != reserved_.end())
            return fail("a number, a name, a function or '('");
        const std::size_t name_at = at_;
        at_ += name.size();
        if (peek() == '(') return call(name, name_at);
        if (std::ranges::find(unary_, name) != unary_.end() || std::ranges::find(binary_, name) != binary_.end()) {
            at_ = name_at;
            return fail(std::format("a name, but '{}' is a function: write {}(...)", name, name));
        }
        if (name == u_) return Field<double>::u();
        if (name == v_) return Field<double>::v();
        return Field<double>::symbol(std::string(name));
    }

    Result<Field<double>> call(std::string_view name, std::size_t name_at) {
        const bool one = std::ranges::find(unary_, name) != unary_.end();
        if (!one && std::ranges::find(binary_, name) == binary_.end()) {
            at_ = name_at;
            return fail("sin, cos, log, exp, sqrt, min, max or less before '('");
        }
        take('(');
        auto a = expression();
        if (!a) return a;
        if (one) {
            if (!take(')')) return fail(std::format("')' after the one argument of {}", name));
            if (name == "sin") return sin(std::move(*a));
            if (name == "cos") return cos(std::move(*a));
            if (name == "log") return log(std::move(*a));
            if (name == "exp") return exp(std::move(*a));
            return sqrt(std::move(*a));
        }
        if (!take(',')) return fail(std::format("',' between the two arguments of {}", name));
        auto b = expression();
        if (!b) return b;
        if (!take(')')) return fail(std::format("')' after the two arguments of {}", name));
        if (name == "min") return min(std::move(*a), *b);
        if (name == "max") return max(std::move(*a), *b);
        return less(std::move(*a), *b);
    }

    static constexpr std::array<std::string_view, 5> unary_{"sin", "cos", "log", "exp", "sqrt"};
    static constexpr std::array<std::string_view, 3> binary_{"min", "max", "less"};

    std::string_view s_;
    std::size_t at_ = 0;
    std::string u_, v_;
    std::vector<std::string> reserved_;
};

inline Result<Field<double>> parse(std::string_view text, std::string u = "u", std::string v = "v") {
    FieldReader r(text, std::move(u), std::move(v));
    auto f = r.expression();
    if (!f) return f;
    if (!r.done()) return r.fail("an operator or the end");
    return f;
}

} // namespace spatium::io::build
