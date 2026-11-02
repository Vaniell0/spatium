#pragma once

// A decimal floating-point number of 16 digits, IEEE 754 decimal64: the arithmetic in which ten
// times 0.1 is 1 and a sum of money rounds the way it is written. It is GCC's
// `std::decimal::decimal64` (ISO/IEC TR 24733, `<decimal/decimal>`), which the C++ standard
// itself does not have; where the compiler has no such header -- Clang, and so Emscripten's
// WebAssembly -- this header compiles to nothing and `SPATIUM_HAS_DECIMAL` is 0, the same way
// core/precision.hpp does for Boost.
//
// The technical report gives the four operations and order, and no function: there is no square
// root or sine here, and an evaluation that needs one is refused by the op's name before it runs
// (`exact_op` in io/field_pod.hpp). A number read from a double is read through its shortest
// decimal form, so `Decimal(0.1)` is the 0.1 a person typed.

#if __has_include(<decimal/decimal>) && defined(__GNUC__) && !defined(__clang__)
#  define SPATIUM_HAS_DECIMAL 1
#else
#  define SPATIUM_HAS_DECIMAL 0
#endif

#if SPATIUM_HAS_DECIMAL

#ifndef SPATIUM_BUILDING_MODULE
#  include <charconv>
#  include <compare>
#  include <cstdint>
#  include <decimal/decimal>
#  include <optional>
#  include <string>
#  include <string_view>
#endif

#include <spatium/_export_macro.hpp>

SPATIUM_EXPORT namespace spatium {

class Decimal {
public:
    using Raw = std::decimal::decimal64;

    Decimal() = default;
    Decimal(int n) : x_(n) {}   // NOLINT: an integer is a decimal
    explicit Decimal(Raw x) : x_(x) {}
    // The double's shortest decimal form; NaN when it is not finite.
    explicit Decimal(double x) {
        char buf[32];
        const auto [end, ec] = std::to_chars(buf, buf + sizeof buf, x);
        const auto r = ec == std::errc{} ? from_decimal(std::string_view(buf, end)) : std::nullopt;
        x_ = r ? r->x_ : nan();
    }

    // Digits with an optional sign, point and exponent, at most 16 significant digits.
    static std::optional<Decimal> from_decimal(std::string_view s) {
        bool negative = false;
        std::size_t i = 0;
        if (i < s.size() && (s[i] == '-' || s[i] == '+')) negative = s[i++] == '-';
        long long digits = 0;
        int scale = 0, count = 0;
        bool point = false, any = false;
        for (; i < s.size() && ((s[i] >= '0' && s[i] <= '9') || s[i] == '.'); ++i) {
            if (s[i] == '.') {
                if (point) return std::nullopt;
                point = true;
                continue;
            }
            any = true;
            if (point) --scale;
            if (count == 0 && s[i] == '0') continue;   // leading zeros carry no digit
            if (++count > 16) return std::nullopt;
            digits = digits * 10 + (s[i] - '0');
        }
        if (!any) return std::nullopt;
        if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
            int e = 0;
            const std::size_t from = i + 1 + (i + 1 < s.size() && s[i + 1] == '+');
            const auto [end, ec] = std::from_chars(s.data() + from, s.data() + s.size(), e);
            if (ec != std::errc{} || end != s.data() + s.size()) return std::nullopt;
            scale += e;
            i = s.size();
        }
        if (i != s.size()) return std::nullopt;
        Raw x = std::decimal::make_decimal64(digits, scale);
        return Decimal(negative ? -x : x);
    }

    Raw raw() const { return x_; }
    bool finite() const { return x_ == x_ && x_ - x_ == x_ - x_; }
    double to_double() const { return std::decimal::decimal64_to_double(x_); }

    // The digits exactly as held, with no exponent for ordinary sizes: `5.25`, `-0.001`, `1E+20`.
    std::string text() const {
        if (x_ != x_) return "nan";
        if (!finite()) return x_ < Raw(0) ? "-inf" : "inf";
        Raw a = x_ < Raw(0) ? -x_ : x_;
        // the coefficient and its power of ten: past 16 digits a tenth is exact, below a whole
        // number ten times is, until the coefficient is whole
        int e = 0, scale = 0;
        while (a >= Raw(10000000000000000LL)) a /= Raw(10), ++e;
        while (Raw(std::decimal::decimal_to_long_long(a)) != a && scale < 398) a *= Raw(10), ++scale;
        long long c = std::decimal::decimal_to_long_long(a);
        while (e > 0 && c % 10 == 0 && c != 0) c /= 10, ++e;
        std::string digits = std::to_string(c), out;
        if (e > 0) {
            out = digits + "E+" + std::to_string(e);
        } else if (scale > 0) {
            if (static_cast<int>(digits.size()) <= scale) digits.insert(0, scale - digits.size() + 1, '0');
            out = digits.substr(0, digits.size() - scale) + "." + digits.substr(digits.size() - scale);
            while (out.back() == '0') out.pop_back();
            if (out.back() == '.') out.pop_back();
        } else {
            out = digits;
        }
        return (x_ < Raw(0) ? "-" : "") + out;
    }

    friend Decimal operator+(Decimal a, Decimal b) { return Decimal(a.x_ + b.x_); }
    friend Decimal operator-(Decimal a, Decimal b) { return Decimal(a.x_ - b.x_); }
    friend Decimal operator*(Decimal a, Decimal b) { return Decimal(a.x_ * b.x_); }
    friend Decimal operator/(Decimal a, Decimal b) { return Decimal(a.x_ / b.x_); }
    friend Decimal operator-(Decimal a) { return Decimal(-a.x_); }
    Decimal& operator+=(Decimal b) { return *this = *this + b; }
    Decimal& operator-=(Decimal b) { return *this = *this - b; }
    Decimal& operator*=(Decimal b) { return *this = *this * b; }
    Decimal& operator/=(Decimal b) { return *this = *this / b; }
    friend Decimal abs(Decimal a) { return a.x_ < Raw(0) ? -a : a; }

    friend bool operator==(Decimal a, Decimal b) { return a.x_ == b.x_; }
    friend std::partial_ordering operator<=>(Decimal a, Decimal b) {
        if (a.x_ < b.x_) return std::partial_ordering::less;
        if (b.x_ < a.x_) return std::partial_ordering::greater;
        if (a.x_ == b.x_) return std::partial_ordering::equivalent;
        return std::partial_ordering::unordered;
    }

private:
    static Raw nan() { return Raw(0) / Raw(0); }
    Raw x_{0};
};

} // namespace spatium

#endif // SPATIUM_HAS_DECIMAL
