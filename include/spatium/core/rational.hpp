#pragma once

// A fraction of two 64-bit integers, kept in lowest terms with a positive denominator: the
// arithmetic in which 1/3 + 1/6 is 1/2 and not 0.49999999999999994.
//
// It is closed under + - * / and order, and under nothing else: there is no square root or sine
// of a fraction in general, so it has none, and an evaluation that needs one is refused by the
// op's name before it runs (`exact_op` in io/field_pod.hpp). A result that does not fit -- a
// numerator or denominator past 64 bits, or a division by zero -- is not wrapped or rounded: the
// fraction is spoiled (`ok()` false) and stays so through everything computed from it, the way a
// NaN does, so one check at the end tells whether the whole computation was exact.
//
// A number read from text is the fraction its decimal digits write (`from_decimal("0.1")` is
// 1/10), and a double is read through its shortest decimal form (`Rational(0.1)` is 1/10 too):
// the literal a person typed, not the binary fraction nearest to it.

#ifndef SPATIUM_BUILDING_MODULE
#  include <cctype>
#  include <charconv>
#  include <cmath>
#  include <compare>
#  include <concepts>
#  include <cstdint>
#  include <optional>
#  include <string>
#  include <string_view>
#  include <utility>
#endif

#include <spatium/_export_macro.hpp>

SPATIUM_EXPORT namespace spatium {

class Rational {
public:
    constexpr Rational() = default;
    template<std::integral I>
    constexpr Rational(I n) : num_(static_cast<long long>(n)) {}   // NOLINT: an integer is a fraction
    Rational(long long n, long long d) { *this = reduce(n, d); }
    // The double's shortest decimal form, as a fraction; spoiled when it is not finite or does
    // not fit.
    explicit Rational(double x) {
        char buf[32];
        const auto [end, ec] = std::to_chars(buf, buf + sizeof buf, x);
        const auto r = ec == std::errc{} && std::isfinite(x) ? from_decimal(std::string_view(buf, end)) : std::nullopt;
        *this = r ? *r : spoiled();
    }

    // Digits with an optional sign, point and exponent (`-12.5e-3`); nothing else.
    static std::optional<Rational> from_decimal(std::string_view s) {
        bool negative = false;
        std::size_t i = 0;
        if (i < s.size() && (s[i] == '-' || s[i] == '+')) negative = s[i++] == '-';
        __int128 digits = 0;
        int scale = 0, count = 0;
        bool point = false;
        for (; i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i])) || s[i] == '.'); ++i) {
            if (s[i] == '.') {
                if (point) return std::nullopt;
                point = true;
                continue;
            }
            if (digits > (__int128(1) << 100)) return std::nullopt;
            digits = digits * 10 + (s[i] - '0');
            ++count;
            if (point) --scale;
        }
        if (count == 0) return std::nullopt;
        if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
            int e = 0;
            const std::size_t from = i + 1 + (i + 1 < s.size() && s[i + 1] == '+');
            const auto [end, ec] = std::from_chars(s.data() + from, s.data() + s.size(), e);
            if (ec != std::errc{} || end != s.data() + s.size()) return std::nullopt;
            scale += e;
            i = s.size();
        }
        if (i != s.size()) return std::nullopt;
        __int128 den = 1;
        for (; scale > 0; --scale) {
            digits *= 10;
            if (digits > max_) return std::nullopt;
        }
        for (; scale < 0; ++scale) {
            den *= 10;
            if (den > (__int128(1) << 100)) return std::nullopt;
        }
        const Rational r = reduce(negative ? -digits : digits, den);
        if (!r.ok()) return std::nullopt;
        return r;
    }

    constexpr long long num() const { return num_; }
    constexpr long long den() const { return den_; }
    // False once a step did not fit or divided by zero.
    constexpr bool ok() const { return den_ != 0; }

    double to_double() const { return ok() ? static_cast<double>(num_) / static_cast<double>(den_) : std::nan(""); }
    // `21/4`, `-3`, or `spoiled`
    std::string text() const {
        if (!ok()) return "spoiled";
        return den_ == 1 ? std::to_string(num_) : std::to_string(num_) + "/" + std::to_string(den_);
    }

    friend Rational operator+(Rational a, Rational b) {
        if (!a.ok() || !b.ok()) return spoiled();
        return reduce(__int128(a.num_) * b.den_ + __int128(b.num_) * a.den_, __int128(a.den_) * b.den_);
    }
    friend Rational operator-(Rational a, Rational b) { return a + (-b); }
    friend Rational operator*(Rational a, Rational b) {
        if (!a.ok() || !b.ok()) return spoiled();
        return reduce(__int128(a.num_) * b.num_, __int128(a.den_) * b.den_);
    }
    friend Rational operator/(Rational a, Rational b) {
        if (!a.ok() || !b.ok() || b.num_ == 0) return spoiled();
        return reduce(__int128(a.num_) * b.den_, __int128(a.den_) * b.num_);
    }
    friend Rational operator-(Rational a) {
        if (!a.ok() || a.num_ == min_) return spoiled();
        a.num_ = -a.num_;
        return a;
    }
    Rational& operator+=(Rational b) { return *this = *this + b; }
    Rational& operator-=(Rational b) { return *this = *this - b; }
    Rational& operator*=(Rational b) { return *this = *this * b; }
    Rational& operator/=(Rational b) { return *this = *this / b; }
    friend Rational abs(Rational a) { return a.num_ < 0 ? -a : a; }

    // Spoiled fractions are equal to each other and below every fraction, so the order stays total.
    friend constexpr bool operator==(const Rational&, const Rational&) = default;
    friend std::strong_ordering operator<=>(Rational a, Rational b) {
        if (!a.ok() || !b.ok()) return a.ok() <=> b.ok();
        return __int128(a.num_) * b.den_ <=> __int128(b.num_) * a.den_;
    }

private:
    static constexpr __int128 max_ = INT64_MAX;
    static constexpr long long min_ = INT64_MIN;
    static constexpr Rational spoiled() {
        Rational r;
        r.num_ = 0, r.den_ = 0;
        return r;
    }
    static Rational reduce(__int128 n, __int128 d) {
        if (d == 0) return spoiled();
        if (d < 0) n = -n, d = -d;
        __int128 a = n < 0 ? -n : n, b = d;
        while (b != 0) a %= b, std::swap(a, b);
        if (a > 1) n /= a, d /= a;
        if (n > max_ || n < -max_ || d > max_) return spoiled();
        Rational r;
        r.num_ = static_cast<long long>(n), r.den_ = static_cast<long long>(d);
        return r;
    }

    long long num_ = 0, den_ = 1;
};

} // namespace spatium
