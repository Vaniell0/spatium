#pragma once

// A fraction of two 64-bit integers, kept in lowest terms with a positive denominator: the
// arithmetic in which 1/3 + 1/6 is 1/2 and not 0.49999999999999994.
//
// It is closed under + - * / and order, and under nothing else: there is no sine of a fraction
// in general, and a square root only where the fraction is a square (`sqrt`, spoiled otherwise),
// so an evaluation that needs one is refused by the op's name before it runs (`exact_op` in
// io/field_pod.hpp). A result that does not fit -- a
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
    // The decimal it is, where that decimal ends -- a denominator of twos and fives only: 391/4 is
    // 97.75, 3/10 is 0.3 -- worked out in integers, never through a double. None for 1/3.
    std::optional<std::string> decimal() const {
        if (!ok()) return std::nullopt;
        long long d = den_;
        int twos = 0, fives = 0;
        while (d % 2 == 0) d /= 2, ++twos;
        while (d % 5 == 0) d /= 5, ++fives;
        if (d != 1) return std::nullopt;
        const int places = twos > fives ? twos : fives;
        if (places > 19) return std::nullopt;   // past 128 bits once scaled; such a decimal is no one's answer
        __int128 scaled = num_ < 0 ? -__int128(num_) : __int128(num_);
        for (int k = 0; k < places; ++k) scaled *= 10;
        scaled /= den_;
        std::string digits;
        do digits.insert(digits.begin(), static_cast<char>('0' + static_cast<int>(scaled % 10))), scaled /= 10;
        while (scaled > 0);
        if (places > 0) {
            if (digits.size() <= static_cast<std::size_t>(places)) digits.insert(0, static_cast<std::size_t>(places) + 1 - digits.size(), '0');
            digits.insert(digits.size() - static_cast<std::size_t>(places), ".");
        }
        return (num_ < 0 ? "-" : "") + digits;
    }

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
    // The square root where it is a fraction -- of 9/4 it is 3/2 -- and spoiled where it is not
    // (of 2, or of a negative): a fraction has no other. A solver written for any scalar
    // (`solve_quadratic_steps`) finds it by ADL and so gives exact roots where they are fractions.
    friend Rational sqrt(Rational a) {
        if (!a.ok() || a.num_ < 0) return spoiled();
        const auto root = [](long long n) -> long long {
            auto r = static_cast<long long>(std::sqrt(static_cast<double>(n)));
            while (r > 0 && r * r > n) --r;
            while (__int128(r + 1) * (r + 1) <= n) ++r;
            return __int128(r) * r == n ? r : -1;
        };
        const long long n = root(a.num_), d = root(a.den_);
        return n < 0 || d < 0 ? spoiled() : Rational(n, d);
    }

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
