#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <format>
#  include <ostream>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// ── operator<< for Vec ─────────────────────────────────────────

template<Scalar T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const Vec<T, N>& v) {
    os << '(';
    for (std::size_t i = 0; i < N; ++i) {
        if (i > 0) os << ", ";
        os << v[i];
    }
    os << ')';
    return os;
}

// ── std::format support for Vec ────────────────────────────────

} // namespace algebra
} // namespace spatium

template<spatium::Scalar T, std::size_t N>
struct std::formatter<spatium::Vec<T, N>> : std::formatter<T> {
    auto format(const spatium::Vec<T, N>& v, auto& ctx) const {
        auto out = ctx.out();
        *out++ = '(';
        for (std::size_t i = 0; i < N; ++i) {
            if (i > 0) {
                *out++ = ',';
                *out++ = ' ';
            }
            ctx.advance_to(out);
            out = std::formatter<T>::format(v[i], ctx);
        }
        *out++ = ')';
        return out;
    }
};
