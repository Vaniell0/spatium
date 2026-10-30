#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/error.hpp>
#  include <cmath>
#  include <string>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// ── A root on a bracket ───────────────────────────────────────
//
// A zero of a continuous f on [lo, hi] where f changes sign, by bisection to the last bit: the
// bracket is halved until its midpoint is one of its ends, so the answer is one of two adjacent
// numbers of T and the root lies between them. Bisection rather than anything faster because it
// is the one method whose answer is guaranteed by the sign change alone -- no derivative, no
// smoothness, no starting guess -- and its cost is about one evaluation per bit of the bracket.
//
// No sign change is not a root: an even-multiplicity zero, or none, or a bracket around a pole,
// are told apart by nothing a sign can see, so they are refused rather than answered. An end at
// which f is exactly zero is that root. A value of f that is not finite is refused where it is
// met: a NaN has no sign to bisect on.

template<Scalar T, typename F>
Result<T> find_root(F&& f, T lo, T hi) {
    using std::isfinite;
    T flo = f(lo), fhi = f(hi);
    if (!isfinite(flo) || !isfinite(fhi))
        return std::unexpected(Error(ErrorCode::InvalidArgument,
            "find_root: f is not finite at an end of the bracket"));
    if (flo == T{0}) return lo;
    if (fhi == T{0}) return hi;
    if ((flo < T{0}) == (fhi < T{0}))
        return std::unexpected(Error(ErrorCode::InvalidArgument,
            "find_root: f does not change sign on the bracket, so it holds no root a sign can find"));
    // Enough halvings for any width down to one unit in the last place, denormals included.
    for (int i = 0; i < 4096; ++i) {
        const T mid = lo + (hi - lo) / T{2};
        if (mid == lo || mid == hi) break;
        const T fm = f(mid);
        if (!isfinite(fm))
            return std::unexpected(Error(ErrorCode::InvalidArgument,
                "find_root: f is not finite inside the bracket"));
        if (fm == T{0}) return mid;
        if ((fm < T{0}) == (flo < T{0})) lo = mid, flo = fm;
        else hi = mid, fhi = fm;
    }
    using std::abs;
    return abs(flo) <= abs(fhi) ? lo : hi;
}

} // inline namespace algebra
} // namespace spatium
