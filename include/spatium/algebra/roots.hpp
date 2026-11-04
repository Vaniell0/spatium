#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/core/interval.hpp>
#  include <spatium/algebra/dual.hpp>
#  include <algorithm>
#  include <cmath>
#  include <limits>
#  include <string>
#  include <vector>
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

// ── Every root on a box, proved ───────────────────────────────
//
// The real roots of f on a box, each in an interval proved to hold exactly one, and the pieces
// where one may still be. `prove_on` asks of each piece whether every root on it is accounted for:
// yes where f's interval leaves zero out, and yes where Newton's interval step
//     N(X) = m - f(m) / f'(X),   m the middle of X,
// with f' enclosed through a dual number of intervals, misses X or lands inside it. Every root of
// X is in N(X), so a step that misses X leaves none, and one inside X holds exactly one there:
// f' keeps one sign on X, so f crosses zero once (Moore's theorem). A root on a cut is the end of
// two pieces and inside neither, so a step that contracts but sticks out is retried on itself
// widened (Rump's epsilon-inflation): every root of X is in it too, and one proved there accounts
// for X's. A proved root is narrowed by the same step until it stops shrinking. A double root
// stays undecided -- f' holds zero there, and no step tells one root from two or none -- and
// adjacent undecided pieces are joined. A piece of one magnitude where f at both its ends is past
// what T holds is left uncut: near the end of the doubles x^2 overflows, and no cut there tells more.
//
// f is called with an Interval<T> and with a Dual<Interval<T>>: a generic lambda, or a lowered field.

template<Scalar T>
struct Roots {
    std::vector<Interval<T>> roots;       // each holds exactly one root of f, in increasing order
    std::vector<Interval<T>> undecided;   // where a root may still be, neither excluded nor proved
    Interval<T> box;                      // what was searched
    std::size_t evaluations = 0;          // of f, as an interval or a dual number
    bool complete() const { return undecided.empty(); }
};

template<Scalar T, class F>
Roots<T> roots_on(F&& f, const Interval<T>& box, const Budget& budget = {}) {
    using I = Interval<T>;
    using D = Dual<I>;
    using std::abs;
    Roots<T> out;
    out.box = box;
    const auto misses_zero = [](const I& y) { return y.is_empty() || T(0) < y.lo || y.hi < T(0); };
    const auto within = [](const I& a, const I& b) { return b.lo <= a.lo && a.hi <= b.hi; };
    const auto at = [&](const I& x) {
        ++out.evaluations;
        return f(D::variable(x));
    };
    // Newton's step from X, given f and f' on X; empty where f' may vanish on X
    const auto step = [&](const I& x, const D& fx) {
        if (!misses_zero(fx.deriv) || fx.deriv.is_empty()) return I::empty();
        const T m = x.mid();
        ++out.evaluations;
        return I{m} - f(I{m}) / fx.deriv;
    };
    const auto unbounded = [](const I& y) { return y.is_empty() || !(abs(y.lo) <= std::numeric_limits<T>::max()) ||
                                                    !(abs(y.hi) <= std::numeric_limits<T>::max()); };
    // one magnitude, as a cut at its middle would split: there f past T at both ends is past it throughout
    const auto one_scale = [](const I& x) {
        return (T(0) < x.lo && x.hi <= T(16) * x.lo) || (x.hi < T(0) && T(16) * x.hi <= x.lo);
    };
    const auto undecided = [&](const I& x) {
        if (!one_scale(x)) return Decided{};
        out.evaluations += 2;
        return Decided{Truth::Undecided, unbounded(f(I{x.lo})) && unbounded(f(I{x.hi}))};
    };
    const auto statement = [&](const Box<T, 1>& b) {
        const I x = b[0];
        // first f alone: a derivative can guess a branch the value does not (abs at its kink)
        ++out.evaluations;
        if (misses_zero(f(x))) return Decided{Truth::Yes};
        const Undecided watch;
        const D fx = at(x);
        I n = step(x, fx);
        if (n.is_empty()) return undecided(x);
        if (n.hi < x.lo || x.hi < n.lo) return Decided{Truth::Yes};
        I y = n;
        if (!within(n, x)) {
            if (!(n.width() < x.width())) return undecided(x);
            bool held = false;
            for (int k = 0; k < 4 && !held; ++k) {
                const T grow = std::max(T(y.width() / T(8)),
                                        T(std::numeric_limits<T>::epsilon() * std::max({T(1), T(abs(y.lo)), T(abs(y.hi))})));
                const I wider{T(y.lo - grow), T(y.hi + grow)};
                const I m = step(wider, at(wider));
                if (m.is_empty()) return undecided(x);
                held = within(m, wider);
                y = held ? m : wider;
            }
            if (!held) return undecided(x);
        }
        // y holds exactly one root and every root of X: narrowed until the step stops shrinking it
        for (int k = 0; k < 64; ++k) {
            const I m = step(y, at(y));
            if (m.is_empty()) break;
            const I z{std::max(m.lo, y.lo), std::min(m.hi, y.hi)};
            if (!(z.width() < y.width())) break;
            y = z;
        }
        if (watch.guessed()) return Decided{};
        out.roots.push_back(y);
        return Decided{Truth::Yes};
    };
    const Proof<T, 1> proof = prove_on(Box<T, 1>{box}, statement, budget);
    const auto by_lo = [](const I& a, const I& b) { return a.lo < b.lo; };
    // a root on a cut, found from both sides: one, where f' keeps its sign across the two
    std::ranges::sort(out.roots, by_lo);
    std::vector<I> roots;
    for (const I& r : out.roots) {
        if (!roots.empty() && !(roots.back().hi < r.lo)) {
            const I hull{roots.back().lo, std::max(roots.back().hi, r.hi)};
            if (misses_zero(at(hull).deriv)) {
                roots.back() = I{std::max(roots.back().lo, r.lo), std::min(roots.back().hi, r.hi)};
                continue;
            }
        }
        roots.push_back(r);
    }
    out.roots = std::move(roots);
    for (const auto& [b, t] : proof.pieces)
        if (t == Truth::Undecided) out.undecided.push_back(b[0]);
    std::ranges::sort(out.undecided, by_lo);
    std::vector<I> joined;
    for (const I& u : out.undecided) {
        if (!joined.empty() && !(joined.back().hi < u.lo)) joined.back().hi = std::max(joined.back().hi, u.hi);
        else joined.push_back(u);
    }
    out.undecided = std::move(joined);
    return out;
}

} // inline namespace algebra
} // namespace spatium
