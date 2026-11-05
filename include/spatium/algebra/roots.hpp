#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/core/interval.hpp>
#  include <spatium/algebra/dual.hpp>
#  include <spatium/algebra/linear_solve.hpp>
#  include <array>
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
    std::size_t pieces = 0;               // the box was cut into, each decided or left undecided
    std::size_t excluded = 0;             // of them proved to hold no root
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
        if (misses_zero(f(x))) return ++out.excluded, Decided{Truth::Yes};
        const Undecided watch;
        const D fx = at(x);
        I n = step(x, fx);
        if (n.is_empty()) return undecided(x);
        if (n.hi < x.lo || x.hi < n.lo) return ++out.excluded, Decided{Truth::Yes};
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
    out.pieces = proof.pieces.size();
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

// ── Every root of a system on a box, proved ─────────────────
//
// The same in N unknowns: f maps N intervals to N, and a root is a point where all N vanish.
// Krawczyk's operator stands for Newton's step,
//     K(X) = m - Y f(m) + (I - Y J(X)) (X - m),
// m the middle of X, J(X) the Jacobian on X (N passes of dual numbers, one per unknown), Y the
// inverse of J's middle. Every root of X is in K(X): K(X) missing X leaves none, and K(X) inside X
// holds exactly one (Krawczyk; Moore and Kioustelidis). A piece where some f_i leaves zero out
// holds none. A root on a cut is retried on K(X) widened, as on a segment, and a root found from
// two pieces is one where their boxes meet. A root where J is singular stays undecided.

template<Scalar T, std::size_t N>
struct RootsN {
    std::vector<Box<T, N>> roots;       // each holds exactly one root of f
    std::vector<Box<T, N>> undecided;   // where a root may still be
    Box<T, N> box;
    std::size_t evaluations = 0;        // of f, in intervals or in dual numbers
    std::size_t pieces = 0;             // the box was cut into, each decided or left undecided
    std::size_t excluded = 0;           // of them proved to hold no root
    bool complete() const { return undecided.empty(); }
};

template<Scalar T, std::size_t N, class F>
RootsN<T, N> roots_on(F&& f, const Box<T, N>& box, const Budget& budget = {}) {
    using I = Interval<T>;
    using D = Dual<I>;
    RootsN<T, N> out;
    out.box = box;
    const auto misses_zero = [](const I& y) { return y.is_empty() || T(0) < y.lo || y.hi < T(0); };
    const auto inside = [](const Box<T, N>& a, const Box<T, N>& b) {
        for (std::size_t i = 0; i < N; ++i)
            if (!(b[i].lo < a[i].lo && a[i].hi < b[i].hi)) return false;
        return true;
    };
    const auto apart = [](const Box<T, N>& a, const Box<T, N>& b) {
        for (std::size_t i = 0; i < N; ++i)
            if (a[i].hi < b[i].lo || b[i].hi < a[i].lo) return true;
        return false;
    };
    const auto values = [&](const Box<T, N>& x) {
        ++out.evaluations;
        return f(x);
    };
    // Krawczyk's K(X), or none where J's middle has no inverse or f has no value
    const auto krawczyk = [&](const Box<T, N>& x) -> std::optional<Box<T, N>> {
        Box<T, N> m;
        for (std::size_t i = 0; i < N; ++i) m[i] = I{x[i].mid()};
        const Box<T, N> fm = values(m);
        std::array<std::array<I, N>, N> J;   // J[i][j] = d f_i / d x_j on X
        for (std::size_t j = 0; j < N; ++j) {
            std::array<D, N> xd;
            for (std::size_t k = 0; k < N; ++k) xd[k] = D{x[k], I{k == j ? T(1) : T(0)}};
            ++out.evaluations;
            const auto fd = f(xd);
            for (std::size_t i = 0; i < N; ++i) J[i][j] = fd[i].deriv;
        }
        Matrix<T, N, N> c{};
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) {
                if (J[i][j].is_empty() || !std::isfinite(J[i][j].lo) || !std::isfinite(J[i][j].hi)) return std::nullopt;
                c(i, j) = J[i][j].mid();
            }
        for (const I& v : fm)
            if (v.is_empty()) return std::nullopt;
        const auto y = invert(c);
        if (!y) return std::nullopt;
        Box<T, N> k;
        for (std::size_t i = 0; i < N; ++i) {
            I sum = m[i];
            for (std::size_t j = 0; j < N; ++j) sum = sum - I{(*y)(i, j)} * fm[j];
            for (std::size_t j = 0; j < N; ++j) {
                I yj{T(i == j ? 1 : 0)};
                for (std::size_t l = 0; l < N; ++l) yj = yj - I{(*y)(i, l)} * J[l][j];
                sum = sum + yj * (x[j] - m[j]);
            }
            k[i] = sum;
        }
        return k;
    };
    const auto statement = [&](const Box<T, N>& x) {
        const Box<T, N> fx = values(x);
        for (const I& v : fx)
            if (misses_zero(v)) return ++out.excluded, Decided{Truth::Yes};
        const Undecided watch;
        const auto k = krawczyk(x);
        if (!k) return Decided{};
        if (apart(*k, x)) return ++out.excluded, Decided{Truth::Yes};
        Box<T, N> y = *k;
        if (!inside(*k, x)) {
            // a root on a cut: K(X) widened, which holds every root of X, proved to hold one
            bool held = false;
            for (int r = 0; r < 4 && !held; ++r) {
                Box<T, N> wider = y;
                for (std::size_t i = 0; i < N; ++i) {
                    const T grow = std::max(T(y[i].width() / T(8)), T(std::numeric_limits<T>::epsilon() *
                                                                      std::max({T(1), T(std::abs(y[i].lo)), T(std::abs(y[i].hi))})));
                    wider[i] = I{T(y[i].lo - grow), T(y[i].hi + grow)};
                }
                const auto kw = krawczyk(wider);
                if (!kw) return Decided{};
                held = inside(*kw, wider);
                y = held ? *kw : wider;
            }
            // K(X) holds X's roots and lies in the box widened from it, whose roots are all in the
            // one K proved there: X's are that one
            if (!held) return Decided{};
        }
        // narrowed while K shrinks it
        for (int r = 0; r < 64; ++r) {
            const auto n = krawczyk(y);
            if (!n) break;
            Box<T, N> z;
            bool smaller = false;
            for (std::size_t i = 0; i < N; ++i) {
                z[i] = I{std::max((*n)[i].lo, y[i].lo), std::min((*n)[i].hi, y[i].hi)};
                smaller |= z[i].width() < y[i].width();
            }
            if (!smaller) break;
            y = z;
        }
        if (watch.guessed()) return Decided{};
        out.roots.push_back(y);
        return Decided{Truth::Yes};
    };
    const Proof<T, N> proof = prove_on(box, statement, budget);
    out.pieces = proof.pieces.size();
    // a root found from two pieces is one where their boxes meet
    std::vector<Box<T, N>> roots;
    for (const auto& r : out.roots) {
        const auto same = std::ranges::find_if(roots, [&](const Box<T, N>& q) { return !apart(q, r); });
        if (same == roots.end()) roots.push_back(r);
        else
            for (std::size_t i = 0; i < N; ++i) (*same)[i] = I{std::max((*same)[i].lo, r[i].lo), std::min((*same)[i].hi, r[i].hi)};
    }
    std::ranges::sort(roots, {}, [](const Box<T, N>& b) { return b[0].lo; });
    out.roots = std::move(roots);
    for (const auto& [b, t] : proof.pieces)
        if (t == Truth::Undecided) out.undecided.push_back(b);
    return out;
}

} // inline namespace algebra
} // namespace spatium
