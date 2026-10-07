#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/quadrature.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/dimension.hpp>
#  include <algorithm>
#  include <cmath>
#  include <cstddef>
#  include <vector>
#endif

// L^2([a, b]): the space of square-integrable functions on an interval -- the first
// space of this library with no finite basis. Its dimension is `Dimension::infinite()`,
// and that is all the library needs to know of it: `Euclidean`'s algorithms read a
// dimension, an algorithm written against an inner product does not.
//
// A point is a truncated Legendre series on the interval, f(x) ~ sum c_n P_n(t) with
// t = (2x - a - b)/(b - a). The Legendre polynomials are orthogonal for the plain L^2 weight,
// so the inner product is Parseval's identity on the coefficients,
//
//   <f, g> = (b - a)/2  sum  2/(2n + 1)  c_n d_n ,
//
// with no quadrature -- and a quadrature of f g is an independent path to the same number,
// which `tests/test_l2.cpp` holds it to. A series carries an estimate of what it neglected:
// the L^2 mass of its last three coefficients, which for a function the series suits is
// the size of the whole tail and for one it does not suit (a kink, |x|) is large, so a
// projection says when it cannot be trusted. It is an estimate, not a bound.
//
// Built from a function by `L2Interval::from_function(f, n)`, the coefficients being
// (2k + 1)/2 times the integrals of f P_k by the door of `integrate`.

SPATIUM_EXPORT namespace spatium {

// P_k(t) by the three-term recurrence (n + 1) P_{n+1} = (2n + 1) t P_n - n P_{n-1}.
template<Scalar T>
T legendre(std::size_t k, const T& t) {
    if (k == 0) return T{1};
    T p0{1}, p1 = t;
    for (std::size_t n = 1; n < k; ++n) {
        const T p2 = T((T(static_cast<double>(2 * n + 1)) * t * p1 - T(static_cast<double>(n)) * p0) / T(static_cast<double>(n + 1)));
        p0 = p1;
        p1 = p2;
    }
    return p1;
}

template<Scalar T = double>
struct LegendreSeries {
    std::vector<T> c;       // coefficients of P_0, P_1, ...
    T tail{0};              // estimate of the L^2 norm neglected beyond c

    friend LegendreSeries operator+(const LegendreSeries& x, const LegendreSeries& y) {
        LegendreSeries r;
        r.c.assign(std::max(x.c.size(), y.c.size()), T{0});
        for (std::size_t i = 0; i < x.c.size(); ++i) r.c[i] = r.c[i] + x.c[i];
        for (std::size_t i = 0; i < y.c.size(); ++i) r.c[i] = r.c[i] + y.c[i];
        r.tail = T(x.tail + y.tail);
        return r;
    }
    friend LegendreSeries operator-(const LegendreSeries& x) {
        LegendreSeries r = x;
        for (auto& v : r.c) v = -v;
        return r;
    }
    friend LegendreSeries operator-(const LegendreSeries& x, const LegendreSeries& y) { return x + (-y); }
    friend LegendreSeries operator*(T s, const LegendreSeries& x) {
        using std::abs;
        LegendreSeries r = x;
        for (auto& v : r.c) v = T(s * v);
        r.tail = T(abs(s) * x.tail);
        return r;
    }
    friend LegendreSeries operator*(const LegendreSeries& x, T s) { return s * x; }
    // Equal as far as the coefficients go (a missing one is zero).
    friend bool operator==(const LegendreSeries& x, const LegendreSeries& y) {
        const std::size_t n = std::max(x.c.size(), y.c.size());
        for (std::size_t i = 0; i < n; ++i)
            if (!((i < x.c.size() ? x.c[i] : T{0}) == (i < y.c.size() ? y.c[i] : T{0}))) return false;
        return true;
    }
};

template<Scalar T = double>
struct L2Interval {
    using ScalarType = T;
    using PointType = LegendreSeries<T>;
    using VectorType = LegendreSeries<T>;
    using TangentVector = LegendreSeries<T>;

    static constexpr Dimension dimension = Dimension::infinite();
    static constexpr bool is_complete = true;

    T a{-1}, b{1};

    constexpr L2Interval() = default;
    constexpr L2Interval(T lo, T hi) : a(lo), b(hi) {}

    constexpr bool contains(const PointType&) const { return true; }

    T inner(const PointType& f, const PointType& g) const {
        T sum{0};
        const std::size_t n = std::min(f.c.size(), g.c.size());
        for (std::size_t k = 0; k < n; ++k)
            sum = T(sum + T{2} / T(static_cast<double>(2 * k + 1)) * f.c[k] * g.c[k]);
        return T((b - a) / T{2} * sum);
    }
    T norm(const PointType& f) const { using std::sqrt; return sqrt(inner(f, f)); }
    T distance(const PointType& f, const PointType& g) const { return norm(f - g); }

    // The value of the series at x in [a, b].
    T evaluate(const PointType& f, T x) const {
        const T t = T((T{2} * x - a - b) / (b - a));
        T p_prev{1}, p = t, sum{0};
        if (!f.c.empty()) sum = T(f.c[0]);
        if (f.c.size() > 1) sum = T(sum + f.c[1] * t);
        for (std::size_t n = 1; n + 1 < f.c.size(); ++n) {
            const T next = T((T(static_cast<double>(2 * n + 1)) * t * p - T(static_cast<double>(n)) * p_prev) / T(static_cast<double>(n + 1)));
            p_prev = p;
            p = next;
            sum = T(sum + f.c[n + 1] * p);
        }
        return sum;
    }

    // The series of a function with n terms, and the L^2 mass of its last three coefficients as
    // the estimate of what was left out.
    template<class F>
        requires Function<F, T, T>
    PointType from_function(F&& f, std::size_t n) const {
        PointType r;
        r.c.assign(n, T{0});
        for (std::size_t k = 0; k < n; ++k) {
            const auto integrand = [&](T t) -> T {
                const T x = T((t * (b - a) + a + b) / T{2});
                return T(f(x) * legendre(k, t));
            };
            const auto integral = integrate(integrand, Finite<T>{T{-1}, T{1}});
            r.c[k] = T(T(static_cast<double>(2 * k + 1)) / T{2} * integral.value);
        }
        T mass{0};
        for (std::size_t k = n > 3 ? n - 3 : 0; k < n; ++k)
            mass = T(mass + T{2} / T(static_cast<double>(2 * k + 1)) * r.c[k] * r.c[k]);
        using std::sqrt;
        r.tail = T(sqrt(T((b - a) / T{2} * mass)));
        return r;
    }

    // The k-th Legendre polynomial as a point: the orthogonal basis of the space.
    PointType basis(std::size_t k) const {
        PointType r;
        r.c.assign(k + 1, T{0});
        r.c[k] = T{1};
        return r;
    }
};

static_assert(HilbertSpace<L2Interval<double>>);
static_assert(!EuclideanSpace<L2Interval<double>>, "no finite basis: not a Euclidean space");

} // namespace spatium
