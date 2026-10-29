#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <cassert>
#  include <cmath>
#  include <cstddef>
#  include <limits>
#  include <numbers>
#  include <span>
#endif

SPATIUM_EXPORT namespace spatium::statistics {

// ── How much a distribution says, and how far two are apart ───────
//
// Over probability vectors given as spans of the same length, in nats (natural logarithms);
// `in_bits` turns nats into bits. A term with p = 0 is 0 whatever q is (the limit of p log p).
//
//   entropy(p)            H(p)    = -sum p log p      how spread p is; log n for n equal chances
//   cross_entropy(p, q)   H(p, q) = -sum p log q      the cost of coding draws of p with a code
//                                                     built for q: in bits, a description length
//   kl_divergence(p, q)   KL(p|q) = sum p log(p / q)  what coding p with q costs beyond H(p)
//
// KL(p|q) = H(p, q) - H(p), it is zero only when p = q, and it is not symmetric. A q that gives
// zero to an outcome p can draw makes the cross-entropy and the divergence infinite, as they are.
//
// Templated on the scalar, so a derivative through Dual is a second route to any gradient
// written from these.

template <class T> T entropy(std::span<const T> p) {
    using std::log;
    T h{0};
    for (const T& x : p)
        if (x > T{0})
            h = h - x * log(x);
    return h;
}

template <class T> T cross_entropy(std::span<const T> p, std::span<const T> q) {
    using std::log;
    assert(p.size() == q.size());
    T h{0};
    for (std::size_t i = 0; i < p.size(); ++i) {
        if (!(p[i] > T{0}))
            continue;
        if (!(q[i] > T{0}))
            return T(std::numeric_limits<double>::infinity());
        h = h - p[i] * log(q[i]);
    }
    return h;
}

template <class T> T kl_divergence(std::span<const T> p, std::span<const T> q) {
    using std::log;
    assert(p.size() == q.size());
    T d{0};
    for (std::size_t i = 0; i < p.size(); ++i) {
        if (!(p[i] > T{0}))
            continue;
        if (!(q[i] > T{0}))
            return T(std::numeric_limits<double>::infinity());
        d = d + p[i] * log(p[i] / q[i]);
    }
    return d;
}

template <class T> T in_bits(T nats) { return nats / T(std::numbers::ln2); }

} // namespace spatium::statistics
