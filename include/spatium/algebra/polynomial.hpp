#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/complex.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <array>
#  include <cmath>
#  include <numbers>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// ── Quadratic solver ──────────────────────────────────────────
// ax² + bx + c = 0 → 2 roots (possibly complex)

template<Scalar T>
std::array<Complex<T>, 2> solve_quadratic(T a, T b, T c) {
    using std::sqrt; // ADL: lets non-std Scalar T (e.g. Real50) provide its own sqrt
    auto disc = b * b - T{4} * a * c;

    if (disc >= T{0}) {
        auto sq = sqrt(disc);
        return {Complex<T>{(-b + sq) / (T{2} * a)},
                Complex<T>{(-b - sq) / (T{2} * a)}};
    } else {
        auto sq = sqrt(-disc);
        return {Complex<T>{-b / (T{2} * a),  sq / (T{2} * a)},
                Complex<T>{-b / (T{2} * a), -sq / (T{2} * a)}};
    }
}

// Convenience: real roots only
template<Scalar T>
std::vector<T> real_roots_quadratic(T a, T b, T c, T eps = epsilon<T>()) {
    auto roots = solve_quadratic(a, b, c);
    std::vector<T> result;
    for (auto& r : roots) {
        if (r.is_real(eps)) result.push_back(r.re);
    }
    return result;
}

// ── Cubic solver (Cardano) ────────────────────────────────────
// ax³ + bx² + cx + d = 0 → 3 roots

template<Scalar T>
std::array<Complex<T>, 3> solve_cubic(T a, T b, T c, T d) {
    // ADL: lets non-std Scalar T (e.g. Real50) provide its own overloads.
    using std::sqrt, std::cbrt, std::abs, std::acos, std::cos, std::atan;

    // Normalize: x³ + px² + qx + r = 0
    auto p = b / a;
    auto q = c / a;
    auto r = d / a;

    // Depressed cubic: t³ + pt + q = 0  where  x = t - p/3
    auto p1 = q - p * p / T{3};
    auto q1 = r - p * q / T{3} + T{2} * p * p * p / T{27};

    auto disc = q1 * q1 / T{4} + p1 * p1 * p1 / T{27};

    auto shift = -p / T{3};
    // std::numbers::pi_v<T> is only defined for the standard floating-point
    // types, not Real50/Real100 -- 4*atan(1) is portable to any Scalar T
    // with an ADL-visible atan, and this function is already non-constexpr
    // (sqrt/cbrt/cos aren't either) so nothing is lost by computing it here.
    auto pi = T{4} * atan(T{1});

    if (disc > epsilon<T>()) {
        // One real root, two complex conjugate
        auto sq = sqrt(disc);
        auto u = cbrt(-q1 / T{2} + sq);
        auto v = cbrt(-q1 / T{2} - sq);
        auto real_root = u + v + shift;

        // Complex roots via Vieta
        auto re_part = -(u + v) / T{2} + shift;
        auto im_part = (u - v) * sqrt(T{3}) / T{2};
        return {Complex<T>{real_root},
                Complex<T>{re_part, im_part},
                Complex<T>{re_part, -im_part}};
    } else if (abs(disc) <= epsilon<T>()) {
        // All real, at least two equal. u = cbrt(-q1/2) in general, but
        // cbrt has an infinite derivative at 0 -- for the triple-root case
        // (q1 itself already within tolerance of zero, t^3 = 0 exactly),
        // computing cbrt of whatever tiny residual arithmetic noise q1
        // actually carries amplifies that residual by roughly a cube-root
        // power law instead of damping it, turning a machine-epsilon-scale
        // q1 error into a much larger (epsilon^(1/3)-scale) error in the
        // returned root. Real, reproduced case: solve_cubic<Real50>(1,-6,
        // 12,-8) -- exact integer coefficients for (x-2)^3 -- passes under
        // one Boost.Multiprecision build (q1 computed near bit-exact) but
        // failed CI under another (a few-ULP q1 residual, amplified via
        // cbrt into a visible ~1e-24 root error against a 1e-25 tolerance;
        // the two returned duplicate roots' errors coming out equal and
        // the third exactly double is the -u/-u/2u structure below showing
        // through). Skip the amplification entirely when it can't matter.
        auto u = (abs(q1) <= epsilon<T>()) ? T{0} : cbrt(-q1 / T{2});
        return {Complex<T>{T{2} * u + shift},
                Complex<T>{-u + shift},
                Complex<T>{-u + shift}};
    } else {
        // Three distinct real roots (casus irreducibilis)
        auto m = T{2} * sqrt(-p1 / T{3});
        auto theta = acos(T{3} * q1 / (p1 * m)) / T{3};
        return {Complex<T>{m * cos(theta) + shift},
                Complex<T>{m * cos(theta - T{2} * pi / T{3}) + shift},
                Complex<T>{m * cos(theta - T{4} * pi / T{3}) + shift}};
    }
}

// ── Quartic solver (Ferrari) ──────────────────────────────────
// ax⁴ + bx³ + cx² + dx + e = 0 → 4 roots

template<Scalar T>
std::array<Complex<T>, 4> solve_quartic(T a, T b, T c, T d, T e) {
    using std::abs, std::sqrt; // ADL: lets non-std Scalar T (e.g. Real50) provide these
    // Normalize: x⁴ + px³ + qx² + rx + s = 0
    auto p = b / a;
    auto q = c / a;
    auto r = d / a;
    auto s = e / a;

    // Depressed quartic: y⁴ + αy² + βy + γ = 0  where  x = y - p/4
    auto shift = -p / T{4};
    auto alpha = q - T{3} * p * p / T{8};
    auto beta  = r - p * q / T{2} + p * p * p / T{8};
    auto gamma = s - p * r / T{4} + p * p * q / T{16} - T{3} * p * p * p * p / T{256};

    if (abs(beta) < epsilon<T>()) {
        // Biquadratic: y⁴ + αy² + γ = 0 → solve as quadratic in y²
        auto sub = solve_quadratic(T{1}, T(alpha), T(gamma));
        auto r0 = sqrt(sub[0]);
        auto r1 = sqrt(sub[1]);
        // T(...) forces eager evaluation: Boost.Multiprecision's number<>
        // arithmetic returns lazy expression-template types, not T itself,
        // and Complex<T>{expr, expr} can't implicitly convert both of a
        // 2-argument brace-init at once -- an explicit cast per argument
        // sidesteps that instead of relying on an implicit conversion.
        return {Complex<T>{T(r0.re + shift), T(r0.im)},
                Complex<T>{T(-r0.re + shift), T(-r0.im)},
                Complex<T>{T(r1.re + shift), T(r1.im)},
                Complex<T>{T(-r1.re + shift), T(-r1.im)}};
    }

    // Resolvent cubic: 8m³ + 8αm² + 2(α²-4γ)m - β² = 0
    //
    // Derivation: completing the square on the depressed quartic gives, as
    // an identity in y for any m (no need to assume the quartic holds),
    //   (y²+α/2+m)² - [2m·y² - βy + (m²+αm+α²/4-γ)] = y⁴+αy²+βy+γ.
    // The bracket is a perfect square in y exactly when its discriminant
    // vanishes: β² - 4·(2m)·(m²+αm+α²/4-γ) = 0, i.e. this cubic. A
    // previous version of this comment (and the code) had -4αm² here --
    // an outright wrong resolvent, not a root-selection nuance -- which
    // meant NO choice of resolvent root reconstructed the original
    // quartic's roots correctly except by accident.
    auto cubic_roots = solve_cubic(T{8}, T(T{8} * alpha),
                                    T(T{2} * (alpha * alpha - T{4} * gamma)),
                                    T(-(beta * beta)));

    // Pick the LARGEST real root of the resolvent, not "the first real
    // root above epsilon, else any real root" (the old logic, which could
    // and did land on a near-zero or negative m -- see torus_quartic's
    // comment for the observed fallout). At m=0 the resolvent evaluates
    // to -β² < 0 (this branch only runs once β is confirmed non-negligible
    // above), and its leading 8m³ term dominates as m→+∞, so by the
    // intermediate value theorem a positive real root is always present
    // when the coefficients are real. Taking the largest real root both
    // guarantees landing on that positive one and keeps 2m as large as
    // the resolvent allows, which keeps sqrt(2m) away from zero -- the
    // same amplify-tiny-residuals concern the solve_cubic repeated-root
    // fix documents, just showing up here as a division by a near-zero
    // sqrt(2m) instead of a cbrt of a near-zero residual.
    T m{};
    bool have_m = false;
    for (auto& cr : cubic_roots) {
        if (cr.is_real() && (!have_m || cr.re > m)) {
            m = cr.re;
            have_m = true;
        }
    }

    // Declared as T, not auto: Boost.Multiprecision's number<> composes
    // abs()/sqrt() into lazy expression-template objects, and letting auto
    // deduce that expression type here (rather than the concrete T it
    // evaluates to) leaves sq_2m holding a dangling reference into a
    // temporary from this statement -- garbage silently, not a compile
    // error. Reproduced directly: solve_quartic<Real50> on an asymmetric-
    // root case landed sq_2m on a nonsense value despite m itself being
    // exactly right, only once a test actually exercised this branch
    // (every prior Real50 quartic test used symmetric roots, which take
    // the biquadratic shortcut above and never reach this line). abs() is
    // a NaN guard against sub-epsilon negative noise on an m that should
    // be exactly >= 0 by the argument above, not a correctness mechanism
    // anymore.
    T sq_2m = sqrt(abs(T{2} * m));

    // Ferrari's difference-of-squares step factors the quartic as
    //   (y² - √(2m)·y + (α/2+m+β/(2√(2m)))) · (y² + √(2m)·y + (α/2+m-β/(2√(2m)))) = 0
    // Both the missing α/2 term and the sign pairing between √(2m) and
    // β/(2√(2m)) were wrong in the old code (verified by reconstructing a
    // quartic with known roots by hand): dropping α/2 and swapping which
    // sign of √(2m) goes with which sign of β/(2√(2m)) both produce a
    // factorization that does NOT multiply back out to the original
    // quartic, for any choice of m.
    // T, not auto, for the same dangling-expression-template reason as
    // sq_2m above.
    T half_beta_over_sq = beta / (T{2} * sq_2m);
    T half_alpha = alpha / T{2};

    // T(...) on every argument: Boost.Multiprecision's number<> arithmetic
    // returns lazy expression-template types, not T -- passed straight into
    // solve_quadratic<T>, its OWN template parameter would need deducing
    // from three differently-typed expressions at once instead of a single
    // consistent T. Forcing each to T here is what makes that deduction
    // (and everything downstream that reads .re/.im off the result) honest.
    auto roots1 = solve_quadratic(T{1}, T(-sq_2m), T(half_alpha + m + half_beta_over_sq));
    auto roots2 = solve_quadratic(T{1}, T(sq_2m), T(half_alpha + m - half_beta_over_sq));

    return {Complex<T>{T(roots1[0].re + shift), T(roots1[0].im)},
            Complex<T>{T(roots1[1].re + shift), T(roots1[1].im)},
            Complex<T>{T(roots2[0].re + shift), T(roots2[0].im)},
            Complex<T>{T(roots2[1].re + shift), T(roots2[1].im)}};
}

// Convenience: real roots only
template<Scalar T>
std::vector<T> real_roots_cubic(T a, T b, T c, T d, T eps = epsilon<T>()) {
    auto roots = solve_cubic(a, b, c, d);
    std::vector<T> result;
    for (auto& r : roots) {
        if (r.is_real(eps)) result.push_back(r.re);
    }
    return result;
}

template<Scalar T>
std::vector<T> real_roots_quartic(T a, T b, T c, T d, T e, T eps = epsilon<T>()) {
    auto roots = solve_quartic(a, b, c, d, e);
    std::vector<T> result;
    for (auto& r : roots) {
        if (r.is_real(eps)) result.push_back(r.re);
    }
    return result;
}

} // namespace algebra
} // namespace spatium
