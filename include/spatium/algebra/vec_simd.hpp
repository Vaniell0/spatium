#pragma once

// SIMD specializations for Vec<float,4> and Vec<double,4>.
// Guarded by ISA macros — no-op if compiler doesn't target SSE/AVX.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <cstddef>

#  if defined(__SSE2__)
#    include <emmintrin.h>  // SSE2
#  endif

#  if defined(__SSE4_1__)
#    include <smmintrin.h>  // SSE4.1 (_mm_dp_ps)
#  endif

#  if defined(__AVX2__)
#    include <immintrin.h>  // AVX2
#  endif
#endif

SPATIUM_EXPORT namespace spatium::simd {

// ── float×4 ───────────────────────────────────────────────────

#if defined(__SSE2__)

inline __m128 load_f(const float* p) { return _mm_loadu_ps(p); }
inline void store_f(float* p, __m128 v) { _mm_storeu_ps(p, v); }

inline void add_f4(const float* a, const float* b, float* out) {
    store_f(out, _mm_add_ps(load_f(a), load_f(b)));
}

inline void sub_f4(const float* a, const float* b, float* out) {
    store_f(out, _mm_sub_ps(load_f(a), load_f(b)));
}

inline void mul_scalar_f4(const float* a, float s, float* out) {
    store_f(out, _mm_mul_ps(load_f(a), _mm_set1_ps(s)));
}

inline void div_scalar_f4(const float* a, float s, float* out) {
    store_f(out, _mm_div_ps(load_f(a), _mm_set1_ps(s)));
}

inline float dot_f4(const float* a, const float* b) {
#if defined(__SSE4_1__)
    // _mm_dp_ps mask 0xFF: multiply all 4 lanes, store to all 4
    __m128 d = _mm_dp_ps(load_f(a), load_f(b), 0xFF);
    return _mm_cvtss_f32(d);
#else
    // SSE2 fallback: mul + hadd manually
    __m128 va = load_f(a);
    __m128 vb = load_f(b);
    __m128 m = _mm_mul_ps(va, vb);
    // m = [m0, m1, m2, m3]
    __m128 shuf = _mm_shuffle_ps(m, m, _MM_SHUFFLE(2, 3, 0, 1)); // [m1, m0, m3, m2]
    __m128 sums = _mm_add_ps(m, shuf);                            // [m0+m1, m1+m0, m2+m3, m3+m2]
    shuf = _mm_movehl_ps(shuf, sums);                             // [m2+m3, m3+m2, ...]
    sums = _mm_add_ss(sums, shuf);                                // [m0+m1+m2+m3, ...]
    return _mm_cvtss_f32(sums);
#endif
}

// Cross product for 3-of-4 lanes (w ignored)
inline void cross_f4(const float* a, const float* b, float* out) {
    __m128 va = load_f(a);
    __m128 vb = load_f(b);
    // a1*b2 - a2*b1, a2*b0 - a0*b2, a0*b1 - a1*b0
    __m128 a_yzx = _mm_shuffle_ps(va, va, _MM_SHUFFLE(3, 0, 2, 1));
    __m128 b_yzx = _mm_shuffle_ps(vb, vb, _MM_SHUFFLE(3, 0, 2, 1));
    __m128 c = _mm_sub_ps(_mm_mul_ps(va, b_yzx), _mm_mul_ps(a_yzx, vb));
    c = _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 0, 2, 1));
    store_f(out, c);
}

#endif // __SSE2__

// ── double×4 ──────────────────────────────────────────────────

#if defined(__AVX2__)

inline __m256d load_d(const double* p) { return _mm256_loadu_pd(p); }
inline void store_d(double* p, __m256d v) { _mm256_storeu_pd(p, v); }

inline void add_d4(const double* a, const double* b, double* out) {
    store_d(out, _mm256_add_pd(load_d(a), load_d(b)));
}

inline void sub_d4(const double* a, const double* b, double* out) {
    store_d(out, _mm256_sub_pd(load_d(a), load_d(b)));
}

inline void mul_scalar_d4(const double* a, double s, double* out) {
    store_d(out, _mm256_mul_pd(load_d(a), _mm256_set1_pd(s)));
}

inline void div_scalar_d4(const double* a, double s, double* out) {
    store_d(out, _mm256_div_pd(load_d(a), _mm256_set1_pd(s)));
}

inline double dot_d4(const double* a, const double* b) {
    __m256d va = load_d(a);
    __m256d vb = load_d(b);
    __m256d m = _mm256_mul_pd(va, vb);
    // m = [m0, m1, m2, m3]
    __m256d sum1 = _mm256_hadd_pd(m, m);          // [m0+m1, m0+m1, m2+m3, m2+m3]
    __m128d lo = _mm256_castpd256_pd128(sum1);     // [m0+m1, m0+m1]
    __m128d hi = _mm256_extractf128_pd(sum1, 1);   // [m2+m3, m2+m3]
    __m128d total = _mm_add_sd(lo, hi);            // [m0+m1+m2+m3, ...]
    return _mm_cvtsd_f64(total);
}

#else // !__AVX2__ — stubs so if-constexpr dead branches compile

inline void add_d4(const double*, const double*, double*) {}
inline void sub_d4(const double*, const double*, double*) {}
inline void mul_scalar_d4(const double*, double, double*) {}
inline void div_scalar_d4(const double*, double, double*) {}
inline double dot_d4(const double*, const double*) { return 0.0; }

#endif // __AVX2__

} // namespace spatium::simd
