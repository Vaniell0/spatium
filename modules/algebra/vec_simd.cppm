module;
#include <cassert>
#if defined(__SSE2__)
#  include <emmintrin.h>
#endif
#if defined(__SSE4_1__)
#  include <smmintrin.h>
#endif
#if defined(__AVX2__)
#  include <immintrin.h>
#endif
export module spatium.algebra:vec_simd;
import std.compat;
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/algebra/vec_simd.hpp>
