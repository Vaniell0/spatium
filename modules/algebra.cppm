// Primary module unit for spatium.algebra. Re-exports all partitions.
// `vec_expr` is folded into `:vector` (see vector.cppm rationale).
export module spatium.algebra;
export import :concepts;
export import :vec_simd;
export import :vector;
export import :matrix;
export import :format;
export import :literals;
export import :complex;
export import :polynomial;
export import :functions;
export import :quaternion;
export import :verify;
export import :groups_so3;
export import :groups_se3;
export import :dual;
export import :ode;
export import :linear_solve;
export import :calculus;
// algebra/eigen_interop.hpp stays header-only in Phase 1 — SSE intrinsics
// pulled by Eigen/Core and by vec_simd's GMF cause conflicting language linkage
// when both BMIs land in the same TU. Revisit once header units (or a SIMD
// rework) lets us share intrinsic decls. Module consumers can still
// `#include <spatium/algebra/eigen_interop.hpp>` directly.
