#pragma once

// ── spatium/core.hpp ─────────────────────────────────────────
// Focused entry point covering the typical "math + geometry"
// usage: algebra (Vec/Matrix/Quat/Complex/Polynomial/literals),
// concepts and verifiers, the canonical spaces (Euclidean,
// Sphere, Hyperbolic, Product, ParametricSurface), all geometry
// primitives and operations, the BVH, and the textual IO helpers
// (Table). Pulls in NO mesh data structure and NO physics models.
// Prefer this over <spatium/spatium.hpp> when you do not need
// mesh subdivision/heat method, SVG/OBJ/STL writers, or atoms.

#include <spatium/core/concepts.hpp>
#include <spatium/core/error.hpp>
#include <spatium/core/verify.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/algebra/vector.hpp>
#include <spatium/algebra/matrix.hpp>
#include <spatium/algebra/quaternion.hpp>
#include <spatium/algebra/complex.hpp>
#include <spatium/algebra/polynomial.hpp>
#include <spatium/algebra/functions.hpp>
#include <spatium/algebra/literals.hpp>
#include <spatium/algebra/format.hpp>
#include <spatium/point.hpp>
#include <spatium/morphism.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/spaces/product.hpp>
#include <spatium/geometry/geometry.hpp>
#include <spatium/spatial/bvh.hpp>
#include <spatium/io/table.hpp>
