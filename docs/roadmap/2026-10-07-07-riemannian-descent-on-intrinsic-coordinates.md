---
date: 2026-10-07
title: Riemannian descent on a space with intrinsic coordinates
stage: 2.1
open:
  - the geodesic and Lie-group algorithms on the customization points
  - riemannian_minimize reads exp and metric as members, not through core/access.hpp
resolves: [2026-10-07-05#1, 2026-10-07-06#2]
---
Stage 2.1. `riemannian_minimize` required `HasNormal`, because it projected the raised gradient onto the tangent space by the normal -- which only a hypersurface of an ambient space has. A space whose points are coordinates of the manifold itself (the kappa family, a chart) has no normal and nothing to project: the raised gradient is already tangent. The projection is now a compile-time branch (`if constexpr (HasNormal<S>)`), and the requirement is `RiemannianManifold` alone; Sphere and Hyperbolic go through exactly the code they did (the existing tests are unchanged). Held by a pair of paths that share nothing: the minimiser of the summed squared distances to three points found from the gradient of that sum through a `Dual` (`riemannian_minimize`) and as the fixed point of the mean of the logs (`frechet_mean`), on `ConstantCurvature<2>` at kappa = -1, 0 and 1, to 1e-6; the test did not compile before the change. Not done: the geodesic and Lie-group algorithms on the customization points, and `riemannian_minimize` reading exp and metric through `core/access.hpp` rather than as members.
