---
date: 2026-10-06
title: The volume element beside the connection
resolves: [2026-10-06-04#4]
---
`volume_element(metric, x)` = sqrt|det g| in `spaces/metric_chart.hpp`, next to `christoffel`, `geodesic_rhs` and `geodesic_step`, which were already N-dimensional -- and `MetricChart::volume_element(p)`. A consequence of the metric, not a separate notion: the review that asked for it was right that it belongs with the connection, and the two are tied by an identity that holds for every metric, Gamma^k_{ki} = d_i ln sqrt|det g|. Held by it, two descriptions of one number with nothing shared -- the contraction of the Christoffel symbols, and a `Dual` through the determinant -- on the sphere, the half-plane, a metric with off-diagonal terms, and Kerr in Boyer-Lindquist coordinates (four dimensions, det g < 0, hence the absolute value), to 1e-10; and the area of the sphere as the integral of the volume element, 4 pi, by `integrate_with_error`.
