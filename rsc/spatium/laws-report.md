# Laws of the spaces

Written by `rsc/spatium/rsc_laws.py discover`; checked by `rsc_laws.py check`. Terms up to size 7 over log, exp, the metric and the distance, 790 kept and 357 merged, on 25 spaces of the connectivity matrix, 732 inputs each.

834 laws in 24 classes. A class is the set of spaces a law holds on; it is set against the library's concepts on the spaces that are nodes of graph.json (E1, E3, E4, S1, S2, S3, S7, H1, H2, H4, SPD2-logE, SPD2-affine, Sphere-implicit, Sphere-implicit-typed, SPD2-affinexS2).

| spaces | laws | the library's concept on those spaces | shortest law |
|---|---|---|---|
| every space | 437 | -- | `0 <= g_p(u,u)` |
| E1, E3, E4, S1, S2, S3, S7, H1, H2, H4, S2xE1, SPD2-logE, SPD2-affine, Torus-typed, Sphere-implicit-typed, S2-derived, H2-derived, Cylinder-ADL, S2xH2, H2xE3, S2xS2, SPD2-affinexS2, CylinderxE3 | 330 | Distanced | `0 = d(p,p)` |
| E1, E3, E4, S1, S2, S3, S7, H1, H2, H4, S2xE1, SPD2-logE, SPD2-affine, Sphere-implicit-typed, S2-derived, H2-derived, Cylinder-ADL, S2xH2, H2xE3, S2xS2, SPD2-affinexS2, CylinderxE3 | 6 | Distanced | `-log(p,exp(p,u)) = log(p,exp(p,-u))` |
| E1, E3, E4, S1, S2, S3, S7, H2, H4, S2xE1, SPD2-logE, SPD2-affine, Torus-typed, Sphere-implicit-typed, S2-derived, H2-derived, Cylinder-ADL, S2xH2, H2xE3, S2xS2, SPD2-affinexS2, CylinderxE3 | 2 | -- | `(d(p,q)+d(q,q)) <= (d(p,r)+d(q,r))` |
| E1, E3, E4, S1, S2, S3, S7, H1, H2, H4, S2xE1, SPD2-logE, SPD2-affine, Sphere-implicit-typed, S2-derived, H2-derived, Cylinder-ADL, S2xH2, H2xE3, SPD2-affinexS2, CylinderxE3 | 2 | Distanced | `d(p,exp(p,u/2)) <= 2*d(p,exp(p,u))` |
| E1, E3, E4, S1, S2, S3, S7, S2xE1, SPD2-logE, SPD2-affine, Torus-typed, Sphere-implicit-typed, S2-derived, H2-derived, Cylinder-ADL, S2xH2, H2xE3, S2xS2, SPD2-affinexS2, CylinderxE3 | 2 | -- | `d(p,q) = (d(p,q)+d(q,q))` |
| E1, E3, E4, S1, S2, S3, S7, S2xE1, SPD2-logE, SPD2-affine, Torus-typed, Sphere-implicit-typed, S2-derived, H2-derived, Cylinder-ADL, H2xE3, S2xS2, SPD2-affinexS2, CylinderxE3 | 1 | -- | `d(p,p) = d(exp(p,w),exp(p,w))` |
| E1, E3, E4, S1, S2, S3, S7, S2xE1, SPD2-logE, SPD2-affine, Torus-typed, Sphere-implicit-typed, S2-derived, H2-derived, Cylinder-ADL, S2xS2, SPD2-affinexS2, CylinderxE3 | 4 | -- | `0 = d(q,q)` |
| E1, E3, E4, H1, H2, H4, SPD2-logE, SPD2-affine, Torus-typed, S2-derived, H2-derived, S2xH2, H2xE3, S2xS2, SPD2-affinexS2, CylinderxE3 | 1 | -- | `1/2*d(p,exp(p,w)) <= d(p,exp(p,2w))` |
| E1, E3, E4, H1, H2, H4, S2xE1, SPD2-logE, SPD2-affine, S2-derived, H2-derived, Cylinder-ADL, S2xH2, H2xE3, SPD2-affinexS2, CylinderxE3 | 1 | -- | `d(p,exp(p,w/2)) <= d(p,exp(p,w))` |
| E1, E3, E4, H1, H2, H4, SPD2-logE, SPD2-affine, Torus-typed, S2-derived, H2-derived, Cylinder-ADL, H2xE3, SPD2-affinexS2, CylinderxE3 | 6 | -- | `u = 2log(p,exp(p,u/2))` |
| E1, E3, E4, H1, H2, H4, SPD2-logE, SPD2-affine, S2-derived, H2-derived, Cylinder-ADL, S2xH2, H2xE3, SPD2-affinexS2, CylinderxE3 | 1 | -- | `d(p,exp(p,u/2)) <= d(p,exp(p,u))` |
| E1, E3, E4, H1, H2, H4, SPD2-logE, SPD2-affine, Torus-typed, S2-derived, H2-derived, S2xH2, H2xE3, SPD2-affinexS2, CylinderxE3 | 1 | -- | `1/2*d(p,exp(p,u)) <= d(p,exp(p,2u))` |
| E1, E3, E4, H1, H2, H4, SPD2-logE, SPD2-affine, S2-derived, H2-derived, S2xH2, H2xE3, SPD2-affinexS2, CylinderxE3 | 2 | -- | `d(p,exp(p,u/2)) <= d(p,exp(p,2u))` |
| E1, E3, E4, H1, H2, H4, SPD2-logE, SPD2-affine, Torus-typed, S2-derived, H2-derived, H2xE3, SPD2-affinexS2 | 1 | -- | `d(p,exp(p,w)) <= d(p,exp(p,2w))` |
| E1, E3, E4, H1, H2, H4, SPD2-logE, SPD2-affine, Torus-typed, S2-derived, H2-derived, S2xH2, H2xE3 | 1 | -- | `d(p,exp(p,u)) <= d(p,exp(p,2u))` |
| E1, E3, E4, H1, H2, H4, SPD2-logE, SPD2-affine, S2-derived, H2-derived, H2xE3 | 14 | -- | `u = log(p,exp(p,u))` |
| E1, E3, E4, H1, H2, H4, SPD2-logE, SPD2-affine, Torus-typed, S2-derived, H2xE3 | 2 | -- | `2*d(p,exp(p,u)) <= d(p,exp(p,2u))` |
| E1, E3, E4, H1, H2, H4, SPD2-logE, SPD2-affine, S2-derived, H2xE3 | 12 | -- | `u = log(p,exp(p,2u))/2` |
| H4, H2xE3, S2xS2, SPD2-affinexS2 | 1 | -- | `1/2*1/2*d(p,q) <= d(q,exp(p,(u+w)))` |
| S7, H2xE3, CylinderxE3 | 2 | -- | `g_p(u,-(u+w)) <= g_p(w,2(u+w))` |
| S7, S2xH2, SPD2-affinexS2 | 1 | -- | `1/2*1/2*d(p,r) <= d(r,exp(p,(u+w)))` |
| E4, S7 | 2 | -- | `g_p(u,-(u+w)) <= g_p(w,(u+w)/2)` |
| S7 | 2 | -- | `1/2*d(p,exp(p,u)) <= d(p,exp(p,(u+w)))` |

## Tasks (5 open, 4 triaged)

- H1: breaks a law of its neighbours (5 laws): `d(p,p) = d(exp(p,w),exp(p,w))`
- H2: breaks a law of its neighbours (3 laws): `d(p,p) = d(exp(p,w),exp(p,w))`
- H4: breaks a law of its neighbours (3 laws): `d(p,p) = d(exp(p,w),exp(p,w))`
- S2xH2: breaks a law of its neighbours (1 laws): `d(p,p) = d(exp(p,w),exp(p,w))`
- SPD2-affine: missing operation: `unary minus on a tangent`

Triaged, no defect:

- S2xS2: `d(p,exp(p,u/2)) <= 2*d(p,exp(p,u))` -- a compact space: past the cut locus exp(p, u) can come back near p while exp(p, u/2) does not, so the law is local; whether a long sample lands there depends on how u splits between the factors
- Sphere-implicit: `p = exp(q,log(q,p))` -- an erased level set: its log is a retraction, the chord projected on the tangent plane, and the library withholds its distance for that reason (core/access.hpp)
- Torus-chart: `p = exp(q,log(q,p))` -- an erased chart: its exp is project(p + v t), a retraction, so exp after log is not the identity
- Torus-typed: `-log(p,exp(p,u)) = log(p,exp(p,-u))` -- the law holds where the geodesic symmetry at p is an isometry -- a symmetric space -- and a torus of revolution is not one; past the cut locus the two sides are different points' logs
