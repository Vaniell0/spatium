# Relativity

## Where it stands

Metric-agnostic geodesic integration (`physics/relativity/`): Schwarzschild, Kerr, a metric as ten IR fields, its GLSL, the accretion disk, and the black-hole renderer (`examples/blackhole_live`).

## Now

- nothing is taken up here now

## Later

- two black holes merging, reporting on itself

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area relativity`.

## Reasoning and earlier history

### General relativity + analytical rendering

- `physics/relativity/` — metric-agnostic geodesic integration: `schwarzschild.hpp` and `kerr.hpp` (both metrics as templated, `Dual<T>`-substitutable callables), `geodesic.hpp` (exact Christoffel symbols via `Dual<T>` partials, Killing-vector conserved quantities), `accretion_disk.hpp` (Schwarzschild thin-disk redshift; Kerr's own BPT-1972 equatorial-orbit/ISCO/photon-orbit/redshift formulas).
- `render/` — `supersample_pixel()`, NxN jittered-grid antialiasing for any CPU raytracer built on Spatium's ray-surface primitives.
- `examples/tumbling_body_demo.cpp` — LGVI-integrated Dzhanibekov-effect (tumbling-body instability) demo.
- `examples/blackhole_gr_demo.cpp`, `wormhole_demo.cpp`, `geodesic_curvature_grid_demo.cpp` — gallery renders built on the above.
