# Lowering

## Where it stands

What a computed field becomes on a device: GLSL for Vulkan compute (std430) and for WebGL2 (GLSL ES 3.00), a symbol read from a uniform slot and a net's weights from a uniform array, so a moved number or a trained net is new data and not a new shader (`io/field_glsl.hpp`). The graph of a function is a fragment shader (`io/plot_glsl.hpp`), checked in a headless browser against the core's fp32 twin. One GPU path: Vulkan compute from the IR (`docs/gpu.md`).

## Now

- nothing is taken up here now

## Later

- slots and a net's weights in Vulkan
- CUDA C from the IR, for those who port the library
- a net wider than three inputs, which needs vectors in the IR

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area lowering`.

## Reasoning and earlier history

### ImGui viewer, Bohr model

- ImGui viewer: element picker, shell toggles, view mode controls
- Bohr model electron orbital visualization
- Dirty buffers: lazy GPU upload on mesh change

### Analytical-render dispatcher, benchmark suite, table fix (2026-04-23)

- **`examples/primitives_demo.cpp`** — multi-scene Vulkan dispatcher: `--scene primitives` (unified sphere/box/torus/cylinder/cone/ellipsoid/triangle with BVH raycast hit-cloud + console report), `--scene torus` (Clifford torus S³→R³), `--scene klein` (animated Klein bottle R⁴→R³). Absorbs and removes the former `viewer_demo`, `torus_demo` and `klein_demo`.
- **`benchmarks/bench_raycast.cpp`** — BVH (60–71 ns) vs brute (8.7 µs – 519 µs) vs analytical `ray_quadric` (23–29 ns) across 320 / 5 120 / 20 480 triangles
- **Table renderer fix** (`include/spatium/io/table.hpp`) — separator width off-by-one (`w+1` → `w+2`); console tables now render with aligned borders. (Path corrected 2026-09-15: this said `src/io/table.cpp`, which has never existed in this tree — `src/` holds three files and none of them is the table renderer, which is header-only like the rest of `io/`. Found by `scripts/check_doc_file_refs.py`, not by hand.)
- **Eigen integration decision** — remain optional via `SPATIUM_EIGEN=ON`; heat method + differential + interop stay isolated. Vec/Matrix not migrated to Eigen; `Eigen::Ref<…>` planned for public APIs only once dense QR/SVD lands (see Backlog → Native math)
- ~~CPU raytracer: Quadric → viewer texture (analytical render without mesh)~~ done: `ray_quadric`/`ray_quadric_proximity` in `geometry/ray_surface.hpp`, exercised by `examples/parametric_analytical_demo.cpp`'s `glow_sphere.png`; ships PNG output rather than a live viewer texture, which covers the same "analytical render without mesh" goal for the gallery use case
- ~~Ray-ParametricSurface: Newton UV iteration for arbitrary f(u,v) surfaces~~ done: `geometry/ray_parametric.hpp` — no tessellation, no BVH, Newton-solves `S(u,v) = o + t·d`; exercised by `examples/parametric_analytical_demo.cpp` producing `klein_analytical.png`/`mobius_analytical.png`/`bumpy_analytical.png`/`parametric_gallery.png`
- ~~Ray-quartic: torus intersection via solve_quartic~~ done: `geometry/ray_surface.hpp`'s `Torus<T>` + `ray_torus()`/`ray_torus_proximity()`, exercised by `primitives_demo`'s analytic-torus raycast row and `torus_analytical.png` — this had drifted stale in Backlog → Analytical rendering (fixed 2026-08-31)

### GPU rendering (CUDA) — kernels built and cross-validated

`gpu/` exists and works, and is genuinely necessary, not optional polish: CPU-only was directly measured and extrapolated to ~20 days for the actual target render (1920×1080, ~750 frames, Kerr flyby) — not a guess, a real number from a real benchmark, after Mojo was investigated and rejected in favor of plain CUDA C++ (T4 fp64 throughput is 1/32 of fp32, which settles the language question regardless of Mojo-vs-CUDA specifics) and an fp32 precision gate was verified to pass cleanly across 22 cross-validated cases. `geodesic_kernel.cu` (batch Schwarzschild/Kerr null-geodesic tracing), `schwarzschild_render_kernel.cu`/`kerr_render_kernel.cu` (full per-ray volumetric disk emission-absorption, ~2K LOC total), `render_4k_frame.cpp` (the actual production render driver) — cross-checked against the CPU path per-ray at 160×90 (`gpu/verify_cuda_render.cpp`, abs/rel diff against the same real library calls the CPU demo uses), not just visually eyeballed.

Closed-form Christoffel symbols were hand-derived for exactly the two metrics the GPU path needs, rather than sharing the generic `Dual<T>`-templated CPU headers directly as `__host__ __device__` code — that was considered and deliberately rejected, since `Dual<T>`'s general autodiff machinery is branch/indirection-heavy in a way that's specifically bad for SIMT throughput. See Backlog → GPU rendering for what's still open (landing the production render, and the one remaining mechanization gap in how those closed-form symbols get into the header).

### Analytical rendering

- **[dare]** Fragment shader raymarcher: GPU-native analytical render (GLSL quadric math)
- **[want]** Dual Quaternion
- **[want]** A real rainbow, and then interference — light given the same wave treatment the sound work gives sound. Proposed 2026-09-17, and the two halves are deliberately one entry because the second is what makes the first worth building. A rainbow alone is dispersion: one more `ray_hit` on a sphere with a wavelength-dependent index, two internal reflections, the 42° caustic falling out rather than being drawn. That is a demo. The step is the second half — treating light as the wave the `wave_string.hpp`/`wave_membrane.hpp` steppers already treat sound as, so that a two-slit pattern is *integrated* rather than plotted from the closed-form fringe spacing. Phase becomes a quantity a ray carries, and the same machinery that gives a membrane its modes gives a field its fringes. The honest caveat up front: this is classical wave optics, and the distance from there to quantum mechanics is a real one — amplitudes that interfere are not yet probability amplitudes, and saying otherwise would be exactly the overclaim this file keeps catching. What it genuinely is: the point where the library's wave equation stops belonging to acoustics and becomes a shared piece.

### GPU rendering

One path, Vulkan compute from the IR ([`gpu.md`](gpu.md)); what crosses out of the library and how it is lowered is [`gpu-abi-design.md`](gpu-abi-design.md). Decided 2026-10-30: no C ABI for the DSL -- the boundary is the document as text or JSON, read by the same reader in the host and in WebAssembly.

- **[course]** A symbol lowered to a read from a uniform or a buffer, so that structure is code and numbers are data: a shader compiled when a document's structure changes, not when a parameter, the time or a weight moves. `lower()` in `io/field_pod.hpp` and `io/field_glsl.hpp` refuses a symbol today; the POD interpreter is the host twin.
- **[course]** A dense layer of a network as one IR operation with its weights in a buffer, lowered like any field, `TanhNet<double>` its host twin; a neural shader (an image filter, a learned field in a scene) is then a field.
- **[course]** GLSL for WebGL2 from the same IR, uniforms and textures in place of std430 storage, for the page's scene windows.
- **[course]** Land the production render (1920x1080, ~750 frames, Kerr flyby) with `blackhole_live --video`; `gallery/blackhole_gr.mp4` is a partial preview, not the full sequence.
