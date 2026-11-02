# Sound

## Where it stands

The synth engines and scores; the reasoning is below. The current work lives on a local branch and is not taken now.

## Now

- nothing is taken up here now

## Later

- the work on its local branch, when it is taken up

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area sound`.

## Reasoning and earlier history

### Sound synthesis

- **[want]** Physically-based sound synthesis — resonant modes of a surface via Laplace-Beltrami eigenvalues, geometric room acoustics via ray-traced reflections on the existing BVH, and time-domain wave-equation integration (string/membrane/percussion) reusing `algebra/ode.hpp`.

  **Shipped:** the wave-equation half. `physics/mechanics/wave_string.hpp` and `wave_membrane.hpp` are explicit central-difference steppers checked against closed-form fundamental frequencies, `io/wav.hpp` is the zero-dependency PCM writer, and `examples/sound_synthesis_demo.cpp` renders with both. Audio output is done, not pending.

  **Blocked:** the modal half, and on a bigger thing than this entry used to claim. Two sentences here were false and are worth recording as such rather than quietly deleting, because both pointed the work at the wrong place.

  The first said modes were *"already computable through the existing heat-method infrastructure"*. Measured, they are not. `HeatSolver` factorises `A + t·L`, which is `t·(L − σ·A)` at `σ = −1/t` with `t = h²` — so the shift is pinned to the mesh spacing rather than to the spectrum. On a 162-vertex sphere that puts `σ ≈ −11.2` against a first non-trivial eigenvalue of `λ₁ ≈ 2.0`, and shift-invert's discrimination `μᵢ = 1/(λᵢ − σ)` comes out at 2.0× across the first ten modes where a spectrum-chosen shift gives 568×. It also gets worse under refinement, not better: `h → 0` drives `σ → −∞` while the physical frequencies stay put. The second cached factorisation is no better — `L_pinned` has no mass term at all and is clamped at vertex 0, which splits the exact threefold `ℓ=1` degeneracy into `(0.152192, 0.152192, 0.193486)`, so it would distort mode *shapes* and not merely shift frequencies.

  The second said *"the only missing piece is audio output"*. The real missing piece is a **sparse generalised eigensolver**, which is a dependency decision rather than an afternoon: Eigen's own eigensolvers are all dense, and its sparse one is a wrapper that links ARPACK, a Fortran library not in the tree. The three routes are ARPACK (a first non-header mandatory dependency), vendoring Spectra header-only under `include/spatium/vendor/` the way `stb_image_write.h` already is, or hand-rolling shift-invert Lanczos over the `SimplicialLDLT` the project already uses — the last being closest to how `json.hpp`, `wav.hpp` and the Jacobi eigendecomposition were done.

  `build_laplacian()` and `build_mass_matrix()` in `mesh/differential.hpp` return `Eigen::SparseMatrix` and need no conversion; what they need is a fresh factorisation of `K − σM` at a small negative shift, since `K` is exactly singular.

- **[want]** Modes of a solid body, for engineering as much as for sound. A surface's Laplace-Beltrami modes are the sound of a membrane stretched over it, and its bending modes are a thin shell's; neither is how a solid part rings. That is linear elasticity on a volume: `K u = ω² M u` with `K` from a tetrahedral discretisation of the strain energy `physics/mechanics/continuum.hpp` already states (Saint-Venant-Kirchhoff, linearised) and `M` the lumped or consistent mass. The same sparse generalised eigensolver the surface modes need serves it, at three unknowns per vertex. The consumer is not a game: natural frequencies and mode shapes of a part are a standard CAD/CAE result, and a stiffness check by modal analysis is the first thing an engineering user would ask of a geometry library that already does meshes and physics. Parked behind the surface path because it needs what does not exist yet: a tetrahedral mesher, and assembled elasticity matrices (`continuum.hpp` is types and energies, no FEM).
