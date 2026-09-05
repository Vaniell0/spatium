spatium_add_example(NAME geometry_demo)
spatium_add_example(NAME showcase)
spatium_add_example(NAME sets_demo)
spatium_add_example(NAME spd_relationship_demo)

# Live terminal animation, header-only (spatium::render::terminal_canvas) --
# -O3 for the tight per-frame theta/phi trig loop, same reasoning as the
# offline raytracers below (real per-frame cost, not meant to be stepped
# through in a debugger).
spatium_add_example(NAME terminal_donut_demo OPTIONS -O3 -march=native)

# CPU-only analytical renderer — no Vulkan required.
#
# -O3 -march=native regardless of the project's overall CMAKE_BUILD_TYPE:
# this is a per-pixel CPU ray tracer meant to produce final output images,
# not to be stepped through in a debugger, and its cost is real (a full
# gallery render was 2m22s before the ray_parametric.hpp fixes and still
# benefits materially from optimization on top of those) -- 2026-08-25.
spatium_add_example(NAME parametric_analytical_demo OPTIONS -O3 -march=native)

# Schwarzschild geodesic ray-tracer — null-geodesic backward integration.
# -O3 -march=native: same reasoning as parametric_analytical_demo above.
spatium_add_example(NAME blackhole_demo OPTIONS -O3 -march=native)

# Rigorous-GR Schwarzschild ray-tracer -- full 4-coordinate geodesic
# integration via physics/relativity/ (Dual<T>-exact Christoffel
# symbols), real tilted/lensed accretion disk, camera choreography
# (frame sequence). Built alongside blackhole_demo.cpp, not replacing
# it -- see this file's header comment for why. -O3 -march=native: same
# reasoning as blackhole_demo above, multiplied by per-ray numerical
# integration instead of a closed-form ODE.
spatium_add_example(NAME blackhole_gr_demo OPTIONS -O3 -march=native)

# Ray-bundle geodesic diagram -- not a shaded render, so no AA/-O3
# tuning concerns on the level of the raytracers above, but still
# benefits from -O3 given the geodesic integration itself.
spatium_add_example(NAME geodesic_curvature_grid_demo OPTIONS -O3 -march=native)

# Ellis/Morris-Thorne wormhole ray-tracer -- camera flythrough (frame
# sequence). -O3 -march=native: same reasoning as blackhole_demo above,
# multiplied by rendering N_FRAMES full images per run.
spatium_add_example(NAME wormhole_demo OPTIONS -O3 -march=native)

# Dzhanibekov-effect rigid-body tumble (LGVI) -- camera-static frame
# sequence. -O3 -march=native: same reasoning as blackhole_demo above.
# Reference example for the render/ engine primitives (Camera,
# parallel_for_rows, supersample_pixel, write_png_rgb) -- see file header.
spatium_add_example(NAME tumbling_body_demo OPTIONS -O3 -march=native)

# Geodesic Voronoi regions on a marching-cubes-tessellated implicit
# "planet" -- BVH<Triangle3> ray_cast against a real mesh, not an
# analytical primitive, so -O3 -march=native matters for the per-pixel
# BVH traversal same as the raytracers above.
spatium_add_example(NAME geodesic_procgen_demo OPTIONS -O3 -march=native)

# Hyperbolic<3> sphere-tracing -- per-pixel exp_map() + distance() calls in
# the marching loop, same per-pixel-cost reasoning as the raytracers above.
spatium_add_example(NAME hyperbolic_tessellation_demo OPTIONS -O3 -march=native)

# Wave equation as cellular-automaton on a (θ, φ) sphere grid.
spatium_add_example(NAME wave_ca_demo)

# Collatz orbits as 4D parametric surface (height view + R^4->R^3 stereographic).
spatium_add_example(NAME collatz_demo)

# Visual cloth demos (cloth_on_sphere/_torus, cloth_visual) stay removed
# pending a real viewer-facing implicit-contact pipeline — explicit
# position-based contact couldn't keep a stiff cloth with heavy
# overhanging corners on an analytical obstacle; tunnelling and
# stretching both surface as failure modes. Sources are in git history
# (commit 1b3ce2c62). cloth_sphere_probe.cpp — the headless diagnostic --
# is back, rebuilt on ipc-toolkit's Newton + log-barrier + CCD-filtered
# line search instead of XPBD; gated on SPATIUM_IPC_TOOLKIT since it
# links ipc::toolkit directly (not through spatium_sdk — see
# CMakeLists.txt's ipc-toolkit section for why).
if(SPATIUM_IPC_TOOLKIT)
    spatium_add_example(NAME cloth_sphere_probe LIBS spatium_core ipc::toolkit)
endif()

if(SPATIUM_HAS_VIEWER)
    spatium_add_example(NAME atom_demo LIBS spatium_viewer)

    # Unified Vulkan demo — dispatches between scenes:
    #   --scene primitives (default)   all primitives + BVH raycast
    #   --scene torus                  Clifford torus S³ → R³
    #   --scene klein                  Klein bottle R⁴ → R³, animated
    spatium_add_example(NAME primitives_demo LIBS spatium_viewer)

    # Burning Ship fractal — Vulkan window with a fragment-shader
    # iterator and Jacobian-based distance estimation.  Deliberately
    # standalone (no spatium_viewer dep) so it can later double as a
    # template for other GLSL fractal demos.
    spatium_add_example(NAME burning_ship_demo
        LIBS spatium_core Vulkan::Vulkan glfw ${SHADERC_LIBRARIES}
        INCLUDE_DIRS ${SHADERC_INCLUDE_DIRS})
endif()
