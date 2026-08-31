// spatium.physics — atom_model, orbital, elements, atom_svg, bohr_model,
// mechanics (units/state/body/force/integrators/Lie-group/variational/LGVI/
// contact), relativity (Schwarzschild/Kerr geodesic integration, accretion
// disk). Single partition for the same reason Spaces/Geometry/Mesh are:
// static_assert and implicit Vec instantiations inside atom_model/bohr_model
// force one TU.
//
// NOTE: atom_model.hpp and bohr_model.hpp include mesh/primitives.hpp, which
// stays HEADER-ONLY (Phase 4 pitfall — mesh cycle with spaces). That's OK:
// primitives' symbols are header-only, this module's use of them is local,
// and the consumer doesn't `import` primitives through spatium.physics.
// Module consumers that need primitives directly continue `#include`ing.
//
// The elements(...) lookup tables live in src/physics/elements.cpp, which
// becomes a module implementation unit (see elements_impl.cpp).
module;
export module spatium.physics;
import std.compat;
import spatium.core;
import spatium.algebra;
import spatium.spaces;
import spatium.geometry;   // Torus / Quadric — needed by narrow_phase.
import spatium.mesh;
import spatium.io;
#define SPATIUM_BUILDING_MODULE 1
// --- spatium::physics::atomic (visualization subsystem) ---
#include <spatium/physics/elements.hpp>
#include <spatium/physics/atomic/orbital.hpp>
#include <spatium/physics/atomic/atom_palette.hpp>
// atom_model and bohr_model pull mesh/primitives.hpp which is header-only.
// When included here with SPATIUM_BUILDING_MODULE defined, primitives' symbols
// end up attached to this module — which is fine because no other module
// attempts to export them.
#include <spatium/mesh/primitives.hpp>
#include <spatium/physics/atomic/atom_model.hpp>
#include <spatium/physics/atomic/atom_svg.hpp>
#include <spatium/physics/atomic/bohr_model.hpp>
// --- spatium::physics::mechanics (foundations) ---
#include <spatium/physics/mechanics/units.hpp>
#include <spatium/physics/mechanics/state.hpp>
#include <spatium/physics/mechanics/body.hpp>
#include <spatium/physics/mechanics/force.hpp>
#include <spatium/physics/mechanics/integrator.hpp>
// Lie-group integrators (RKMK precursor: Lie-Euler / Lie-midpoint).
// Note: mesh/dec.hpp stays header-only because <Eigen/Sparse> would clash with
// the SSE intrinsics already in spatium.algebra :vec_simd BMI.
#include <spatium/physics/mechanics/lie_integrator.hpp>
#include <spatium/physics/mechanics/manifold_body.hpp>
// Variational integrator (Marsden-West) for separable Lagrangians on flat
// Euclidean, then SymplecticManifold concept + CotangentBundle<M>.
#include <spatium/physics/mechanics/variational.hpp>
#include <spatium/physics/mechanics/symplectic.hpp>
// LGVI for SO(3) free rigid body + continuum scaffolding.
#include <spatium/physics/mechanics/lgvi.hpp>
#include <spatium/physics/mechanics/continuum.hpp>
// IPC barrier potential + XPBD constraint solver.
#include <spatium/physics/mechanics/contact.hpp>
#include <spatium/physics/mechanics/xpbd.hpp>
// narrow_phase wraps geometry/ray_* as point-to-surface ContactQuery;
// it lives in physics because it bridges the IPC barrier with the
// existing analytical surface stack.
#include <spatium/physics/mechanics/narrow_phase.hpp>
// --- spatium::physics::relativity (metric-agnostic geodesic integration) ---
// Independent headers (no cross-includes among the four); schwarzschild.hpp
// and kerr.hpp supply concrete metrics, geodesic.hpp integrates any
// Dual<T>-substitutable metric callable via algebra's :dual/:ode/:linear_solve
// (both already re-exported through the `import spatium.algebra;` above).
#include <spatium/physics/relativity/schwarzschild.hpp>
#include <spatium/physics/relativity/kerr.hpp>
#include <spatium/physics/relativity/geodesic.hpp>
#include <spatium/physics/relativity/accretion_disk.hpp>
