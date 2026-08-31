#pragma once

// ── spatium/physics.hpp ──────────────────────────────────────
// Focused entry point for the physics modules: orbital solvers,
// chemistry/element data, the atom + Bohr models, atom SVG
// rendering, the mechanics submodule (body/state/integrator/
// contact/continuum/manifold/variational/XPBD/LGVI/Lie integrator
// /narrow-phase), and the relativity submodule (metric-agnostic
// geodesic integration via Dual<T>-exact Christoffel symbols,
// Schwarzschild as the one metric shipped today). <spatium/core.hpp>
// is pulled in for the algebra and space primitives the physics
// layer is built on; mesh-bearing pieces (atom_svg uses io/svg.hpp
// which needs Mesh) come from <spatium/mesh.hpp>.

#include <spatium/core.hpp>
#include <spatium/mesh.hpp>

#include <spatium/physics/atomic/orbital.hpp>
#include <spatium/physics/elements.hpp>
#include <spatium/physics/atomic/atom_model.hpp>
#include <spatium/physics/atomic/atom_svg.hpp>
#include <spatium/physics/atomic/bohr_model.hpp>
#include <spatium/physics/mechanics/body.hpp>
#include <spatium/physics/mechanics/state.hpp>
#include <spatium/physics/mechanics/force.hpp>
#include <spatium/physics/mechanics/units.hpp>
#include <spatium/physics/mechanics/integrator.hpp>
#include <spatium/physics/mechanics/symplectic.hpp>
#include <spatium/physics/mechanics/lie_integrator.hpp>
#include <spatium/physics/mechanics/lgvi.hpp>
#include <spatium/physics/mechanics/variational.hpp>
#include <spatium/physics/mechanics/manifold_body.hpp>
#include <spatium/physics/mechanics/contact.hpp>
#include <spatium/physics/mechanics/continuum.hpp>
#include <spatium/physics/mechanics/narrow_phase.hpp>
#include <spatium/physics/mechanics/xpbd.hpp>
#include <spatium/physics/relativity/schwarzschild.hpp>
#include <spatium/physics/relativity/kerr.hpp>
#include <spatium/physics/relativity/geodesic.hpp>
#include <spatium/physics/relativity/accretion_disk.hpp>
