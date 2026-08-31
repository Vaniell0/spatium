#pragma once

// ── spatium/spatium.hpp ──────────────────────────────────────
// Convenience umbrella that pulls in every public header.
//
// DEPRECATED for new code: prefer the focused entry points
//
//     #include <spatium/core.hpp>     // ~80% of usage
//     #include <spatium/mesh.hpp>     // mesh + IO writers + heat method
//     #include <spatium/physics.hpp>  // atoms, mechanics, integrators
//
// They keep compile times tighter and document intent. This file
// re-exports the three so existing translation units that include
// <spatium/spatium.hpp> keep compiling unchanged.
//
// NOTE: by default this header still injects sub-namespaces
// (geometry, mesh, io, spatial, physics) into ::spatium so that
// `using namespace spatium;` exposes free functions like tri(),
// seg(), distance(), intersect(), etc. without further qualifiers.
// Define SPATIUM_NO_USING before including to opt out and write
// fully-qualified names — spatium::geometry::tri(...) etc. The
// focused headers (core.hpp / mesh.hpp / physics.hpp) never inject
// using-directives; this convenience belongs to the umbrella only
// and is itself a candidate for removal in a future major release.

#include <spatium/core.hpp>
#include <spatium/mesh.hpp>
#include <spatium/physics.hpp>

#ifndef SPATIUM_NO_USING
namespace spatium {
    using namespace geometry;
    using namespace mesh;
    using namespace io;
    using namespace spatial;
    using namespace physics;
}
#endif
