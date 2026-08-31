#pragma once

// ── spatium/mesh.hpp ─────────────────────────────────────────
// Focused entry point for mesh-bearing code: the Mesh<Surface>
// data structure, subdivision and LOD, primitive icosahedron/
// tetrahedron generators, topology + geodesic + transport +
// voronoi, mesh quality, the heat geodesic solver (which
// requires Eigen at compile time when SPATIUM_EIGEN is ON),
// SVG / OBJ / STL writers (each pulls Mesh internally), and
// the surface flavours that are themselves built on Mesh
// (ImplicitSurface marching cubes, ParametricSurface tessellator).
//
// Including this header is heavier than <spatium/core.hpp>;
// reach for it only when you actually need mesh structures.
// <spatium/core.hpp> is included for convenience.

#include <spatium/core.hpp>

#include <spatium/spaces/parametric.hpp>
#include <spatium/spaces/implicit.hpp>
#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/mesh/lod.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/topology.hpp>
#include <spatium/mesh/geodesic.hpp>
#include <spatium/mesh/transport.hpp>
#include <spatium/mesh/voronoi.hpp>
#include <spatium/mesh/operations.hpp>
#include <spatium/mesh/quality.hpp>
#include <spatium/io/svg.hpp>
#include <spatium/io/obj.hpp>
#include <spatium/io/stl.hpp>
