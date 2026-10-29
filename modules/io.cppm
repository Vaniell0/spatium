// spatium.io — table, svg, obj, stl, chess. svg.hpp pulls mesh; obj/stl pull spaces;
// chess.hpp pulls discrete.
// All implementations are inline; no .cpp companion.
module;
export module spatium.io;
import std.compat;
import spatium.core;
import spatium.algebra;
import spatium.spaces;
import spatium.mesh;
import spatium.discrete;
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/io/table.hpp>
#include <spatium/io/svg.hpp>
#include <spatium/io/obj.hpp>
#include <spatium/io/stl.hpp>
#include <spatium/io/chess.hpp>
