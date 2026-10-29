// spatium.discrete — finite sets (pull geometry for set ops) and the rules of chess on the
// 8x8 board and on boxes of any number of axes (pull algebra's groups).
module;
export module spatium.discrete;
import std.compat;
import spatium.core;
import spatium.algebra;
import spatium.geometry;
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/discrete/finite_set.hpp>
#include <spatium/discrete/bitboard.hpp>
#include <spatium/discrete/attacks.hpp>
#include <spatium/discrete/lattice.hpp>
#include <spatium/discrete/move_features.hpp>
#include <spatium/discrete/position.hpp>
