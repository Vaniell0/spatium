// spatium.statistics — fitting parameters to a sample, growing a tree from one, and entropy,
// cross-entropy and divergence. The tree it grows is algebra's.
module;
export module spatium.statistics;
import std.compat;
import spatium.core;
import spatium.algebra;
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/statistics/fit.hpp>
#include <spatium/statistics/information.hpp>
#include <spatium/statistics/tree_fit.hpp>
