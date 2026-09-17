module;
#include <cassert>  // assert is a macro — must come from a real header, not import std
export module spatium.core:up_to;
import std.compat;
import :error;
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/core/up_to.hpp>
