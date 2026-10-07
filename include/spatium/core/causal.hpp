#pragma once

#include <spatium/_export_macro.hpp>

SPATIUM_EXPORT namespace spatium {

// How two events of a Lorentzian space are related: the sign of the squared
// interval between them, with the signature (-, +, ..., +) the library uses
// (the convention of `Hyperbolic`'s Minkowski form and of `lorentzian_at`).
//   Timelike   a signal slower than light can join them (interval < 0)
//   Null       only light can (interval = 0), up to a tolerance
//   Spacelike  nothing can (interval > 0)
enum class Causal { Timelike, Null, Spacelike };

} // namespace spatium
