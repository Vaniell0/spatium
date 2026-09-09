#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/spaces/parametric.hpp>
#  include <functional>
#endif

SPATIUM_EXPORT namespace spatium {

// offset_surface: compose a new analytic ParametricSurface from `base`,
// pushed out along base's own normal by `thickness` (constant, or a field
// over the parameter domain) -- e.g. icing draped over dough, expressed as
// a function composition, not a tessellated/projected mesh. Same domain
// and periodicity as `base`; every Surface method (project/normal/
// exp_map/geodesics/...) keeps working on the result, since it's still a
// genuine ParametricSurface, not a fixed point sample of one.

template<Scalar T = double>
ParametricSurface<T> offset_surface(const ParametricSurface<T>& base,
                                     std::function<T(T, T)> thickness) {
    return ParametricSurface<T>(
        [=](T u, T v) -> Vec<T, 3> {
            return base.evaluate(u, v) + base.normal_at(u, v) * thickness(u, v);
        },
        base.domain(), base.periodic_u(), base.periodic_v());
}

template<Scalar T = double>
ParametricSurface<T> offset_surface(const ParametricSurface<T>& base, T thickness = T{0}) {
    return offset_surface<T>(base, std::function<T(T, T)>{
        [thickness](T, T) { return thickness; }});
}

} // namespace spatium
