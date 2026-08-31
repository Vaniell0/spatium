#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/topology.hpp>  // no_vertex, Surface
#  include <cstdint>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

template<Surface S>
struct DistanceField {
    using ScalarT = typename S::ScalarType;

    std::vector<ScalarT> distances;
    std::vector<uint32_t> predecessors;  // no_vertex = source or unreachable
};

template<Surface S>
struct GeodesicPath {
    using ScalarT = typename S::ScalarType;

    std::vector<uint32_t> vertices;
    ScalarT total_length{};
};

} // namespace spatium::mesh
