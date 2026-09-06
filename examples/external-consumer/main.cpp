// Minimal real consumer of Spatium, built entirely from outside the
// Spatium repository via CMake FetchContent (see this directory's
// CMakeLists.txt and README.md).

#include <spatium/core.hpp>

#include <numbers>
#include <print>

int main() {
    using namespace spatium;

    // Unit 2-sphere (radius 1) embedded in R^3.
    Sphere<2> sphere;

    // North pole and a point on the equator.
    Vec<double, 3> north{0.0, 0.0, 1.0};
    Vec<double, 3> equator_point{1.0, 0.0, 0.0};

    // Great-circle (geodesic) distance between the two points.
    double d = sphere.distance(north, equator_point);

    std::println("Sphere<2> geodesic distance (north pole -> equator): {:.6f}", d);
    std::println("Expected: pi/2 = {:.6f}", std::numbers::pi / 2.0);

    return 0;
}
