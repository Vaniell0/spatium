#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/io/scene.hpp>

#include <filesystem>
#include <fstream>
#include <numbers>

using namespace spatium;
using namespace spatium::io;
using Catch::Matchers::WithinAbs;

namespace {

Scene<double> make_sample_scene() {
    Scene<double> scene;

    SceneObject<double> sphere;
    sphere.name = "ball";
    sphere.shape_kind = "sphere";
    sphere.params.set("radius", 2.5);
    sphere.position = Vec<double, 3>{1.0, 2.0, 3.0};
    sphere.orientation = Vec<double, 3>{0.0, 0.0, 0.0};
    sphere.scale = Vec<double, 3>{1.0, 1.0, 1.0};
    sphere.material.base_color = Vec<double, 3>{0.9, 0.1, 0.1};
    scene.objects.push_back(sphere);

    SceneObject<double> box;
    box.name = "crate";
    box.shape_kind = "box";
    box.params.set("half_extents", [] {
        auto a = JsonValue::array();
        a.push_back(0.5); a.push_back(1.0); a.push_back(1.5);
        return a;
    }());
    box.position = Vec<double, 3>{-3.0, 0.0, 0.5};
    box.orientation = Vec<double, 3>{0.0, std::numbers::pi / 4.0, 0.0};
    box.material.base_color = Vec<double, 3>{0.2, 0.6, 0.9};
    scene.objects.push_back(box);

    SceneObject<double> torus;
    torus.name = "ring";
    torus.shape_kind = "torus";
    torus.params.set("major_radius", 1.5);
    torus.params.set("minor_radius", 0.3);
    torus.position = Vec<double, 3>{0.0, 5.0, 0.0};
    scene.objects.push_back(torus);

    scene.camera = render::Camera<double>{
        .position = {0.0, -10.0, 2.0}, .target = {0.0, 0.0, 0.0},
        .up = {0.0, 0.0, 1.0}, .fov_deg = 40.0};

    return scene;
}

} // namespace

// ── Load/save round-trip ────────────────────────────────────────

TEST_CASE("Scene save/load round-trip preserves objects and camera", "[scene]") {
    auto scene = make_sample_scene();
    auto path = std::filesystem::temp_directory_path() / "spatium_test_scene.json";

    auto saved = save_scene(path, scene);
    REQUIRE(saved.has_value());

    auto loaded = load_scene<double>(path);
    REQUIRE(loaded.has_value());

    auto& s = *loaded;
    REQUIRE(s.objects.size() == 3);

    auto& sphere = s.objects[0];
    CHECK(sphere.name == "ball");
    CHECK(sphere.shape_kind == "sphere");
    CHECK_THAT(sphere.params.number_or("radius", -1.0), WithinAbs(2.5, 1e-12));
    CHECK_THAT(sphere.position[0], WithinAbs(1.0, 1e-12));
    CHECK_THAT(sphere.position[1], WithinAbs(2.0, 1e-12));
    CHECK_THAT(sphere.position[2], WithinAbs(3.0, 1e-12));
    CHECK_THAT(sphere.material.base_color[0], WithinAbs(0.9, 1e-12));
    CHECK_THAT(sphere.material.base_color[1], WithinAbs(0.1, 1e-12));

    auto& box = s.objects[1];
    CHECK(box.name == "crate");
    CHECK(box.shape_kind == "box");
    auto* half = box.params.find("half_extents");
    REQUIRE(half);
    REQUIRE(half->as_array().size() == 3);
    CHECK_THAT(half->as_array()[1].as_number(), WithinAbs(1.0, 1e-12));
    CHECK_THAT(box.orientation[1], WithinAbs(std::numbers::pi / 4.0, 1e-12));

    auto& torus = s.objects[2];
    CHECK(torus.name == "ring");
    CHECK(torus.shape_kind == "torus");
    CHECK_THAT(torus.params.number_or("major_radius", -1.0), WithinAbs(1.5, 1e-12));
    CHECK_THAT(torus.params.number_or("minor_radius", -1.0), WithinAbs(0.3, 1e-12));

    REQUIRE(s.camera.has_value());
    CHECK_THAT(s.camera->position[1], WithinAbs(-10.0, 1e-12));
    CHECK_THAT(s.camera->fov_deg, WithinAbs(40.0, 1e-12));

    std::filesystem::remove(path);
}

TEST_CASE("Scene load fails cleanly on a missing file", "[scene]") {
    auto result = load_scene<double>("/nonexistent/path/spatium_no_such_scene.json");
    CHECK(!result.has_value());
    CHECK(result.error().code == ErrorCode::InvalidArgument);
}

TEST_CASE("Scene load fails on malformed JSON", "[scene]") {
    auto path = std::filesystem::temp_directory_path() / "spatium_bad_scene.json";
    {
        std::ofstream bad(path);
        bad << "{ not valid json ";
    }
    auto result = load_scene<double>(path);
    CHECK(!result.has_value());
    std::filesystem::remove(path);
}

TEST_CASE("Scene object missing a shape kind fails to parse", "[scene]") {
    auto j = parse(R"({"objects": [{"name": "oops"}]})");
    REQUIRE(j);
    auto scene = scene_from_json<double>(*j);
    CHECK(!scene.has_value());
    CHECK(scene.error().code == ErrorCode::ParseError);
}

// ── Registry: built-ins ──────────────────────────────────────────

TEST_CASE("resolve_shape hits a registered sphere", "[scene]") {
    SceneObject<double> obj;
    obj.shape_kind = "sphere";
    obj.params.set("radius", 2.0);
    obj.position = Vec<double, 3>{5.0, 0.0, 0.0};

    auto resolved = resolve_shape(obj);
    REQUIRE(resolved.has_value());

    // Ray straight along +X through the sphere's center.
    geometry::Ray<3, double> ray{Vec<double, 3>{-5.0, 0.0, 0.0}, Vec<double, 3>{1.0, 0.0, 0.0}};
    auto hits = resolved->ray_hits(ray);
    REQUIRE(hits.size() == 2);
    CHECK_THAT(hits[0].point[0], WithinAbs(3.0, 1e-9));  // 5 - 2
    CHECK_THAT(hits[1].point[0], WithinAbs(7.0, 1e-9));  // 5 + 2

    // A ray that passes well clear of the sphere misses entirely.
    geometry::Ray<3, double> miss{Vec<double, 3>{-5.0, 10.0, 0.0}, Vec<double, 3>{1.0, 0.0, 0.0}};
    CHECK(resolved->ray_hits(miss).empty());
}

TEST_CASE("resolve_shape hits a registered box at identity orientation", "[scene]") {
    SceneObject<double> obj;
    obj.shape_kind = "box";
    obj.params.set("half_extents", [] {
        auto a = JsonValue::array();
        a.push_back(1.0); a.push_back(1.0); a.push_back(1.0);
        return a;
    }());
    obj.position = Vec<double, 3>{0.0, 0.0, 0.0};

    auto resolved = resolve_shape(obj);
    REQUIRE(resolved.has_value());

    geometry::Ray<3, double> ray{Vec<double, 3>{-5.0, 0.0, 0.0}, Vec<double, 3>{1.0, 0.0, 0.0}};
    auto hits = resolved->ray_hits(ray);
    REQUIRE(hits.size() == 1);
    CHECK_THAT(hits[0].point[0], WithinAbs(-1.0, 1e-9));
    CHECK_THAT(hits[0].normal[0], WithinAbs(-1.0, 1e-9));
}

TEST_CASE("resolve_shape's box hit-test actually applies the object's orientation", "[scene]") {
    SceneObject<double> obj;
    obj.shape_kind = "box";
    obj.params.set("half_extents", [] {
        auto a = JsonValue::array();
        a.push_back(1.0); a.push_back(2.0); a.push_back(3.0);
        return a;
    }());
    obj.position = Vec<double, 3>{0.0, 0.0, 0.0};
    // 90-degree rotation about Z: swaps the box's local X/Y axes into
    // world Y/-X, so a ray along -X now reaches the half_extents[1]=2
    // face (originally the Y face) instead of the half_extents[0]=1 one.
    obj.orientation = Vec<double, 3>{0.0, 0.0, std::numbers::pi / 2.0};

    auto resolved = resolve_shape(obj);
    REQUIRE(resolved.has_value());

    geometry::Ray<3, double> ray{Vec<double, 3>{10.0, 0.0, 0.0}, Vec<double, 3>{-1.0, 0.0, 0.0}};
    auto hits = resolved->ray_hits(ray);
    REQUIRE(hits.size() == 1);
    CHECK_THAT(hits[0].point[0], WithinAbs(2.0, 1e-9));
    CHECK_THAT(hits[0].normal[0], WithinAbs(1.0, 1e-9));
    CHECK_THAT(hits[0].normal[1], WithinAbs(0.0, 1e-9));
}

// geometry::ray_torus()'s quartic solver turned out to be unreliable for
// most oblique rays -- confirmed by calling it directly, independent of
// this registry entirely (see io/scene.hpp's make_torus() comment and
// this PR's description for the repro: a sweep at a range of viewing
// angles found essentially every ray except the exactly-symmetric one
// below either misses a genuine hit or returns a root that fails the
// torus's own implicit equation). This test deliberately stays on that
// one reliable, symmetric configuration -- a ray in the ring's own
// plane through its center -- so it exercises real, currently-working
// behavior rather than asserting on a known upstream limitation.
TEST_CASE("resolve_shape hits a registered torus on its symmetric axis", "[scene]") {
    SceneObject<double> obj;
    obj.shape_kind = "torus";
    obj.params.set("major_radius", 1.4);
    obj.params.set("minor_radius", 0.35);
    obj.position = Vec<double, 3>{0.0, 0.0, 0.0};

    auto resolved = resolve_shape(obj);
    REQUIRE(resolved.has_value());

    geometry::Ray<3, double> ray{Vec<double, 3>{0.0, -6.0, 0.0}, Vec<double, 3>{0.0, 1.0, 0.0}};
    auto hits = resolved->ray_hits(ray);
    REQUIRE(hits.size() == 4);
    CHECK_THAT(hits[0].point[1], WithinAbs(-1.75, 1e-9)); // -(R+r)
    CHECK_THAT(hits[1].point[1], WithinAbs(-1.05, 1e-9)); // -(R-r)
    CHECK_THAT(hits[2].point[1], WithinAbs(1.05, 1e-9));  //  (R-r)
    CHECK_THAT(hits[3].point[1], WithinAbs(1.75, 1e-9));  //  (R+r)
}

TEST_CASE("resolve_shape returns an error for an unknown shape kind", "[scene]") {
    SceneObject<double> obj;
    obj.shape_kind = "does_not_exist_kind";

    auto resolved = resolve_shape(obj);
    REQUIRE(!resolved.has_value());
    CHECK(resolved.error().code == ErrorCode::InvalidArgument);
}

// ── Registry: extensibility with a custom kind ────────────────────

namespace {

// A toy shape kind not built into scene.hpp at all: a flat disc lying
// in the object's local XY plane. Proves the registry lets an entirely
// new kind plug in from outside this file without touching scene.hpp's
// parsing or the built-in factories.
Result<ResolvedShape<double>> make_toy_disc(const SceneObject<double>& obj) {
    double radius = obj.params.number_or("radius", 1.0);
    if (radius <= 0.0)
        return std::unexpected(Error{ErrorCode::InvalidArgument, "toy_disc radius must be positive"});

    Vec<double, 3> center = obj.position;
    ResolvedShape<double> shape;
    shape.ray_hits = [center, radius](const geometry::Ray<3, double>& ray) {
        std::vector<geometry::RayHit<double>> hits;
        // Intersect with the z = center.z plane, then check the radius.
        if (std::abs(ray.direction[2]) < 1e-12) return hits;
        double t = (center[2] - ray.origin[2]) / ray.direction[2];
        if (t < 0.0) return hits;
        auto p = ray.origin + ray.direction * t;
        auto d = p - center;
        if (d[0] * d[0] + d[1] * d[1] > radius * radius) return hits;
        hits.push_back({t, p, Vec<double, 3>{0.0, 0.0, 1.0}});
        return hits;
    };
    return shape;
}

} // namespace

TEST_CASE("A custom-registered shape kind round-trips through save/load and resolves",
          "[scene]") {
    register_shape_kind<double>("toy_disc", &make_toy_disc);

    Scene<double> scene;
    SceneObject<double> obj;
    obj.name = "coaster";
    obj.shape_kind = "toy_disc";
    obj.params.set("radius", 1.25);
    obj.position = Vec<double, 3>{2.0, 2.0, 4.0};
    scene.objects.push_back(obj);

    auto path = std::filesystem::temp_directory_path() / "spatium_test_custom_scene.json";
    REQUIRE(save_scene(path, scene).has_value());

    auto loaded = load_scene<double>(path);
    REQUIRE(loaded.has_value());
    REQUIRE(loaded->objects.size() == 1);

    auto& loaded_obj = loaded->objects[0];
    CHECK(loaded_obj.name == "coaster");
    CHECK(loaded_obj.shape_kind == "toy_disc");
    CHECK_THAT(loaded_obj.params.number_or("radius", -1.0), WithinAbs(1.25, 1e-12));

    // scene.hpp never knew "toy_disc" existed while parsing the file
    // above -- resolving it now proves the registered factory (from
    // this test file, not scene.hpp) is what actually answers for it.
    auto resolved = resolve_shape(loaded_obj);
    REQUIRE(resolved.has_value());

    geometry::Ray<3, double> hit_ray{Vec<double, 3>{2.0, 2.0, 0.0}, Vec<double, 3>{0.0, 0.0, 1.0}};
    auto hits = resolved->ray_hits(hit_ray);
    REQUIRE(hits.size() == 1);
    CHECK_THAT(hits[0].point[2], WithinAbs(4.0, 1e-9));

    geometry::Ray<3, double> miss_ray{Vec<double, 3>{10.0, 10.0, 0.0}, Vec<double, 3>{0.0, 0.0, 1.0}};
    CHECK(resolved->ray_hits(miss_ray).empty());

    std::filesystem::remove(path);
}
