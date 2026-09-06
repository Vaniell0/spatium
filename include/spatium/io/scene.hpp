#pragma once

// Scene description format: a small JSON file listing placed geometric
// objects (a shape-kind tag + JSON params + transform + material) and
// an optional camera, built on io/json.hpp.
//
// Open by design, not a closed enum of shape kinds -- the project wants
// to grow in independent add-on chunks rather than by bloating a
// monolithic core. A ShapeRegistry maps a shape_kind string to a
// factory that turns one object's JSON params + placement into a
// ResolvedShape a renderer can ray-test uniformly; a new kind is added
// by registering a factory, anywhere, without editing this file's
// parsing logic. load_scene()/save_scene() never consult the registry
// at all -- an object with an unregistered shape_kind still loads and
// saves correctly, carrying its shape_kind string and raw JSON params
// through untouched. The registry only matters once something (a
// renderer, a validator) calls resolve_shape() on an object.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/groups/so3.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/intersection.hpp>
#  include <spatium/geometry/ray_surface.hpp>
#  include <spatium/io/json.hpp>
#  include <spatium/render/camera.hpp>
#  include <algorithm>
#  include <filesystem>
#  include <fstream>
#  include <functional>
#  include <optional>
#  include <sstream>
#  include <string>
#  include <unordered_map>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io {

// ── Material ─────────────────────────────────────────────────────

template<Scalar T = double>
struct Material {
    // Linear RGB in [0,1] -- matches the convention render/viewer
    // already use for a per-object color (viewer::App::add_mesh's
    // Vec4f color{1,1,1,1}), not the [0,255] convention
    // render::hsv_to_rgb255()/blackbody_to_rgb255() use for a final
    // pixel value.
    Vec<T, 3> base_color{T{0.8}, T{0.8}, T{0.8}};
};

// ── Scene object ─────────────────────────────────────────────────

// One placed object: a shape-kind tag + opaque JSON params (meaningful
// only to whatever factory is registered for that kind), a rigid+scale
// transform, and a material. `orientation` is an axis-angle vector --
// algebra::SO3<T>::AlgebraType, the same representation
// examples/tumbling_body_demo.cpp and primitives_demo.cpp's rotavg
// scene already build rotations from via SO3<T>::exp()/log() -- rather
// than a Quaternion, which this codebase's examples don't otherwise use
// for object orientation.
template<Scalar T = double>
struct SceneObject {
    std::string name;
    std::string shape_kind;
    JsonValue params = JsonValue::object();

    Vec<T, 3> position{};
    Vec<T, 3> orientation{};              // axis-angle; SO3<T>{}.exp(orientation) = rotation matrix
    Vec<T, 3> scale{T{1}, T{1}, T{1}};

    Material<T> material{};
};

// ── Scene ────────────────────────────────────────────────────────

template<Scalar T = double>
struct Scene {
    std::vector<SceneObject<T>> objects;
    std::optional<render::Camera<T>> camera;
};

// ── Shape registry ───────────────────────────────────────────────

// A shape resolved from an object's (shape_kind, params) plus its
// placement, ready to ray-test in world space. Type-erased via
// std::function rather than a closed set of concrete geometry types:
// geometry::RayHit<T> is already the common currency ray_quadric()/
// ray_torus() both return, so any factory -- built-in or
// user-registered -- can hand back a closure in that same shape and a
// renderer never needs to know which concrete shape it is calling.
template<Scalar T = double>
struct ResolvedShape {
    std::function<std::vector<geometry::RayHit<T>>(const geometry::Ray<3, T>&)> ray_hits;
};

template<Scalar T = double>
using ShapeFactory = std::function<Result<ResolvedShape<T>>(const SceneObject<T>&)>;

namespace detail {

// World <-> object-local frame for ray-hit testing. Rotation +
// translation only (an isometry): every built-in factory below bakes
// its object's `scale` directly into the shape's own local parameters
// (radius, half-extents, torus radii) rather than applying it as a ray
// transform, so a hit's t and the rotation of its normal carry over
// unchanged with no inverse-transpose normal correction needed.
template<Scalar T>
geometry::Ray<3, T> to_local_ray(const geometry::Ray<3, T>& ray,
                                  const Vec<T, 3>& position,
                                  const Matrix<T, 3, 3>& rotation) {
    auto Rt = rotation.transpose();
    return {Rt * (ray.origin - position), Rt * ray.direction};
}

template<Scalar T>
geometry::RayHit<T> to_world_hit(const geometry::RayHit<T>& local,
                                  const Vec<T, 3>& position,
                                  const Matrix<T, 3, 3>& rotation) {
    return {local.t, position + rotation * local.point, rotation * local.normal};
}

} // namespace detail

// internal -- backs register_shape_kind()/resolve_shape() below. A
// per-scalar-type singleton (ShapeRegistry<float> and
// ShapeRegistry<double> are independent registries), since factories
// are themselves templated on T.
template<Scalar T = double>
class ShapeRegistry {
public:
    static ShapeRegistry& instance() {
        static ShapeRegistry registry;
        return registry;
    }

    void register_kind(std::string name, ShapeFactory<T> factory) {
        factories_[std::move(name)] = std::move(factory);
    }

    Result<ResolvedShape<T>> resolve(const SceneObject<T>& obj) const {
        auto it = factories_.find(obj.shape_kind);
        if (it == factories_.end())
            return std::unexpected(Error{ErrorCode::InvalidArgument,
                "unknown shape kind: " + obj.shape_kind});
        return it->second(obj);
    }

private:
    std::unordered_map<std::string, ShapeFactory<T>> factories_;
};

// Register a factory for a shape_kind string. New object kinds are
// added purely by calling this -- anywhere, at any point before the
// first resolve_shape() call for that kind -- never by editing this
// file's parsing logic or the registry's internals.
template<Scalar T = double>
void register_shape_kind(std::string name, ShapeFactory<T> factory) {
    ShapeRegistry<T>::instance().register_kind(std::move(name), std::move(factory));
}

// Resolve one scene object into a ray-testable shape, via whatever
// factory is registered for its shape_kind. Result<T>, not assert: an
// unregistered kind (typo, or a scene file written against a build with
// more factories registered than this one) is a real, recoverable
// failure the caller should be able to observe, not an internal
// invariant.
template<Scalar T = double>
Result<ResolvedShape<T>> resolve_shape(const SceneObject<T>& obj) {
    return ShapeRegistry<T>::instance().resolve(obj);
}

// ── Built-in shape kinds ─────────────────────────────────────────

namespace detail {

// "sphere": {"radius": r}. Scale is taken as isotropic (obj.scale[0]
// only) -- a non-uniformly-scaled sphere is really an ellipsoid, which
// this factory doesn't attempt; round-sphere is the case it targets.
template<Scalar T>
Result<ResolvedShape<T>> make_sphere(const SceneObject<T>& obj) {
    T radius = static_cast<T>(obj.params.number_or("radius", T{1}));
    if (radius <= T{0})
        return std::unexpected(Error{ErrorCode::InvalidArgument, "sphere radius must be positive"});

    auto quadric = geometry::Quadric<T>::sphere(radius * obj.scale[0]);
    Vec<T, 3> position = obj.position;
    auto rotation = spatium::algebra::SO3<T>{}.exp(obj.orientation);

    ResolvedShape<T> shape;
    shape.ray_hits = [quadric, position, rotation](const geometry::Ray<3, T>& ray) {
        auto local_ray = to_local_ray(ray, position, rotation);
        auto hits = geometry::ray_quadric(local_ray, quadric);
        for (auto& h : hits) h = to_world_hit(h, position, rotation);
        return hits;
    };
    return shape;
}

// "box": {"half_extents": [hx, hy, hz]}, default [0.5, 0.5, 0.5]. Scale
// multiplies each half-extent component directly -- full anisotropic
// support, unlike sphere/torus: a box has no round part that would make
// non-uniform scale ill-defined.
template<Scalar T>
Result<ResolvedShape<T>> make_box(const SceneObject<T>& obj) {
    Vec<T, 3> half{T{0.5}, T{0.5}, T{0.5}};
    if (auto* he = obj.params.find("half_extents")) {
        auto& arr = he->as_array();
        if (arr.size() != 3)
            return std::unexpected(Error{ErrorCode::InvalidArgument,
                "box half_extents must be a 3-element array"});
        half = Vec<T, 3>{static_cast<T>(arr[0].as_number()),
                         static_cast<T>(arr[1].as_number()),
                         static_cast<T>(arr[2].as_number())};
    }
    Vec<T, 3> scaled_half{half[0] * obj.scale[0], half[1] * obj.scale[1], half[2] * obj.scale[2]};
    if (scaled_half[0] <= T{0} || scaled_half[1] <= T{0} || scaled_half[2] <= T{0})
        return std::unexpected(Error{ErrorCode::InvalidArgument, "box half_extents must be positive"});

    auto local_box = geometry::Box<3, T>::from_center_half_extents(Vec<T, 3>{}, scaled_half);
    Vec<T, 3> position = obj.position;
    auto rotation = spatium::algebra::SO3<T>{}.exp(obj.orientation);

    ResolvedShape<T> shape;
    shape.ray_hits = [local_box, position, rotation](const geometry::Ray<3, T>& ray)
        -> std::vector<geometry::RayHit<T>> {
        auto local_ray = to_local_ray(ray, position, rotation);
        auto params = geometry::intersect_parameters(local_ray, local_box);
        if (!params) return {};
        T t = params->first;

        // geometry::intersect()/intersect_parameters() report the hit
        // interval but not which face it lies on -- recover that here
        // by re-running the same per-axis slab test and keeping the
        // axis whose entry-t matches the winning tmin, rather than
        // reimplementing ray-box intersection itself.
        Vec<T, 3> normal{};
        bool found = false;
        for (std::size_t i = 0; i < 3 && !found; ++i) {
            if (std::abs(local_ray.direction[i]) < epsilon<T>()) continue;
            auto inv_d = T{1} / local_ray.direction[i];
            auto t1 = (local_box.min_corner[i] - local_ray.origin[i]) * inv_d;
            auto t2 = (local_box.max_corner[i] - local_ray.origin[i]) * inv_d;
            T entry = std::min(t1, t2);
            if (std::abs(entry - t) < epsilon<T>() * (std::abs(t) + T{1})) {
                normal[i] = (t1 < t2) ? T{-1} : T{1};
                found = true;
            }
        }
        if (!found) return {};

        geometry::RayHit<T> local_hit{t, local_ray.origin + local_ray.direction * t, normal};
        return {to_world_hit(local_hit, position, rotation)};
    };
    return shape;
}

// "torus": {"major_radius": R, "minor_radius": r}, default {1, 0.25}.
// Scale is taken as isotropic (obj.scale[0]) -- ray_torus()'s quartic
// assumes a circular tube, so a non-uniform scale would silently render
// the wrong (ellipse-tube) shape rather than fail loudly.
template<Scalar T>
Result<ResolvedShape<T>> make_torus(const SceneObject<T>& obj) {
    T major = static_cast<T>(obj.params.number_or("major_radius", T{1}));
    T minor = static_cast<T>(obj.params.number_or("minor_radius", T{0.25}));
    if (major <= T{0} || minor <= T{0} || minor >= major)
        return std::unexpected(Error{ErrorCode::InvalidArgument,
            "torus needs 0 < minor_radius < major_radius"});

    T s = obj.scale[0];
    geometry::Torus<T> torus{.center = {}, .axis = {T{0}, T{0}, T{1}},
                              .major_radius = major * s, .minor_radius = minor * s};
    Vec<T, 3> position = obj.position;
    auto rotation = spatium::algebra::SO3<T>{}.exp(obj.orientation);

    ResolvedShape<T> shape;
    shape.ray_hits = [torus, position, rotation](const geometry::Ray<3, T>& ray) {
        auto local_ray = to_local_ray(ray, position, rotation);
        auto hits = geometry::ray_torus(local_ray, torus);

        // ray_torus() sits on algebra::solve_quartic() (Ferrari's
        // method), which has a real, root-caused robustness bug: its
        // resolvent-cubic root selection (algebra/polynomial.hpp,
        // "Pick a real root of the resolvent" and the fallback right
        // after it) takes the first real root above epsilon, or -- if
        // none clears that bar -- silently falls back to ANY real root,
        // including a negative or near-zero one. `sq_2m =
        // sqrt(abs(2*m))` and the later `beta/(2*sq_2m)` divide don't
        // care whether that m actually corresponds to a valid
        // factorization: a near-zero m blows the division up to +-inf,
        // and a wrong-sign m reconstructs two quadratics that don't
        // multiply back out to the original quartic at all. Confirmed
        // independent of this file: geometry::ray_torus() called
        // directly against a torus at the origin (major_radius 1.4,
        // minor_radius 0.35), varying ONLY the camera's distance along
        // the ring's own symmetric axis (o=(0,-dist,0), d=(0,1,0),
        // always exactly through the tube's true hit points at
        // y = +-1.05/+-1.75) -- distances 3,4,6,8,10,11,11.6 came back
        // exactly right, while 5, 9 and 15 returned a single (t, y) =
        // (inf, inf) "hit" and 7, 12, 20 returned two hits sitting near
        // y=0, nowhere close to the real intersection points. Moving off
        // that one symmetric axis makes it worse, not better: a sweep
        // of viewing elevations away from it found zero mathematically
        // valid hits at every angle except the exactly-symmetric one.
        // test_scene.cpp's own torus test deliberately sits on a
        // verified-good distance for exactly this reason. None of this
        // is introduced by this factory or by the to_local_ray()/
        // to_world_hit() transform above -- see this PR's description
        // for the full repro. The filter below at least stops the
        // spurious (but finite) roots from rendering as visibly wrong
        // hits; it cannot recover the genuine hits solve_quartic() fails
        // to find in the first place, so a torus resolved through this
        // factory can legitimately render sparser than it should for
        // most camera angles. A real fix -- e.g. always taking the
        // largest real resolvent root, or polishing with a Newton step
        // -- belongs in solve_quartic()/ray_torus() itself, out of scope
        // here; solve_quartic() has no other caller in this codebase
        // besides ray_torus(), so the blast radius is contained to
        // torus rendering.
        T R = torus.major_radius, r = torus.minor_radius;
        std::vector<geometry::RayHit<T>> valid;
        valid.reserve(hits.size());
        for (auto& h : hits) {
            T lp2 = h.point.dot(h.point);
            T s = lp2 + R * R - r * r;
            T lhs = s * s;
            T rhs = T{4} * R * R * (h.point[0] * h.point[0] + h.point[1] * h.point[1]);
            T scale = std::max({std::abs(lhs), std::abs(rhs), T{1}});
            if (std::abs(lhs - rhs) < T{1e-6} * scale)
                valid.push_back(to_world_hit(h, position, rotation));
        }
        return valid;
    };
    return shape;
}

// Runs once per translation unit at static-init time. An `inline`
// variable is header-safe (every TU that includes this header
// initializes the same one; ODR keeps them in sync) and registers the
// three built-in kinds into ShapeRegistry<double> -- the scalar type
// every example and test in this codebase actually instantiates Scene
// with. A program using Scene<float> would need its own registrations;
// not done automatically since the registry is genuinely per-T.
inline bool register_builtin_shapes = [] {
    register_shape_kind<double>("sphere", &make_sphere<double>);
    register_shape_kind<double>("box", &make_box<double>);
    register_shape_kind<double>("torus", &make_torus<double>);
    return true;
}();

} // namespace detail

// ── JSON (de)serialization ───────────────────────────────────────

namespace detail {

template<Scalar T>
JsonValue vec3_to_json(const Vec<T, 3>& v) {
    auto a = JsonValue::array();
    a.push_back(static_cast<double>(v[0]));
    a.push_back(static_cast<double>(v[1]));
    a.push_back(static_cast<double>(v[2]));
    return a;
}

template<Scalar T>
Result<Vec<T, 3>> vec3_from_json(const JsonValue& v, const char* field) {
    if (!v.is_array() || v.as_array().size() != 3)
        return std::unexpected(Error{ErrorCode::ParseError,
            std::string(field) + " must be a 3-element array"});
    auto& a = v.as_array();
    return Vec<T, 3>{static_cast<T>(a[0].as_number()),
                     static_cast<T>(a[1].as_number()),
                     static_cast<T>(a[2].as_number())};
}

template<Scalar T>
JsonValue object_to_json(const SceneObject<T>& obj) {
    auto j = JsonValue::object();
    j.set("name", obj.name);
    j.set("shape", obj.shape_kind);
    j.set("params", obj.params);
    j.set("position", vec3_to_json(obj.position));
    j.set("orientation", vec3_to_json(obj.orientation));
    j.set("scale", vec3_to_json(obj.scale));
    j.set("base_color", vec3_to_json(obj.material.base_color));
    return j;
}

template<Scalar T>
Result<SceneObject<T>> object_from_json(const JsonValue& j) {
    if (!j.is_object())
        return std::unexpected(Error{ErrorCode::ParseError, "scene object must be a JSON object"});

    SceneObject<T> obj;
    obj.name = j.string_or("name", "");
    obj.shape_kind = j.string_or("shape", "");
    if (obj.shape_kind.empty())
        return std::unexpected(Error{ErrorCode::ParseError, "scene object missing \"shape\""});
    if (auto* p = j.find("params")) obj.params = *p;

    obj.scale = Vec<T, 3>{T{1}, T{1}, T{1}};
    if (auto* p = j.find("position")) {
        auto v = vec3_from_json<T>(*p, "position");
        if (!v) return std::unexpected(v.error());
        obj.position = *v;
    }
    if (auto* p = j.find("orientation")) {
        auto v = vec3_from_json<T>(*p, "orientation");
        if (!v) return std::unexpected(v.error());
        obj.orientation = *v;
    }
    if (auto* p = j.find("scale")) {
        auto v = vec3_from_json<T>(*p, "scale");
        if (!v) return std::unexpected(v.error());
        obj.scale = *v;
    }
    if (auto* p = j.find("base_color")) {
        auto v = vec3_from_json<T>(*p, "base_color");
        if (!v) return std::unexpected(v.error());
        obj.material.base_color = *v;
    }
    return obj;
}

template<Scalar T>
JsonValue camera_to_json(const render::Camera<T>& cam) {
    auto j = JsonValue::object();
    j.set("position", vec3_to_json(cam.position));
    j.set("target", vec3_to_json(cam.target));
    j.set("up", vec3_to_json(cam.up));
    j.set("fov_deg", static_cast<double>(cam.fov_deg));
    return j;
}

template<Scalar T>
Result<render::Camera<T>> camera_from_json(const JsonValue& j) {
    render::Camera<T> cam;
    if (auto* p = j.find("position")) {
        auto v = vec3_from_json<T>(*p, "camera.position");
        if (!v) return std::unexpected(v.error());
        cam.position = *v;
    }
    if (auto* p = j.find("target")) {
        auto v = vec3_from_json<T>(*p, "camera.target");
        if (!v) return std::unexpected(v.error());
        cam.target = *v;
    }
    if (auto* p = j.find("up")) {
        auto v = vec3_from_json<T>(*p, "camera.up");
        if (!v) return std::unexpected(v.error());
        cam.up = *v;
    }
    cam.fov_deg = static_cast<T>(j.number_or("fov_deg", static_cast<double>(cam.fov_deg)));
    return cam;
}

} // namespace detail

template<Scalar T = double>
JsonValue scene_to_json(const Scene<T>& scene) {
    auto j = JsonValue::object();
    auto objects = JsonValue::array();
    for (auto& obj : scene.objects) objects.push_back(detail::object_to_json(obj));
    j.set("objects", std::move(objects));
    if (scene.camera) j.set("camera", detail::camera_to_json<T>(*scene.camera));
    return j;
}

template<Scalar T = double>
Result<Scene<T>> scene_from_json(const JsonValue& j) {
    if (!j.is_object())
        return std::unexpected(Error{ErrorCode::ParseError, "scene must be a JSON object"});

    Scene<T> scene;
    if (auto* objects = j.find("objects")) {
        for (auto& item : objects->as_array()) {
            auto obj = detail::object_from_json<T>(item);
            if (!obj) return std::unexpected(obj.error());
            scene.objects.push_back(std::move(*obj));
        }
    }
    if (auto* cam = j.find("camera")) {
        auto c = detail::camera_from_json<T>(*cam);
        if (!c) return std::unexpected(c.error());
        scene.camera = *c;
    }
    return scene;
}

// ── File I/O ──────────────────────────────────────────────────────

template<Scalar T = double>
Result<Scene<T>> load_scene(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file)
        return std::unexpected(Error{ErrorCode::InvalidArgument, "cannot open file: " + path.string()});

    std::ostringstream ss;
    ss << file.rdbuf();

    auto parsed = parse(ss.str());
    if (!parsed) return std::unexpected(parsed.error());
    return scene_from_json<T>(*parsed);
}

template<Scalar T = double>
Result<void> save_scene(const std::filesystem::path& path, const Scene<T>& scene) {
    std::ofstream file(path);
    if (!file)
        return std::unexpected(Error{ErrorCode::InvalidArgument,
            "cannot open file for writing: " + path.string()});

    file << scene_to_json(scene).dump();
    return {};
}

} // namespace spatium::io
