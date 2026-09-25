#include <catch2/catch_test_macros.hpp>
#include <spatium/algebra/noise.hpp>
#include <spatium/io/build.hpp>
#include <spatium/io/field_glsl.hpp>
#include <spatium/io/field_pod.hpp>

#include "donut_scene.hpp"

#include <bit>
#include <cmath>
#include <cstdint>
#include <format>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

using namespace spatium;
namespace bd = spatium::io::build;

namespace {

std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

bool same(const Vec<double, 3>& a, const Vec<double, 3>& b) {
    return bits(a[0]) == bits(b[0]) && bits(a[1]) == bits(b[1]) && bits(a[2]) == bits(b[2]);
}

bool mentions(const std::string& s, const std::string& what) {
    return s.find(what) != std::string::npos;
}

} // namespace

// Lowering has exactly one way to fail, and it has to say where. A
// message that only said "not lowerable" would send whoever reads it back
// to walking the field by hand, which is the thing the report exists to
// save them.
TEST_CASE("lower() refuses a closure and names the op it found",
          "[field][pod]") {
    using F = bd::ScalarField<double>;
    using V = bd::VecField<double>;

    auto leaf = F{1.0} + F{[](double u, double) { return u; }};
    auto r = bd::lower(leaf);
    REQUIRE_FALSE(r.has_value());
    CHECK(mentions(r.error().message, "Opaque"));

    auto opaque_motion = V::opaque_of_time([](double t) { return Vec<double, 3>{t, 0.0, 0.0}; });
    auto rv = bd::lower(opaque_motion);
    REQUIRE_FALSE(rv.has_value());
    CHECK(mentions(rv.error().message, "Opaque"));

    // A factor closure is the blind spot the report once had, so it is
    // the one lowering must not share.
    auto closure_factor = scaled(V::point(), [](double t) { return 1.0 + t; });
    auto rs = bd::lower(closure_factor);
    REQUIRE_FALSE(rs.has_value());
    CHECK(mentions(rs.error().message, "Scale"));

    // And an opaque leaf that came in through the expression door is
    // named with the factor that carried it.
    auto sneaky = rotated(V::point(), F{0.0}, F{[](double t, double) { return t; }}, F{0.0});
    auto rr = bd::lower(sneaky);
    REQUIRE_FALSE(rr.has_value());
    CHECK(mentions(rr.error().message, "Rotate"));
    CHECK(mentions(rr.error().message, "factor 1"));
    CHECK(mentions(rr.error().message, "Opaque"));
}

// Every op the vocabulary has, in one field, so that an op added to `Op`
// and forgotten in the interpreter fails here by name rather than waiting
// for a scene to use it.
TEST_CASE("The interpreter covers every op, bit for bit", "[field][pod]") {
    using F = bd::ScalarField<double>;
    algebra::PerlinNoise n(5);

    auto table = std::make_shared<const std::vector<Vec<double, 3>>>(
        std::vector<Vec<double, 3>>{{1.0, 2.0, 3.0}, {-4.0, 5.5, 0.25}, {7.0, -8.0, 9.5}});
    auto f = min(F::u(), F::v()) + max(F::u(), F{0.5}) * sin(F::u()) -
             cos(F::v()) / F{3.0} + less(F::u(), F::v()) * noise(n, F::u(), F::v(), F{0.25}) +
             F::id() * F{1e-3} + F::origin(1) - F::point(2) * F::hash(7u) +
             sqrt(F::hash(11u) + F{0.5}) + gather(table, F::hash(3u) * F{3.0}, 2) +
             gather(table, F::u(), 0);

    auto lowered = bd::lower(f);
    REQUIRE(lowered.has_value());
    for (auto op = std::uint8_t{0}; op < static_cast<std::uint8_t>(bd::Op::Opaque); ++op) {
        bool present = false;
        for (const auto& p : lowered->ops) present = present || p.code == op;
        INFO("op " << bd::op_name(static_cast<bd::Op>(op)) << " is not exercised");
        CHECK(present);
    }

    // Past both ends of the gather's table on purpose (u from -2 to 3.7
    // over three points), where the clamp is what is being tested.
    for (double u : {-2.0, -0.5, 0.0, 0.25, 1.0, 3.7})
        for (double v : {-1.0, 0.0, 0.25, 2.0, 5.5})
            for (std::uint32_t id : {0u, 1u, 41u, 1999999u}) {
                bd::FieldInputs<double> in{u, v, id, {0.1 * u, -v, 2.0}, {v, u, 0.5 * id}};
                CHECK(bits(bd::interpret(*lowered, in)) == bits(f(in)));
            }
}

// A vector made of three scalar fields lowers and interprets like a
// Rotate's three factors do -- in component order, from one shared pool.
TEST_CASE("A Make of three scalar fields interprets bit for bit", "[field][pod]") {
    using F = bd::ScalarField<double>;
    using V = bd::VecField<double>;
    auto m = V::make(F::hash(1u) * F::origin(0), F::t() + F::hash(2u), F::point(1) * F{2.0}) +
             V::point();
    CHECK(m.is_structural());
    CHECK_FALSE(m.is_placement());   // one component reads the point

    auto pod = bd::lower(m);
    REQUIRE(pod.has_value());
    for (double t : {-1.0, 0.0, 0.7, 3.2})
        for (std::uint32_t id : {0u, 5u, 123456u}) {
            bd::MotionEnv<double> env{{0.3, -1.2, 2.0}, t, {0.5, 0.25, -0.75}, id};
            CHECK(same(bd::interpret(*pod, env), m(env)));
        }

    // Without the point it is a translation per instance, and stays a
    // placement -- which is what lets a particle cloud built from it
    // instance.
    auto moved = V::point() + V::make(F::hash(1u), F::hash(2u), F::t());
    CHECK(moved.is_placement());
}

// The test `docs/gpu-abi-design.md` calls the first step, run over the
// real scene rather than a sample of it: every field on every node the
// report walks, lowered and interpreted beside `eval_into`, compared on
// the bits.
//
// No tolerance, deliberately. A tolerance would pass an interpreter that
// computed `a * b * c` as `a * (b * c)`, and a reordering is exactly what
// an export introduces without anything else noticing. A mismatch here is
// a finding -- it names an op whose meaning is not what the interpreter
// thinks -- not noise to be absorbed.
TEST_CASE("Every structural field in the donut scene interprets bit for bit",
          "[field][pod][donut]") {
    auto boom = donut::load_boom_points(SPATIUM_SOURCE_DIR "/examples/data/boom_points.txt");
    REQUIRE_FALSE(boom.empty());

    bd::Trace<double> scene;
    donut::build_scene(scene, 11000, 35200, boom);

    // Inputs for the scalar fields, which on this scene are all surface
    // thicknesses over (u, v): past both ends of the domain on purpose,
    // where a clamp is the only thing doing any work.
    std::vector<double> us, vs;
    for (int i = 0; i < 48; ++i) us.push_back(-1.0 + 8.5 * i / 47.0);
    for (int j = 0; j < 24; ++j) vs.push_back(-1.0 + 5.5 * j / 23.0);

    // Times for the point fields, dense across the build-up and exact at
    // every edge the scene's factors switch on -- the cube's step and one
    // ulp either side of it included.
    std::vector<double> ts;
    for (int k = 0; k < 64; ++k) ts.push_back(-1.0 + 6.0 * k / 63.0);
    for (double t : {0.12, std::nextafter(0.12, 0.0), std::nextafter(0.12, 1.0),
                     donut::T_EXPLODE, donut::T_HOLD, donut::T_DISSOLVE,
                     donut::T_DONUT_START, donut::T_DONUT_END})
        ts.push_back(t);
    const std::vector<Vec<double, 3>> points{
        {0.0, 0.0, 0.0}, {1.0, -2.0, 0.5}, {-0.3, 0.7, 2.2}, {3.0, 3.0, -1.0}};
    const std::vector<Vec<double, 3>> origins{{0.0, 0.0, 0.0}, {0.0006, -0.0003, 0.0007}};

    std::size_t fields = 0, lowered = 0, evaluations = 0, mismatches = 0;
    std::vector<std::string> refused;

    auto scalar = [&](const bd::ScalarField<double>& f, std::size_t node) {
        ++fields;
        auto pod = bd::lower(f);
        CHECK(pod.has_value() == f.is_structural());
        if (!pod) { refused.push_back(std::format("node {}: {}", node, pod.error().message)); return; }
        ++lowered;
        for (double u : us)
            for (double v : vs) {
                ++evaluations;
                if (bits(bd::interpret(*pod, u, v)) != bits(f(u, v))) ++mismatches;
            }
    };

    auto point = [&](const bd::PointField<double>& f, std::size_t node) {
        ++fields;
        auto pod = bd::lower(f);
        CHECK(pod.has_value() == f.is_structural());
        if (!pod) { refused.push_back(std::format("node {}: {}", node, pod.error().message)); return; }
        ++lowered;
        for (double t : ts)
            for (const auto& p : points)
                for (const auto& o : origins)
                    for (std::uint32_t id : {0u, 1u, 35199u}) {
                        bd::MotionEnv<double> env{p, t, o, id};
                        ++evaluations;
                        if (!same(bd::interpret(*pod, env), f(env))) ++mismatches;
                    }
    };

    // The same slots `field_report()` walks, in the same order. A slot
    // added there and not here would make the count below disagree, which
    // is the point of comparing them.
    for (std::size_t i = 0; i < scene.size(); ++i) {
        const auto& n = scene.node(i);
        if (n.kind == bd::Kind::Offset) scalar(n.thickness, i);
        point(n.transform, i);
        point(n.color_fn, i);
        point(n.emissive_fn, i);
    }

    const auto report = bd::field_report(scene);
    for (const auto& why : refused) UNSCOPED_INFO(why);
    INFO(lowered << " of " << fields << " fields lowered, " << evaluations
                 << " evaluations, " << mismatches << " mismatches");

    CHECK(fields == report.fields);
    CHECK(lowered == report.structural_fields);
    CHECK(refused.size() == report.opaque_fields);
    CHECK(evaluations > 0);
    CHECK(mismatches == 0);
}

// The GLSL emitter refuses what lowering refuses, and a deformation has no
// placement to emit. What it produces is compiled and checked on a device
// by `donut_live --dust-check`; here only the refusals and the data it
// carries, which need no GPU.
TEST_CASE("GLSL emission refuses a closure and a deformation, and carries the tables",
          "[field][pod][glsl]") {
    using F = bd::ScalarField<double>;
    using V = bd::VecField<double>;
    bd::GlslModule m;
    CHECK_FALSE(bd::emit_vector(m, V::opaque_of_time([](double t) { return Vec<double, 3>{t, 0, 0}; }),
                                "a").has_value());
    CHECK_FALSE(bd::emit_placement(m, V::make(F::point(0), F{0.0}, F{0.0}) + V::point(), "b")
                    .has_value());

    algebra::PerlinNoise n(3);
    auto table = std::make_shared<const std::vector<Vec<double, 3>>>(
        std::vector<Vec<double, 3>>{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
    auto motion = V::point() + V::make(noise(n, F::hash(1u), F::t(), F{0.0}),
                                       gather(table, F::hash(2u) * F{2.0}, 1), F{0.0});
    REQUIRE(motion.is_placement());
    REQUIRE(bd::emit_placement(m, motion, "c").has_value());
    CHECK(m.perm.size() == 512);
    CHECK(m.points.size() == 6);
    CHECK(m.code.find("void c_place(") != std::string::npos);
    CHECK(m.code.find("field_noise(0u") != std::string::npos);
}
