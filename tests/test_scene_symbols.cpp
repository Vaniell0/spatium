// Named numbers in a scene: given as they are or as the root of an equation, read by name in its
// fields, and the scene written as LaTeX.
//
// Each claim is held against a second route: the thickness an equation defines against its
// closed form, and the surface built with it against the geometry (every vertex of the shell at
// the tube's radius plus the thickness). Each check gets a case it must fail: a field that reads
// a name nothing defined is refused at the call that hands it over -- for a thickness, a motion
// and a colour alike -- so no symbol reaches cook() or materialize() to become a NaN there.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/roots.hpp>
#include <spatium/io/build.hpp>
#include <spatium/io/trace_text.hpp>

#include <cmath>
#include <format>
#include <stdexcept>
#include <string>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinULP;
namespace bd = spatium::io::build;
using F = bd::Field<double>;
using V = bd::PointField<double>;

namespace {

// the torus's tube as a scene: R, r given, the shell's thickness h the root of an equation
struct ShellScene {
    bd::Trace<double> trace;
    bd::Named<double> big, tube, h;
    bd::Handle<double> torus, shell;

    ShellScene() {
        big = trace.let("R", 2.0);
        tube = trace.let("r", 0.5);
        // the shell holds as much as the tube: pi (r + h)^2 - pi r^2 = pi r^2
        const F r = tube.field(), t = F::symbol("h");
        h = trace.let("h", (r + t) * (r + t) - F(2.0) * r * r, 0.0, 1.0);
        torus = trace.torus(big, tube);
        shell = trace.offset(torus, h.field());
    }
};

bool contains(const std::string& text, const std::string& part) { return text.find(part) != std::string::npos; }

template<typename Call>
std::string refusal(Call&& call) {
    try {
        call();
    } catch (const std::invalid_argument& e) {
        return e.what();
    }
    return {};
}

} // namespace

TEST_CASE("find_root: the root a sign change holds, to the last bit, and no root without one", "[algebra][roots]") {
    const auto two = spatium::find_root([](double x) { return x * x - 2; }, 0.0, 2.0);
    REQUIRE(two.has_value());
    CHECK_THAT(*two, WithinULP(std::sqrt(2.0), 1));
    CHECK(*spatium::find_root([](double x) { return x - 1; }, 1.0, 3.0) == 1.0);   // an end that is the root
    // the cases it must fail: no change of sign (a double root, and none), and a NaN met inside
    CHECK_FALSE(spatium::find_root([](double x) { return (x - 1) * (x - 1); }, 0.0, 2.0).has_value());
    CHECK_FALSE(spatium::find_root([](double x) { return x * x + 1; }, -1.0, 1.0).has_value());
    CHECK_FALSE(spatium::find_root([](double x) { return x < 0.5 ? -1.0 : std::sqrt(-x); }, 0.0, 1.0).has_value());
}

TEST_CASE("scene symbols: a thickness given by an equation is its closed form, and the shell is where it says",
          "[build_dsl][symbol]") {
    ShellScene s;
    CHECK_THAT(s.h.value, WithinULP(0.5 * (std::sqrt(2.0) - 1), 4));
    REQUIRE(s.trace.definitions().size() == 3);
    CHECK(s.trace.definitions()[2].equation.has_value());

    // the node computes the bound thickness and keeps the written one
    const auto& node = s.trace.node(s.shell.index);
    CHECK(node.thickness.free_symbols().empty());
    CHECK(node.thickness_written.free_symbols() == std::vector<std::string>{"h"});

    // every vertex of the shell at the tube's radius plus the thickness, and none of them NaN
    const auto mesh = bd::materialize_mesh(s.trace, s.shell.index);
    REQUIRE(mesh.vertex_count() > 0);
    for (const auto& p : mesh.vertices) {
        REQUIRE((std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2])));
        const double ring = std::hypot(p[0], p[1]) - s.big.value;
        CHECK_THAT(std::hypot(ring, p[2]), WithinAbs(s.tube.value + s.h.value, 1e-12));
    }
}

TEST_CASE("scene symbols: a name nothing defined is refused where it is handed over, never a NaN later",
          "[build_dsl][symbol]") {
    ShellScene s;
    CHECK(contains(refusal([&] { s.trace.offset(s.torus, F::symbol("q")); }), "'q'"));
    CHECK(contains(refusal([&] { s.trace.let("k", F::symbol("k") - F::symbol("w"), 0.0, 1.0); }), "'w'"));
    CHECK(contains(refusal([&] { s.trace.let("k", F::symbol("R") - F(1.0), 0.0, 1.0); }), "no open name"));
    CHECK(contains(refusal([&] { s.trace.let("k", F::symbol("k") * F::symbol("k") + F(1.0), -1.0, 1.0); }),
                   "change sign"));
    CHECK(contains(refusal([&] { s.trace.let("r", 1.0); }), "already defined"));

    // a motion and a colour: their factors are fields too, bound or refused the same way
    auto ball = s.trace.sphere(1.0);
    CHECK(contains(refusal([&] { ball.moving(scaled(V::point(), F::symbol("grow"))); }), "'grow'"));
    CHECK(contains(refusal([&] { ball.colored(V::make(F::symbol("red"), F(0.0), F(0.0))); }), "'red'"));

    // defined, the motion reaches cook() bound: the object is scaled by the number the name has
    const auto grow = s.trace.let("grow", 3.0);
    auto moved = s.trace.sphere(1.0).moving(scaled(V::point(), grow.field()));
    CHECK(s.trace.node(moved.index).transform.free_symbols().empty());
    const auto cooked = bd::cook(s.trace, moved.index);
    REQUIRE(cooked.object_count() == 1);
    CHECK(cooked.objects().front().scale == 3.0);
    for (const auto& p : bd::materialize_mesh(s.trace, moved.index).vertices)
        CHECK_THAT(p.norm(), WithinAbs(3.0, 1e-12));
}

TEST_CASE("scene symbols: the scene as LaTeX says what was asked and what it came to", "[build_dsl][latex]") {
    ShellScene s;
    const auto tex = bd::to_latex(s.trace, s.shell.index);
    REQUIRE(tex.has_value());
    const std::string want = std::format(
        "\\begin{{aligned}}\n"
        "R &= 2 \\\\\n"
        "r &= 0.5 \\\\\n"
        "& \\left(r + h\\right) \\cdot \\left(r + h\\right) - 2 \\cdot r \\cdot r = 0,\\quad h \\in \\left[0, 1\\right] "
        "\\;\\Rightarrow\\; h = {} \\\\\n"
        "S_{{0}} &= \\operatorname{{torus}}\\left(2, 0.5\\right) \\\\\n"
        "S_{{1}} &= \\operatorname{{offset}}\\left(S_{{0}}, h\\right)\n"
        "\\end{{aligned}}",
        s.h.value);
    CHECK(*tex == want);
    // the nodes the root does not reach are not in it
    s.trace.sphere(4.0);
    CHECK(*bd::to_latex(s.trace, s.shell.index) == want);
}
