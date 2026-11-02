// A net of one hidden tanh layer as one op of the field IR, its weights numbers the op reads.
//
// Each claim is held against a second route: the field computes what `TanhNet` computes, to the
// bit; lowered, the plain-data interpreter computes the same to the bit, and its fp32 twin to a
// float's rounding. The claim the op exists for -- the weights are data, so a net trained
// further is new numbers for the same program -- is held as two nets with different weights
// emitting the same shader text and different weight buffers. The cases it must fail: weights
// of the wrong count are refused, a target with no buffer for them refuses the net, and nets
// that differ in one weight are not merged.

#include <catch2/catch_test_macros.hpp>
#include <spatium/algebra/models/network.hpp>
#include <spatium/io/field.hpp>
#include <spatium/io/field_glsl.hpp>
#include <spatium/io/field_pod.hpp>
#include <spatium/io/field_text.hpp>

#include <bit>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace bd = spatium::io::build;
using F = bd::Field<double>;

namespace {
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }
}

TEST_CASE("net: the field computes what the net computes, to the bit, and so does its lowered form", "[field][net]") {
    const auto net = spatium::TanhNet<2, 5>::random(7);
    const F f = bd::net<2, 5>(net.weights(), {F::u(), F::v() * F(0.5)});
    const auto pod = bd::lower(f);
    REQUIRE(pod.has_value());
    CHECK(pod->weights.size() == spatium::TanhNet<2, 5>::weight_count);
    const auto fp32 = bd::narrowed(*pod);
    double worst = 0;
    for (double u = -2; u <= 2; u += 0.25)
        for (double v = -2; v <= 2; v += 0.5) {
            const double want = net.value({u, v * 0.5});
            CHECK(bits(f(u, v)) == bits(want));
            CHECK(bits(bd::interpret(*pod, u, v)) == bits(want));
            worst = std::max(worst, std::abs(double(bd::interpret(fp32, float(u), float(v))) - want));
        }
    CHECK(worst < 1e-6);
    CHECK(worst > 0);   // and fp32 is not the double: the twin is a second route, not the same one
    // one input, and three
    const auto one = spatium::TanhNet<1, 3>::random(1);
    CHECK(bits(bd::net<1, 3>(one.weights(), {F::u()})(0.3, 0.0)) == bits(one.value({0.3})));
    const auto three = spatium::TanhNet<3, 4>::random(2);
    const F g = bd::net<3, 4>(three.weights(), {F::u(), F::v(), F::u() * F::v()});
    CHECK(bits(g(0.7, -1.1)) == bits(three.value({0.7, -1.1, 0.7 * -1.1})));
    CHECK(*bd::to_text(g) == "net[4](u, v, u * v)");
    CHECK(*bd::to_latex(bd::net<1, 3>(one.weights(), {F::u()})) == "\\operatorname{net}_{3}\\left(u\\right)");
}

TEST_CASE("net: the weights are data -- two nets of one shape are one program", "[field][net][glsl]") {
    const auto a = spatium::TanhNet<1, 8>::random(11), b = spatium::TanhNet<1, 8>::random(12);
    const auto emit = [](const spatium::TanhNet<1, 8>& n) {
        auto m = bd::glsl_es_module();
        REQUIRE(bd::emit_scalar_es(m, bd::net<1, 8>(n.weights(), {F::u()}), "f").has_value());
        return m;
    };
    const auto ma = emit(a), mb = emit(b);
    CHECK(ma.code == mb.code);
    CHECK(ma.code.find("field_net(0, 1, 8, vec3(") != std::string::npos);
    CHECK(ma.weights != mb.weights);
    CHECK(bd::glsl_es_prelude(ma).find("uniform vec4 field_weights[7];") != std::string::npos);   // 25 weights
    CHECK(bd::glsl_es_weights(ma).size() == 28);   // padded to whole vec4s
    // a net used twice is one block of weights; nets that differ are two
    const F twice = bd::net<1, 8>(a.weights(), {F::u()}) + bd::net<1, 8>(a.weights(), {F::v()});
    CHECK(bd::lower(twice)->weights.size() == 25);
    const F both = bd::net<1, 8>(a.weights(), {F::u()}) + bd::net<1, 8>(b.weights(), {F::u()});
    CHECK(bd::lower(both)->weights.size() == 50);
    // the cases it must fail: one weight changed is another net, not merged by consing
    auto c = a;
    c.weights()[3] += 1e-12;
    const F same = bd::net<1, 8>(a.weights(), {F::u()}), other = bd::net<1, 8>(c.weights(), {F::u()});
    std::vector<std::uint64_t> ka, kb;
    REQUIRE(same.canonical_key(ka));
    REQUIRE(other.canonical_key(kb));
    CHECK(ka != kb);
}

TEST_CASE("net: what is refused, and nothing folds it", "[field][net]") {
    // a count that is not F H + 2 H + 1
    const auto wrong = F::net_of(std::vector<double>(9, 0.1), 3, {F::u()});
    REQUIRE_FALSE(wrong.has_value());
    CHECK_FALSE(F::net_of(std::vector<double>(9, 0.1), 2, {F::u(), F::u(), F::u(), F::u()}).has_value());
    CHECK(wrong.error().message ==
          "net: 1 inputs through 3 units take 10 weights, and 9 were given; a net reads one to three inputs");
    CHECK(F::net_of(std::vector<double>(10, 0.1), 3, {F::u()}).has_value());
    // a Vulkan module binds no buffer for the weights
    const auto n = spatium::TanhNet<1, 2>::random(3);
    bd::GlslModule vk;
    const auto r = bd::emit_scalar(vk, bd::net<1, 2>(n.weights(), {F::u()}), "f");
    REQUIRE_FALSE(r.has_value());
    CHECK(r.error().message.find("binds no buffer to read its weights from") != std::string::npos);
    // a net of constants stays a net: folding it would write its weights into the program
    const F k = bd::net<1, 2>(n.weights(), {F(0.5)}).simplified();
    bool has_net = false;
    for (std::size_t i = 0; i < k.size(); ++i) has_net = has_net || k.op(i).op == bd::Op::Net;
    CHECK(has_net);
}
