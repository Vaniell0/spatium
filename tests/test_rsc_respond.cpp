// One door for a question typed as text, and the answer as the JSON a page reads.
//
// Each claim is held against the route it dispatches to: an equation answered through `respond`
// carries the steps `solve_closed_form` gives, a question with a parameter the answer `answer`
// gives, with the certificate of its pieces. The bound written for the error is above every
// piece's error and within a factor of ten of the largest. Each check gets a case it must fail: a
// text that is neither kind of question is refused, with the library's own reason.

#include <catch2/catch_test_macros.hpp>
#include <rsc/spatium/respond.hpp>

#include <cmath>
#include <string>

namespace {

using spatium::io::JsonValue;

const JsonValue& at(const JsonValue& v, std::string_view key) {
    for (const auto& [k, x] : v.as_object())
        if (k == key) return x;
    static const JsonValue none;
    return none;
}

} // namespace

TEST_CASE("respond: an expression with no unknown is evaluated, as a calculator does", "[rsc][respond]") {
    const auto four = rsc::respond("2+2");
    CHECK(at(four, "route").as_string() == "evaluate");
    CHECK(at(four, "value").as_number() == 4);
    CHECK(at(four, "latex").as_string() == "2 + 2 = 4");
    CHECK(at(rsc::respond("2 + 2 ="), "value").as_number() == 4);
    // the value is the field's, evaluated once: the same bits as the expression in code
    CHECK(at(rsc::respond("sqrt(2) * 3 - 1/3"), "value").as_number() == std::sqrt(2.0) * 3 - 1.0 / 3);
    CHECK(at(rsc::respond("2(3 + 4)^2"), "value").as_number() == 98);
    // the case it must fail: a name with no value is not evaluated, and the reason says what to write
    CHECK(at(rsc::respond("x + 1"), "refused").as_string().find("write an equation (x + 1 = 0)") != std::string::npos);
    CHECK(at(rsc::respond("2 +"), "refused").as_string().find("column") != std::string::npos);
}

TEST_CASE("respond: an equation is answered in closed form with its steps", "[rsc][respond]") {
    const auto r = rsc::respond("2x^2 - 3x - 5 = 0");
    CHECK(at(r, "route").as_string() == "closed form");
    CHECK(at(r, "unknown").as_string() == "x");
    const auto& steps = at(r, "steps").as_array();
    const auto direct = rsc::solve_closed_form("2x^2 - 3x - 5 = 0");
    REQUIRE(direct.has_value());
    REQUIRE(steps.size() == direct->steps.size());
    for (std::size_t i = 0; i < steps.size(); ++i) {
        CHECK(at(steps[i], "label").as_string() == direct->steps[i].label);
        CHECK(at(steps[i], "latex").as_string() == direct->steps[i].latex);
    }
    // the word solve in front changes nothing but the question as typed
    const auto worded = rsc::respond("  solve 2x^2 - 3x - 5 = 0 ");
    CHECK(at(worded, "question").as_string() == "solve 2x^2 - 3x - 5 = 0");
    CHECK(at(worded, "steps").dump() == at(r, "steps").dump());
}

TEST_CASE("respond: a question with a parameter is answered by synthesis, with its certificate", "[rsc][respond]") {
    const std::string q = "solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]";
    const auto r = rsc::respond(q, 1e-10);
    CHECK(at(r, "route").as_string() == "synthesis");
    CHECK(at(r, "answered").as_bool());
    const auto& c = at(r, "certificate");
    CHECK(at(c, "strength").as_string() == "sampled");
    CHECK(at(c, "points").as_number() == 2001);
    CHECK(at(c, "scalar").as_string() == "double");
    CHECK(at(c, "computed").as_string() == "host");
    CHECK(at(at(r, "answer").as_array().front(), "latex").as_string() == "\\sqrt{p}");

    // the bound is above every piece's error, and not ten times the largest
    auto asked = rsc::ask(q);
    REQUIRE(asked.has_value());
    asked->tolerance = 1e-10;
    const auto a = rsc::answer(*asked);
    double worst = 0;
    for (const auto& s : a.pieces) worst = std::max(worst, s.error);
    const double below = at(c, "error_below").as_number();
    CHECK(worst <= below);
    CHECK(below <= 10 * worst);
    CHECK(at(r, "pieces").as_number() == double(a.pieces.size()));
}

TEST_CASE("respond: a text that is no question is refused with the library's reason", "[rsc][respond]") {
    CHECK(at(rsc::respond("x^3 - 1 = 0"), "refused").as_string().find("degree at most two") != std::string::npos);
    CHECK(at(rsc::respond("hello"), "refused").as_string().find("reads 'hello'") != std::string::npos);
    CHECK(at(rsc::respond("solve x - p for x in [0, 1]"), "refused").as_string().find("parameter") != std::string::npos);
    CHECK(at(rsc::respond("solve x^2 + p = 0 for x in [0, 10], p in [0.5, 4]"), "refused").as_string().find("no value")
          != std::string::npos);
    CHECK(at(rsc::respond("2x^2 - 3x - 5 = 0"), "refused").is_null());
}
