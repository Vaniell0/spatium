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
    // a fourfold root: no closed form, and no step can prove it, so the refusal says both and where it may be
    const auto fourfold = at(rsc::respond("x^4 = 0"), "refused").as_string();
    CHECK(fourfold.find("degree at most two") != std::string::npos);
    CHECK(fourfold.find("one may lie in") != std::string::npos);
    CHECK(at(rsc::respond("hello"), "refused").as_string().find("reads 'hello'") != std::string::npos);
    CHECK(at(rsc::respond("solve x - p for x in [0, 1]"), "refused").as_string().find("parameter") != std::string::npos);
    CHECK(at(rsc::respond("solve x^2 + p = 0 for x in [0, 10], p in [0.5, 4]"), "refused").as_string().find("no value")
          != std::string::npos);
    CHECK(at(rsc::respond("2x^2 - 3x - 5 = 0"), "refused").is_null());
}

TEST_CASE("respond: several lines are a document and its question", "[rsc][respond]") {
    // a name defined above is no unknown: 3a with a = 2 + 2 is a value, not an equation to write
    const auto twelve = rsc::respond("a = 2 + 2\n3a");
    CHECK(at(twelve, "route").as_string() == "evaluate");
    CHECK(at(twelve, "value").as_number() == 12);
    // one line that is a statement, not a shorthand: answered as the same document asks for it
    const auto root = rsc::respond("x: x^2 - 2 = 0 on [0, 2]");
    const auto by_document = rsc::evaluate_request(R"({"document": "x: x^2 - 2 = 0 on [0, 2]", "name": "x"})");
    CHECK(at(root, "route").as_string() == "closed form");
    CHECK(at(root, "steps").dump() == at(by_document, "steps").dump());
    // the cases it must fail: a line above that does not read, and a question reading an undefined name
    const auto bad = rsc::respond("p in [1, 0.5\nx^2 - p = 0");
    CHECK(at(bad, "route").as_string() == "document");
    CHECK(at(bad, "refused").as_string().find("line 1") != std::string::npos);
    CHECK(at(rsc::respond("a = 2\nb + a"), "refused").as_string().find("reads 'b'") != std::string::npos);
}

TEST_CASE("respond: the rows of a LaTeX editor answer as the lines they read as", "[rsc][respond]") {
    // as the editor's parser writes them (latex-syntax 0.151.0): p \in [0.5, 4] above x^2 - p = 0, x \in [0, 10]
    const std::string rows = R"({"rows": [["Element", "p", ["Interval", {"num": "0.5"}, 4]],
        ["Delimiter", ["Sequence", ["Equal", ["Subtract", ["Power", "x", 2], "p"], {"num": "0"}],
         ["Element", "x", ["Interval", {"num": "0"}, 10]]], {"str": ","}]]})";
    const auto by_rows = rsc::respond_rows(rows);
    const auto by_text = rsc::respond("p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]");
    CHECK(at(by_rows, "route").as_string() == "synthesis");
    CHECK(by_rows.dump() == by_text.dump());
    // 2x^2 - 3x - 5 = 0 as typed in the editor, and as typed in the language: the same steps
    const auto typed = rsc::respond_rows(R"({"rows": [["Equal", ["Subtract", ["Subtract", ["InvisibleOperator", 2, ["Power", "x", 2]],
        ["InvisibleOperator", 3, "x"]], 5], {"num": "0"}]]})");
    CHECK(at(typed, "route").as_string() == "closed form");
    CHECK(at(typed, "steps").dump() == at(rsc::respond("2x^2 - 3x - 5 = 0"), "steps").dump());
    // a row may be a line already in the language beside rows of MathJSON: the same document
    const auto mixed = rsc::respond_rows(R"({"rows": [["Element", "p", ["Interval", {"num": "0.5"}, 4]],
        {"text": "x: x^2 - p = 0 on [0, 10]"}]})");
    CHECK(mixed.dump() == by_text.dump());
    // the cases it must fail: a row the language cannot say is refused with its number, and no rows
    const auto integral_inside = rsc::respond_rows(R"({"rows": [["Element", "p", ["Interval", 0, 1]], ["Add", "p", ["D", "p"]]]})");
    CHECK(at(integral_inside, "row").as_number() == 2);
    CHECK(at(integral_inside, "refused").as_string() == "row 2: at [2]: 'D' is not in the language");
    CHECK(at(rsc::respond_rows(R"({"rows": []})"), "refused").as_string() == "the request has no rows");
    const auto bad_text = rsc::respond_rows(R"({"rows": [{"text": "x: x^2 - p = 0 on [0, 10]"}, {"text": "1 +"}]})");
    CHECK(at(bad_text, "refused").is_string());
}

TEST_CASE("respond: an answer carries the graph a page draws of it", "[rsc][respond]") {
    // an equation: its left side less its right along the unknown, over the roots 5/2 and -1 and a margin
    const auto quadratic = rsc::respond("2x^2 - 3x - 5 = 0");
    const auto& g = at(quadratic, "graph");
    CHECK(at(g, "document").as_string() == "x in [-2.75, 4.25]\ngraph = (2x^2 - 3x - 5) - (0)");
    // the second route: the graph drawn changes sign over its interval, as a function with those roots does
    const auto drawn = rsc::plot_request(g.dump());
    REQUIRE(at(drawn, "refused").is_null());
    const auto& range = at(drawn, "range").as_array();
    CHECK(range[0].as_number() < 0);
    CHECK(range[1].as_number() > 0);
    // a synthesised answer: its expression along its parameter, reaching sqrt(4) at the end
    const auto root = rsc::respond("solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]");
    CHECK(at(at(root, "graph"), "document").as_string() == "p in [0.5, 4]\nx = sqrt(p)");
    const auto root_drawn = rsc::plot_request(at(root, "graph").dump());
    const auto& r = at(root_drawn, "range").as_array();
    REQUIRE(r.size() == 2);
    CHECK(r[0].as_number() >= std::sqrt(0.5f) - 1e-6);
    CHECK(r[1].as_number() <= 2.0 + 1e-6);
    CHECK(r[1].as_number() > 1.99);
    // the cases it must fail: a number has no graph, and neither has a refusal
    CHECK(at(rsc::respond("2 + 2"), "graph").is_null());
    CHECK(at(rsc::respond("x^4 = 0"), "graph").is_null());
    // a root proved where no closed form is: drawn about it, as a closed form's are
    const auto cubic = rsc::respond("x^3 - 1 = 0");
    CHECK(at(cubic, "route").as_string() == "every root, proved");
    const std::string cubic_drawn = at(at(cubic, "graph"), "document").as_string();
    CHECK(cubic_drawn.starts_with("x in ["));
    CHECK(cubic_drawn.ends_with("graph = (x^3 - 1) - (0)"));
}

TEST_CASE("respond: a new name on the left is a definition, not an equation in it", "[rsc][respond]") {
    const auto four = rsc::respond("a = 2 + 2");
    CHECK(at(four, "route").as_string() == "evaluate");
    CHECK(at(four, "value").as_number() == 4);
    // a value of one parameter, as a person writes it under the parameter, has a graph along it
    const auto f = rsc::respond("a in [0.5, 2]\nf = a^2 + 1");
    CHECK(at(f, "refused").is_null());
    CHECK(at(at(f, "graph"), "along").as_string() == "a");
    // the case it must fail: a name the right side reads is an equation in it
    CHECK(at(rsc::respond("x = x^2 - 2"), "route").as_string() == "closed form");
}
