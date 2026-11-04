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
    // the question back as the row of a formula editor, written by the library
    CHECK(at(r, "question_latex").as_array().front().as_string() == "2x^{2} - 3x - 5 = 0");
}

TEST_CASE("respond: a question with a parameter is answered by synthesis, with its certificate", "[rsc][respond]") {
    const std::string q = "p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]";
    const auto r = rsc::respond(q, {.tolerance = 1e-10});
    CHECK(at(r, "route").as_string() == "synthesis");
    CHECK(at(r, "answered").as_bool());
    const auto& c = at(r, "certificate");
    CHECK(at(c, "strength").as_string() == "sampled");
    CHECK(at(c, "points").as_number() == 2001);
    CHECK(at(c, "scalar").as_string() == "double");
    CHECK(at(c, "computed").as_string() == "host");
    CHECK(at(at(r, "answer").as_array().front(), "latex").as_string() == "\\sqrt{p}");

    // the bound is above every piece's error, and not ten times the largest
    const auto doc = spatium::io::build::read_document(q);
    REQUIRE(doc.has_value());
    auto asked = rsc::query_of(*doc, "x");
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
    CHECK(at(rsc::respond("x: x - p = 0 on [0, 1]"), "refused").as_string().find("'p'") != std::string::npos);
    CHECK(at(rsc::respond("p in [0.5, 4]\nx: x^2 + p = 0 on [0, 10]"), "refused").as_string().find("no value")
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

TEST_CASE("respond: rows that are equations in as many unknowns are one system", "[rsc][respond]") {
    // in fractions: x = 1, y = 2 exactly, and each equation holds at it
    const auto r = rsc::respond("x + 2y = 5\n3x - y = 1");
    CHECK(at(r, "route").as_string() == "linear system");
    CHECK(at(at(r, "exact"), "y").as_string() == "2");
    const auto& s = at(r, "sheet");
    CHECK(at(s, "outcome").as_string() == "system");
    CHECK(at(s, "guarantee").as_string() == "exact");
    CHECK(at(s, "result").as_string() == "x = 1,\\quad y = 2");
    CHECK(at(r, "question_latex").as_array().size() == 2);
    // the same rows from the editor, and the same system written as a document's line
    const auto rows = rsc::respond_rows(R"({"rows": [["Equal", ["Add", "x", ["Multiply", 2, "y"]], 5], ["Equal", ["Subtract", ["Multiply", 3, "x"], "y"], 1]]})");
    CHECK(at(rows, "sheet").dump() == s.dump());
    CHECK(at(rsc::respond("x, y: x + 2y = 5, 3x - y = 1"), "sheet").dump() == s.dump());
    // fractions where they are: 1/3 is a third, not its double
    CHECK(at(at(rsc::respond("3x + y = 1\nx - y = 0"), "exact"), "x").as_string() == "1/4");
    // a coefficient with no fraction, sqrt(2): in doubles, checked by putting the solution back
    const auto in_double = rsc::respond("sqrt(2) x + y = 1\nx - y = 0");
    CHECK(at(at(in_double, "sheet"), "guarantee").as_string() == "computed");
    CHECK(std::abs(at(at(in_double, "solution"), "x").as_number() - 1 / (std::sqrt(2.0) + 1)) < 1e-15);
    // the cases it must fail: an equation that is not linear, and equations that are not independent
    CHECK(at(rsc::respond("x^2 + y = 1\nx - y = 0"), "refused").as_string() == "equation 1 is not linear in the unknowns");
    CHECK(at(rsc::respond("x + y = 1\n2x + 2y = 2"), "refused").as_string().find("not independent") != std::string::npos);
    // rows in more names than equations stay rows of a document
    CHECK(at(rsc::respond("x + y = 1\nx - y = z"), "route").as_string() != "linear system");
}

TEST_CASE("respond: rows are computed in the arithmetic their engine names", "[rsc][respond]") {
    // x^2 - 2 = 0 in fractions: the roots as the closed form writes them, the steps evaluate_request gives
    const std::string row = R"(["Equal", ["Subtract", ["Power", "x", 2], 2], 0])";
    const auto exact = rsc::respond_rows(R"({"rows": [)" + row + R"(], "engine": {"arithmetic": "rational"}})");
    const auto in_double = rsc::respond_rows(R"({"rows": [)" + row + "]}");
    const auto by_document = rsc::evaluate_request(R"({"document": "x: x^2 - 2 = 0", "name": "x", "engine": {"arithmetic": "rational"}})");
    CHECK(at(at(exact, "sheet"), "result").as_string() == "x_{1,2} = \\pm \\sqrt{2}");
    CHECK(at(exact, "steps").dump() == at(by_document, "steps").dump());
    CHECK(at(at(in_double, "sheet"), "result").as_string().find("1.414213562373095") != std::string::npos);
    // the case it must fail: an arithmetic the engine does not have is refused, naming those it has
    const auto unknown = rsc::respond_rows(R"({"rows": [)" + row + R"(], "engine": {"arithmetic": "quaternion"}})");
    CHECK(at(unknown, "refused").as_string().starts_with("the engine has no arithmetic 'quaternion'; it has double"));
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
    const auto root = rsc::respond("p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]");
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

TEST_CASE("respond: an exercise as a page of problems writes it, exact, with its steps", "[rsc][respond][steps]") {
    // ((3,2 - 1,7) : 0,003 / ((29/35 - 3/7) * 4 : 0,2) - (1 13/20 - 1,5) * 1,5 / ((2,44 + 1 14/25) * 1/8))
    //   : 62 1/20 + 1,364 : 0,124, as the formula editor's parser writes it once the page has put its
    // colons and decimal commas in that parser's notation (web/src/answers.ts, `notation`)
    const std::string row = R"(["Add",["Divide",["Delimiter",["Subtract",["Divide",["Divide",["Delimiter",["Subtract",{"num":"3.2"},{"num":"1.7"}]],{"num":"0.003"}],["Multiply",["Delimiter",["Subtract",["Divide",29,35],["Divide",3,7]]],["Divide",4,{"num":"0.2"}]]],["Divide",["Multiply",["Delimiter",["Subtract",["InvisibleOperator",1,["Divide",13,20]],{"num":"1.5"}]],{"num":"1.5"}],["Multiply",["Delimiter",["Add",{"num":"2.44"},["InvisibleOperator",1,["Divide",14,25]]]],["Divide",1,8]]]]],["InvisibleOperator",62,["Divide",1,20]]],["Divide",{"num":"1.364"},{"num":"0.124"}]])";
    const auto r = rsc::respond_rows(R"({"rows": [)" + row + "]}");
    CHECK(at(r, "route").as_string() == "evaluate");
    // twelve, exactly: the double's 12.000000000000002 is not the answer
    CHECK(at(at(r, "proved"), "exact").as_bool());
    CHECK(at(at(r, "proved"), "text").as_string() == "12");
    CHECK(at(r, "latex").as_string().ends_with(" = 12"));
    // a mixed number is a sum, 1 13/20 = 33/20, not the product a bare juxtaposition is
    const auto& steps = at(r, "steps").as_array();
    REQUIRE(steps.size() == 18);
    CHECK(at(steps[6], "latex").as_string() == "1 + \\frac{13}{20} = \\frac{33}{20}");
    CHECK(at(steps[15], "latex").as_string() == "\\frac{\\frac{1241}{20}}{\\frac{1241}{20}} = 1");
    CHECK(at(steps[16], "latex").as_string() == "\\frac{1.364}{0.124} = 11");
    CHECK(at(steps.back(), "latex").as_string() == "1 + 11 = 12");
}

TEST_CASE("respond: every answer is one sheet, its sections in one order", "[rsc][respond][sheet]") {
    struct Want { const char* question; const char* outcome; const char* guarantee; const char* label; bool ok; };
    const std::vector<std::string> order{"Input", "Steps", "Roots", "Certificate", "Check", "Pieces", "Properties", "How it was found"};
    for (const Want& w : std::initializer_list<Want>{
             {"2(3 + 4)^2 - 1/4", "number", "exact", "Exact", true},
             {"sqrt(2) + 1", "number", "proved", "Proved to 16 digits", true},
             {"10^20 + 1 - 10^20", "number", "none", "Evaluated; no digit proved, within 10⁴ of zero", false},
             {"2x^2 - 3x - 5 = 0", "roots", "proved", "Closed form, checked by Vieta's formulas", true},
             {"2^x = 8", "roots", "complete", "1 root proved, and no other on the box", true},
             {"x: sin(x) = 0.5 on [0, 10]", "roots", "complete", "4 roots proved, and no other on the box", true},
             {"x + 2y = 5\n3x - y = 1", "system", "exact", "Exact, by elimination in fractions", true},
             {"p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]", "function", "sampled", "Certified, relative error below 10⁻¹⁵", true},
             {"I = integral exp(-t^2) dt on [0, 1]", "number", "estimated", "An integral, error below 10⁻¹⁶ (estimated)", true},
             {"p in [1, 2]", "parameter", "exact", "Parameter", false},
             {"x^4 = 0", "refused", "none", "Not answered", false},
         }) {
        INFO(w.question);
        const auto r = rsc::respond(w.question);
        const auto& s = at(r, "sheet");
        CHECK(at(s, "outcome").as_string() == w.outcome);
        CHECK(at(s, "guarantee").as_string() == w.guarantee);
        CHECK(at(s, "label").as_string() == w.label);
        CHECK(at(s, "ok").as_bool() == w.ok);
        // the sections in the one order, each with what it shows
        std::size_t last = 0;
        for (const auto& sec : at(s, "sections").as_array()) {
            // a number's steps are titled by their count, "18 steps"
            const auto title = at(sec, "title").as_string().ends_with(" steps") ? std::string("Steps") : at(sec, "title").as_string();
            const auto k = std::ranges::find(order, title) - order.begin();
            REQUIRE(static_cast<std::size_t>(k) < order.size());
            CHECK(static_cast<std::size_t>(k) >= last);
            last = static_cast<std::size_t>(k);
            CHECK((sec.find("rows") || sec.find("steps") || sec.find("facts") || sec.find("table")));
        }
    }
    // the closed form's roots, each inside an interval the second route proved to hold one
    const auto q = at(rsc::respond("2x^2 - 3x - 5 = 0"), "sheet");
    bool held = false;
    for (const auto& sec : at(q, "sections").as_array())
        if (at(sec, "title").as_string() == "Check")
            for (const auto& f : at(sec, "facts").as_array())
                if (at(f, "label").as_string() == "Every real root, proved") held = at(f, "ok").as_bool();
    CHECK(held);
}

TEST_CASE("respond: rows above the last are rows of the page, an expression among them", "[rsc][respond][rows]") {
    // a row below an expression is answered: the expression is a statement of its own, answer_1
    CHECK(at(at(rsc::respond("2(3 + 4)^2 - 1/4\n1 + 2"), "sheet"), "result").as_string() == "3");
    CHECK(at(at(rsc::respond("2(3 + 4)^2 - 1/4\nanswer_1 + 1"), "sheet"), "result").as_string() == "\\frac{395}{4} = 98.75");
    // an equation above, then a parameter and a root in a document's own words
    const auto below = rsc::respond("2x^2 - 3x - 5 = 0\np in [0.5, 4]\ny: y^2 - p = 0 on [0, 10]");
    CHECK(at(at(below, "sheet"), "result").as_string() == "y = \\sqrt{p}");
    // the cases it must refuse, by the row they are on
    CHECK(at(rsc::respond("1 +\n2 + 2"), "refused").as_string().starts_with("line 1: parse: at column 4"));
    CHECK(at(rsc::respond("p in [0.5, 4]\nx: x^2 - q = 0 on [0, 10]"), "refused").as_string() ==
          "line 2: x reads 'q', which no line above defines");
}
