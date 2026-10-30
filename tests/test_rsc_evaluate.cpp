// A definition of a document given its value, by the route it takes.
//
// What a definition is worth follows from what it still reads once the chosen parameters are put
// in: nothing, a number; one parameter, a function of it; more, a refusal. Each mode is held
// against another path to the same answer: the function synthesised with nothing chosen against
// the number with the parameter chosen, a closed form against the root on the same bracket, a
// value written as a function against the same value chosen, the calculator's shorthand against
// the document line it stands for, and a document sent as text against the same document as JSON.

#include <catch2/catch_test_macros.hpp>
#include <rsc/spatium/evaluate.hpp>
#include <rsc/spatium/respond.hpp>
#include <spatium/io/document.hpp>
#include <spatium/io/document_trace.hpp>

#include <cmath>
#include <numbers>
#include <string>

namespace {

namespace bd = spatium::io::build;
using spatium::io::JsonValue;

bd::Document doc(std::string_view text) {
    auto d = bd::read_document(text);
    REQUIRE(d);
    return *d;
}

bool contains(const std::string& s, std::string_view part) { return s.find(part) != std::string::npos; }

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)); }

// a candidate is taken exactly when it is the first that gives a value
void taken_is_first_with_a_value(const rsc::Evaluation& e) {
    std::size_t taken = 0;
    for (std::size_t i = 0; i < e.candidates.size(); ++i) {
        if (e.taken == i) {
            ++taken;
            CHECK(e.candidates[i].refused.empty());
            for (std::size_t j = 0; j < i; ++j) CHECK_FALSE(e.candidates[j].refused.empty());
        }
    }
    CHECK(taken == (e.taken ? 1u : 0u));
    CHECK(e.refused.empty() == e.taken.has_value());
}

} // namespace

TEST_CASE("evaluate: a root is a function of its parameter with nothing chosen, and a number with it chosen", "[rsc][evaluate]") {
    const auto d = doc("p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]");
    const auto f = rsc::evaluate(d, "x");
    REQUIRE(f.route() == rsc::Route::Synthesis);
    REQUIRE(f.taken);
    CHECK(f.free == std::vector<std::string>{"p"});
    REQUIRE(f.answer);
    CHECK(f.answer->answered);
    REQUIRE(f.answer->pieces.size() == 1);
    CHECK(*bd::to_text(f.answer->pieces.front().expr, "p") == "sqrt(p)");
    CHECK(f.candidates.front().cost > 0);
    taken_is_first_with_a_value(f);

    const auto two = rsc::evaluate(d, "x", {{"p", 2.0}});
    REQUIRE(two.route() == rsc::Route::ClosedForm);
    CHECK(two.free.empty());
    REQUIRE(two.values.size() == 1);   // the negative root is off the bracket
    CHECK(close(two.values.front(), std::sqrt(2.0), 1e-15));
    // the closed form's steps are counted: three values for the coefficients, five to check the degree
    CHECK(two.candidates[0].cost == 8);
    // the root on the same bracket is the second route, and agrees
    REQUIRE(two.candidates.size() == 2);
    CHECK(two.candidates[1].route == rsc::Route::Bracket);
    CHECK(two.candidates[1].refused.empty());
    CHECK(two.candidates[1].cost > 0);
    REQUIRE(two.candidates[1].values.size() == 1);
    CHECK(close(two.candidates[1].values.front(), two.values.front(), 1e-15));
    taken_is_first_with_a_value(two);
    // the function at the chosen value is the number, to the synthesis's tolerance
    for (double p : {0.5, 1.3, 2.0, 3.7, 4.0}) {
        const auto at = rsc::evaluate(d, "x", {{"p", p}});
        REQUIRE(at.values.size() == 1);
        CHECK(close(f.answer->pieces.front().expr(p, 0.0), at.values.front(), 1e-8));
    }
}

TEST_CASE("evaluate: a parameter left free makes the value a function of it, two make a refusal", "[rsc][evaluate]") {
    const auto d = doc("p in [0.5, 4]\nq in [1, 2]\nx: x^2 - p * q = 0 on [0, 10]");
    const auto both = rsc::evaluate(d, "x");
    CHECK_FALSE(both.taken);
    CHECK(contains(both.refused, "reads the parameters p, q"));
    CHECK(both.candidates.empty());

    const auto one = rsc::evaluate(d, "x", {{"q", 1.0}});
    REQUIRE(one.route() == rsc::Route::Synthesis);
    CHECK(one.free == std::vector<std::string>{"p"});
    REQUIRE(one.answer);
    CHECK(*bd::to_text(one.answer->pieces.front().expr, "p") == "sqrt(p)");

    const auto none = rsc::evaluate(d, "x", {{"p", 2.0}, {"q", 2.0}});
    REQUIRE(none.route() == rsc::Route::ClosedForm);
    REQUIRE(none.values.size() == 1);
    CHECK(close(none.values.front(), 2.0, 1e-15));
}

TEST_CASE("evaluate: a value is a number, or its own expression when it reads a free parameter", "[rsc][evaluate]") {
    const auto d = doc("r = 0.5\np in [1, 2]\na = 2 * r * p");
    const auto f = rsc::evaluate(d, "a");
    REQUIRE(f.route() == rsc::Route::Expression);
    CHECK(f.free == std::vector<std::string>{"p"});
    CHECK(f.written.free_symbols() == std::vector<std::string>{"p"});
    for (double p : {1.0, 1.25, 2.0}) {
        const auto at = rsc::evaluate(d, "a", {{"p", p}});
        REQUIRE(at.route() == rsc::Route::Evaluate);
        REQUIRE(at.values.size() == 1);
        CHECK(at.values.front() == f.written.bind("p", bd::Field<double>(p))(0.0, 0.0));
    }
    // a parameter is worth what it is given
    const auto p = rsc::evaluate(d, "p", {{"p", 1.5}});
    REQUIRE(p.values.size() == 1);
    CHECK(p.values.front() == 1.5);
    CHECK(rsc::evaluate(d, "p").route() == rsc::Route::Expression);
    CHECK(rsc::evaluate(d, "r").values == std::vector<double>{0.5});
}

TEST_CASE("evaluate: a root with no bracket is solved in closed form, or refused for want of one", "[rsc][evaluate]") {
    const auto d = doc("x: 2x^2 - 3x - 5 = 0\ny: y^3 - 1 = 0\nz: z^3 - 1 = 0 on [0, 2]");
    // the line without a bracket reads and prints back as written, as text and as JSON
    CHECK(bd::to_text(d) == "x: 2x^2 - 3x - 5 = 0\ny: y^3 - 1 = 0\nz: z^3 - 1 = 0 on [0, 2]\n");
    const auto again = bd::document_from_json(bd::to_json(d));
    REQUIRE(again);
    CHECK(bd::same_document(d, *again));

    const auto x = rsc::evaluate(d, "x");
    REQUIRE(x.route() == rsc::Route::ClosedForm);
    CHECK(x.values == std::vector<double>{2.5, -1.0});
    CHECK(x.candidates.size() == 1);

    const auto y = rsc::evaluate(d, "y");
    CHECK_FALSE(y.taken);
    CHECK(contains(y.refused, "degree at most two"));
    CHECK(contains(y.refused, "y needs a bracket"));

    const auto z = rsc::evaluate(d, "z");
    REQUIRE(z.route() == rsc::Route::Bracket);
    CHECK_FALSE(z.candidates.front().refused.empty());
    REQUIRE(z.values.size() == 1);
    CHECK(close(z.values.front(), 1.0, 1e-15));
    taken_is_first_with_a_value(z);

    // what needs a bracket says so: a question, and a scene
    CHECK(contains(rsc::query_of(d, "x").error().message, "x has no bracket"));
    const auto scene = bd::read_document("x: x^2 - 2 = 0\nS = sphere(x)");
    REQUIRE(scene);
    CHECK(contains(bd::to_trace(*scene).error().message, "a root with no bracket"));
}

TEST_CASE("evaluate: an integral with nothing free is a number with its error", "[rsc][evaluate]") {
    const auto d = doc("I = integral exp(-t^2) dt on [0, 1]\np in [0.5, 2]\nJ = integral exp(-p * t^2) dt on [0, 1]");
    const auto i = rsc::evaluate(d, "I");
    REQUIRE(i.route() == rsc::Route::Integral);
    REQUIRE(i.values.size() == 1);
    CHECK(close(i.values.front(), std::sqrt(std::numbers::pi) / 2 * std::erf(1.0), 1e-10));
    REQUIRE(i.error);
    CHECK(i.candidates.front().cost > 0);
    const auto j = rsc::evaluate(d, "J", {{"p", 1.0}});
    REQUIRE(j.values.size() == 1);
    CHECK(j.values.front() == i.values.front());
    CHECK(rsc::evaluate(d, "J").route() == rsc::Route::Synthesis);
}

TEST_CASE("evaluate: what is not asked well is refused, by name", "[rsc][evaluate]") {
    const auto d = doc("p in [0.5, 4]\nr = 1\nx: x^2 - p = 0 on [0, 10]\nS = sphere(r)");
    CHECK(contains(rsc::evaluate(d, "x", {{"p", 5.0}}).refused, "p = 5 is not on its interval [0.5, 4]"));
    CHECK(contains(rsc::evaluate(d, "x", {{"r", 1.0}}).refused, "the document has no parameter r"));
    CHECK(contains(rsc::evaluate(d, "w").refused, "defines no 'w'"));
    CHECK(contains(rsc::evaluate(d, "S").refused, "S is a shape"));
    // a refusal has no value and no route taken; a value has no refusal
    for (const char* name : {"x", "r", "w", "S"}) taken_is_first_with_a_value(rsc::evaluate(d, name));
}

TEST_CASE("evaluate: the calculator's shorthand is the document line it stands for", "[rsc][evaluate]") {
    const auto drop = [](JsonValue v, std::initializer_list<std::string_view> keys) {
        JsonValue::ObjectType kept;
        for (const auto& [k, x] : v.as_object())
            if (std::ranges::find(keys, k) == keys.end()) kept.emplace_back(k, x);
        return JsonValue(kept).dump();
    };
    const auto line = [](std::string_view text, std::string_view name) {
        return JsonValue(rsc::evaluation_json(rsc::evaluate(doc(text), name))).dump();
    };
    CHECK(drop(rsc::respond("2x^2 - 3x - 5 = 0"), {"question", "document"}) == line("x: 2x^2 - 3x - 5 = 0", "x"));
    CHECK(drop(rsc::respond("2 + 2"), {"question", "document"}) == line("answer = 2 + 2", "answer"));
    CHECK(drop(rsc::respond("solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]"), {"question", "document"}) ==
          line("p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]", "x"));
    // the answer carries the document it was read as
    const auto r = rsc::respond("(y - 1)(y + 2) = 4 - y");
    const JsonValue* written = r.find("document");
    REQUIRE(written);
    const auto back = bd::document_from_json(*written);
    REQUIRE(back);
    CHECK(bd::to_text(*back) == "y: (y - 1)(y + 2) = 4 - y\n");
}

TEST_CASE("evaluate: a document sent as text and as JSON is answered the same", "[rsc][evaluate]") {
    const std::string text = "p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]";
    JsonValue::ObjectType as_text{{"document", text}, {"name", "x"}, {"chosen", JsonValue::ObjectType{{"p", 2.0}}}};
    JsonValue::ObjectType as_json{{"document", bd::to_json(doc(text))}, {"name", "x"}, {"chosen", JsonValue::ObjectType{{"p", 2.0}}}};
    const auto a = rsc::evaluate_request(JsonValue(as_text).dump()), b = rsc::evaluate_request(JsonValue(as_json).dump());
    CHECK(a.dump() == b.dump());
    REQUIRE(a.find("values"));
    CHECK(close(a.find("values")->as_array().front().as_number(), std::sqrt(2.0), 1e-15));
    CHECK(a.find("candidates")->as_array().size() == 2);
    // and a request that is not one says why
    CHECK(contains(rsc::evaluate_request("[1]").string_or("refused", ""), "not a JSON object"));
    CHECK(contains(rsc::evaluate_request(R"({"name": "x"})").string_or("refused", ""), "no document"));
    CHECK(contains(rsc::evaluate_request(R"({"document": "x: = 1", "name": "x"})").string_or("refused", ""), "line 1"));
    CHECK(contains(rsc::evaluate_request(R"({"document": "p in [0, 1]", "name": "p", "chosen": {"p": "a"}})").string_or("refused", ""),
                   "not a number"));
}
