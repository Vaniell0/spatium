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

#include <algorithm>
#include <array>
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
    // the root on the same bracket by two more routes, proved and bisected, and they agree
    REQUIRE(two.candidates.size() == 3);
    for (std::size_t k : {1u, 2u}) {
        CHECK(two.candidates[k].refused.empty());
        CHECK(two.candidates[k].cost > 0);
        REQUIRE(two.candidates[k].values.size() == 1);
        CHECK(close(two.candidates[k].values.front(), two.values.front(), 1e-15));
    }
    CHECK(two.candidates[1].route == rsc::Route::Enclosed);
    CHECK(two.candidates[2].route == rsc::Route::Bracket);
    REQUIRE(two.enclosed);
    CHECK(two.enclosed->roots.front().contains(two.values.front()));
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

TEST_CASE("evaluate: a root with no bracket is solved in closed form, or proved on the doubles' whole range", "[rsc][evaluate]") {
    const auto d = doc("x: 2x^2 - 3x - 5 = 0\ny: y^3 - 1 = 0\nz: z^3 - 1 = 0 on [0, 2]");
    // the line without a bracket reads and prints back as written, as text and as JSON
    CHECK(bd::to_text(d) == "x: 2x^2 - 3x - 5 = 0\ny: y^3 - 1 = 0\nz: z^3 - 1 = 0 on [0, 2]\n");
    const auto again = bd::document_from_json(bd::to_json(d));
    REQUIRE(again);
    CHECK(bd::same_document(d, *again));

    const auto x = rsc::evaluate(d, "x");
    REQUIRE(x.route() == rsc::Route::ClosedForm);
    CHECK(x.values == std::vector<double>{2.5, -1.0});
    // the second path: the same two roots, each proved alone in its interval, on the whole line
    REQUIRE(x.candidates.size() == 2);
    REQUIRE(x.enclosed);
    REQUIRE(x.enclosed->roots.size() == 2);
    // where 2x^2 overflows a double there is no telling: only there, past a quarter of the largest double, nothing is proved
    REQUIRE(x.enclosed->undecided.size() == 1);
    CHECK(x.enclosed->undecided.front().lo > 4e307);
    CHECK(x.enclosed->roots[0].contains(-1.0));
    CHECK(x.enclosed->roots[1].contains(2.5));

    // a cubic has no closed form here: its one real root is proved, and the line holds no other
    const auto y = rsc::evaluate(d, "y");
    REQUIRE(y.route() == rsc::Route::Enclosed);
    CHECK(contains(y.candidates.front().refused, "degree at most two"));
    REQUIRE(y.values.size() == 1);
    REQUIRE(y.enclosed);
    CHECK(y.enclosed->complete());
    CHECK(y.enclosed->roots.front().contains(1.0));
    CHECK(spatium::proved(y.enclosed->roots.front()).digits >= 15);

    // on a bracket, the proof and the bisection agree
    const auto z = rsc::evaluate(d, "z");
    REQUIRE(z.route() == rsc::Route::Enclosed);
    CHECK_FALSE(z.candidates.front().refused.empty());
    REQUIRE(z.values.size() == 1);
    CHECK(close(z.values.front(), 1.0, 1e-15));
    REQUIRE(z.candidates.back().route == rsc::Route::Bracket);
    CHECK(z.enclosed->roots.front().contains(z.candidates.back().values.front()));
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
    // the shorthand adds what it was read as and the graph a page draws of it
    const auto drop = [](JsonValue v) {
        constexpr std::array<std::string_view, 3> added{"question", "document", "graph"};
        JsonValue::ObjectType kept;
        for (const auto& [k, x] : v.as_object())
            if (!std::ranges::contains(added, k)) kept.emplace_back(k, x);
        return JsonValue(kept).dump();
    };
    const auto line = [](std::string_view text, std::string_view name) {
        return JsonValue(rsc::evaluation_json(rsc::evaluate(doc(text), name))).dump();
    };
    CHECK(drop(rsc::respond("2x^2 - 3x - 5 = 0")) == line("x: 2x^2 - 3x - 5 = 0", "x"));
    CHECK(drop(rsc::respond("2 + 2")) == line("answer = 2 + 2", "answer"));
    CHECK(drop(rsc::respond("solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]")) ==
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
    CHECK(a.find("candidates")->as_array().size() == 3);
    // and a request that is not one says why
    CHECK(contains(rsc::evaluate_request("[1]").string_or("refused", ""), "not a JSON object"));
    CHECK(contains(rsc::evaluate_request(R"({"name": "x"})").string_or("refused", ""), "no document"));
    CHECK(contains(rsc::evaluate_request(R"({"document": "x: = 1", "name": "x"})").string_or("refused", ""), "line 1"));
    CHECK(contains(rsc::evaluate_request(R"({"document": "p in [0, 1]", "name": "p", "chosen": {"p": "a"}})").string_or("refused", ""),
                   "not a number"));
}

TEST_CASE("plot: a value is drawn along one parameter, the others its shader's slots", "[rsc][plot]") {
    const std::string text = "a in [0.5, 2]\nr = 0.1\nx in [-4, 4]\nf = a * sin(x) + r * x^2";
    const auto request = [&](std::string extra) {
        return rsc::plot_request(R"({"document": ")" + std::string("a in [0.5, 2]\\nr = 0.1\\nx in [-4, 4]\\nf = a * sin(x) + r * x^2") +
                                 R"(", "name": "f", "along": "x")" + extra + "}");
    };
    const auto p = request("");
    REQUIRE_FALSE(p.find("refused"));
    const auto& symbols = p.find("symbols")->as_array();
    REQUIRE(symbols.size() == 1);   // r is a value, put in; a is a parameter, left a slot
    CHECK(symbols[0].string_or("name", "") == "a");
    CHECK(symbols[0].find("interval")->as_array()[1].as_number() == 2);
    CHECK(p.find("interval")->as_array()[0].as_number() == -4);
    const std::string fragment = p.string_or("fragment", "");
    CHECK(fragment.starts_with("#version 300 es"));
    CHECK(contains(fragment, "uniform float field_symbols[1];"));
    CHECK(contains(fragment, "fwidth(F)"));
    // with a grid, F at each pixel's centre as the fp32 twin computes it: f(x) - y
    const auto g = request(R"(, "grid": {"width": 4, "height": 2, "view": [-4, -2, 4, 2], "slots": [1.5]})");
    const auto& values = g.find("values")->as_array();
    REQUIRE(values.size() == 8);
    for (int py = 0; py < 2; ++py)
        for (int px = 0; px < 4; ++px) {
            const double x = -4 + 8 * (px + 0.5) / 4, y = -2 + 4 * (py + 0.5) / 2;
            CHECK(close(values[static_cast<std::size_t>(py * 4 + px)].as_number(), 1.5 * std::sin(x) + 0.1 * x * x - y, 1e-5));
        }
    // the cases it must fail: a root is not drawn here, nor a parameter the value does not read
    CHECK(contains(rsc::plot_request(R"({"document": "p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]", "name": "x", "along": "p"})")
                       .string_or("refused", ""), "a graph is drawn of a value"));
    CHECK(contains(request("").string_or("refused", "none"), "none"));
    CHECK(contains(rsc::plot_request(R"({"document": "x in [0, 1]\nq in [0, 1]\nf = x", "name": "f", "along": "q"})")
                       .string_or("refused", ""), "does not read the parameter 'q'"));
    (void)text;
}

TEST_CASE("evaluate: a value is computed in the engine's arithmetic, each held against a double", "[rsc][evaluate][engine]") {
    const auto d = doc("a = 1/3 + 1/6\nb = 0.1 * 3\nc = abs(2 - 7) + min(1/4, 2) * 2^3\ns = sqrt(2) + 1\n"
                       "big = 2^62 * 4\nnone = 1 / (2 - 2)");
    const auto in = [&](const char* name, rsc::Arithmetic a, bool strict = false) {
        return rsc::evaluate(d, name, {}, rsc::Engine{.arithmetic = a, .strict = strict});
    };
    const auto& taken = [](const rsc::Evaluation& e) -> const rsc::Candidate& { return e.candidates[*e.taken]; };

    // a fraction is exact where a double is not, and says so
    CHECK(taken(in("a", rsc::Arithmetic::Rational)).text == "1/2");
    CHECK(taken(in("b", rsc::Arithmetic::Rational)).text == "3/10");
    CHECK(taken(in("c", rsc::Arithmetic::Rational)).text == "7");
    // a float writes its own shortest digits: 0.1 * 3 in a float is the float nearest 0.3
    CHECK(taken(in("b", rsc::Arithmetic::Float)).text == "0.3");
    // the case a double must fail: 0.1 * 3 is not 3/10 in binary
    CHECK(in("b", rsc::Arithmetic::Double).values.front() != 0.3);
    CHECK(in("b", rsc::Arithmetic::Rational).values.front() == 0.3);

    // the pair: every arithmetic's value, as a double, is the double's to its own precision
    for (const char* name : {"a", "b", "c"}) {
        INFO(name);
        const double want = in(name, rsc::Arithmetic::Double).values.front();
        CHECK(std::abs(in(name, rsc::Arithmetic::Rational).values.front() - want) <= 4e-16 * std::abs(want));
        const auto f = in(name, rsc::Arithmetic::Float);
        CHECK(taken(f).arithmetic == rsc::Arithmetic::Float);
        CHECK(std::abs(f.values.front() - want) <= 2e-7 * std::abs(want));
        if (std::string_view(name) == "b") CHECK(f.values.front() != want);   // and 0.1 * 3 rounds otherwise in a float
#if SPATIUM_HAS_DECIMAL
        CHECK(std::abs(in(name, rsc::Arithmetic::Decimal).values.front() - want) <= 4e-16 * std::abs(want));
#endif
    }

    // what an arithmetic cannot compute: a double after it, or with `strict`, a refusal
    const auto s = in("s", rsc::Arithmetic::Rational);
    REQUIRE(s.candidates.size() == 2);
    CHECK(s.candidates[0].refused == "rational has no Sqrt in it");
    CHECK(taken(s).arithmetic == rsc::Arithmetic::Double);
    CHECK(in("s", rsc::Arithmetic::Rational, true).refused == "rational has no Sqrt in it");
    CHECK(contains(in("big", rsc::Arithmetic::Rational, true).refused, "does not fit in 64 bits"));
    CHECK(contains(in("none", rsc::Arithmetic::Rational, true).refused, "divides by zero"));
    CHECK(in("big", rsc::Arithmetic::Rational).values.front() == std::ldexp(1.0, 64));

#if SPATIUM_HAS_DECIMAL
    // a decimal: ten digits of 0.1 are 1, and a third has sixteen threes
    CHECK(taken(in("b", rsc::Arithmetic::Decimal)).text == "0.3");
    CHECK(taken(in("a", rsc::Arithmetic::Decimal)).text == "0.5");
    CHECK(rsc::evaluate(doc("t = 1/3"), "t", {}, rsc::Engine{.arithmetic = rsc::Arithmetic::Decimal}).candidates[0].text ==
          "0.3333333333333333");
#else
    CHECK(contains(in("b", rsc::Arithmetic::Decimal, true).refused, "no decimal"));
#endif
#if defined(SPATIUM_HAS_BOOST_MULTIPRECISION) && SPATIUM_HAS_BOOST_MULTIPRECISION
    // Real50 has a square root, so it computes what a fraction cannot
    const auto r = in("s", rsc::Arithmetic::Real50, true);
    REQUIRE(r.taken);
    CHECK(taken(r).text.starts_with("2.41421356237309504880168872420969807856967187537"));
    CHECK(r.values.front() == std::sqrt(2.0) + 1);
#endif
}

TEST_CASE("evaluate: a number's digits are proved over the reals, and held against fifty digits", "[rsc][evaluate][proved]") {
    const auto d = doc("a = (3 + 7/12 - (2 + 11/18)) / (19/84)\nb = 0.1 * 3\ns = sqrt(2) + exp(1) * sin(1/3)\n"
                       "z = sin(3.141592653589793)\nfour = 2 + 2\nr = sqrt(0.1 * 3 - 0.3)");
    const auto in = [&](const char* name, rsc::Arithmetic a = rsc::Arithmetic::Double) {
        return rsc::evaluate(d, name, {}, rsc::Engine{.arithmetic = a});
    };
    for (const char* name : {"a", "b", "s"}) {
        INFO(name);
        const auto e = in(name);
        REQUIRE(e.proved);
        CHECK(e.proved->digits >= 14);
        // the pair: the proved digits are fifty digits' own, to half a unit of the last
#if defined(SPATIUM_HAS_BOOST_MULTIPRECISION) && SPATIUM_HAS_BOOST_MULTIPRECISION
        const auto deep = in(name, rsc::Arithmetic::Real50);
        REQUIRE(deep.proved);
        CHECK(deep.proved->digits >= 45);
        const spatium::Real50 proved(e.proved->text), truth(deep.proved->text);
        CHECK(abs(proved - truth) <= spatium::Real50(5) * pow(spatium::Real50(10), -e.proved->digits) * abs(truth));
#endif
    }
    // the case a double must fail: sin of the decimal written is 2.38e-16, a double says 1.22e-16,
    // and the interval proves no digit of it, only that it is within 6e-16 of zero
    const auto z = in("z");
    CHECK(close(z.values.front(), 1.2246467991473532e-16, 1e-15));
    REQUIRE(z.proved);
    CHECK(z.proved->digits == 0);
    CHECK(z.proved->bound < 6e-16);
#if defined(SPATIUM_HAS_BOOST_MULTIPRECISION) && SPATIUM_HAS_BOOST_MULTIPRECISION
    const auto z50 = in("z", rsc::Arithmetic::Real50);
    REQUIRE(z50.proved);
    // pi less its first sixteen digits; fifty digits of the argument leave about thirty-three of it
    CHECK(z50.proved->text.starts_with("2.384626433832795028841971693"));
    CHECK(z50.proved->text.ends_with("e-16"));
    CHECK(z50.proved->digits >= 30);
#endif
    // an exact sum is one number; a root of what is zero over the reals but not in a double is not proved
    CHECK(in("four").proved->exact);
    const auto r = in("r");
    CHECK(r.values.front() > 0);
    CHECK_FALSE(r.proved);
    CHECK(contains(r.unproved, "below zero"));
    // and the page reads it
    const auto json = rsc::respond("2 + 2");
    CHECK(contains(json.dump(), "\"proved\""));
}

TEST_CASE("evaluate: every root of an equation on its box, each proved alone, from the question line", "[rsc][evaluate]") {
    const auto d = doc("x: 2^x = 8\nw: sin(w) = 0.5 on [0, 10]\nn: n^2 + 1 = 0\nt: t^4 = 0 on [-1, 1]\ns: less(s, 0.5) - 0.5 = 0 on [0, 1]");
    // no closed form for a power: its one root on the whole line, proved to the double's digits
    const auto x = rsc::evaluate(d, "x");
    REQUIRE(x.route() == rsc::Route::Enclosed);
    REQUIRE(x.enclosed);
    CHECK(x.enclosed->complete());
    REQUIRE(x.enclosed->roots.size() == 1);
    CHECK(spatium::proved(x.enclosed->roots.front()).text == "3.000000000000000");
    CHECK(x.enclosed->box.hi == std::numeric_limits<double>::max());
    // four on the bracket, against their closed forms pi/6 + 2 pi k and 5 pi/6 + 2 pi k
    const auto w = rsc::evaluate(d, "w");
    REQUIRE(w.route() == rsc::Route::Enclosed);
    REQUIRE(w.values.size() == 4);
    const double pi = std::numbers::pi, at[4] = {pi / 6, 5 * pi / 6, 13 * pi / 6, 17 * pi / 6};
    for (int k = 0; k < 4; ++k) CHECK(close(w.values[k], at[k], 1e-14));
    // none real: the closed form's two are complex, and the proof agrees that the line holds none
    const auto n = rsc::evaluate(d, "n");
    REQUIRE(n.route() == rsc::Route::ClosedForm);
    CHECK(n.values.empty());
    REQUIRE(n.candidates.size() == 2);
    CHECK(contains(n.candidates[1].refused, "every piece of it is proved to hold none"));
    // a fourfold root is not proved, nor found by a sign it does not change, and the refusal says where it may be
    const auto t = rsc::evaluate(d, "t");
    CHECK_FALSE(t.taken);
    CHECK(contains(t.refused, "no root proved on [-1, 1]; one may lie in"));
    // a jump across zero is no root: the bisection lands on it, the proof leaves it undecided
    const auto s = rsc::evaluate(d, "s");
    REQUIRE(s.route() == rsc::Route::Bracket);
    CHECK(close(s.values.front(), 0.5, 1e-15));
    REQUIRE(s.enclosed);
    CHECK(s.enclosed->roots.empty());
    REQUIRE(s.enclosed->undecided.size() == 1);
    CHECK(s.enclosed->undecided.front().contains(0.5));
    CHECK(s.enclosed->undecided.front().width() < 1e-8);
    // and the page reads it, from the question line with its box
    const auto json = rsc::respond("sin(w) = 0.5 on [0, 10]").dump();
    CHECK(contains(json, "every root, proved"));
    CHECK(contains(json, "\"enclosed\""));
    CHECK(contains(json, "\"complete\""));
}

TEST_CASE("evaluate: the engine as JSON, beside the document", "[rsc][evaluate][engine]") {
    const auto ask = [](std::string_view engine) {
        return rsc::evaluate_request(std::string(R"({"document": "a = 1/3 + 1/6 - 1", "name": "a", "engine": )") +
                                     std::string(engine) + "}");
    };
    const auto r = ask(R"({"arithmetic": "rational"})");
    CHECK(r.string_or("text", "") == "-1/2");
    CHECK(r.string_or("latex", "") == "\\frac{1}{3} + \\frac{1}{6} - 1 = -\\frac{1}{2}");
    CHECK(r.find("candidates")->as_array()[0].string_or("arithmetic", "") == "rational");
    // with no engine, a double, and every candidate says so
    const auto plain = rsc::evaluate_request(R"({"document": "a = 1/3 + 1/6 - 1", "name": "a"})");
    CHECK_FALSE(plain.find("text"));
    CHECK(plain.find("candidates")->as_array()[0].string_or("arithmetic", "") == "double");
    // the case it must fail: an arithmetic the engine does not have, refused with those it has
    CHECK(contains(ask(R"({"arithmetic": "ternary"})").string_or("refused", ""),
                   "no arithmetic 'ternary'; it has double, float, rational, decimal, real50"));
}

TEST_CASE("evaluate: in fractions, a closed form's roots are what they are", "[rsc][evaluate][engine]") {
    using spatium::Rational;
    // the fraction's own square root: where it is one, and spoiled where it is not
    CHECK(sqrt(Rational(9, 4)) == Rational(3, 2));
    CHECK_FALSE(sqrt(Rational(2)).ok());
    CHECK_FALSE(sqrt(Rational(-4)).ok());

    const auto roots = [](std::string_view equation, bool strict = false) {
        return rsc::evaluate(doc(std::string("x: ") + std::string(equation)), "x", {},
                             rsc::Engine{.arithmetic = rsc::Arithmetic::Rational, .strict = strict});
    };
    const auto& taken = [](const rsc::Evaluation& e) -> const rsc::Candidate& { return e.candidates[*e.taken]; };

    const auto two = roots("2x^2 - 3x - 5 = 0");
    CHECK(taken(two).arithmetic == rsc::Arithmetic::Rational);
    CHECK(taken(two).text == "5/2, -1");
    CHECK(two.values == std::vector<double>{2.5, -1.0});
    CHECK(two.closed->steps[4].latex == "x_{1} = \\frac{5}{2},\\quad x_{2} = -1");
    CHECK(roots("x^2 - 2x + 1 = 0").closed->exact->text == "1, 1");
    CHECK(roots("0.5x^2 - 0.125 = 0").closed->exact->text == "1/2, -1/2");
    // a constant is read as written, not as the double it folds to: 1/3 is a third
    CHECK(roots("3x^2 - 1/3 = 0").closed->exact->text == "1/3, -1/3");
    CHECK(roots("x - 1/3 = 0").closed->exact->text == "1/3");

    // where the discriminant is not a square: (p ± q√m)/d, m free of square factors
    const auto surd = roots("x^2 - 3x - 2 = 0");
    CHECK(taken(surd).text == "(3 ± √17)/2");
    CHECK(surd.closed->exact->latex == "\\frac{3 \\pm \\sqrt{17}}{2}");
    CHECK(roots("x^2 - 8 = 0").closed->exact->text == "±2√2");
    CHECK(roots("x^2 + 2x + 5 = 0").closed->exact->text == "-1 ± 2i");
    CHECK(roots("3x^2 - 2x - 2 = 0").closed->exact->text == "(1 ± √7)/3");
    // the pair: p ± q√m as doubles are the double solver's roots
    for (const char* eq : {"x^2 - 3x - 2 = 0", "x^2 - 8 = 0", "3x^2 - 2x - 2 = 0", "x^2 - 0.1 = 0"}) {
        INFO(eq);
        const auto e = roots(eq);
        REQUIRE(e.closed->exact);
        const auto& x = *e.closed->exact;
        const double p = x.p.to_double(), q = x.q.to_double() * std::sqrt(static_cast<double>(x.m));
        const double r0 = e.closed->solution.roots[0].re, r1 = e.closed->solution.roots[1].re;
        CHECK(std::abs(std::max(r0, r1) - (p + q)) <= 4e-16 * std::abs(p + q) + 1e-300);
        CHECK(std::abs(std::min(r0, r1) - (p - q)) <= 4e-16 * std::abs(p - q) + 1e-300);
    }

    // the case it must fail: coefficients that are not fractions, and the double's route stands
    const auto irrational = roots("sqrt(2) x^2 - 1 = 0");
    REQUIRE(irrational.candidates.size() >= 2);
    CHECK(irrational.candidates[0].refused == "the equation reads Sqrt, which has no fraction");
    CHECK(taken(irrational).arithmetic == rsc::Arithmetic::Double);
    CHECK(roots("sqrt(2) x^2 - 1 = 0", true).refused == "the equation reads Sqrt, which has no fraction");
    // and with no engine, the closed form is the double's, with nothing exact
    CHECK_FALSE(rsc::evaluate(doc("x: x^2 - 3x - 2 = 0"), "x").closed->exact);
}
