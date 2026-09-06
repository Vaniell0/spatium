#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/io/json.hpp>

using namespace spatium::io;
using Catch::Matchers::WithinAbs;

// ── Scalars ──────────────────────────────────────────────────────

TEST_CASE("JSON parses null/true/false", "[json]") {
    auto n = parse("null");
    REQUIRE(n);
    CHECK(n->is_null());

    auto t = parse("true");
    REQUIRE(t);
    CHECK(t->is_bool());
    CHECK(t->as_bool() == true);

    auto f = parse("false");
    REQUIRE(f);
    CHECK(f->as_bool(true) == false);
}

TEST_CASE("JSON parses numbers", "[json]") {
    auto a = parse("42");
    REQUIRE(a);
    CHECK(a->as_number() == 42.0);

    auto b = parse("-17");
    REQUIRE(b);
    CHECK(b->as_number() == -17.0);

    auto c = parse("3.14");
    REQUIRE(c);
    CHECK_THAT(c->as_number(), WithinAbs(3.14, 1e-12));

    auto d = parse("-2.5e3");
    REQUIRE(d);
    CHECK_THAT(d->as_number(), WithinAbs(-2500.0, 1e-9));

    auto e = parse("1E2");
    REQUIRE(e);
    CHECK_THAT(e->as_number(), WithinAbs(100.0, 1e-9));

    auto z = parse("0");
    REQUIRE(z);
    CHECK(z->as_number() == 0.0);
}

TEST_CASE("JSON parses strings with escapes", "[json]") {
    auto s = parse(R"("hello\nworld")");
    REQUIRE(s);
    CHECK(s->as_string() == "hello\nworld");

    auto q = parse(R"("say \"hi\"")");
    REQUIRE(q);
    CHECK(q->as_string() == "say \"hi\"");

    auto u = parse(R"("AB")");
    REQUIRE(u);
    CHECK(u->as_string() == "AB");

    auto slash = parse(R"("a\/b")");
    REQUIRE(slash);
    CHECK(slash->as_string() == "a/b");

    auto plain = parse(R"("no escapes here")");
    REQUIRE(plain);
    CHECK(plain->as_string() == "no escapes here");
}

// ── Containers ───────────────────────────────────────────────────

TEST_CASE("JSON parses arrays", "[json]") {
    auto a = parse("[1, 2, 3]");
    REQUIRE(a);
    REQUIRE(a->is_array());
    REQUIRE(a->as_array().size() == 3);
    CHECK(a->as_array()[0].as_number() == 1.0);
    CHECK(a->as_array()[2].as_number() == 3.0);

    auto empty = parse("[]");
    REQUIRE(empty);
    CHECK(empty->is_array());
    CHECK(empty->as_array().empty());

    auto mixed = parse(R"([1, "two", true, null, [3]])");
    REQUIRE(mixed);
    REQUIRE(mixed->as_array().size() == 5);
    CHECK(mixed->as_array()[1].as_string() == "two");
    CHECK(mixed->as_array()[3].is_null());
    CHECK(mixed->as_array()[4].as_array()[0].as_number() == 3.0);
}

TEST_CASE("JSON parses nested objects and arrays", "[json]") {
    auto j = parse(
        R"({"name": "sphere", "pos": [1.5, -2.0, 0.0],
            "meta": {"visible": true, "tags": ["a", "b"]}})");
    REQUIRE(j);
    REQUIRE(j->is_object());
    CHECK(j->string_or("name", "") == "sphere");

    auto* pos = j->find("pos");
    REQUIRE(pos);
    REQUIRE(pos->as_array().size() == 3);
    CHECK_THAT(pos->as_array()[1].as_number(), WithinAbs(-2.0, 1e-12));

    auto* meta = j->find("meta");
    REQUIRE(meta);
    REQUIRE(meta->find("visible"));
    CHECK(meta->find("visible")->as_bool() == true);

    auto* tags = meta->find("tags");
    REQUIRE(tags);
    REQUIRE(tags->as_array().size() == 2);
    CHECK(tags->as_array()[0].as_string() == "a");
    CHECK(tags->as_array()[1].as_string() == "b");
}

TEST_CASE("JSON empty object parses", "[json]") {
    auto j = parse("{}");
    REQUIRE(j);
    CHECK(j->is_object());
    CHECK(j->as_object().empty());
}

// ── Malformed input ──────────────────────────────────────────────

TEST_CASE("JSON rejects malformed input", "[json]") {
    CHECK(!parse(""));
    CHECK(!parse("{"));
    CHECK(!parse("[1, 2,]"));           // trailing comma
    CHECK(!parse(R"({"a":1)"));          // unterminated object
    CHECK(!parse(R"("unterminated)"));   // unterminated string
    CHECK(!parse("{} extra"));           // trailing content after value
    CHECK(!parse("nul"));                // truncated literal
    CHECK(!parse("-"));                  // bare sign, no digits
    CHECK(!parse(R"({"a" 1})"));          // missing colon
    CHECK(!parse("[1 2]"));              // missing comma
}

// ── Round-trip via the JsonValue builder API ──────────────────────

TEST_CASE("JSON round-trips through dump()/parse()", "[json]") {
    auto obj = JsonValue::object();
    obj.set("name", std::string("test"));
    obj.set("count", 3);
    obj.set("ratio", 0.5);
    obj.set("active", true);
    obj.set("nothing", nullptr);

    auto arr = JsonValue::array();
    arr.push_back(1);
    arr.push_back(2);
    arr.push_back(3);
    obj.set("values", arr);

    auto nested = JsonValue::object();
    nested.set("inner", std::string("deep\"quote"));
    obj.set("nested", nested);

    std::string text = obj.dump();
    auto reparsed = parse(text);
    REQUIRE(reparsed);

    CHECK(reparsed->string_or("name", "") == "test");
    CHECK(reparsed->number_or("count", -1) == 3.0);
    CHECK_THAT(reparsed->number_or("ratio", -1), WithinAbs(0.5, 1e-12));
    REQUIRE(reparsed->find("active"));
    CHECK(reparsed->find("active")->as_bool() == true);
    REQUIRE(reparsed->find("nothing"));
    CHECK(reparsed->find("nothing")->is_null());

    auto* values = reparsed->find("values");
    REQUIRE(values);
    REQUIRE(values->as_array().size() == 3);
    CHECK(values->as_array()[2].as_number() == 3.0);

    auto* nested2 = reparsed->find("nested");
    REQUIRE(nested2);
    REQUIRE(nested2->find("inner"));
    CHECK(nested2->find("inner")->as_string() == "deep\"quote");
}

TEST_CASE("JSON round-trips a string with control characters", "[json]") {
    JsonValue v(std::string("tab\there\x01end"));
    auto reparsed = parse(v.dump());
    REQUIRE(reparsed);
    CHECK(reparsed->as_string() == "tab\there\x01end");
}

TEST_CASE("JSON integral doubles serialize without a decimal point", "[json]") {
    CHECK(JsonValue(5.0).dump() == "5");
    CHECK(JsonValue(-3.0).dump() == "-3");
    CHECK(JsonValue(0.0).dump() == "0");
    CHECK(JsonValue(2.5).dump() != "2");
}

TEST_CASE("JSON find() and _or() accessors on non-object/missing keys", "[json]") {
    JsonValue arr = JsonValue::array();
    CHECK(arr.find("x") == nullptr);
    CHECK(arr.number_or("x", 7.0) == 7.0);
    CHECK(arr.string_or("x", "fallback") == "fallback");

    auto obj = JsonValue::object();
    obj.set("present", 1);
    CHECK(obj.number_or("absent", 42.0) == 42.0);
    CHECK(obj.number_or("present", 42.0) == 1.0);
}
