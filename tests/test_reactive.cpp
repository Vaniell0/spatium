// Pull-based reactive cells and memos (core/reactive.hpp).

#include <catch2/catch_test_macros.hpp>
#include <spatium/core/reactive.hpp>

using namespace spatium;

TEST_CASE("an unread memo never runs", "[reactive]") {
    Cell<int> a(1);
    auto m = computed([&] { return a.get() * 2; });
    a.set(5);
    REQUIRE(m.runs() == 0);
}

TEST_CASE("a second read is cached", "[reactive]") {
    Cell<int> a(3);
    auto m = computed([&] { return a.get() + 1; });
    REQUIRE(m.get() == 4);
    REQUIRE(m.get() == 4);
    REQUIRE(m.runs() == 1);
}

TEST_CASE("setting a cell to its current value does nothing", "[reactive]") {
    Cell<int> a(3);
    auto m = computed([&] { return a.get() + 1; });
    (void)m.get();
    const auto v = a.version();
    a.set(3);
    REQUIRE(a.version() == v);
    (void)m.get();
    REQUIRE(m.runs() == 1);
    a.set(4);
    REQUIRE(m.get() == 5);
    REQUIRE(m.runs() == 2);
}

TEST_CASE("a chain recomputes only what changed", "[reactive]") {
    Cell<int> a(1), x(10);
    auto b = computed([&] { return a.get() + 1; });
    auto c = computed([&] { return b.get() * 2; });
    auto d = computed([&] { return x.get(); });
    REQUIRE(c.get() == 4);
    (void)d.get();
    a.set(2);
    REQUIRE(c.get() == 6);
    REQUIRE(b.runs() == 2);
    REQUIRE(c.runs() == 2);
    REQUIRE(d.runs() == 1);
}

TEST_CASE("an unchanged memo does not invalidate its dependents", "[reactive]") {
    Cell<int> a(4);
    auto parity = computed([&] { return a.get() % 2; });
    auto c = computed([&] { return parity.get() + 100; });
    REQUIRE(c.get() == 100);
    a.set(6);
    REQUIRE(c.get() == 100);
    REQUIRE(parity.runs() == 2);
    REQUIRE(c.runs() == 1);
    a.set(7);
    REQUIRE(c.get() == 101);
    REQUIRE(c.runs() == 2);
}

TEST_CASE("a diamond runs its join once", "[reactive]") {
    Cell<int> a(1);
    auto l = computed([&] { return a.get() + 1; });
    auto r = computed([&] { return a.get() * 10; });
    auto join = computed([&] { return l.get() + r.get(); });
    REQUIRE(join.get() == 12);
    a.set(2);
    REQUIRE(join.get() == 23);
    REQUIRE(l.runs() == 2);
    REQUIRE(r.runs() == 2);
    REQUIRE(join.runs() == 2);
}

TEST_CASE("dependencies are re-tracked on every run", "[reactive]") {
    Cell<bool> use_a(true);
    Cell<int> a(1), b(2);
    auto m = computed([&] { return use_a.get() ? a.get() : b.get(); });
    REQUIRE(m.get() == 1);
    b.set(20);  // not read yet: no rerun
    REQUIRE(m.get() == 1);
    REQUIRE(m.runs() == 1);
    use_a.set(false);
    REQUIRE(m.get() == 20);
    REQUIRE(m.runs() == 2);
    a.set(99);  // no longer read: no rerun
    REQUIRE(m.get() == 20);
    REQUIRE(m.runs() == 2);
    b.set(21);
    REQUIRE(m.get() == 21);
    REQUIRE(m.runs() == 3);
}
