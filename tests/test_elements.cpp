#include <catch2/catch_test_macros.hpp>
#include <spatium/physics/elements.hpp>
#include <spatium/physics/atomic/atom_model.hpp>
#include <spatium/physics/atomic/atom_svg.hpp>
#include <filesystem>

using namespace spatium::physics::atomic;

// ── Element lookup ────────────────────────────────────────────

TEST_CASE("Element by Z hydrogen", "[elements]") {
    auto& h = element(1);
    CHECK(h.Z == 1);
    CHECK(std::string_view(h.symbol) == "H");
    CHECK(std::string_view(h.name) == "Hydrogen");
}

TEST_CASE("Element by symbol", "[elements]") {
    auto& fe = element("Fe");
    CHECK(fe.Z == 26);
    CHECK(fe.period == 4);
}

TEST_CASE("Element electron count matches Z", "[elements]") {
    for (int Z : {1, 2, 6, 8, 10, 18, 26, 29, 36}) {
        auto& e = element(Z);
        CHECK(e.total_electrons() == Z);
    }
}

TEST_CASE("Element config subshells", "[elements]") {
    auto& c = element(6);  // Carbon: 1s² 2s² 2p²
    auto config = c.electron_config();
    CHECK(config.size() == 3);
    CHECK(config[0].n == 1);
    CHECK(config[0].l == 0);
    CHECK(config[0].electrons == 2);
    CHECK(config[2].l == 1);  // 2p
    CHECK(config[2].electrons == 2);
}

TEST_CASE("Element invalid Z throws", "[elements]") {
    CHECK_THROWS(element(999));
}

// ── Atom model ────────────────────────────────────────────────

TEST_CASE("AtomModel hydrogen", "[atom_model]") {
    auto model = AtomModel<>::build(1, 12);
    CHECK(model.Z == 1);
    CHECK(model.nucleus.vertices.size() > 0);
    // Hydrogen has 1s¹ → 1 orbital (l=0, m=0)
    // Marching cubes might not find it at low resolution, so just check it runs
}

TEST_CASE("AtomModel carbon", "[atom_model]") {
    auto model = AtomModel<>::build(6, 10);
    CHECK(model.Z == 6);
    CHECK(model.nucleus.faces.size() > 0);
}

// ── SVG export ────────────────────────────────────────────────

TEST_CASE("Orbital cross-section SVG", "[atom_svg]") {
    auto svg = orbital_cross_section_svg(2, 1, 0, 400, 50);
    CHECK(svg.content.size() > 100);
    auto path = std::filesystem::temp_directory_path() / "spatium_orbital_2p0.svg";
    CHECK(svg.save(path.string()));
    CHECK(std::filesystem::file_size(path) > 100);
    std::filesystem::remove(path);
}

TEST_CASE("Periodic table SVG", "[atom_svg]") {
    auto svg = periodic_table_svg();
    CHECK(svg.content.size() > 1000);
    auto path = std::filesystem::temp_directory_path() / "spatium_periodic.svg";
    CHECK(svg.save(path.string()));
    CHECK(std::filesystem::file_size(path) > 1000);
    std::filesystem::remove(path);
}

TEST_CASE("1s orbital cross-section SVG", "[atom_svg]") {
    auto svg = orbital_cross_section_svg(1, 0, 0, 400, 50);
    CHECK(svg.content.size() > 100);
}
