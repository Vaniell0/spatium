#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <array>
#  include <cstdint>
#  include <span>
#  include <string_view>
#endif

SPATIUM_EXPORT namespace spatium::physics::atomic {

struct Subshell {
    uint8_t n, l, electrons;  // e.g. {1, 0, 2} = 1s²
};

struct Element {
    uint8_t Z;
    const char* symbol;
    const char* name;
    double atomic_mass;
    uint8_t period;
    uint8_t group;  // 0 for lanthanides/actinides
    uint8_t subshell_count;
    std::array<Subshell, 24> config;  // max 24 subshells (enough for all elements)

    std::span<const Subshell> electron_config() const {
        return {config.data(), subshell_count};
    }

    int total_electrons() const {
        int sum = 0;
        for (uint8_t i = 0; i < subshell_count; ++i)
            sum += config[i].electrons;
        return sum;
    }
};

// Access by atomic number (1-118) or symbol
const Element& element(int Z);
const Element& element(std::string_view symbol);

constexpr int ELEMENT_COUNT = 118;

// Subshell letter
constexpr char subshell_letter(int l) {
    constexpr char letters[] = "spdfghiklmno";
    return (l >= 0 && l < 12) ? letters[l] : '?';
}

} // namespace spatium::physics::atomic
