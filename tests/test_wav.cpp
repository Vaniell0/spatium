#include <catch2/catch_test_macros.hpp>
#include <spatium/io/wav.hpp>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace spatium;

namespace {

std::vector<unsigned char> read_all_bytes(const std::filesystem::path& p) {
    std::ifstream f(p, std::ios::binary);
    return std::vector<unsigned char>(std::istreambuf_iterator<char>(f),
                                       std::istreambuf_iterator<char>());
}

std::string tag_at(const std::vector<unsigned char>& b, std::size_t off) {
    return std::string(b.begin() + static_cast<long>(off), b.begin() + static_cast<long>(off) + 4);
}

std::uint32_t read_u32(const std::vector<unsigned char>& b, std::size_t off) {
    std::uint32_t v;
    std::memcpy(&v, b.data() + off, 4);
    return v;
}

std::uint16_t read_u16(const std::vector<unsigned char>& b, std::size_t off) {
    std::uint16_t v;
    std::memcpy(&v, b.data() + off, 2);
    return v;
}

std::int16_t read_i16(const std::vector<unsigned char>& b, std::size_t off) {
    std::int16_t v;
    std::memcpy(&v, b.data() + off, 2);
    return v;
}

} // namespace

TEST_CASE("save_wav writes a byte-exact 44-byte canonical PCM header", "[wav]") {
    std::vector<double> samples = {0.0, 1.0, -1.0, 0.5};
    auto path = std::filesystem::temp_directory_path() / "spatium_test_wav_header.wav";

    io::WavFormat fmt;
    fmt.sample_rate  = 44100;
    fmt.num_channels = 1;

    auto result = io::save_wav(samples, path, fmt);
    REQUIRE(result.has_value());

    auto bytes = read_all_bytes(path);
    REQUIRE(bytes.size() == 44 + samples.size() * 2); // header + int16 samples

    CHECK(tag_at(bytes, 0) == "RIFF");
    CHECK(read_u32(bytes, 4) == 36 + samples.size() * 2); // ChunkSize
    CHECK(tag_at(bytes, 8) == "WAVE");

    CHECK(tag_at(bytes, 12) == "fmt ");
    CHECK(read_u32(bytes, 16) == 16); // Subchunk1Size (PCM, no extension)
    CHECK(read_u16(bytes, 20) == 1);  // AudioFormat = PCM integer
    CHECK(read_u16(bytes, 22) == 1);  // NumChannels
    CHECK(read_u32(bytes, 24) == 44100);      // SampleRate
    CHECK(read_u32(bytes, 28) == 44100 * 2);  // ByteRate = rate * channels * bytes/sample
    CHECK(read_u16(bytes, 32) == 2);          // BlockAlign = channels * bytes/sample
    CHECK(read_u16(bytes, 34) == 16);         // BitsPerSample

    CHECK(tag_at(bytes, 36) == "data");
    CHECK(read_u32(bytes, 40) == samples.size() * 2); // Subchunk2Size

    // Sample payload, decoded by hand against the documented [-1,1] -> int16
    // mapping (symmetric *32767 both directions, round-half-away-from-zero).
    CHECK(read_i16(bytes, 44) == 0);
    CHECK(read_i16(bytes, 46) == 32767);
    CHECK(read_i16(bytes, 48) == -32767);
    CHECK(read_i16(bytes, 50) == 16384); // round(0.5 * 32767) = round(16383.5) = 16384
}

TEST_CASE("save_wav header reflects a non-default sample rate and stereo channel count", "[wav]") {
    std::vector<double> samples = {0.0, 0.0, 0.25, -0.25}; // 2 stereo frames
    auto path = std::filesystem::temp_directory_path() / "spatium_test_wav_stereo.wav";

    io::WavFormat fmt;
    fmt.sample_rate  = 22050;
    fmt.num_channels = 2;

    REQUIRE(io::save_wav(samples, path, fmt).has_value());
    auto bytes = read_all_bytes(path);

    REQUIRE(bytes.size() == 44 + samples.size() * 2);
    CHECK(read_u16(bytes, 22) == 2);          // NumChannels
    CHECK(read_u32(bytes, 24) == 22050);      // SampleRate
    CHECK(read_u32(bytes, 28) == 22050 * 2 * 2); // ByteRate
    CHECK(read_u16(bytes, 32) == 4);          // BlockAlign = 2 channels * 2 bytes
    CHECK(read_u32(bytes, 4) == 36 + samples.size() * 2);
    CHECK(read_u32(bytes, 40) == samples.size() * 2);
}

TEST_CASE("save_wav clamps out-of-range samples instead of wrapping", "[wav]") {
    std::vector<double> samples = {2.0, -3.0, 1.0000001, -1.0000001};
    auto path = std::filesystem::temp_directory_path() / "spatium_test_wav_clamp.wav";

    REQUIRE(io::save_wav(samples, path).has_value()); // default format
    auto bytes = read_all_bytes(path);

    CHECK(read_i16(bytes, 44) == 32767);
    CHECK(read_i16(bytes, 46) == -32767);
    CHECK(read_i16(bytes, 48) == 32767);
    CHECK(read_i16(bytes, 50) == -32767);
}

TEST_CASE("save_wav rejects zero channels and a sample count that doesn't divide evenly", "[wav]") {
    auto path = std::filesystem::temp_directory_path() / "spatium_test_wav_invalid.wav";
    std::vector<double> samples = {0.0, 0.0, 0.0};

    io::WavFormat zero_channels;
    zero_channels.num_channels = 0;
    auto r1 = io::save_wav(samples, path, zero_channels);
    CHECK_FALSE(r1.has_value());

    io::WavFormat stereo;
    stereo.num_channels = 2;
    auto r2 = io::save_wav(samples, path, stereo); // 3 samples, not a multiple of 2
    CHECK_FALSE(r2.has_value());
}

TEST_CASE("save_wav default format is 44100 Hz mono", "[wav]") {
    std::vector<double> samples = {0.0};
    auto path = std::filesystem::temp_directory_path() / "spatium_test_wav_default.wav";
    REQUIRE(io::save_wav(samples, path).has_value());

    auto bytes = read_all_bytes(path);
    CHECK(read_u32(bytes, 24) == 44100);
    CHECK(read_u16(bytes, 22) == 1);
}
