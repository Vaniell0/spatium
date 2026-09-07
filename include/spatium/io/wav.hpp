#pragma once

// Hand-rolled WAV (RIFF/WAVE) writer -- zero external audio dependency,
// same house style as this directory's obj.hpp/stl.hpp (a raw binary
// chunk format written field-by-field via std::ofstream, no library).
//
// Chunk layout (canonical 44-byte PCM header, per the Microsoft
// Multimedia Programming Interface WAVE spec):
//
//   offset  size  field
//   0       4     "RIFF"
//   4       4     ChunkSize        = 36 + data_bytes
//   8       4     "WAVE"
//   12      4     "fmt "
//   16      4     Subchunk1Size    = 16 (PCM, no extension)
//   20      2     AudioFormat      = 1  (PCM integer)
//   22      2     NumChannels
//   24      4     SampleRate
//   28      4     ByteRate         = SampleRate * NumChannels * BytesPerSample
//   32      2     BlockAlign       = NumChannels * BytesPerSample
//   34      2     BitsPerSample    = 16
//   36      4     "data"
//   40      4     Subchunk2Size    = NumSamples * NumChannels * BytesPerSample
//   44      ...   interleaved little-endian int16 samples
//
// Sample format: 16-bit signed PCM, chosen over 32-bit float (WAV format
// tag 3) because it's the one every WAV-capable tool/player handles
// unconditionally, and its fixed-width integer layout makes the file
// byte-exact and trivially diffable/testable -- no rounding ambiguity
// from a float encoding to account for. The differential-equation
// displacement fields this is meant to export (see
// physics/mechanics/wave_string.hpp) don't need more than 16-bit dynamic
// range to sound clean.
//
// Like save_stl's raw `float data[12]`/reinterpret_cast writes, this
// assumes a little-endian host (true for every platform this project
// targets) rather than carrying its own byte-swap path -- no existing
// io/ writer swaps bytes either.
//
// Scope: writer only, no load_wav. Nothing in this codebase needs to
// read a WAV file back in (the OBJ/STL pairing exists because meshes
// round-trip through both directions in real workflows); correctness is
// instead verified by parsing the written bytes directly in
// tests/test_wav.cpp. A reader can be added later if something needs it.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <algorithm>
#  include <cmath>
#  include <cstdint>
#  include <filesystem>
#  include <fstream>
#  include <span>
#endif

SPATIUM_EXPORT namespace spatium::io {

// ── WAV format parameters ───────────────────────────────────────

struct WavFormat {
    std::uint32_t sample_rate  = 44100;
    std::uint16_t num_channels = 1;
};

// ── Save WAV (16-bit PCM) ────────────────────────────────────────
//
// `samples` are interleaved per channel (frame0-ch0, frame0-ch1, ...,
// frame1-ch0, ...) and expected in [-1, 1]; values outside that range are
// clamped rather than left to wrap around on the int16 cast, since a
// silent overflow-wraparound would sound far worse than a hard clip.
inline Result<void> save_wav(std::span<const double> samples,
                              const std::filesystem::path& path,
                              WavFormat format = {}) {
    if (format.num_channels == 0)
        return std::unexpected(Error{ErrorCode::InvalidArgument, "num_channels must be >= 1"});
    if (samples.size() % format.num_channels != 0)
        return std::unexpected(Error{ErrorCode::InvalidArgument,
                                      "sample count is not a multiple of num_channels"});

    std::ofstream file(path, std::ios::binary);
    if (!file)
        return std::unexpected(Error{ErrorCode::InvalidArgument, "cannot open file for writing"});

    constexpr std::uint16_t bits_per_sample  = 16;
    constexpr std::uint16_t bytes_per_sample = bits_per_sample / 8;
    constexpr std::uint16_t audio_format_pcm = 1;
    constexpr std::uint32_t fmt_chunk_size   = 16;

    const auto data_bytes  = static_cast<std::uint32_t>(samples.size() * bytes_per_sample);
    const auto byte_rate   = static_cast<std::uint32_t>(format.sample_rate) * format.num_channels * bytes_per_sample;
    const auto block_align = static_cast<std::uint16_t>(format.num_channels * bytes_per_sample);
    const std::uint32_t riff_chunk_size = 36 + data_bytes;

    auto write_tag = [&](const char* tag) { file.write(tag, 4); };
    auto write_u32 = [&](std::uint32_t v) { file.write(reinterpret_cast<const char*>(&v), 4); };
    auto write_u16 = [&](std::uint16_t v) { file.write(reinterpret_cast<const char*>(&v), 2); };

    write_tag("RIFF");
    write_u32(riff_chunk_size);
    write_tag("WAVE");

    write_tag("fmt ");
    write_u32(fmt_chunk_size);
    write_u16(audio_format_pcm);
    write_u16(format.num_channels);
    write_u32(format.sample_rate);
    write_u32(byte_rate);
    write_u16(block_align);
    write_u16(bits_per_sample);

    write_tag("data");
    write_u32(data_bytes);
    for (double s : samples) {
        double clamped = std::clamp(s, -1.0, 1.0);
        auto quantized = static_cast<std::int16_t>(std::lround(clamped * 32767.0));
        write_u16(static_cast<std::uint16_t>(quantized));
    }

    if (!file)
        return std::unexpected(Error{ErrorCode::InvalidArgument, "write failed"});

    return {};
}

} // namespace spatium::io
