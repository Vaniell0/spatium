#pragma once

// Minimal hand-rolled JSON value type + parser/serializer.
//
// General-purpose utility -- deliberately scene-agnostic, so anything
// else that later wants a small config/data format can reuse it without
// pulling in spatium::io::scene. Written as a tokenizer/recursive-
// descent parser in the same hand-rolled style as io/obj.hpp rather than
// vendoring a JSON library: this project avoids new dependencies (see
// CONTRIBUTING.md), and a full JSON grammar is small enough to own.
//
// Supports the standard JSON value kinds (null, bool, number, string,
// array, object) with the usual string escapes (\" \\ \/ \b \f \n \r \t
// \uXXXX -- BMP code points only, no surrogate-pair decoding for values
// above U+FFFF). Object keys keep insertion order in a
// std::vector<pair<string, JsonValue>> rather than a map: JSON documents
// in this project are small hand-written files (scene descriptions,
// configs), not large data dumps that would need O(log n) key lookup or
// care about key ordering being sorted.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <charconv>
#  include <cmath>
#  include <string>
#  include <string_view>
#  include <utility>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io {

// A JSON value: null, bool, number (always double -- JSON has no
// separate integer type), string, array, or object.
class JsonValue {
public:
    enum class Type { Null, Bool, Number, String, Array, Object };
    using ArrayType  = std::vector<JsonValue>;
    using ObjectType = std::vector<std::pair<std::string, JsonValue>>;

    JsonValue() = default;
    JsonValue(std::nullptr_t) {}
    JsonValue(bool b) : type_(Type::Bool), bool_(b) {}
    JsonValue(double n) : type_(Type::Number), number_(n) {}
    JsonValue(int n) : type_(Type::Number), number_(static_cast<double>(n)) {}
    JsonValue(std::size_t n) : type_(Type::Number), number_(static_cast<double>(n)) {}
    JsonValue(std::string s) : type_(Type::String), string_(std::move(s)) {}
    JsonValue(const char* s) : type_(Type::String), string_(s) {}
    JsonValue(ArrayType a) : type_(Type::Array), array_(std::move(a)) {}
    JsonValue(ObjectType o) : type_(Type::Object), object_(std::move(o)) {}

    static JsonValue object() { return JsonValue(ObjectType{}); }
    static JsonValue array() { return JsonValue(ArrayType{}); }

    Type type() const { return type_; }
    bool is_null() const   { return type_ == Type::Null; }
    bool is_bool() const   { return type_ == Type::Bool; }
    bool is_number() const { return type_ == Type::Number; }
    bool is_string() const { return type_ == Type::String; }
    bool is_array() const  { return type_ == Type::Array; }
    bool is_object() const { return type_ == Type::Object; }

    bool as_bool(bool fallback = false) const { return type_ == Type::Bool ? bool_ : fallback; }
    double as_number(double fallback = 0.0) const { return type_ == Type::Number ? number_ : fallback; }

    const std::string& as_string() const {
        static const std::string empty;
        return type_ == Type::String ? string_ : empty;
    }

    const ArrayType& as_array() const {
        static const ArrayType empty;
        return type_ == Type::Array ? array_ : empty;
    }

    const ObjectType& as_object() const {
        static const ObjectType empty;
        return type_ == Type::Object ? object_ : empty;
    }

    // Field lookup on an object value -- nullptr if this is not an
    // object or the key is absent.
    const JsonValue* find(std::string_view key) const {
        if (type_ != Type::Object) return nullptr;
        for (auto& [k, v] : object_)
            if (k == key) return &v;
        return nullptr;
    }

    // Convenience accessors for optional scene/config fields: fall back
    // silently instead of forcing every call site to null-check find().
    double number_or(std::string_view key, double fallback) const {
        auto* v = find(key);
        return v ? v->as_number(fallback) : fallback;
    }

    std::string string_or(std::string_view key, std::string fallback) const {
        auto* v = find(key);
        return (v && v->is_string()) ? v->as_string() : fallback;
    }

    // Insert or overwrite a field. Turns a still-Null value into an
    // empty object first, so `JsonValue{}.set("x", 1).set("y", 2)` can
    // build an object from scratch without a separate object() call.
    JsonValue& set(std::string key, JsonValue value) {
        if (type_ == Type::Null) type_ = Type::Object;
        for (auto& [k, v] : object_) {
            if (k == key) { v = std::move(value); return v; }
        }
        object_.emplace_back(std::move(key), std::move(value));
        return object_.back().second;
    }

    // Append to an array. Turns a still-Null value into an empty array
    // first, mirroring set()'s auto-upgrade for objects.
    void push_back(JsonValue value) {
        if (type_ == Type::Null) type_ = Type::Array;
        array_.push_back(std::move(value));
    }

    // Serialize to compact JSON text (no pretty-printing -- scene files
    // in this project are small enough that a human editing them by
    // hand cares more about round-trip correctness than indentation).
    std::string dump() const {
        std::string out;
        write(out);
        return out;
    }

private:
    Type type_ = Type::Null;
    bool bool_ = false;
    double number_ = 0.0;
    std::string string_;
    ArrayType array_;
    ObjectType object_;

    static void append_hex4(std::string& out, unsigned code) {
        static constexpr char hex[] = "0123456789abcdef";
        out += "\\u";
        out += hex[(code >> 12) & 0xF];
        out += hex[(code >> 8) & 0xF];
        out += hex[(code >> 4) & 0xF];
        out += hex[code & 0xF];
    }

    static void write_escaped(std::string& out, std::string_view s) {
        out += '"';
        for (unsigned char c : s) {
            switch (c) {
                case '"':  out += "\\\""; break;
                case '\\': out += "\\\\"; break;
                case '\n': out += "\\n"; break;
                case '\r': out += "\\r"; break;
                case '\t': out += "\\t"; break;
                default:
                    if (c < 0x20) append_hex4(out, c);
                    else out += static_cast<char>(c);
            }
        }
        out += '"';
    }

    void write(std::string& out) const {
        switch (type_) {
            case Type::Null: out += "null"; break;
            case Type::Bool: out += bool_ ? "true" : "false"; break;
            case Type::Number: {
                // Integral-valued doubles print without a trailing
                // ".0" -- most scene fields (indices, small parameter
                // counts) are conceptually integers even though
                // JsonValue only stores double. Anything else uses
                // to_chars's shortest round-tripping representation.
                // The magnitude check runs BEFORE truncating to a
                // long long -- doing it the other way risks undefined
                // behavior converting an out-of-range double.
                if (std::abs(number_) < 1e15 && number_ == std::floor(number_)) {
                    out += std::to_string(static_cast<long long>(number_));
                } else {
                    char buf[64];
                    auto res = std::to_chars(buf, buf + sizeof(buf), number_);
                    out.append(buf, res.ptr);
                }
                break;
            }
            case Type::String: write_escaped(out, string_); break;
            case Type::Array: {
                out += '[';
                for (std::size_t i = 0; i < array_.size(); ++i) {
                    if (i) out += ',';
                    array_[i].write(out);
                }
                out += ']';
                break;
            }
            case Type::Object: {
                out += '{';
                for (std::size_t i = 0; i < object_.size(); ++i) {
                    if (i) out += ',';
                    write_escaped(out, object_[i].first);
                    out += ':';
                    object_[i].second.write(out);
                }
                out += '}';
                break;
            }
        }
    }
};

// ── Parse ──────────────────────────────────────────────────────

// internal -- do not use, no API stability. JsonParser is the
// recursive-descent implementation shared by every parse_*() helper
// below; call parse() instead and treat this as private, same
// convention as ray_surface.hpp's `detail` namespace.
namespace detail {

class JsonParser {
public:
    explicit JsonParser(std::string_view text) : s_(text) {}

    Result<JsonValue> parse() {
        auto v = parse_value();
        if (!v) return v;
        skip_ws();
        if (pos_ != s_.size())
            return std::unexpected(Error{ErrorCode::ParseError, "trailing content after JSON value"});
        return v;
    }

private:
    std::string_view s_;
    std::size_t pos_ = 0;

    bool eof() const { return pos_ >= s_.size(); }
    char peek() const { return s_[pos_]; }

    void skip_ws() {
        while (!eof() && (peek() == ' ' || peek() == '\t' || peek() == '\n' || peek() == '\r'))
            ++pos_;
    }

    bool consume(char c) {
        if (!eof() && peek() == c) { ++pos_; return true; }
        return false;
    }

    bool consume_literal(std::string_view lit) {
        if (s_.substr(pos_, lit.size()) == lit) { pos_ += lit.size(); return true; }
        return false;
    }

    Result<JsonValue> parse_value() {
        skip_ws();
        if (eof())
            return std::unexpected(Error{ErrorCode::ParseError, "unexpected end of input"});
        char c = peek();
        if (c == '{') return parse_object();
        if (c == '[') return parse_array();
        if (c == '"') return parse_string_value();
        if (c == 't' || c == 'f') return parse_bool();
        if (c == 'n') return parse_null();
        if (c == '-' || (c >= '0' && c <= '9')) return parse_number();
        return std::unexpected(Error{ErrorCode::ParseError,
            std::string("unexpected character '") + c + "'"});
    }

    Result<JsonValue> parse_null() {
        if (!consume_literal("null"))
            return std::unexpected(Error{ErrorCode::ParseError, "expected 'null'"});
        return JsonValue{};
    }

    Result<JsonValue> parse_bool() {
        if (consume_literal("true")) return JsonValue(true);
        if (consume_literal("false")) return JsonValue(false);
        return std::unexpected(Error{ErrorCode::ParseError, "expected 'true' or 'false'"});
    }

    Result<JsonValue> parse_number() {
        std::size_t start = pos_;
        if (!eof() && peek() == '-') ++pos_;
        while (!eof() && peek() >= '0' && peek() <= '9') ++pos_;
        if (!eof() && peek() == '.') {
            ++pos_;
            while (!eof() && peek() >= '0' && peek() <= '9') ++pos_;
        }
        if (!eof() && (peek() == 'e' || peek() == 'E')) {
            ++pos_;
            if (!eof() && (peek() == '+' || peek() == '-')) ++pos_;
            while (!eof() && peek() >= '0' && peek() <= '9') ++pos_;
        }
        std::string_view tok = s_.substr(start, pos_ - start);
        if (tok.empty() || tok == "-")
            return std::unexpected(Error{ErrorCode::ParseError, "bad number literal"});
        double value = 0.0;
        auto res = std::from_chars(tok.data(), tok.data() + tok.size(), value);
        if (res.ec != std::errc{})
            return std::unexpected(Error{ErrorCode::ParseError, "bad number literal: " + std::string(tok)});
        return JsonValue(value);
    }

    Result<std::string> parse_string_raw() {
        if (!consume('"'))
            return std::unexpected(Error{ErrorCode::ParseError, "expected '\"'"});
        std::string out;
        while (true) {
            if (eof())
                return std::unexpected(Error{ErrorCode::ParseError, "unterminated string"});
            char c = s_[pos_++];
            if (c == '"') break;
            if (c != '\\') { out += c; continue; }

            if (eof())
                return std::unexpected(Error{ErrorCode::ParseError, "unterminated escape"});
            char e = s_[pos_++];
            switch (e) {
                case '"':  out += '"';  break;
                case '\\': out += '\\'; break;
                case '/':  out += '/';  break;
                case 'b':  out += '\b'; break;
                case 'f':  out += '\f'; break;
                case 'n':  out += '\n'; break;
                case 'r':  out += '\r'; break;
                case 't':  out += '\t'; break;
                case 'u': {
                    if (pos_ + 4 > s_.size())
                        return std::unexpected(Error{ErrorCode::ParseError, "bad \\u escape"});
                    unsigned code = 0;
                    auto res = std::from_chars(s_.data() + pos_, s_.data() + pos_ + 4, code, 16);
                    if (res.ec != std::errc{})
                        return std::unexpected(Error{ErrorCode::ParseError, "bad \\u escape"});
                    pos_ += 4;
                    // BMP-only UTF-8 encoding -- no surrogate-pair
                    // decoding for \uD800-\uDFFF pairs above U+FFFF.
                    if (code < 0x80) {
                        out += static_cast<char>(code);
                    } else if (code < 0x800) {
                        out += static_cast<char>(0xC0 | (code >> 6));
                        out += static_cast<char>(0x80 | (code & 0x3F));
                    } else {
                        out += static_cast<char>(0xE0 | (code >> 12));
                        out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
                        out += static_cast<char>(0x80 | (code & 0x3F));
                    }
                    break;
                }
                default:
                    return std::unexpected(Error{ErrorCode::ParseError, "bad escape character"});
            }
        }
        return out;
    }

    Result<JsonValue> parse_string_value() {
        auto s = parse_string_raw();
        if (!s) return std::unexpected(s.error());
        return JsonValue(std::move(*s));
    }

    Result<JsonValue> parse_array() {
        if (!consume('[')) return std::unexpected(Error{ErrorCode::ParseError, "expected '['"});
        auto result = JsonValue::array();
        skip_ws();
        if (consume(']')) return result;
        while (true) {
            auto v = parse_value();
            if (!v) return std::unexpected(v.error());
            result.push_back(std::move(*v));
            skip_ws();
            if (consume(',')) continue;
            if (consume(']')) break;
            return std::unexpected(Error{ErrorCode::ParseError, "expected ',' or ']' in array"});
        }
        return result;
    }

    Result<JsonValue> parse_object() {
        if (!consume('{')) return std::unexpected(Error{ErrorCode::ParseError, "expected '{'"});
        auto result = JsonValue::object();
        skip_ws();
        if (consume('}')) return result;
        while (true) {
            skip_ws();
            auto key = parse_string_raw();
            if (!key) return std::unexpected(key.error());
            skip_ws();
            if (!consume(':'))
                return std::unexpected(Error{ErrorCode::ParseError, "expected ':' in object"});
            auto v = parse_value();
            if (!v) return std::unexpected(v.error());
            result.set(std::move(*key), std::move(*v));
            skip_ws();
            if (consume(',')) continue;
            if (consume('}')) break;
            return std::unexpected(Error{ErrorCode::ParseError, "expected ',' or '}' in object"});
        }
        return result;
    }
};

} // namespace detail

// Parse a JSON document from text. Rejects trailing content after the
// top-level value (a lone extra `}` or a second document concatenated
// on), same strictness as a standard JSON parser.
inline Result<JsonValue> parse(std::string_view text) {
    detail::JsonParser parser(text);
    return parser.parse();
}

} // namespace spatium::io
