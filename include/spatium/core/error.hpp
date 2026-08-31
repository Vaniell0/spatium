#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <concepts>
#  include <cstdint>
#  include <exception>
#  include <expected>
#  include <format>
#  include <functional>
#  include <string>
#  include <type_traits>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium {

enum class ErrorCode : std::uint32_t {
    Ok              = 0,
    InvalidArgument = 100,
    DegenerateInput = 101,
    OutOfDomain     = 102,
    NoIntersection  = 103,
    DegenerateShape = 104,
    SingularMatrix  = 105,
    ZeroNorm        = 106,
    ParseError      = 107,
    NotConverged    = 200,
    DimensionMismatch = 300,
    NotImplemented  = 900,
};

struct Error {
    ErrorCode   code;
    std::string message;

    inline explicit Error(ErrorCode c, std::string msg = {})
        : code(c), message(std::move(msg)) {}

    [[nodiscard]] inline std::string to_string() const {
        if (message.empty())
            return std::format("Error({})", static_cast<std::uint32_t>(code));
        return std::format("Error({}): {}",
                           static_cast<std::uint32_t>(code), message);
    }
};

template<typename T>
using Result = std::expected<T, Error>;

// ── BadResult exception ───────────────────────────────────────
// Thrown by unwrap()/expect() when the Result is an unexpected.
// Carries the original Error so callers can introspect code +
// message instead of getting a flat string-only std::runtime_error.

class BadResult : public std::exception {
public:
    explicit BadResult(Error e)
        : err_(std::move(e)), what_(err_.to_string()) {}

    [[nodiscard]] const Error& error() const noexcept { return err_; }
    [[nodiscard]] ErrorCode    code()  const noexcept { return err_.code; }
    [[nodiscard]] const char*  what()  const noexcept override { return what_.c_str(); }

private:
    Error       err_;
    std::string what_;
};

// Ergonomic Result unwrapping
template<typename T>
T unwrap(Result<T>&& r) {
    if (!r) throw BadResult(std::move(r.error()));
    return std::move(*r);
}

template<typename T>
T unwrap(const Result<T>& r) {
    if (!r) throw BadResult(r.error());
    return *r;
}

template<typename T>
T expect(Result<T>&& r, const char* msg) {
    if (!r) throw BadResult(Error{r.error().code,
                                  std::string{msg} + ": " + r.error().message});
    return std::move(*r);
}

// ── Generic Result<T> | Fn pipe ───────────────────────────────
// Any callable that takes T (by value, ref, or const-ref) and
// returns U or Result<U> can be piped onto a Result<T>:
//
//     Result<int> r = ...;
//     auto s = r | [](int x) { return x + 1; };       // Result<int>
//     auto v = r | [](int x) -> void { use(x); };     // Result<void>
//     auto t = r | [](int x) -> Result<int> { ... };  // Result<int>, flattened
//
// Failure short-circuits: if r holds an Error, the closure is not
// invoked and the Error propagates. The flatten overload mirrors
// std::expected::and_then so users do not nest Result<Result<U>>.
//
// To avoid clashing with the type-specific operator| pipes already
// shipped (Morphism, AffineTransform, lazy TransformExpr, geometry
// intersect chains), the generic pipe is gated on PlainPipeFn<F>:
// F must be a real callable but must NOT carry the structural tags
// of the type-specific overloads.

namespace detail {

template<typename F>
concept LooksLikeMorphism = requires { typename F::DomainPoint; typename F::CodomainPoint; };

template<typename F>
concept LooksLikeTransformLike = requires {
    typename F::scalar_type;
    { F::dim } -> std::convertible_to<std::size_t>;
};

template<typename F>
concept LooksLikeShape = requires {
    typename F::ScalarType;
    typename F::PointType;
    { F::ambient_dimension } -> std::convertible_to<std::size_t>;
};

template<typename F>
concept PlainPipeFn = !LooksLikeMorphism<std::remove_cvref_t<F>>
                  && !LooksLikeTransformLike<std::remove_cvref_t<F>>
                  && !LooksLikeShape<std::remove_cvref_t<F>>;

template<typename U>
struct unwrap_result { using type = U; };
template<typename U>
struct unwrap_result<Result<U>> { using type = U; };
template<typename U>
using unwrap_result_t = typename unwrap_result<U>::type;

} // namespace detail

template<typename T, typename F>
    requires std::invocable<F, T> && detail::PlainPipeFn<F>
constexpr auto operator|(Result<T>&& r, F&& f)
    -> Result<detail::unwrap_result_t<std::invoke_result_t<F, T>>>
{
    using Raw = std::invoke_result_t<F, T>;
    using Out = detail::unwrap_result_t<Raw>;
    if (!r) return std::unexpected(std::move(r.error()));
    if constexpr (std::is_void_v<Raw>) {
        std::invoke(std::forward<F>(f), *std::move(r));
        return Result<void>{};
    } else if constexpr (std::is_same_v<Raw, Result<Out>>) {
        return std::invoke(std::forward<F>(f), *std::move(r));
    } else {
        return std::invoke(std::forward<F>(f), *std::move(r));
    }
}

template<typename T, typename F>
    requires std::invocable<F, const T&> && detail::PlainPipeFn<F>
constexpr auto operator|(const Result<T>& r, F&& f)
    -> Result<detail::unwrap_result_t<std::invoke_result_t<F, const T&>>>
{
    using Raw = std::invoke_result_t<F, const T&>;
    using Out = detail::unwrap_result_t<Raw>;
    if (!r) return std::unexpected(r.error());
    if constexpr (std::is_void_v<Raw>) {
        std::invoke(std::forward<F>(f), *r);
        return Result<void>{};
    } else if constexpr (std::is_same_v<Raw, Result<Out>>) {
        return std::invoke(std::forward<F>(f), *r);
    } else {
        return std::invoke(std::forward<F>(f), *r);
    }
}

} // namespace spatium
