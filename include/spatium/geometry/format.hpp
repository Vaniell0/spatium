#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/format.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/circle.hpp>
#  include <spatium/geometry/hyperplane.hpp>
#  include <spatium/geometry/line.hpp>
#  include <spatium/geometry/polygon.hpp>
#  include <spatium/geometry/ray_surface.hpp>
#  include <spatium/geometry/simplex.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <spatium/point.hpp>
#  include <format>
#  include <ostream>
#endif

// ── operator<< for geometry types ──────────────────────────────

SPATIUM_EXPORT namespace spatium::geometry {

template<std::size_t N, Scalar T>
std::ostream& operator<<(std::ostream& os, const Triangle<N, T>& t) {
    return os << "\u25b3[" << t[0] << ", " << t[1] << ", " << t[2] << ']';
}

template<std::size_t N, Scalar T>
std::ostream& operator<<(std::ostream& os, const Segment<N, T>& s) {
    return os << "[\u2014" << s.a << "\u2014" << s.b << ']';
}

template<std::size_t N, Scalar T>
std::ostream& operator<<(std::ostream& os, const Box<N, T>& b) {
    return os << "\u25a1[" << b.min_corner << "\u2194" << b.max_corner << ']';
}

template<std::size_t N, Scalar T>
std::ostream& operator<<(std::ostream& os, const Line<N, T>& l) {
    return os << "\u2190" << l.origin << " dir" << l.direction << "\u2192";
}

template<std::size_t N, Scalar T>
std::ostream& operator<<(std::ostream& os, const Ray<N, T>& r) {
    return os << r.origin << " \u2192" << r.direction;
}

template<std::size_t N, Scalar T>
std::ostream& operator<<(std::ostream& os, const Hyperplane<N, T>& h) {
    return os << "\u27c2 n=" << h.normal << " d=" << h.offset;
}

template<std::size_t N, Scalar T>
std::ostream& operator<<(std::ostream& os, const Circle<N, T>& c) {
    return os << "\u25ef c=" << c.center << " r=" << c.radius;
}

template<std::size_t N, Scalar T>
std::ostream& operator<<(std::ostream& os, const Disk<N, T>& d) {
    return os << "\u25cf c=" << d.boundary.center << " r=" << d.boundary.radius;
}

template<std::size_t N, Scalar T>
std::ostream& operator<<(std::ostream& os, const Polygon<N, T>& p) {
    os << "\u2b21[";
    for (std::size_t i = 0; i < p.vertices.size(); ++i) {
        if (i > 0) os << ", ";
        os << p.vertices[i];
    }
    return os << ']';
}

template<std::size_t N, std::size_t K, Scalar T>
    requires (K <= N)
std::ostream& operator<<(std::ostream& os, const Simplex<N, K, T>& s) {
    os << "\u0394" << K << '[';
    for (std::size_t i = 0; i <= K; ++i) {
        if (i > 0) os << ", ";
        os << s[i];
    }
    return os << ']';
}

} // namespace spatium::geometry

// ── std::formatter specializations ─────────────────────────────

template<std::size_t N, spatium::Scalar T>
struct std::formatter<spatium::geometry::Triangle<N, T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::Triangle<N, T>& t, auto& ctx) const {
        return std::format_to(ctx.out(), "\u25b3[{}, {}, {}]", t[0], t[1], t[2]);
    }
};

template<std::size_t N, spatium::Scalar T>
struct std::formatter<spatium::geometry::Segment<N, T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::Segment<N, T>& s, auto& ctx) const {
        return std::format_to(ctx.out(), "[\u2014{}\u2014{}]", s.a, s.b);
    }
};

template<std::size_t N, spatium::Scalar T>
struct std::formatter<spatium::geometry::Box<N, T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::Box<N, T>& b, auto& ctx) const {
        return std::format_to(ctx.out(), "\u25a1[{}\u2194{}]", b.min_corner, b.max_corner);
    }
};

template<std::size_t N, spatium::Scalar T>
struct std::formatter<spatium::geometry::Line<N, T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::Line<N, T>& l, auto& ctx) const {
        return std::format_to(ctx.out(), "\u2190{} dir{}\u2192", l.origin, l.direction);
    }
};

template<std::size_t N, spatium::Scalar T>
struct std::formatter<spatium::geometry::Ray<N, T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::Ray<N, T>& r, auto& ctx) const {
        return std::format_to(ctx.out(), "{} \u2192{}", r.origin, r.direction);
    }
};

template<std::size_t N, spatium::Scalar T>
struct std::formatter<spatium::geometry::Hyperplane<N, T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::Hyperplane<N, T>& h, auto& ctx) const {
        return std::format_to(ctx.out(), "\u27c2 n={} d={}", h.normal, h.offset);
    }
};

template<std::size_t N, spatium::Scalar T>
struct std::formatter<spatium::geometry::Circle<N, T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::Circle<N, T>& c, auto& ctx) const {
        return std::format_to(ctx.out(), "\u25ef c={} r={}", c.center, c.radius);
    }
};

template<std::size_t N, spatium::Scalar T>
struct std::formatter<spatium::geometry::Polygon<N, T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::Polygon<N, T>& p, auto& ctx) const {
        auto out = std::format_to(ctx.out(), "Polygon{{");
        for (std::size_t i = 0; i < p.vertices.size(); ++i) {
            if (i > 0) out = std::format_to(out, ", ");
            out = std::format_to(out, "{}", p.vertices[i]);
        }
        return std::format_to(out, "}}");
    }
};

template<std::size_t N, spatium::Scalar T>
struct std::formatter<spatium::geometry::Disk<N, T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::Disk<N, T>& d, auto& ctx) const {
        return std::format_to(ctx.out(), "Disk{{center={}, r={}}}",
                              d.boundary.center, d.boundary.radius);
    }
};

template<std::size_t N, std::size_t K, spatium::Scalar T>
    requires (K <= N)
struct std::formatter<spatium::geometry::Simplex<N, K, T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::Simplex<N, K, T>& s, auto& ctx) const {
        auto out = std::format_to(ctx.out(), "Simplex{}{{", K);
        for (std::size_t i = 0; i <= K; ++i) {
            if (i > 0) out = std::format_to(out, ", ");
            out = std::format_to(out, "{}", s[i]);
        }
        return std::format_to(out, "}}");
    }
};

template<spatium::Scalar T>
struct std::formatter<spatium::geometry::Quadric<T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::Quadric<T>& q, auto& ctx) const {
        // Quadric Q is the symmetric 4x4 matrix p^T Q p = 0.
        // Print the upper triangle in row-major order — six Q*xx
        // entries fully describe the surface.
        return std::format_to(ctx.out(),
                              "Quadric{{Qxx={}, Qyy={}, Qzz={}, "
                              "Qww={}, Qxy={}, Qxz={}, Qxw={}, "
                              "Qyz={}, Qyw={}, Qzw={}}}",
                              q.Q(0, 0), q.Q(1, 1), q.Q(2, 2), q.Q(3, 3),
                              q.Q(0, 1), q.Q(0, 2), q.Q(0, 3),
                              q.Q(1, 2), q.Q(1, 3), q.Q(2, 3));
    }
};

template<spatium::Scalar T>
struct std::formatter<spatium::geometry::Torus<T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::Torus<T>& t, auto& ctx) const {
        return std::format_to(ctx.out(),
                              "Torus{{center={}, axis={}, R={}, r={}}}",
                              t.center, t.axis, t.major_radius, t.minor_radius);
    }
};

template<spatium::Set S>
struct std::formatter<spatium::Point<S>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::Point<S>& p, auto& ctx) const {
        return std::format_to(ctx.out(), "P{}", p.raw());
    }
};
