#pragma once
// A cooked scene laid out for a ray tracer: world-space triangles in one
// tree, instanced closed forms in another.
//
// This lived inside `donut_demo.cpp` until a second renderer needed it.
// The GPU path has to see exactly the scene the CPU path sees, and two
// copies of "how a Cooked becomes triangles" would disagree the first
// time one of them changed -- the same way `materialize_mesh()` and
// `cook()` once disagreed about every scattered object's orientation
// because nothing compared them.
//
// Not part of the `spatium.render` module, for the reason `io/build.hpp`
// is not part of `spatium.io`: it sits on the scene DSL, which is
// header-only outside the module build.
#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/geometry/ray_hit.hpp>
#  include <spatium/geometry/ray_surface.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <spatium/io/build.hpp>
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/spaces/euclidean.hpp>
#  include <spatium/spatial/bvh.hpp>
#  include <algorithm>
#  include <any>
#  include <array>
#  include <cstddef>
#  include <cstdint>
#  include <cstring>
#  include <type_traits>
#  include <unordered_map>
#  include <utility>
#  include <vector>
#endif

namespace spatium::render {

// Smooth normals per face *corner*, not per vertex, so an edge can stay
// sharp. Flat per-triangle normals on a coarse tessellation read as
// faceted plastic, and averaging every face around a vertex -- what this
// did until the cube's diagonal bands were noticed -- smooths across edges
// that are meant to be there: a cube's corner normal came out pointing
// along the diagonal, and each face shaded as a gradient between three.
//
// So a corner averages only the faces around its vertex that lie within
// `crease_cos` of its own face (60 degrees by default). A smooth surface's
// neighbours are all well inside that, and for it the result is the
// area-weighted average it always was, summed in the same face order and
// therefore to the same bits; a cube's are at 90 degrees and drop out.
template<Scalar T = double>
std::vector<std::array<Vec<T, 3>, 3>> corner_normals(const mesh::Mesh<Euclidean<3, T>>& m,
                                                     T crease_cos = T{0.5}) {
    const std::size_t F = m.faces.size(), V = m.vertex_count();
    std::vector<Vec<T, 3>> area_n(F), unit_n(F);
    for (std::size_t f = 0; f < F; ++f) {
        const auto& a = m.vertices[m.faces[f][0]];
        const auto& b = m.vertices[m.faces[f][1]];
        const auto& c = m.vertices[m.faces[f][2]];
        area_n[f] = Vec<T, 3>{Vec<T, 3>{b - a}.cross(Vec<T, 3>{c - a})};  // area-weighted
        const T len = area_n[f].norm();
        unit_n[f] = len > T{0} ? Vec<T, 3>{area_n[f] / len} : Vec<T, 3>{};
    }

    // Faces around each vertex, as one flat array with offsets, filled in
    // face order so the sums below run in the order the old per-vertex
    // accumulation did.
    std::vector<std::uint32_t> start(V + 1, 0), around;
    for (const auto& face : m.faces)
        for (auto v : face) ++start[v + 1];
    for (std::size_t v = 0; v < V; ++v) start[v + 1] += start[v];
    around.resize(start[V]);
    std::vector<std::uint32_t> fill(start.begin(), start.end() - 1);
    for (std::size_t f = 0; f < F; ++f)
        for (auto v : m.faces[f]) around[fill[v]++] = static_cast<std::uint32_t>(f);

    std::vector<std::array<Vec<T, 3>, 3>> out(F);
    for (std::size_t f = 0; f < F; ++f)
        for (std::size_t k = 0; k < 3; ++k) {
            const auto v = m.faces[f][k];
            Vec<T, 3> sum{0, 0, 0};
            for (auto i = start[v]; i < start[v + 1]; ++i) {
                const auto g = around[i];
                if (unit_n[g].dot(unit_n[f]) >= crease_cos) sum = Vec<T, 3>{sum + area_n[g]};
            }
            const T len = sum.norm();
            out[f][k] = len > T{1e-12} ? Vec<T, 3>{sum / len} : Vec<T, 3>{0, 0, 1};
        }
    return out;
}

// A triangle as the shader sees it. No opacity, and that is the current
// behaviour rather than an omission made here: only instances are ever
// see-through.
template<Scalar T = double>
struct ShadedTriangle {
    geometry::Triangle<3, T> tri;
    std::array<Vec<T, 3>, 3> vertex_normals;
    Vec<T, 3> color;
    T roughness;
    Vec<T, 3> emissive;
};

// What an instance looks like, parallel to the instance tree's input.
template<Scalar T = double>
struct InstanceLook {
    Vec<T, 3> color;
    T roughness;
    Vec<T, 3> emissive;
    T opacity;
};

// Where an instance came from, when it is known: the generator that made it
// (an index into `CookedScene::generators`) and its site within it. An
// instance that is one object of its own, or whose source has no version,
// has no origin, and is compared to the next scene's instances one by one.
struct InstanceOrigin {
    static constexpr std::uint32_t none = ~std::uint32_t{0};
    std::uint32_t generator = none;
    std::uint32_t site = 0;
};

// A scatter's instances as laid out, kept so the next scene can be compared
// to this one by what generated them instead of by the instances: the
// generator's identity (see `io::build::NodeVersion`), the shape every one
// of its instances points at, and how many there are.
struct LaidGenerator {
    std::uint64_t id = 0;
    std::vector<std::uint64_t> key;
    std::size_t shape = 0;
    std::size_t count = 0;
};

// Two trees, because the scene has two kinds of object in it.
//
// Most of it is ordinary geometry, each piece different, and a triangle
// tree is the right home for it. The rest is many copies of one closed
// form, and an instance holds a *reference* to the shape plus where this
// copy sits, so the tree holds one leaf per copy rather than one per
// triangle of every copy.
//
// What makes an object eligible is not that it looks small. It is that
// its motion is a *placement* -- affine in the point, so it moves the
// object without reshaping it (see VecField::is_placement) -- and that
// its shape carries a closed form.
//
// Move-only: every instance points into `quadrics`. Moving a vector keeps
// its buffer, copying does not, so a copy would leave the instance tree
// pointing into the original.
template<Scalar T = double>
struct CookedScene {
    std::vector<ShadedTriangle<T>> prims;          // parallel to the triangle tree's input
    std::vector<geometry::BoundedQuadric<T>> quadrics;  // one slot per shape, never grown
    std::vector<InstanceLook<T>> looks;            // parallel to the instance tree's input
    std::vector<InstanceOrigin> origins;           // parallel to the instance tree's input
    std::vector<LaidGenerator> generators;         // what `origins` refer to
    spatial::BVH<geometry::Triangle<3, T>> triangles;
    spatial::BVH<geometry::Instanced<geometry::BoundedQuadric<T>>> instances;

    CookedScene() = default;
    CookedScene(CookedScene&&) noexcept = default;
    CookedScene& operator=(CookedScene&&) noexcept = default;
    CookedScene(const CookedScene&) = delete;
    CookedScene& operator=(const CookedScene&) = delete;
};

// The triangle tree of everything that is not an instance, kept between
// lay-outs. That tree costs about 35 ms (-O3) for a 240x120 torus and depends on
// nothing but the triangles that go into it: the number of scattered
// copies, their seed and their placement only change the instance tree.
// So a scene regenerated with the same static part can reuse it.
//
// The key is the triangles themselves, not a description of how they were
// made. A recipe, its steps, a transform and a placement all end up in
// the vertex bytes, so one hash of those bytes covers every input that
// can change the tree, including ones added later -- and it costs about a
// hundredth of the build it replaces. The hash only picks the entry; a
// hit is confirmed by comparing the bytes, so a collision is a miss, never
// a wrong tree. (-0.0 and 0.0 differ in bytes: that is a miss, not an
// error.)
//
// The caller owns it and decides how long it lives; nothing is global.
// It keeps at most `capacity` trees, dropping the oldest, and each holds
// a copy of its triangles.
template<Scalar T = double>
class StaticCache {
public:
    using Tri = geometry::Triangle<3, T>;
    using Tree = spatial::BVH<Tri>;

    explicit StaticCache(std::size_t capacity = 4) : capacity_(capacity) {}

    std::size_t hits() const { return hits_; }
    std::size_t misses() const { return misses_; }
    std::size_t size() const { return entries_.size(); }
    void clear() { entries_.clear(); }

    // The tree over `tris`, built on a miss and copied out of the cache on
    // a hit. Takes the triangles by value for the same reason `build` does.
    Tree tree(std::vector<Tri> tris) {
        static_assert(std::is_trivially_copyable_v<Tri> && sizeof(Tri) == 9 * sizeof(T),
                      "hashing the bytes of a triangle needs it to have no padding");
        const std::uint64_t key = hash(tris);
        for (const auto& e : entries_)
            if (e.key == key && same(e.tree.shapes(), tris)) {
                ++hits_;
                return e.tree;
            }
        ++misses_;
        Tree built = Tree::build(std::move(tris));
        if (capacity_ > 0) {
            if (entries_.size() >= capacity_) entries_.erase(entries_.begin());
            entries_.push_back({key, built});
        }
        return built;
    }

private:
    struct Entry {
        std::uint64_t key;
        Tree tree;
    };

    // A word at a time, multiply-xorshift: a byte at a time would cost
    // as much as the thing being skipped is worth on a small tree.
    static std::uint64_t hash(const std::vector<Tri>& tris) {
        const auto* p = reinterpret_cast<const unsigned char*>(tris.data());
        const std::size_t n = tris.size() * sizeof(Tri);
        std::uint64_t h = 0x9e3779b97f4a7c15ull ^ n;
        std::size_t i = 0;
        for (; i + 8 <= n; i += 8) {
            std::uint64_t w;
            std::memcpy(&w, p + i, 8);
            h = (h ^ w) * 0xff51afd7ed558ccdull;
            h ^= h >> 32;
        }
        std::uint64_t w = 0;
        if (i < n) std::memcpy(&w, p + i, n - i);
        h = (h ^ w) * 0xff51afd7ed558ccdull;
        return h ^ (h >> 32);
    }

    static bool same(const std::vector<Tri>& a, const std::vector<Tri>& b) {
        return a.size() == b.size() &&
               (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(Tri)) == 0);
    }

    std::size_t capacity_;
    std::size_t hits_ = 0, misses_ = 0;
    std::vector<Entry> entries_;
};

// Lay a cooked scene out for tracing. `extra` is geometry that is drawn
// but is not in the scene -- the build-up frames' gizmos are the editor,
// not the subject, and arrive as Placed views onto their own Trace.
//
// `skip_nodes`, when given, is indexed by trace node: objects from a node
// marked in it are left out -- the instances a device kernel moves itself,
// which would otherwise be laid out twice.
//
// `cache`, when given, supplies the triangle tree if the same triangles
// were laid out before (see StaticCache); the result is the same scene
// either way, and without it nothing is kept.
template<Scalar T = double>
CookedScene<T> lay_out(const io::build::Cooked<T>& cooked,
                       const std::vector<io::build::Placed<T>>& extra = {},
                       const std::vector<char>* skip_nodes = nullptr,
                       StaticCache<T>* cache = nullptr);

// What an update of the instance tree did, for a caller that wants to
// report it (the example's HUD, the benchmark).
struct InstanceUpdate {
    bool reused = false;         // the previous tree was the starting point
    std::size_t kept = 0;        // instances found unchanged in it
    std::size_t removed = 0;
    std::size_t inserted = 0;
    bool rebuilt = false;        // too much changed (or no previous): built from scratch

    // How the change was found. Instances of a generator that was in the
    // previous scene are decided by comparing the generators; the rest are
    // hashed and compared one by one.
    std::size_t generators_matched = 0;
    std::size_t hashed = 0;      // instances that went through the per-instance comparison
};

// Above this fraction of the instances changed, the tree is built again
// rather than updated: see the figures in docs/ROADMAP.md.
inline constexpr double instance_update_max_fraction = 0.25;

namespace detail {

template<Scalar T>
bool same_bits(T a, T b) { return std::memcmp(&a, &b, sizeof(T)) == 0; }

template<Scalar T>
bool same_bits(const Vec<T, 3>& a, const Vec<T, 3>& b) {
    return same_bits(a[0], b[0]) && same_bits(a[1], b[1]) && same_bits(a[2], b[2]);
}

// Two instances are the same when everything that reaches a pixel is the
// same to the bit: the shape, the placement, and what it looks like.
template<Scalar T>
bool same_instance(const geometry::Instanced<geometry::BoundedQuadric<T>>& a, const InstanceLook<T>& la,
                   const geometry::Instanced<geometry::BoundedQuadric<T>>& b, const InstanceLook<T>& lb) {
    if (a.shape != b.shape || !same_bits(a.translation, b.translation) || !same_bits(a.scale, b.scale))
        return false;
    for (std::size_t r = 0; r < 3; ++r)
        for (std::size_t c = 0; c < 3; ++c)
            if (!same_bits(a.rotation(r, c), b.rotation(r, c))) return false;
    return same_bits(la.color, lb.color) && same_bits(la.roughness, lb.roughness) &&
           same_bits(la.emissive, lb.emissive) && same_bits(la.opacity, lb.opacity);
}

template<Scalar T>
std::uint64_t instance_hash(const geometry::Instanced<geometry::BoundedQuadric<T>>& a, const InstanceLook<T>& l) {
    std::uint64_t h = 0x9e3779b97f4a7c15ull ^ reinterpret_cast<std::uintptr_t>(a.shape);
    auto mix = [&h](T v) {
        std::uint64_t w = 0;
        std::memcpy(&w, &v, std::min(sizeof(T), sizeof(w)));
        h = (h ^ w) * 0xff51afd7ed558ccdull;
        h ^= h >> 32;
    };
    for (std::size_t i = 0; i < 3; ++i) mix(a.translation[i]);
    mix(a.scale);
    for (std::size_t r = 0; r < 3; ++r)
        for (std::size_t c = 0; c < 3; ++c) mix(a.rotation(r, c));
    for (std::size_t i = 0; i < 3; ++i) mix(l.color[i]);
    mix(l.roughness);
    for (std::size_t i = 0; i < 3; ++i) mix(l.emissive[i]);
    mix(l.opacity);
    return h;
}

template<Scalar T>
bool same_quadrics(const std::vector<geometry::BoundedQuadric<T>>& a,
                   const std::vector<geometry::BoundedQuadric<T>>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].closed != b[i].closed) return false;
        for (std::size_t r = 0; r < 4; ++r)
            for (std::size_t c = 0; c < 4; ++c)
                if (!same_bits(a[i].surface.Q(r, c), b[i].surface.Q(r, c))) return false;
        if (!same_bits(a[i].clip.min_corner, b[i].clip.min_corner) ||
            !same_bits(a[i].clip.max_corner, b[i].clip.max_corner))
            return false;
    }
    return true;
}

// Which of the old instances equal one of the new ones, bit for bit. The
// old instances are the tree's shapes at `old_slots` (all of them when it
// is null), the new ones are `insts`; sets `old_kept` / `new_kept` for each
// pair found and returns how many. The old ones go into an open-addressing
// table keyed by hash and each new one probes for an old one equal to it,
// so a hash collision is two unmatched instances and never a wrong match.
template<Scalar T>
std::size_t match_instances(const std::vector<geometry::Instanced<geometry::BoundedQuadric<T>>>& old,
                            const std::vector<InstanceLook<T>>& looks, const std::vector<std::uint32_t>* old_slots,
                            const std::vector<geometry::Instanced<geometry::BoundedQuadric<T>>>& insts,
                            const std::vector<InstanceLook<T>>& new_looks, std::vector<char>& old_kept,
                            std::vector<char>& new_kept) {
    struct Slot { std::uint64_t hash; std::uint32_t index; };
    constexpr std::uint32_t none = ~std::uint32_t{0};
    const std::size_t n_old = old_slots ? old_slots->size() : old.size();
    std::size_t cap = 16;
    while (cap < 2 * n_old) cap <<= 1;
    const std::size_t mask = cap - 1;
    std::vector<Slot> table(cap, Slot{0, none});
    for (std::size_t i = 0; i < n_old; ++i) {
        const std::uint32_t o = old_slots ? (*old_slots)[i] : std::uint32_t(i);
        const std::uint64_t h = instance_hash(old[o], looks[o]);
        std::size_t p = h & mask;
        while (table[p].index != none) p = (p + 1) & mask;
        table[p] = {h, o};
    }

    std::size_t kept = 0;
    for (std::size_t n = 0; n < insts.size(); ++n) {
        const std::uint64_t h = instance_hash(insts[n], new_looks[n]);
        for (std::size_t p = h & mask; table[p].index != none; p = (p + 1) & mask) {
            const auto o = table[p].index;
            if (table[p].hash != h || old_kept[o] || !same_instance(old[o], looks[o], insts[n], new_looks[n]))
                continue;
            old_kept[o] = 1;
            new_kept[n] = 1;
            ++kept;
            break;
        }
    }
    return kept;
}

template<Scalar T>
using InstanceTree = spatial::BVH<geometry::Instanced<geometry::BoundedQuadric<T>>>;

// Take out the old instances not marked kept and put in the new ones not
// marked kept, mirroring every removal on `looks` and `origins`: the tree
// moves its last shape into a freed slot, and the data parallel to it must
// do the same. Highest index first: the shape that moves into a freed slot
// comes from the end, and every unmatched index still to come is below it.
template<Scalar T>
void edit_instance_tree(InstanceTree<T>& tree, std::vector<InstanceLook<T>>& looks,
                        std::vector<InstanceOrigin>& origins, const std::vector<char>& old_kept) {
    for (std::size_t k = old_kept.size(); k-- > 0;) {
        if (old_kept[k]) continue;
        tree.remove(k);
        looks[k] = looks.back();
        looks.pop_back();
        origins[k] = origins.back();
        origins.pop_back();
    }
}

// Bring `tree` (with `looks` parallel to its shapes) to hold `insts` with
// `new_looks`, finding what changed by hashing every instance on both sides:
// the instances present in both stay where they are, the others are removed
// and inserted, and when more than `max_fraction` of them differ the tree is
// built again instead. Either way the tree answers as a build over `insts`
// would; only the order of its shapes differs.
//
// This is the way to update a set of instances nothing is known about, and
// the one `lay_out` falls back on for instances without an origin; with an
// origin the generators are compared instead (`update_laid_instances`).
template<Scalar T>
InstanceUpdate update_instance_tree(InstanceTree<T>& tree, std::vector<InstanceLook<T>>& looks,
                                    std::vector<geometry::Instanced<geometry::BoundedQuadric<T>>>&& insts,
                                    std::vector<InstanceLook<T>>&& new_looks, double max_fraction) {
    InstanceUpdate u;
    u.reused = true;
    const std::size_t n_old = tree.size(), n_new = insts.size();
    std::vector<char> old_kept(n_old, 0), new_kept(n_new, 0);
    u.kept = match_instances<T>(tree.shapes(), looks, nullptr, insts, new_looks, old_kept, new_kept);
    u.removed = n_old - u.kept;
    u.inserted = n_new - u.kept;
    u.hashed = n_old + n_new;

    if (double(u.removed + u.inserted) > max_fraction * double(std::max(n_old, n_new))) {
        tree = InstanceTree<T>::build(std::move(insts));
        looks = std::move(new_looks);
        u.rebuilt = true;
        return u;
    }
    std::vector<InstanceOrigin> none(n_old);
    edit_instance_tree<T>(tree, looks, none, old_kept);
    for (std::size_t k = 0; k < n_new; ++k) {
        if (new_kept[k]) continue;
        tree.insert(insts[k]);
        looks.push_back(new_looks[k]);
    }
    return u;
}

template<Scalar T>
CookedScene<T> lay_out_impl(CookedScene<T>* previous, InstanceUpdate* report,
                            const io::build::Cooked<T>& cooked,
                            const std::vector<io::build::Placed<T>>& extra,
                            const std::vector<char>* skip_nodes, StaticCache<T>* cache) {
    using Tri = geometry::Triangle<3, T>;
    using Quadric = geometry::BoundedQuadric<T>;
    using Inst = geometry::Instanced<Quadric>;

    CookedScene<T> out;
    std::vector<Tri> tris;
    std::vector<Inst> insts;
    std::vector<InstanceLook<T>> looks;
    std::vector<InstanceOrigin> origins;

    // One quadric per distinct *shape*, not per object -- which is the
    // whole point of reading a cooked scene rather than a list of
    // materialized nodes. Sized up front and never grown, because every
    // instance below holds a pointer into it.
    out.quadrics.resize(cooked.shape_count());
    std::vector<char> has_quadric(cooked.shape_count(), 0);
    std::vector<std::vector<std::array<Vec<T, 3>, 3>>> shape_normals(cooked.shape_count());
    for (std::size_t i = 0; i < cooked.shape_count(); ++i) {
        const auto& sh = cooked.shapes()[i];
        if (const auto* q = std::any_cast<Quadric>(&sh.exact)) {
            out.quadrics[i] = *q;
            has_quadric[i] = 1;
        }
    }

    // The previous scene's instance tree points into its quadric table. If
    // the table is the same, this scene points into that same buffer, so
    // the tree stays valid and can be updated instead of built.
    const bool reuse = previous && detail::same_quadrics(previous->quadrics, out.quadrics);
    if (reuse) out.quadrics = std::move(previous->quadrics);

    const auto& objs = cooked.objects();
    auto skipped = [&](std::size_t node) { return skip_nodes && node < skip_nodes->size() && (*skip_nodes)[node]; };

    // The generators whose instances can be known by their generator: a
    // scatter with a version, all of whose objects are laid out as
    // instances. An object is not laid out when its scale is zero, and then
    // the instances would no longer be "the first count sites", so such a
    // generator is treated like the objects of a node with no version.
    const auto& gens = cooked.generators();
    constexpr std::uint32_t no_gen = InstanceOrigin::none;
    std::vector<std::uint32_t> laid_index(gens.size(), no_gen);   // gens[i] -> out.generators
    for (std::size_t i = 0; i < gens.size(); ++i) {
        const auto& g = gens[i];
        if (skipped(g.node)) continue;
        if (g.count > 0 && !(objs[g.first_object].instanceable && has_quadric[g.shape])) continue;
        bool laid = true;
        for (std::size_t o = g.first_object; o < g.first_object + g.count && laid; ++o)
            laid = objs[o].scale != T{0};
        if (!laid) continue;
        laid_index[i] = static_cast<std::uint32_t>(out.generators.size());
        out.generators.push_back({g.version.generator, g.version.key, g.shape, g.count});
    }

    // A shape's geometry is its *rest* form when the object placing it is
    // instanceable, and its *placed* form when the object deforms -- in
    // which case cook() hands back an identity transform. So one formula
    // covers both, and the deforming case is not a special case here.
    // Takes the rotation already expanded, rather than reaching into the
    // object for it: this runs per *vertex*, and a quaternion unpacked
    // here would spend the arithmetic the compact storage was supposed to
    // be paying for. The expansion happens once per object, below.
    auto to_world = [](const Vec<T, 3>& v, const Matrix<T, 3, 3>& R,
                       const io::build::Object<T>& o) {
        return Vec<T, 3>{R * Vec<T, 3>{v * o.scale} + o.translation};
    };

    auto emit_triangles = [&](const mesh::Mesh<Euclidean<3, T>>& m,
                              const std::vector<std::array<Vec<T, 3>, 3>>& cn,
                              const io::build::Object<T>& o, const Matrix<T, 3, 3>& R,
                              const io::Material<T>& mat) {
        for (std::size_t i = 0; i < m.faces.size(); ++i) {
            const auto& f = m.faces[i];
            Tri t(to_world(m.vertices[f[0]], R, o), to_world(m.vertices[f[1]], R, o),
                  to_world(m.vertices[f[2]], R, o));
            tris.push_back(t);
            // Normals turn by the rotation alone: the scale is uniform, so
            // it divides out of the inverse transpose.
            out.prims.push_back(ShadedTriangle<T>{t,
                                                  {Vec<T, 3>{R * cn[i][0]},
                                                   Vec<T, 3>{R * cn[i][1]},
                                                   Vec<T, 3>{R * cn[i][2]}},
                                                  mat.base_color, mat.roughness, mat.emissive});
        }
    };

    // Once per object. Everything below uses this, including the
    // instance, whose rotation stays a matrix precisely so the
    // traversal never has to unpack anything.
    auto instance_of = [&](const io::build::Object<T>& obj) {
        const auto& mat = obj.material;
        return std::pair<Inst, InstanceLook<T>>{
            Inst{&out.quadrics[obj.shape], obj.translation, obj.scale, obj.rotation_q.to_matrix()},
            InstanceLook<T>{mat.base_color, mat.roughness, mat.emissive, mat.opacity}};
    };

    // When the previous scene can be updated, the instances of a laid
    // generator are decided below without visiting them here.
    std::size_t gi = 0;
    for (std::size_t o = 0; o < objs.size(); ++o) {
        const auto& obj = objs[o];
        while (gi < gens.size() && gens[gi].first_object + gens[gi].count <= o) ++gi;
        const bool in_gen = gi < gens.size() && o >= gens[gi].first_object;
        if (in_gen && reuse && laid_index[gi] != no_gen) {
            o = gens[gi].first_object + gens[gi].count - 1;
            continue;
        }

        if (skipped(obj.source_node)) continue;
        const auto& mat = obj.material;
        const Matrix<T, 3, 3> R = obj.rotation_q.to_matrix();

        // Scaled to nothing is not there. The donut's dust dissolves by
        // scaling to zero, so at t=3.9 35 200 of its 46 192 objects were
        // still in the instance tree as points: never hit -- `ray_hit`
        // refuses a zero scale -- and still walked, 31 nodes a pixel
        // against 4.6 for a scene without them.
        if (obj.scale == T{0}) continue;

        // Two conditions, and neither is a guess about what the object
        // looks like: the shape must carry a closed form, and the object's
        // motion must be a placement so copies differ only by where they
        // are.
        if (obj.instanceable && has_quadric[obj.shape]) {
            insts.push_back({&out.quadrics[obj.shape], obj.translation, obj.scale, R});
            looks.push_back({mat.base_color, mat.roughness, mat.emissive, mat.opacity});
            origins.push_back(in_gen && laid_index[gi] != no_gen
                                  ? InstanceOrigin{laid_index[gi], obj.site}
                                  : InstanceOrigin{});
            continue;
        }

        // Smooth normals once per shape rather than once per object: a
        // thousand instances of one sprinkle used to recompute the same
        // normals a thousand times, through a merged mesh that had to be
        // built first.
        auto& vn = shape_normals[obj.shape];
        if (vn.empty()) vn = corner_normals(cooked.shapes()[obj.shape].geometry);
        emit_triangles(cooked.shapes()[obj.shape].geometry, vn, obj, R, mat);
    }

    for (const auto& g : extra) {
        auto m = g.mesh();
        auto cn = corner_normals(m);
        auto mat = g.material();
        for (std::size_t i = 0; i < m.faces.size(); ++i) {
            const auto& f = m.faces[i];
            Tri t(m.vertices[f[0]], m.vertices[f[1]], m.vertices[f[2]]);
            tris.push_back(t);
            out.prims.push_back(ShadedTriangle<T>{t, cn[i], mat.base_color, mat.roughness,
                                                  mat.emissive});
        }
    }

    // Moved, not copied. `build` takes its shapes by value and keeps
    // them, so handing it an lvalue leaves two copies of the array alive
    // for the rest of the frame -- 112 bytes per instance, which is 214 MB
    // at two million and was simply being spent.
    out.triangles = cache ? cache->tree(std::move(tris)) : spatial::BVH<Tri>::build(std::move(tris));
    if (!reuse) {
        out.instances = InstanceTree<T>::build(std::move(insts));
        out.looks = std::move(looks);
        out.origins = std::move(origins);
        if (report) {
            *report = {};
            report->rebuilt = true;
            report->inserted = out.instances.size();
            report->hashed = 0;
        }
        return out;
    }

    // ── Updating the previous scene's instance tree ──────────────────
    //
    // `insts` holds only the instances with no generator to compare. For a
    // laid generator the previous scene may hold one with the same
    // identity; if so its first min(old, new) instances are the new ones,
    // by construction, and only the rest is removed or added. A generator
    // with no counterpart replaces its old instances whole.
    InstanceUpdate u;
    u.reused = true;
    out.instances = std::move(previous->instances);
    out.looks = std::move(previous->looks);
    out.origins = std::move(previous->origins);
    const auto& old_gens = previous->generators;

    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> old_by_id;
    for (std::size_t i = 0; i < old_gens.size(); ++i)
        old_by_id[old_gens[i].id].push_back(static_cast<std::uint32_t>(i));
    std::vector<std::uint32_t> old_to_new(old_gens.size(), no_gen);
    std::vector<std::size_t> new_prefix(out.generators.size(), 0);   // known sites of each new generator
    std::vector<char> claimed(old_gens.size(), 0);
    std::vector<std::size_t> keep_below(old_gens.size(), 0);         // old sites to keep, per old generator
    for (std::size_t j = 0; j < out.generators.size(); ++j) {
        const auto& ng = out.generators[j];
        const auto it = old_by_id.find(ng.id);
        if (it == old_by_id.end()) continue;
        for (auto i : it->second) {
            // The key decides, not the id: equal ids are a candidate.
            if (claimed[i] || old_gens[i].shape != ng.shape || old_gens[i].key != ng.key) continue;
            claimed[i] = 1;
            old_to_new[i] = static_cast<std::uint32_t>(j);
            new_prefix[j] = std::min(old_gens[i].count, ng.count);
            keep_below[i] = new_prefix[j];
            ++u.generators_matched;
            break;
        }
    }

    const std::size_t n_old = out.instances.size();
    std::vector<char> old_kept(n_old, 0);
    std::vector<std::uint32_t> old_plain;
    for (std::size_t s = 0; s < n_old; ++s) {
        auto& og = out.origins[s];
        if (og.generator == no_gen) { old_plain.push_back(static_cast<std::uint32_t>(s)); continue; }
        if (og.site < keep_below[og.generator]) {
            old_kept[s] = 1;
            og.generator = old_to_new[og.generator];
        }
    }

    std::vector<char> new_kept(insts.size(), 0);
    std::size_t plain_kept = 0;
    if (!old_plain.empty() && !insts.empty())
        plain_kept = match_instances<T>(out.instances.shapes(), out.looks, &old_plain, insts, looks, old_kept, new_kept);
    u.hashed = old_plain.size() + insts.size();

    std::size_t tails = 0, new_total = insts.size();
    for (std::size_t j = 0; j < out.generators.size(); ++j) {
        tails += out.generators[j].count - new_prefix[j];
        new_total += out.generators[j].count;
    }
    std::size_t kept_old = 0;
    for (char c : old_kept) kept_old += c;
    u.kept = kept_old;
    u.removed = n_old - kept_old;
    u.inserted = (insts.size() - plain_kept) + tails;

    if (double(u.removed + u.inserted) > instance_update_max_fraction * double(std::max(n_old, new_total))) {
        // Built again, from every instance of the new scene.
        for (std::size_t i = 0; i < gens.size(); ++i) {
            if (laid_index[i] == no_gen) continue;
            const auto& g = gens[i];
            for (std::size_t site = 0; site < g.count; ++site) {
                auto [inst, look] = instance_of(objs[g.first_object + site]);
                insts.push_back(inst);
                looks.push_back(look);
                origins.push_back({laid_index[i], static_cast<std::uint32_t>(site)});
            }
        }
        out.instances = InstanceTree<T>::build(std::move(insts));
        out.looks = std::move(looks);
        out.origins = std::move(origins);
        u.rebuilt = true;
        if (report) *report = u;
        return out;
    }

    edit_instance_tree<T>(out.instances, out.looks, out.origins, old_kept);
    for (std::size_t k = 0; k < insts.size(); ++k) {
        if (new_kept[k]) continue;
        out.instances.insert(insts[k]);
        out.looks.push_back(looks[k]);
        out.origins.push_back(origins[k]);
    }
    for (std::size_t i = 0; i < gens.size(); ++i) {
        if (laid_index[i] == no_gen) continue;
        const auto& g = gens[i];
        for (std::size_t site = new_prefix[laid_index[i]]; site < g.count; ++site) {
            auto [inst, look] = instance_of(objs[g.first_object + site]);
            out.instances.insert(inst);
            out.looks.push_back(look);
            out.origins.push_back({laid_index[i], static_cast<std::uint32_t>(site)});
        }
    }
    if (report) *report = u;
    return out;
}

}  // namespace detail

template<Scalar T>
CookedScene<T> lay_out(const io::build::Cooked<T>& cooked, const std::vector<io::build::Placed<T>>& extra,
                       const std::vector<char>* skip_nodes, StaticCache<T>* cache) {
    return detail::lay_out_impl<T>(nullptr, nullptr, cooked, extra, skip_nodes, cache);
}

// The same scene as `lay_out` makes, starting from the scene laid out
// before it: when the instance set changed by a few, the instance tree is
// updated (see BVH::insert / remove) instead of built, which is the part of
// laying out that grows with the number of instances. `previous` is
// consumed; it is the caller's last scene, replaced by this one.
//
// The result gives the same ray hits as a fresh `lay_out`. The instances
// are in a different order (removal moves the last one into the freed
// slot), with `looks` following; nothing else about the scene differs.
// `report`, when given, says what was done.
template<Scalar T>
CookedScene<T> lay_out(CookedScene<T>&& previous, const io::build::Cooked<T>& cooked,
                       const std::vector<io::build::Placed<T>>& extra = {},
                       const std::vector<char>* skip_nodes = nullptr, StaticCache<T>* cache = nullptr,
                       InstanceUpdate* report = nullptr) {
    return detail::lay_out_impl<T>(&previous, report, cooked, extra, skip_nodes, cache);
}

}  // namespace spatium::render
