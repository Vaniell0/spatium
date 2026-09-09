# Getting Started: The Declarative Scene DSL

Every 3D beginner eventually makes a donut — it's the subject of the most
famous introductory tutorial in the field. This is the same idea for
Spatium's `spatium::io::build` DSL: three lines of code, no loops, and by
the end you have dough, icing, and sprinkles. If you've never written C++
before, you can still follow this — read each step, run it, see what
changed.

The full, runnable version of this file is `examples/donut_demo.cpp`.

## Step 0: what you're building with

`torus()`, `offset()`, and `scatter()` don't build anything by themselves.
They add one entry to a `Trace` — a plain list, sitting in memory, that
you can look at. Nothing gets computed until you ask for it. This matters
for one reason: it means the *description* of your scene and the
*picture* of your scene are two different things, and you can inspect the
first one before you ever pay for the second.

```cpp
#include <spatium/io/build.hpp>
namespace bd = spatium::io::build;

bd::Trace<double> scene;
```

`scene` is empty right now — zero nodes.

## Step 0: delete the default cube

Every Blender tutorial starts here: a cube sits in the scene, and the very
first thing you do is get rid of it. Same idea, done declaratively:

```cpp
auto cube = scene.cube({0.9, 0.9, 0.9})
                .moving([](const Vec<double, 3>& p, double t) {
                    double e = /* smoothstep(t) */;
                    return p * (1.0 - e); // shrinks to nothing as e -> 1
                });
```

`.moving(f)` is the one motion/mutation slot every node has: a function
from `(point, time)` to a new point. Here it's a plain shrink, so at
`t=0` the cube is there and at `t=1` it's gone — no per-frame state, no
loop stepping a simulation forward; `materialize()` can ask for any `t`
directly and get the right answer. (`examples/donut_demo.cpp`'s version
adds a little `spatium::algebra::PerlinNoise` turbulence on the way out,
so it wobbles instead of shrinking in a perfectly straight line — same
idea, one more term.)

## Step 1: the dough is a torus

```cpp
auto dough = scene.torus(2.0, 1.0)
                 .colored({.base_color = {0.80, 0.55, 0.32}});
```

`2.0` is the donut's overall radius, `1.0` is the tube's radius.
`.colored(...)` sets the material. `dough` is not a shape — it's a
*handle*, a name pointing at node `0` in `scene`. Nothing has been drawn.

## Step 2: the icing is the dough's own surface, pushed outward

```cpp
auto icing = scene.offset(dough, 0.035)
                 .colored({.base_color = {0.95, 0.72, 0.86}});
```

`offset()` doesn't build a second shape and glue it onto the first. It
reads the dough's surface formula and adds a little in the direction the
dough's own surface is already facing (its normal) at every point. The
icing is a *function of the dough* — if you change the dough later, the
icing definition doesn't need to change with it.

(An earlier draft of this demo built icing as an entirely separate torus,
just slightly bigger, and projected it onto the dough. That's redundant —
two copies of the same formula doing the work of one — and it's the
mistake `offset()` exists to make impossible.)

## Step 3: sprinkles scatter across the icing

```cpp
auto sprinkle  = scene.cylinder(0.025, 0.12, 6, 2);
auto sprinkles = scene.scatter(sprinkle, icing, 200)
                     .colored({.base_color = {0.95, 0.20, 0.25}});
```

`scatter()` places 200 copies of the tiny cylinder across the icing's
surface. "Across the surface" is measured honestly — by actual surface
area, not by counting evenly in some coordinate grid, which would bunch
sprinkles up wherever the parametrization happens to compress space (the
donut's inner rim, for instance). No mesh gets built to figure this out;
it's read directly off the icing's own geometry.

## Step 4: group them, then look

```cpp
auto donut = scene.compose({dough, icing, sprinkles});
```

`scene` now has 5 nodes: dough, icing, the sprinkle template, the
scattered sprinkles, and this group. You can print them:

```cpp
for (std::size_t i = 0; i < scene.size(); ++i)
    std::println("[{}] {}", i, kind_name(scene.node(i).kind));
```

```
[0] Space
[1] Offset
[2] Space
[3] Scatter
[4] Compose
```

That's the whole scene, as data, before a single triangle exists.

## Step 5: get a picture (optional)

Everything above works with no GPU and no window — that's the point of
staying declarative this long. When you actually want pixels:

```cpp
auto placed = bd::materialize(scene, donut.index);
```

`placed` is a list of real, triangulated meshes — this is the one place
in the whole pipeline where triangles appear, and only because a screen
needs them. `examples/donut_demo.cpp` feeds `placed` into the same
CPU raytracer every offline demo in this repository uses. Build and run
it yourself:

```bash
./build/examples/donut_demo              # prints the trace, no picture
./build/examples/donut_demo --photo      # also writes donut.png
```

## Why bother being declarative

Nothing above required you to write a loop over 200 sprinkles yourself —
`scatter(sprinkle, icing, 200)` is the whole instruction, and the
`Trace` remembers exactly what you asked for. A `for` loop that pushed
200 cylinders into a vector would get you the same picture, but the
*intent* — "200 things, spread evenly, on that surface" — would be gone,
baked into one-off arithmetic nobody can look at afterward and reuse. The
declarative form keeps that intent around as data, which is what makes it
possible to inspect, replay, or (later) compile.
