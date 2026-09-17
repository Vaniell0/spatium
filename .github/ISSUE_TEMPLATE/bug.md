---
name: Something is wrong
about: A wrong answer, a crash, a build failure, or a picture that is not what it should be
labels: bug
---

## What happened

<!-- The observed behaviour. Exact error text if there is any and it is short. -->

## What you expected

<!--
For a numerical result, the expected value and where it comes from -- a
closed form, a reference implementation, a paper -- is worth more than
"it looks wrong". Spatium's own bugs have mostly been found by someone
knowing what the answer should have been.
-->

## The smallest thing that reproduces it

```cpp
// A `main()` that shows it, ideally one that only includes <spatium/...>.
```

## Build

- Compiler and version:
- CMake options (`-DSPATIUM_EIGEN=ON`, `-DSPATIUM_BOOST=ON`, modules on/off):
- Nix or not:

<!--
Optional dependencies matter more here than in most libraries: the heat
method, DEC and geodesics need `SPATIUM_EIGEN=ON` and simply are not
present without it, and `Real50`/`Real100` need `SPATIUM_BOOST=ON`.
-->
