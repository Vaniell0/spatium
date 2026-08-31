#pragma once
// SPATIUM_EXPORT — single-source-of-truth macro for the modules migration.
//
// In legacy (header) builds, SPATIUM_BUILDING_MODULE is undefined and
// SPATIUM_EXPORT expands to nothing — the header behaves as a normal #include.
//
// Inside a module unit (.cppm), the file defines SPATIUM_BUILDING_MODULE
// before #including the corresponding .hpp; SPATIUM_EXPORT then expands to
// `export`, marking the wrapped namespace block as part of the module's
// exported interface.
//
// This lets every header serve as both a legacy include and the source of
// truth for its module partition — no code duplication during migration.

#ifndef SPATIUM_EXPORT
#  ifdef SPATIUM_BUILDING_MODULE
#    define SPATIUM_EXPORT export
#  else
#    define SPATIUM_EXPORT
#  endif
#endif
