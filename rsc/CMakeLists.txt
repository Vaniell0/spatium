# rsc_route: the routes to an operation on a type of the library, read off rsc/spatium/graph.json.
# Not part of ALL: build it by name (`ninja -C build rsc_route`). It includes rsc/core only; the
# library is not linked, because the graph is data.
add_executable(rsc_route EXCLUDE_FROM_ALL spatium/rsc_route.cpp)
target_include_directories(rsc_route PRIVATE ${PROJECT_SOURCE_DIR})
target_compile_features(rsc_route PRIVATE cxx_std_23)
target_compile_definitions(rsc_route PRIVATE
    RSC_GRAPH_JSON="${PROJECT_SOURCE_DIR}/rsc/spatium/graph.json")
