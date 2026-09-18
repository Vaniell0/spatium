add_executable(train_base train_base.cpp)
target_include_directories(train_base PRIVATE
    ${CMAKE_SOURCE_DIR}/rsc/include
    ${CMAKE_BINARY_DIR}/generated
)
target_link_libraries(train_base PRIVATE spatium_core)
target_compile_options(train_base PRIVATE -Wall -Wextra -Wpedantic)

add_executable(search_capacity search_capacity.cpp)
target_include_directories(search_capacity PRIVATE
    ${CMAKE_SOURCE_DIR}/rsc/include
    ${CMAKE_BINARY_DIR}/generated
)
target_link_libraries(search_capacity PRIVATE spatium_core)
target_compile_options(search_capacity PRIVATE -Wall -Wextra -Wpedantic)

add_executable(dispatch_cost dispatch_cost.cpp)
target_include_directories(dispatch_cost PRIVATE
    ${CMAKE_SOURCE_DIR}/rsc/include
    ${CMAKE_BINARY_DIR}/generated
)
target_link_libraries(dispatch_cost PRIVATE spatium_core)
target_compile_options(dispatch_cost PRIVATE -Wall -Wextra -Wpedantic)

add_executable(distill_tree distill_tree.cpp)
target_include_directories(distill_tree PRIVATE
    ${CMAKE_SOURCE_DIR}/rsc/include
    ${CMAKE_BINARY_DIR}/generated
)
target_link_libraries(distill_tree PRIVATE spatium_core)
target_compile_options(distill_tree PRIVATE -Wall -Wextra -Wpedantic)

add_executable(search_heuristic search_heuristic.cpp)
target_include_directories(search_heuristic PRIVATE
    ${CMAKE_SOURCE_DIR}/rsc/include
    ${CMAKE_BINARY_DIR}/generated
)
target_link_libraries(search_heuristic PRIVATE spatium_core)
target_compile_options(search_heuristic PRIVATE -Wall -Wextra -Wpedantic)
