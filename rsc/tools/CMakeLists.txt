add_executable(train_base train_base.cpp)
target_include_directories(train_base PRIVATE
    ${CMAKE_SOURCE_DIR}/rsc/include
    ${CMAKE_BINARY_DIR}/generated
)
target_link_libraries(train_base PRIVATE spatium_core)
target_compile_options(train_base PRIVATE -Wall -Wextra -Wpedantic)
