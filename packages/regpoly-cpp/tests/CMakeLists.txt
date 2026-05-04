include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

# regpoly_core is defined at the parent CMakeLists.txt scope; we just
# link against it here. The test target additionally needs access to
# private internal headers (concrete generator classes like
# TauswortheGen) to construct fixtures — those headers are NOT part of
# the public API surface.
add_executable(regpoly_cpp_tests
    test_smoke.cpp
    test_combination.cpp
    test_combined_generator.cpp
    test_primitive_search.cpp
    test_primitivity.cpp
    test_random_samplers.cpp
    test_resolution_sets.cpp
    test_seek_search.cpp
    test_tempering_optimizer.cpp
)
target_include_directories(regpoly_cpp_tests PRIVATE
    ../src/include/algebra
    ../src/include/core
    ../src/include/generators
    ../src/include/transforms
    ../src/include/lattice
    ../src/include/search)
target_link_libraries(regpoly_cpp_tests PRIVATE regpoly_core GTest::gtest_main)

gtest_discover_tests(regpoly_cpp_tests)
