include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

# Phase 0: only the slice of pure-C++ source that has no pybind11 dependency.
# `BitVect` lives in src/algebra/bitvect.cpp. Phase 1 promotes a real
# `regpoly_core` library and grows this list.
add_executable(regpoly_cpp_tests
    test_smoke.cpp
    ../src/algebra/bitvect.cpp
)
target_include_directories(regpoly_cpp_tests PRIVATE ../src/include/algebra)
target_link_libraries(regpoly_cpp_tests PRIVATE GTest::gtest_main)

gtest_discover_tests(regpoly_cpp_tests)
