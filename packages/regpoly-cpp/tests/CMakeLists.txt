include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

# Phase 1: promote a real `regpoly_core` static library that the test
# binary, the Python extension, and the future regpoly-cli all link
# against. The static library has no pybind11 dependency at the source
# level — factory.cpp's `register_generator_types` is the only file that
# imports pybind11; we exclude it from regpoly_core and keep it in the
# bindings translation unit.

set(REGPOLY_CORE_SOURCES "")
file(GLOB_RECURSE REGPOLY_CORE_SOURCES_ALL
    CONFIGURE_DEPENDS
    ../src/algebra/*.cpp
    ../src/core/*.cpp
    ../src/generators/*.cpp
    ../src/transforms/*.cpp
    ../src/lattice/*.cpp)
foreach(src IN LISTS REGPOLY_CORE_SOURCES_ALL)
    if(NOT src MATCHES "src/core/factory\\.cpp$")
        list(APPEND REGPOLY_CORE_SOURCES "${src}")
    endif()
endforeach()

add_library(regpoly_core STATIC ${REGPOLY_CORE_SOURCES})
target_include_directories(regpoly_core PUBLIC
    ../src/include/algebra
    ../src/include/core
    ../src/include/generators
    ../src/include/transforms
    ../src/include/lattice)
target_link_libraries(regpoly_core PUBLIC ntl pthread)
target_compile_options(regpoly_core PRIVATE -O3 -march=native)

add_executable(regpoly_cpp_tests
    test_smoke.cpp
    test_combined_generator.cpp
)
target_link_libraries(regpoly_cpp_tests PRIVATE regpoly_core GTest::gtest_main)

gtest_discover_tests(regpoly_cpp_tests)
