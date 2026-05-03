cmake_minimum_required(VERSION 3.20)
project(regpoly_cpp VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(REGPOLY_BUILD_TESTS "Build C++ unit tests via GoogleTest" OFF)
option(REGPOLY_INSTALL_CMAKE_PACKAGE
       "Install regpolyConfig.cmake for find_package(regpoly)" ON)
option(REGPOLY_BUILD_CLI
       "Build the regpoly-cli executable (Phase 4 fills in subcommands)" ON)

find_package(pybind11 CONFIG REQUIRED)

# All C++ source files except factory.cpp (which carries pybind11
# bindings registration). regpoly_core is the pure-C++ static library
# that downstream C++ users link against.
file(GLOB_RECURSE REGPOLY_ALL_SOURCES
     CONFIGURE_DEPENDS
     src/algebra/*.cpp
     src/analyses/*.cpp
     src/core/*.cpp
     src/generators/*.cpp
     src/transforms/*.cpp
     src/lattice/*.cpp)

set(REGPOLY_CORE_SOURCES "")
foreach(src IN LISTS REGPOLY_ALL_SOURCES)
    if(NOT src MATCHES "src/core/factory\\.cpp$")
        list(APPEND REGPOLY_CORE_SOURCES "${src}")
    endif()
endforeach()

set(REGPOLY_CPP_INCLUDE_DIRS
    src/include/core
    src/include/algebra
    src/include/analyses
    src/include/generators
    src/include/transforms
    src/include/lattice)

# Public-API include directory for downstream consumers. The source
# tree ships `include/regpoly/regpoly.h` directly; the per-class
# headers (bitvect.h, generator.h, …) are mirrored from src/include/
# at config time so the BUILD_INTERFACE include path matches the
# INSTALL_INTERFACE layout.
set(REGPOLY_PUBLIC_INCLUDE_DIR include)
set(REGPOLY_PUBLIC_HEADER_BASES
    src/include/algebra/bitvect.h
    src/include/algebra/primitivity.h
    src/include/core/generator.h
    src/include/core/transformation.h
    src/include/core/params.h
    src/include/algebra/gauss.h
    src/include/analyses/resolution_sets.h
    src/include/analyses/equidistribution_runner.h
    src/include/analyses/tuplets_runner.h
    src/include/generators/combined.h
    src/include/lattice/me_helpers.h
    src/include/lattice/me_harase.h
    src/include/lattice/me_notprimitive.h
    src/include/lattice/me_notprimitive_simd.h
    src/include/lattice/temper_optimizer.h
    src/include/lattice/dual_lattice.h
)
foreach(hdr IN LISTS REGPOLY_PUBLIC_HEADER_BASES)
    get_filename_component(name "${hdr}" NAME)
    configure_file("${hdr}"
                   "${CMAKE_CURRENT_BINARY_DIR}/include/regpoly/${name}"
                   COPYONLY)
endforeach()

# ── Pure-C++ static library: regpoly_core ───────────────────────────────
# Linked by the Python extension, the regpoly-cli binary, and any
# downstream C++ project via find_package(regpoly).

add_library(regpoly_core STATIC ${REGPOLY_CORE_SOURCES})
target_include_directories(regpoly_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${REGPOLY_PUBLIC_INCLUDE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${REGPOLY_CPP_INCLUDE_DIRS})
target_link_libraries(regpoly_core PUBLIC ntl pthread)
target_compile_options(regpoly_core PRIVATE -O3 -march=native)
add_library(regpoly::regpoly ALIAS regpoly_core)

# ── Python extension: _regpoly_cpp ─────────────────────────────────────
# Adds factory.cpp (pybind11 binding registration) and bindings.cpp on
# top of regpoly_core.

pybind11_add_module(_regpoly_cpp
    src/bindings/bindings.cpp
    src/core/factory.cpp)
target_include_directories(_regpoly_cpp PRIVATE ${REGPOLY_CPP_INCLUDE_DIRS})
target_link_libraries(_regpoly_cpp PRIVATE regpoly_core)
target_compile_options(_regpoly_cpp PRIVATE -O3 -march=native)

install(TARGETS _regpoly_cpp DESTINATION regpoly_cpp)

# ── regpoly-cli executable (Phase 4 fills in subcommands) ───────────────

if(REGPOLY_BUILD_CLI)
    add_executable(regpoly_cli src/cli/main.cpp)
    set_target_properties(regpoly_cli PROPERTIES OUTPUT_NAME regpoly-cli)
    target_link_libraries(regpoly_cli PRIVATE regpoly_core)
    install(TARGETS regpoly_cli RUNTIME DESTINATION bin)
endif()

# ── Public-API headers ──────────────────────────────────────────────────
# Installed under ${CMAKE_INSTALL_PREFIX}/include/regpoly/.

set(REGPOLY_PUBLIC_HEADERS
    include/regpoly/regpoly.h
    src/include/algebra/bitvect.h
    src/include/algebra/primitivity.h
    src/include/core/generator.h
    src/include/core/transformation.h
    src/include/core/params.h
    src/include/algebra/gauss.h
    src/include/analyses/resolution_sets.h
    src/include/analyses/equidistribution_runner.h
    src/include/analyses/tuplets_runner.h
    src/include/generators/combined.h
    src/include/lattice/me_helpers.h
    src/include/lattice/me_harase.h
    src/include/lattice/me_notprimitive.h
    src/include/lattice/me_notprimitive_simd.h
    src/include/lattice/temper_optimizer.h
    src/include/lattice/dual_lattice.h
)
install(FILES ${REGPOLY_PUBLIC_HEADERS}
        DESTINATION include/regpoly)

# ── CMake config package: find_package(regpoly) ────────────────────────

if(REGPOLY_INSTALL_CMAKE_PACKAGE)
    install(TARGETS regpoly_core
            EXPORT regpolyTargets
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib
            RUNTIME DESTINATION bin
            INCLUDES DESTINATION include)

    install(EXPORT regpolyTargets
            FILE regpolyTargets.cmake
            NAMESPACE regpoly::
            DESTINATION lib/cmake/regpoly)

    include(CMakePackageConfigHelpers)
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/regpolyConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)

    configure_package_config_file(
        cmake/regpolyConfig.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/regpolyConfig.cmake"
        INSTALL_DESTINATION lib/cmake/regpoly)

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/regpolyConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/regpolyConfigVersion.cmake"
        DESTINATION lib/cmake/regpoly)
endif()

# ── Tests ──────────────────────────────────────────────────────────────

if(REGPOLY_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
