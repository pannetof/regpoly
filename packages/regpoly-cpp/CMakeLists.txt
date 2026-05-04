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

# ── yaml-cpp (Phase 3): vendored via FetchContent ─────────────────────
# Used by the C++ Catalog and YAML config parsers. Tests are off; we
# don't need yaml-cpp's own test suite. Built static so it can be
# linked into regpoly_core without an extra runtime dependency.
include(FetchContent)
set(YAML_CPP_BUILD_TESTS  OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS  OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
set(YAML_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG        0.8.0
)
FetchContent_MakeAvailable(yaml-cpp)

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
     src/lattice/*.cpp
     src/library/*.cpp
     src/io/*.cpp
     src/search/*.cpp)

# factory.cpp moved into regpoly_core in Phase 2.4 (pybind11 wiring
# now lives in src/bindings/factory_bindings.cpp). Every .cpp under
# src/{algebra,analyses,core,generators,transforms,lattice,search}/
# now contributes to regpoly_core.
set(REGPOLY_CORE_SOURCES ${REGPOLY_ALL_SOURCES})

set(REGPOLY_CPP_INCLUDE_DIRS
    src/include/core
    src/include/algebra
    src/include/analyses
    src/include/generators
    src/include/transforms
    src/include/lattice
    src/include/library
    src/include/io
    src/include/search)

# Public-API include directory for downstream consumers. The source
# tree ships `include/regpoly/regpoly.h` directly; the per-class
# headers (bitvect.h, generator.h, …) are mirrored from src/include/
# at config time so the BUILD_INTERFACE include path matches the
# INSTALL_INTERFACE layout.
set(REGPOLY_PUBLIC_INCLUDE_DIR include)
set(REGPOLY_PUBLIC_HEADER_BASES
    src/include/algebra/bitvect.h
    src/include/algebra/primitivity.h
    src/include/core/combination.h
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
    src/include/search/search_types.h
    src/include/search/primitive_search.h
    src/include/search/seek_search.h
    src/include/search/tempering_optimizer.h
    src/include/search/tempering_search.h
    src/include/library/catalog.h
    src/include/io/legacy_reader.h
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
# yaml-cpp is an implementation detail of the catalog / config readers;
# wrapped in $<BUILD_INTERFACE:> so it is linked into the regpoly_core
# static archive at build time but is NOT re-exported as an
# INTERFACE_LINK_LIBRARIES entry that downstream find_package(regpoly)
# consumers would need to satisfy. (Static archives bake the .o files
# from yaml-cpp in directly, so consumers see no missing symbols.)
target_link_libraries(regpoly_core PRIVATE
    $<BUILD_INTERFACE:yaml-cpp::yaml-cpp>)
target_compile_options(regpoly_core PRIVATE -O3 -march=native)
add_library(regpoly::regpoly ALIAS regpoly_core)

# ── Python extension: _regpoly_cpp ─────────────────────────────────────
# Adds factory.cpp (pybind11 binding registration) and bindings.cpp on
# top of regpoly_core.

pybind11_add_module(_regpoly_cpp
    src/bindings/bindings.cpp
    src/bindings/factory_bindings.cpp)
target_include_directories(_regpoly_cpp PRIVATE ${REGPOLY_CPP_INCLUDE_DIRS})
target_link_libraries(_regpoly_cpp PRIVATE regpoly_core)
target_compile_options(_regpoly_cpp PRIVATE -O3 -march=native)

install(TARGETS _regpoly_cpp DESTINATION regpoly_cpp)

# ── regpoly-cli executable (Phase 4 fills in subcommands) ───────────────

if(REGPOLY_BUILD_CLI)
    add_executable(regpoly_cli src/cli/main.cpp)
    set_target_properties(regpoly_cli PROPERTIES OUTPUT_NAME regpoly-cli)
    target_include_directories(regpoly_cli PRIVATE ${REGPOLY_CPP_INCLUDE_DIRS})
    target_link_libraries(regpoly_cli PRIVATE regpoly_core)
    install(TARGETS regpoly_cli RUNTIME DESTINATION bin)
endif()

# ── Public-API headers ──────────────────────────────────────────────────
# Installed under ${CMAKE_INSTALL_PREFIX}/include/regpoly/.

set(REGPOLY_PUBLIC_HEADERS
    include/regpoly/regpoly.h
    src/include/algebra/bitvect.h
    src/include/algebra/primitivity.h
    src/include/core/combination.h
    src/include/search/search_types.h
    src/include/search/primitive_search.h
    src/include/search/seek_search.h
    src/include/search/tempering_optimizer.h
    src/include/search/tempering_search.h
    src/include/library/catalog.h
    src/include/io/legacy_reader.h
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
