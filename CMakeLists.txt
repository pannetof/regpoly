cmake_minimum_required(VERSION 3.20)
project(regpoly_cpp VERSION 2.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(REGPOLY_BUILD_TESTS "Build C++ unit tests via GoogleTest" OFF)
option(REGPOLY_INSTALL_CMAKE_PACKAGE
       "Install regpolyConfig.cmake for find_package(regpoly)" ON)
option(REGPOLY_BUILD_CLI
       "Build the regpoly-cli executable (Phase 4 fills in subcommands)" ON)
# Toggle the pybind11 extension. ON for the scikit-build-core / uv path
# (the Python wheel ships `_regpoly_cpp.so`). OFF for pure-C++ consumers
# who only want the static library + `regpoly-cli` and don't have pybind11
# (or Python) available. When OFF, `find_package(pybind11)` is not called.
option(REGPOLY_BUILD_PYTHON_EXTENSION
       "Build the _regpoly_cpp pybind11 extension module" ON)

if(REGPOLY_BUILD_PYTHON_EXTENSION)
    find_package(pybind11 CONFIG REQUIRED)
endif()

# ── yaml-cpp: fetched + pinned at configure time. Was a git submodule
# (third_party/yaml-cpp) until the two-repo split; FetchContent keeps the
# single Library repo submodule-free. Used by the C++ Catalog and YAML
# config parsers. Tests/tools off; built static so it links into
# regpoly_core with no extra runtime dependency. Needs network at
# configure time (CI/docker have it); set FETCHCONTENT_SOURCE_DIR_YAML-CPP
# to a local checkout for offline/air-gapped builds.
include(FetchContent)
set(YAML_CPP_BUILD_TESTS  OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS  OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
set(YAML_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG 0.8.0
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL)
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
     src/search/*.cpp
     src/yaml_config/*.cpp)

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
    src/include/search
    src/include/yaml_config)

# Public-API include directory for downstream consumers. The source
# tree ships `include/regpoly/regpoly.h` directly; the per-class
# headers (bitvect.h, generator.h, …) are mirrored from src/include/
# at config time so the BUILD_INTERFACE include path matches the
# INSTALL_INTERFACE layout.
set(REGPOLY_PUBLIC_INCLUDE_DIR include)
set(REGPOLY_PUBLIC_HEADER_BASES
    src/include/algebra/bitvect.h
    src/include/algebra/primitivity.h
    src/include/core/combo_enumerator.h
    src/include/core/f2_linear_source.h
    src/include/core/generator.h
    src/include/core/i_testable.h
    src/include/core/tempering_chain.h
    src/include/core/transformation.h
    src/include/core/params.h
    src/include/algebra/gauss.h
    src/include/analyses/resolution_sets.h
    src/include/analyses/equidistribution_runner.h
    src/include/analyses/tuplets_runner.h
    src/include/generators/combined_f2_linear_source.h
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
    src/include/yaml_config/seek_config.h
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

# CPU tuning (dockerize-plan Phase 3). Default ON preserves legacy
# `-march=native` for dev builds. The Docker build passes
# -DREGPOLY_NATIVE_ARCH=OFF -DREGPOLY_TUNE_ARCH=neoverse-n1 to produce a
# binary tuned for the Hetzner CAX21's Ampere Altra CPU. Without
# REGPOLY_TUNE_ARCH (and with REGPOLY_NATIVE_ARCH=OFF), only -O3 is set
# — portable but slower.
option(REGPOLY_NATIVE_ARCH "Compile with -march=native (dev builds)" ON)
set(REGPOLY_TUNE_ARCH "" CACHE STRING
    "ARM -mcpu / x86 -march override when REGPOLY_NATIVE_ARCH=OFF")
if(REGPOLY_NATIVE_ARCH)
    set(_REGPOLY_ARCH_FLAGS -O3 -march=native)
elseif(REGPOLY_TUNE_ARCH)
    set(_REGPOLY_ARCH_FLAGS -O3 -mcpu=${REGPOLY_TUNE_ARCH})
else()
    set(_REGPOLY_ARCH_FLAGS -O3)
endif()
target_compile_options(regpoly_core PRIVATE ${_REGPOLY_ARCH_FLAGS})
add_library(regpoly::regpoly ALIAS regpoly_core)

# ── Python extension: _regpoly_cpp ─────────────────────────────────────
# Adds factory.cpp (pybind11 binding registration) and bindings.cpp on
# top of regpoly_core. Gated by REGPOLY_BUILD_PYTHON_EXTENSION so a
# pure-C++ build can skip pybind11 entirely.

if(REGPOLY_BUILD_PYTHON_EXTENSION)
    pybind11_add_module(_regpoly_cpp
        src/bindings/bindings.cpp
        src/bindings/factory_bindings.cpp)
    target_include_directories(_regpoly_cpp PRIVATE ${REGPOLY_CPP_INCLUDE_DIRS})
    target_link_libraries(_regpoly_cpp PRIVATE regpoly_core)
    target_compile_options(_regpoly_cpp PRIVATE ${_REGPOLY_ARCH_FLAGS})

    install(TARGETS _regpoly_cpp DESTINATION regpoly_cpp)
endif()

# ── regpoly-cli executable (Phase 4 fills in subcommands) ───────────────

if(REGPOLY_BUILD_CLI)
    add_executable(regpoly_cli src/cli/main.cpp)
    set_target_properties(regpoly_cli PROPERTIES OUTPUT_NAME regpoly-cli)
    target_include_directories(regpoly_cli PRIVATE ${REGPOLY_CPP_INCLUDE_DIRS})
    # Bake the in-source catalog path in so a binary run from its build
    # tree finds the bundled catalog with no --library flag (the resolver
    # also honours $REGPOLY_CATALOG_DIR). The catalog now lives inside the
    # Python package data so the same files serve the wheel and the CLI.
    target_compile_definitions(regpoly_cli PRIVATE
        REGPOLY_SOURCE_CATALOG="${CMAKE_CURRENT_SOURCE_DIR}/src/regpoly/_data/library")
    # Phase 4.3: cmd_show parses tested-generator YAMLs directly via
    # yaml-cpp (the structure is too free-form to push into a typed
    # loader). Link yaml-cpp into the CLI.
    target_link_libraries(regpoly_cli PRIVATE regpoly_core yaml-cpp::yaml-cpp)
    install(TARGETS regpoly_cli RUNTIME DESTINATION bin)
endif()

# ── Catalog data for installed C++ consumers ───────────────────────────
# The catalog (paper YAMLs) + reference PDFs ship as Python package data
# under src/regpoly/_data/{library,papers}. For a pure-C++ install we copy
# them so the pdf: field ("papers/<file>") still resolves as a sibling of
# the catalog dir: library -> share/regpoly/catalog, papers ->
# share/regpoly/papers. Installed CLI users point --library / $REGPOLY_CATALOG_DIR
# at share/regpoly/catalog.
install(DIRECTORY src/regpoly/_data/library/
        DESTINATION share/regpoly/catalog)
install(DIRECTORY src/regpoly/_data/papers/
        DESTINATION share/regpoly/papers)

# ── Public-API headers ──────────────────────────────────────────────────
# Installed under ${CMAKE_INSTALL_PREFIX}/include/regpoly/.
# Single source of truth: REGPOLY_PUBLIC_HEADER_BASES (defined above)
# plus the umbrella header. The configure_file loop above stages copies
# under ${BINARY_DIR}/include/regpoly/; this install rule pulls them
# from their source paths, so we reuse the same list.

set(REGPOLY_PUBLIC_HEADERS
    include/regpoly/regpoly.h
    ${REGPOLY_PUBLIC_HEADER_BASES}
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
    add_subdirectory(cpp-tests)
endif()
