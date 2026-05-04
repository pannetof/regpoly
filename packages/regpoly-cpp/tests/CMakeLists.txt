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
    test_seek_config.cpp
    test_catalog.cpp
    test_legacy_reader.cpp
    test_tempering_optimizer.cpp
    test_tempering_search.cpp
    test_yaml_cpp_smoke.cpp
)
target_include_directories(regpoly_cpp_tests PRIVATE
    ../src/include/algebra
    ../src/include/core
    ../src/include/generators
    ../src/include/transforms
    ../src/include/lattice
    ../src/include/library
    ../src/include/io
    ../src/include/search
    ../src/include/yaml_config)
target_link_libraries(regpoly_cpp_tests PRIVATE
    regpoly_core
    yaml-cpp::yaml-cpp
    GTest::gtest_main)

gtest_discover_tests(regpoly_cpp_tests)

# ── regpoly-cli subprocess tests ──────────────────────────────────────
# Exercise the Phase 4.1 CLI subcommands by spawning the binary. These
# don't go through GoogleTest — they just check exit code + stdout
# pattern via `add_test` + `set_tests_properties(... PASS_REGULAR_EXPRESSION ...)`.
# The binary lives under the build dir; the docs/library and
# shared/legacy_parameters fixtures live in the source tree.

if(TARGET regpoly_cli)
    set(REGPOLY_CLI $<TARGET_FILE:regpoly_cli>)
    set(REPO_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../..)

    add_test(NAME RegpolyCli.Help
             COMMAND ${REGPOLY_CLI} --help)
    set_tests_properties(RegpolyCli.Help PROPERTIES
        PASS_REGULAR_EXPRESSION "Usage: regpoly-cli")

    add_test(NAME RegpolyCli.Version
             COMMAND ${REGPOLY_CLI} --version)
    set_tests_properties(RegpolyCli.Version PROPERTIES
        PASS_REGULAR_EXPRESSION "regpoly-cli 2\\.0\\.0")

    add_test(NAME RegpolyCli.CatalogList
             COMMAND ${REGPOLY_CLI} catalog list
                     --library ${REPO_ROOT}/docs/library)
    set_tests_properties(RegpolyCli.CatalogList PROPERTIES
        PASS_REGULAR_EXPRESSION "matsumoto-nishimura-1998")

    add_test(NAME RegpolyCli.CatalogShowMT
             COMMAND ${REGPOLY_CLI} catalog show matsumoto-nishimura-1998
                     --library ${REPO_ROOT}/docs/library)
    set_tests_properties(RegpolyCli.CatalogShowMT PROPERTIES
        PASS_REGULAR_EXPRESSION "Mersenne Twister")

    add_test(NAME RegpolyCli.CatalogGenMT
             COMMAND ${REGPOLY_CLI} catalog gen mt19937
                     --library ${REPO_ROOT}/docs/library)
    set_tests_properties(RegpolyCli.CatalogGenMT PROPERTIES
        PASS_REGULAR_EXPRESSION "MersenneTwister")

    add_test(NAME RegpolyCli.LegacyInfoMT
             COMMAND ${REGPOLY_CLI} legacy-info
                     ${REPO_ROOT}/shared/legacy_parameters/MT19937.dat)
    set_tests_properties(RegpolyCli.LegacyInfoMT PROPERTIES
        PASS_REGULAR_EXPRESSION "MersenneTwister")

    add_test(NAME RegpolyCli.LegacyTrans
             COMMAND ${REGPOLY_CLI} legacy-trans
                     ${REPO_ROOT}/shared/legacy_parameters/trans32.dat)
    set_tests_properties(RegpolyCli.LegacyTrans PROPERTIES
        PASS_REGULAR_EXPRESSION "tempMK")

    # Phase 4.2: end-to-end search via the YAML config loader.
    add_test(NAME RegpolyCli.SearchCarry32
             COMMAND ${REGPOLY_CLI} search
                     ${REPO_ROOT}/shared/yaml/equidist/carry32.yaml)
    set_tests_properties(RegpolyCli.SearchCarry32 PROPERTIES
        PASS_REGULAR_EXPRESSION "combos")

    # Negative cases — exit non-zero on bad input.
    add_test(NAME RegpolyCli.UnknownCommand
             COMMAND ${REGPOLY_CLI} bogus-command)
    set_tests_properties(RegpolyCli.UnknownCommand PROPERTIES
        WILL_FAIL TRUE)

    add_test(NAME RegpolyCli.MissingPaper
             COMMAND ${REGPOLY_CLI} catalog show no-such-paper
                     --library ${REPO_ROOT}/docs/library)
    set_tests_properties(RegpolyCli.MissingPaper PROPERTIES
        WILL_FAIL TRUE)
endif()
