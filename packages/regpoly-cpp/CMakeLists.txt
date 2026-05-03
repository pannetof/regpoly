cmake_minimum_required(VERSION 3.20)
project(regpoly_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(pybind11 CONFIG REQUIRED)

file(GLOB_RECURSE REGPOLY_CPP_SOURCES
     CONFIGURE_DEPENDS
     src/algebra/*.cpp
     src/core/*.cpp
     src/generators/*.cpp
     src/transforms/*.cpp
     src/lattice/*.cpp)

set(REGPOLY_CPP_INCLUDE_DIRS
    src/include/core
    src/include/algebra
    src/include/generators
    src/include/transforms
    src/include/lattice)

pybind11_add_module(_regpoly_cpp
    src/bindings/bindings.cpp
    ${REGPOLY_CPP_SOURCES})

target_include_directories(_regpoly_cpp PRIVATE ${REGPOLY_CPP_INCLUDE_DIRS})
target_link_libraries(_regpoly_cpp PRIVATE ntl pthread)
target_compile_options(_regpoly_cpp PRIVATE -O3 -march=native)

install(TARGETS _regpoly_cpp DESTINATION regpoly_cpp)
