cmake_minimum_required(VERSION 3.20)
project(regpoly_example LANGUAGES CXX)

# Standalone consumer of an *installed* regpoly. Build it against an install
# prefix with:
#   cmake -S cpp/examples -B build-ex -DCMAKE_PREFIX_PATH=/path/to/install
# This is the acceptance gate: it proves the public header compiles, the
# imported target links, and NTL + yaml-cpp resolve via find_dependency —
# with no Python anywhere.
find_package(regpoly REQUIRED)

add_executable(regpoly_example main.cpp)
target_link_libraries(regpoly_example PRIVATE regpoly::regpoly)
