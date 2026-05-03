# Using REGPOLY from C++

> **Status**: Phase 0 placeholder. Final content authored in Phase 4–6 once `regpoly-cli` and the public C++ headers are in place.

Pure-C++ usage requires Phase 1+ to land. After that, downstream C++ projects will use `find_package(regpoly)`:

```cmake
find_package(regpoly REQUIRED)
target_link_libraries(my_project PRIVATE regpoly::regpoly)
```

A standalone CLI (`regpoly-cli`) will accept the same YAML configs as the Python `regpoly` CLI.
