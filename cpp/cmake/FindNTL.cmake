# FindNTL.cmake — locate the NTL library and its GMP dependency, and expose
# an imported target NTL::NTL (with GMP propagated). NTL ships no CMake
# config of its own, so consumers rely on this module.
#
# Honours the standard hints NTL_ROOT / GMP_ROOT and CMAKE_PREFIX_PATH.

include(FindPackageHandleStandardArgs)

find_path(NTL_INCLUDE_DIR NAMES NTL/ZZ.h)
find_library(NTL_LIBRARY NAMES ntl libntl)

# NTL is built on GMP; locate it so downstream links resolve.
find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIBRARY NAMES gmp libgmp)

find_package_handle_standard_args(NTL
    REQUIRED_VARS NTL_LIBRARY NTL_INCLUDE_DIR)

if(NTL_FOUND AND NOT TARGET NTL::NTL)
    add_library(NTL::NTL UNKNOWN IMPORTED)
    set_target_properties(NTL::NTL PROPERTIES
        IMPORTED_LOCATION "${NTL_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${NTL_INCLUDE_DIR}")
    if(GMP_LIBRARY)
        set_property(TARGET NTL::NTL APPEND PROPERTY
            INTERFACE_LINK_LIBRARIES "${GMP_LIBRARY}")
    endif()
    if(GMP_INCLUDE_DIR)
        set_property(TARGET NTL::NTL APPEND PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")
    endif()
endif()

mark_as_advanced(NTL_INCLUDE_DIR NTL_LIBRARY GMP_INCLUDE_DIR GMP_LIBRARY)
