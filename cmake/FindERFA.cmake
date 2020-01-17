if (ERFA_INCLUDE_DIR AND Arb_LIBRARIES)
    # Already in cache, be silent
    set(ERFA_FIND_QUIETLY TRUE)
endif ()

find_path(ERFA_INCLUDE_DIR erfa.h)
find_library(ERFA_LIBRARIES NAMES erfa)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(ERFA DEFAULT_MSG ERFA_INCLUDE_DIR ERFA_LIBRARIES)

mark_as_advanced(ERFA_INCLUDE_DIR ERFA_LIBRARIES)