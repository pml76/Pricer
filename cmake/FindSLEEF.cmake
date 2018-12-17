if (SLEEF_INCLUDE_DIR AND Arb_LIBRARIES)
    # Already in cache, be silent
    set(SLEEF_FIND_QUIETLY TRUE)
endif ()

find_path(SLEEF_INCLUDE_DIR sleef.h)
find_library(SLEEF_LIBRARIES NAMES sleef)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(SLEEF DEFAULT_MSG SLEEF_INCLUDE_DIR SLEEF_LIBRARIES)

mark_as_advanced(SLEEF_INCLUDE_DIR SLEEF_LIBRARIES)