cmake_minimum_required(VERSION 3.12)
project(Pricer)


# create the dipatcher generator
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/dispatcher-generator/)
execute_process(
        COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR}/dispatcher-generator/ -G "Unix Makefiles"
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/dispatcher-generator/)

execute_process(
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}/dispatcher-generator/
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/dispatcher-generator/)

find_file(DISPATCHER_GENERATOR dispgen ${CMAKE_BINARY_DIR}/dispatcher-generator/)

include(cmake/Configure.cmake)


OPTION(OPTION_BUILD_TESTS "Build the test-suite" ON)
OPTION(OPTION_BUILD_BENCHMARKS "Build the benchmark-suite" ON)

# prep cmake to accept our local .cmake files
set(CMAKE_PREFIX_PATH "S{CMAKE_PREFIX_PATH} S{CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})


include(cmake/CCache.cmake)


#
# builds libpricer
set(TARGET_LIBPRICER "pricer")
set(PRICER_VERSION 0)
set(PRICER_SOVERSION 1)

set(TARGET_LIBPRICER_MKL "pricer_mkl")

# builds the test-executable
set(TARGET_PRICER_TEST "pricer-test")

# builds the benchmark executable
set(TARGET_PRICER_BENCHMARK "pricer-benchmark")


# add the third party projects ....
add_subdirectory("3rdParty" "3rdParty")


set(CMAKE_CXX_STANDARD 17)

# this is for the package-manager spack:
# enable @rpath in the install name for any shared library being built
# note: it is planned that a future version of CMake will enable this by default
set(CMAKE_MACOSX_RPATH 1)

# Always use full RPATH
# http://www.cmake.org/Wiki/CMake_RPATH_handling
# http://www.kitware.com/blog/home/post/510

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF ("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF ("${isSystemDir}" STREQUAL "-1")


add_subdirectory(src)


if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
    if (${OPTION_BUILD_TESTS})
        add_subdirectory(tests)
    endif ()

    if (${OPTION_BUILD_BENCHMARKS})
        add_subdirectory(benchmarks)
    endif ()
endif()