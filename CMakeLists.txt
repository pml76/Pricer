cmake_minimum_required(VERSION 3.12)
project(Pricer)

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


set(CMAKE_PREFIX_PATH "S{CMAKE_PREFIX_PATH} S{CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -march=ivybridge -fno-omit-frame-pointer -fopt-info-vec-missed -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=ivybridge -fno-omit-frame-pointer -fopt-info-vec-missed -fPIC")


set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

find_package(GoogleBenchmark REQUIRED)
find_package(Catch2 REQUIRED)
find_package(MKL REQUIRED)

find_package(ARB REQUIRED)
find_package(FLINT REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

find_package(SLEEF REQUIRED)


set(CMAKE_VERBOSE_MAKEFILE ON)

set(PROJECT_SOURCES src/math/pricers/mkl_pricer.c)
set(PROJECT_HEADERS src/math/pricers/mkl_pricer.c)

set(PROJECT_TEST_SOURCES tests/math/pricers/Pricer.h tests/math/pricers/mkl-pricer-test.cc)

set(PROJECT_BENCHMARK_SOURCES benchmarks/main.cc benchmarks/math/exp-benchmark.cc benchmarks/math/erf-benchmark.cc benchmarks/math/log-benchmark.cc benchmarks/math/sqrt-benchmark.cc benchmarks/math/mkl-pricer-benchmark.cc src/math/include/erfc.h tests/math/Erf-Tests.c)

include_directories(./ $ENV{SPACK_TRANSITIVE_INCLUDE_PATH} ${MKL_INCLUDE_DIRS})

add_executable(Pricer ${PROJECT_SOURCES} ${PROJECT_HEADERS} src/main.cpp)
target_link_libraries(Pricer pthread dl ${MKL_LIBRARIES})

add_executable(Pricer-Tests tests/main.cc ${PROJECT_TEST_SOURCES} ${PROJECT_SOURCES} ${PROJECT_HEADERS})
target_link_libraries(Pricer-Tests pthread dl ${MKL_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${FLINT_LIBRARIES} ${Arb_LIBRARIES})

add_executable(Pricer-Benchmark benchmarks/main.cc ${PROJECT_BENCHMARK_SOURCES} ${PROJECT_SOURCES} ${PROJECT_HEADERS})
target_link_libraries(Pricer-Benchmark pthread dl ${SLEEF_LIBRARIES} ${MKL_LIBRARIES} ${benchmark_LIBRARIES})

if ("${CMAKE_SYSTEM_NAME}" MATCHES "Windows")
    target_link_libraries(Pricer-Benchmark Shlwapi)
endif ()