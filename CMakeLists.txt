cmake_minimum_required(VERSION 3.12)
project(Pricer)


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


set(CMAKE_PREFIX_PATH "S{CMAKE_PREFIX_PATH} S{CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -fopenmp -Ofast -march=native -fno-omit-frame-pointer -ffp-contract=off -fno-math-errno -fno-trapping-math -O3 -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fopenmp -Ofast -march=native -fno-omit-frame-pointer -ffp-contract=off -fno-math-errno -fno-trapping-math -O3 -fPIC")


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

set(PROJECT_SOURCES src/math/pricers/mkl_pricer.c
        3rdParty/sleef/src/libm/sleefsimddp.c
        3rdParty/sleef/src/libm/rempitab.c
        src/math/pricers/sleef_pricer.c)
set(PROJECT_HEADERS src/math/pricers/mkl_pricer.h
        3rdParty/sleef/src/arch/helperadvsimd.h
        3rdParty/sleef/src/arch/helperavx.h
        3rdParty/sleef/src/arch/helperavx2.h
        3rdParty/sleef/src/arch/helperavx2_128.h
        3rdParty/sleef/src/arch/helperavx512f.h
        3rdParty/sleef/src/arch/helperneon32.h
        3rdParty/sleef/src/arch/helperpower_128.h
        3rdParty/sleef/src/arch/helperpurec.h
        3rdParty/sleef/src/arch/helpersse2.h
        3rdParty/sleef/src/arch/helpersve.h
        3rdParty/sleef/src/arch/helpervecext.h
        3rdParty/sleef/src/common/misc.h
        3rdParty/sleef/src/libm/dd.h
        src/math/pricers/sleef_pricer.h)

add_compile_definitions(ENABLE_AVX=1)

set(PROJECT_TEST_SOURCES tests/math/pricers/Pricer.h tests/math/pricers/mkl-pricer-test.cc tests/math/pricers/sleef-pricer-test.cc)

set(PROJECT_BENCHMARK_SOURCES benchmarks/main.cc benchmarks/math/mkl-pricer-benchmark.cc benchmarks/math/sleef-pricer-benchmark.cc)

include_directories(./ $ENV{SPACK_TRANSITIVE_INCLUDE_PATH} ${MKL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/3rdParty/sleef/src/arch ${CMAKE_SOURCE_DIR}/3rdParty/sleef/src/common ${CMAKE_SOURCE_DIR}/3rdParty/sleef/src/libm)

add_executable(Pricer ${PROJECT_SOURCES} ${PROJECT_HEADERS} src/main.cpp)
target_link_libraries(Pricer pthread dl ${MKL_LIBRARIES})

add_executable(Pricer-Tests tests/main.cc ${PROJECT_TEST_SOURCES} ${PROJECT_SOURCES} ${PROJECT_HEADERS})
target_link_libraries(Pricer-Tests pthread dl ${MKL_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${FLINT_LIBRARIES} ${Arb_LIBRARIES})

add_executable(Pricer-Benchmark benchmarks/main.cc ${PROJECT_BENCHMARK_SOURCES} ${PROJECT_SOURCES} ${PROJECT_HEADERS})
target_link_libraries(Pricer-Benchmark pthread dl ${SLEEF_LIBRARIES} ${MKL_LIBRARIES} ${benchmark_LIBRARIES})

if ("${CMAKE_SYSTEM_NAME}" MATCHES "Windows")
    target_link_libraries(Pricer-Benchmark Shlwapi)
endif ()