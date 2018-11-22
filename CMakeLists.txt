cmake_minimum_required(VERSION 3.12)
project(Pricer)

set(CMAKE_CXX_STANDARD 17)

include(ExternalProject)

# build google benchmark (target: benchmark)
# do not build tests of benchmarking lib
# set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's tests" FORCE)
# add_subdirectory(googlebenchmark)

# build tests (targets: gtest_main, gtest)
# add_subdirectory(googletest)


set(CMAKE_PREFIX_PATH "S{CMAKE_PREFIX_PATH} S{CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=ivybridge -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=ivybridge -fno-omit-frame-pointer")

set(CMAKE_VERBOSE_MAKEFILE ON)

# Download and unpack googlebenchmark at configure time
configure_file(cmake/googlebenchmark.cmake.in googlebenchmark-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googlebenchmark-download)
if (result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif ()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googlebenchmark-download)
if (result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif ()

find_library(GOOGLE_BENCHMARK_LIB NAMES benchmark PATHS "${CMAKE_BINARY_DIR}/googlebenchmark-install/lib")


# Download and unpack batch2 at configure time
configure_file(cmake/catch2.cmake.in catch2-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/catch2-download)
if (result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif ()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/catch2-download)
if (result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif ()


set(PROJECT_SOURCES src/math/erf_wrapper.c src/math/exp_wrapper.c src/math/log_wrapper.c src/math/sqrt_wrapper.c src/math/erf_wrapper.c)
set(PROJECT_HEADERS src/math/erf_wrapper.h src/math/exp_wrapper.h src/math/log_wrapper.h src/math/sqrt_wrapper.h src/math/erf_wrapper.h)

set(PROJECT_TEST_SOURCES tests/math/exp-test.cc tests/math/log-test.cc tests/math/erf-test.cc tests/math/sqrt-test.cc)

set(PROJECT_BENCHMARK_SOURCES benchmarks/main.cc benchmarks/math/exp-benchmark.cc benchmarks/math/erf-benchmark.cc benchmarks/math/log-benchmark.cc benchmarks/math/sqrt-benchmark.cc)

include_directories(./ ${CMAKE_BINARY_DIR}/catch2-install/include ${CMAKE_BINARY_DIR}/googlebenchmark-install/include)
add_executable(Pricer ${PROJECT_SOURCES} ${PROJECT_HEADERS} src/main.cpp)

add_executable(Pricer-Tests tests/main.cc ${PROJECT_TEST_SOURCES} ${PROJECT_SOURCES} ${PROJECT_HEADERS})


add_executable(Pricer-Benchmark benchmarks/main.cc ${PROJECT_BENCHMARK_SOURCES} ${PROJECT_SOURCES} ${PROJECT_HEADERS})
target_link_libraries(Pricer-Benchmark pthread ${GOOGLE_BENCHMARK_LIB})