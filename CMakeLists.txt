cmake_minimum_required(VERSION 3.12)
project(Pricer)

set(CMAKE_CXX_STANDARD 17)


include_directories(./)

add_executable(Pricer src/main.cpp src/math/erf_wrapper.c src/math/exp_wrapper.c src/math/exp_wrapper.h src/math/log_wrapper.c src/math/log_wrapper.h src/math/sqrt_wrapper.c src/math/erf_wrapper.c src/math/erf_wrapper.h src/math/sqrt_wrapper.h)

add_executable(Pricer-Tests tests/test-main.cc tests/math/exp-test.cc src/math/exp_wrapper.h src/math/exp_wrapper.c src/math/log_wrapper.c src/math/log_wrapper.h tests/math/log-test.cc src/math/sqrt_wrapper.c src/math/erf_wrapper.c tests/math/erf-test.cc tests/math/sqrt-test.cc src/math/erf_wrapper.h src/math/sqrt_wrapper.h)

target_link_libraries(Pricer m)