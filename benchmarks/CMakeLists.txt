
find_package(GoogleBenchmark REQUIRED)
find_package(ARB REQUIRED)
find_package(FLINT REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(MKL REQUIRED)

find_package(CLI11 CONFIG REQUIRED)

set(TARGET_PRICER_BENCHMARK_SOURCES
        math/mkl-pricer-benchmark.cc
        math/sleef-pricer-benchmark.cc
        main.cc)

include_directories(
        ./
        ../src/include/
        ../src/
        ${CMAKE_BINARY_DIR}
        $ENV{SPACK_TRANSITIVE_INCLUDE_PATH}
        ${MKL_INCLUDE_DIRS})


add_executable(${TARGET_PRICER_BENCHMARK} ${TARGET_PRICER_BENCHMARK_SOURCES})
target_link_libraries(${TARGET_PRICER_BENCHMARK}
        pthread
        dl
        ${MKL_LIBRARIES}
        ${GMP_LIBRARIES}
        ${MPFR_LIBRARIES}
        ${FLINT_LIBRARIES}
        ${Arb_LIBRARIES}
        sleef
        ${TARGET_LIBPRICER}
        ${benchmark_LIBRARIES})

if ("${CMAKE_SYSTEM_NAME}" MATCHES "Windows")
    target_link_libraries(Pricer-Benchmark Shlwapi)
endif ()