

find_package(GoogleBenchmark REQUIRED)
find_package(ARB REQUIRED)
find_package(FLINT REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(MKL REQUIRED)

set(TARGET_PRICER_BENCHMARKS_SOURCES
        math/mkl-pricer-benchmark.cc
        math/pricers/sleef-pricer-benchmark.cc
        main.cc)


add_executable(${TARGET_PRICER_BENCHMARK} ${TARGET_PRICER_SOURCES})
target_link_libraries(${TARGET_PRICER_BENCHMARK}
        pthread
        dl
        ${MKL_LIBRARIES}
        ${GMP_LIBRARIES}
        ${MPFR_LIBRARIES}
        ${FLINT_LIBRARIES}
        ${Arb_LIBRARIES}
        sleef
        ${TARGET_LIBPRICER})

if ("${CMAKE_SYSTEM_NAME}" MATCHES "Windows")
    target_link_libraries(Pricer-Benchmark Shlwapi)
endif ()