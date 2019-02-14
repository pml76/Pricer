

find_package(GoogleBenchmark REQUIRED)
find_package(Catch2 REQUIRED)
find_package(ARB REQUIRED)
find_package(FLINT REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(MKL REQUIRED)

set(TARGET_PRICER_TESTS_SOURCES
        math/pricers/mkl-pricer-test.cc
        math/pricers/sleef-pricer-test.cc)

set(TARGET_PRICER_Tests_HEADERS
        math/pricers/Pricer.h)

add_executable(${TARGET_PRICER_TEST} ${TARGET_PRICER_SOURCES} ${TARGET_PRICER_HEADERS})
target_link_libraries(${TARGET_PRICER_TEST}
        pthread
        dl
        ${MKL_LIBRARIES}
        ${GMP_LIBRARIES}
        ${MPFR_LIBRARIES}
        ${FLINT_LIBRARIES}
        ${Arb_LIBRARIES}
        sleef
        ${TARGET_LIBPRICER})