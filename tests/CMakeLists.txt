

find_package(GoogleBenchmark REQUIRED)
find_package(Catch2 REQUIRED)
find_package(ARB REQUIRED)
find_package(FLINT REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(MKL REQUIRED)
find_package(ERFA REQUIRED)

set(TARGET_PRICER_TEST_SOURCES
        math/pricers/mkl-pricer-test.cc
        math/pricers/sleef-pricer-test.cc
        main.cc IO/io_tests.cc)

set(TARGET_PRICER_TTEST_HEADERS
        math/pricers/Pricer.h)

include_directories(
        ./
        ../src/
        ../src/include/
        ${CMAKE_BINARY_DIR}
        $ENV{SPACK_TRANSITIVE_INCLUDE_PATH}
        ${MKL_INCLUDE_DIRS}
        ${ERFA_INCLUDE_DIR})


add_executable(${TARGET_PRICER_TEST} ${TARGET_PRICER_TEST_SOURCES} ${TARGET_PRICER_TEST_HEADERS})
target_link_libraries(${TARGET_PRICER_TEST}
        pthread
        dl
        ${MKL_LIBRARIES}
        ${GMP_LIBRARIES}
        ${MPFR_LIBRARIES}
        ${FLINT_LIBRARIES}
        ${Arb_LIBRARIES}
        sleef
        ${TARGET_LIBPRICER}
        ${TARGET_ASMLIB}
        ${ERFA_LIBRARIES}
        Catch2::Catch2)