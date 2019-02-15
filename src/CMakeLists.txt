

find_package(MKL REQUIRED)


set(TARGET_PRICER_SOURCES
        math/pricers/mkl_pricer.c
        math/pricers/sleef_pricer.c)

set(TARGET_PRICER_HEADERS
        math/pricers/mkl_pricer.h
        math/pricers/sleef_pricer.h)


add_compile_definitions(ENABLE_AVX=1)
add_compile_definitions(DORENAME=1)
include_directories(
        ./
        $ENV{SPACK_TRANSITIVE_INCLUDE_PATH}
        ${MKL_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/3rdParty/sleef/src/arch
        ${CMAKE_SOURCE_DIR}/3rdParty/sleef/src/common
        ${CMAKE_SOURCE_DIR}/3rdParty/sleef/src/libm
        ${CMAKE_BINARY_DIR}/3rdParty/sleef/src/libm/include)

add_library(${TARGET_LIBPRICER} ${TARGET_PRICER_SOURCES} ${TARGET_PRICER_HEADERS})
set_target_properties(${TARGET_LIBPRICER} PROPERTIES
        VERSION ${PRICER_VERSION}
        SOVERSION ${PRICER_SOVERSION}
        ${COMMON_TARGET_PROPERTIES}
        )
target_link_libraries(${TARGET_LIBPRICER} sleef)

