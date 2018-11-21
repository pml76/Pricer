//
// Created by peter on 11/20/18.
//

#include <catch2/catch.hpp>

#include <math.h>
#include <src/math/erf_wrapper.h>
#include <string>
#include <iostream>


TEST_CASE("ieee754_erf(x) == erf(x)", "[erf]") {

    for (double x = 0.; x <= 10.; x += 0.01) {
        REQUIRE(fabs(erf(x) - ieee754_erf(x)) < 1.0e-12);
    }
}

TEST_CASE("ieee754_erfc(x) == erfc(x)", "[erf]") {

    for (double x = 0.; x <= 10.; x += 0.01) {
        REQUIRE(fabs(erfc(x) - ieee754_erfc(x)) < 1.0e-12);
    }
}