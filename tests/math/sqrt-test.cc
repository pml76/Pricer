//
// Created by peter on 11/20/18.
//

#include <catch2/catch.hpp>

#include <math.h>
#include <src/math/sqrt_wrapper.h>
#include <string>
#include <iostream>


TEST_CASE("ieee754_sqrt(x) == sqrt(x)", "[sqrt]") {

    for (double x = 0.; x <= 10.; x += 0.01) {
        REQUIRE(fabs(sqrt(x) - ieee754_sqrt(x)) < 1.0e-12);
    }
}


TEST_CASE("ieee754_sqrt(x)*ieee754_sqrt(x) == x", "[sqrt]") {

    for (double x = 0.; x <= 10.; x += 0.01) {
        REQUIRE(fabs(ieee754_sqrt(x) * ieee754_sqrt(x) - x) < 1.0e-12);
    }
}

TEST_CASE("ieee754_sqrt(x*x) == x", "[sqrt]") {

    for (double x = 0.; x <= 10.; x += 0.01) {
        REQUIRE(fabs(ieee754_sqrt(x * x) - x) < 1.0e-12);
    }
}
