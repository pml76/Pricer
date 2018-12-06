//
// Created by peter on 11/20/18.
//

#include <catch2/catch.hpp>

#include <math.h>
#include <src/math/glibc-2.28_log.h>
#include <string>
#include <iostream>


TEST_CASE("log(x) == ieee754_log(x)", "[log]") {

    for (double x = 0.01; x <= 10.; x += 0.01) {
        REQUIRE(fabs(log(x) - ieee754_log(x)) < 1.0e-12);
    }
}

TEST_CASE("ieee754_log(x*x) == 2*ieee754_log(x)", "[log]") {

    for (double x = 0.01; x <= 10.; x += 0.01) {
        double y = fabs(ieee754_log(x * x) - 2 * ieee754_log(x));
        REQUIRE(y < 1.0e-6);
    }
}


TEST_CASE("ieee754_log(exp(x)) == x", "[log]") {

    for (double x = -10.; x <= 10.; x += 0.01) {
        double y = fabs(ieee754_log(exp(x)) - x);
        REQUIRE(y < 1.0e-6);
    }
}

TEST_CASE("exp(ieee754_log(x) == x", "[log]") {

    for (double x = 0.01; x <= 10.; x += 0.01) {
        double y = fabs(exp(ieee754_log(x)) - x);
        REQUIRE(y < 1.0e-6);
    }
}
