//
// Created by peter on 11/20/18.
//

#include <catch2/catch.hpp>

#include <math.h>
#include <src/math/exp_wrapper.h>
#include <string>
#include <iostream>


TEST_CASE("exp(x) == ieee754_exp(x)", "[exp]") {

    for (double x = -10.; x <= 10.; x += 0.01) {
        REQUIRE(fabs(exp(x) - ieee754_exp(x)) < 1.0e-12);
    }
}


TEST_CASE("ieee754_exp(x)*ieee754_exp(x) == ieee754_exp(2*x)", "[exp]") {

    for (double x = -10.; x <= 10.; x += 0.01) {
        double y = fabs(ieee754_exp(x) * ieee754_exp(x) - ieee754_exp(2 * x));
        REQUIRE(y <= 1.0e-6);
    }
}


TEST_CASE("log(ieee_exp(x)) == x", "[exp]") {

    for (double x = -10.; x <= 10.; x += 0.01) {
        double y = fabs(log(ieee754_exp(x)) - x);
        REQUIRE(y <= 1.0e-6);
    }
}


TEST_CASE("ieee_exp(log(x)) == x", "[exp]") {

    for (double x = 0.01; x <= 10.; x += 0.01) {
        double y = fabs(ieee754_exp(log(x)) - x);
        REQUIRE(y <= 1.0e-6);
    }
}
