//
// Created by peter on 11/20/18.
//

#include <catch2/catch.hpp>

#include <math.h>
#include <src/math/glibc-2.28_sqrt.h>
#include <string>
#include <iostream>


TEST_CASE("glibc_sqrt(x) == sqrt(x)", "[sqrt]") {

    for (double x = 0.; x <= 10.; x += 0.01) {
        REQUIRE(fabs(sqrt(x) - glibc_sqrt(x)) < 1.0e-12);
    }
}


TEST_CASE("glibc_sqrt(x)*ieee754_sqrt(x) == x", "[sqrt]") {

    for (double x = 0.; x <= 10.; x += 0.01) {
        REQUIRE(fabs(glibc_sqrt(x) * glibc_sqrt(x) - x) < 1.0e-12);
    }
}

TEST_CASE("glibc_sqrt(x*x) == x", "[sqrt]") {

    for (double x = 0.; x <= 10.; x += 0.01) {
        REQUIRE(fabs(glibc_sqrt(x * x) - x) < 1.0e-12);
    }
}
