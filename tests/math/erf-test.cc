//
// Created by peter on 11/20/18.
//

#include <catch2/catch.hpp>

#include <math.h>
#include <src/math/erf_wrapper.h>
#include <string>
#include <iostream>

#include <mkl.h>


TEST_CASE("ieee754_erf(x) == erf(x)", "[erf]") {

    for (double x = -10.; x <= 10.; x += 0.01) {
        REQUIRE(fabs(erf(x) - ieee754_erf(x)) < 1.0e-12);
    }
}

TEST_CASE("ieee754_erfc(x) == erfc(x)", "[erf]") {

    for (double x = -10.; x <= 10.; x += 0.01) {
        REQUIRE(fabs(erfc(x) - ieee754_erfc(x)) < 1.0e-12);
    }
}


TEST_CASE("MKL_erf(x) == erf(x)", "[erf]") {


    vmlSetMode(VML_LA);
    vmlSetMode(VML_FTZDAZ_OFF);
    vmlSetMode(VML_ERRMODE_NOERR);

    double y;
    for (double x = -10.; x <= 10.; x += 0.01) {
        vdErf(1, &x, &y);
        REQUIRE(fabs(erf(x) - y) < 1.0e-10);
    }
}

TEST_CASE("MKL_erfc(x) == erfc(x)", "[erf]") {

    vmlSetMode(VML_LA);
    vmlSetMode(VML_FTZDAZ_OFF);
    vmlSetMode(VML_ERRMODE_NOERR);

    double y;
    for (double x = -10.; x <= 10.; x += 0.01) {
        vdErfc(1, &x, &y);
        REQUIRE(fabs(erfc(x) - y) < 1.0e-10);
    }
}