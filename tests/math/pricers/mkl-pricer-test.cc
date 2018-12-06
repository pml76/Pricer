//
// Created by peter on 12/5/18.
//

#include <catch2/catch.hpp>
#include <iostream>
#include <src/math/pricers/mkl_pricer.h>
#include <tests/math/pricers/Pricer.h>




TEST_CASE("pricer-class equals mkl-pricer", "[pricer]") {

    double r     = 0.01;    // interest rate
    double s     = 70.;     // stock price
    double t     = 1.2;     // time to maturity
    double tau   = 1./12.;  // length of averaging time-window
    double sigma = 0.3;     // vola
    double x     = 72;      // strike

    double tmp1, tmp2, tmp3, tmp4, tmp5, sigmaA, sigmaA2T2,sigmaAsqrtT,emrt;
    double d1,d2;
    double price;
    MKL_INT64 flags = 0;

    Pricer p;


    p.set_market_data(sigma,t,tau,r,s);

    init_mkl_pricer();
    prepare_mkl_pricer(1,&s,&sigma,&t,&tau,&r,&tmp1,&tmp2,&tmp3,&tmp4,&tmp5,&sigmaA,&sigmaA2T2,&sigmaAsqrtT,&emrt);


    // test long call price
    mkl_pricer(1,&flags,&s,&x,&sigmaA2T2,&sigmaAsqrtT,&emrt,&tmp1,&tmp2,&tmp3,&tmp4,&d1,&d2,&price);
    double reference_pricer_value = p.compute_call_price(x);
    REQUIRE(abs(reference_pricer_value - price) < 1.0e-8 );

    // test short call price
    flags = 2;
    mkl_pricer(1,&flags,&s,&x,&sigmaA2T2,&sigmaAsqrtT,&emrt,&tmp1,&tmp2,&tmp3,&tmp4,&d1,&d2,&price);
    reference_pricer_value = p.compute_call_price(x);
    REQUIRE(abs(reference_pricer_value + price) < 1.0e-8 );

    // test long put price
    flags = 1;
    mkl_pricer(1,&flags,&s,&x,&sigmaA2T2,&sigmaAsqrtT,&emrt,&tmp1,&tmp2,&tmp3,&tmp4,&d1,&d2,&price);
    reference_pricer_value = p.compute_put_price(x);
    REQUIRE(abs(reference_pricer_value - price) < 1.0e-8 );

    // test short put price
    flags = 3;
    mkl_pricer(1,&flags,&s,&x,&sigmaA2T2,&sigmaAsqrtT,&emrt,&tmp1,&tmp2,&tmp3,&tmp4,&d1,&d2,&price);
    reference_pricer_value = p.compute_put_price(x);
    REQUIRE(abs(reference_pricer_value + price) < 1.0e-8 );

}

TEST_CASE("pricer-class equals ddx-mkl-pricer", "[pricer]") {

    double r     = 0.00;    // interest rate
    double s     = 70.;     // stock price
    double t     = 1.2;     // time to maturity
    double tau   = 1./12.;  // length of averaging time-window
    double sigma = 0.3;     // vola
    double x     = 72;      // strike

    double tmp1, tmp2, tmp3, tmp4, tmp5, sigmaA, sigmaA2T2,sigmaAsqrtT,emrt;
    double d1,d2;
    double price, ddx_price;
    const double eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma,t,tau,r,s);

    init_mkl_pricer();
    prepare_mkl_pricer(1,&s,&sigma,&t,&tau,&r,&tmp1,&tmp2,&tmp3,&tmp4,&tmp5,&sigmaA,&sigmaA2T2,&sigmaAsqrtT,&emrt);


    // test long call price
    MKL_INT64 flags = 0;
    mkl_pricer(1,&flags,&s,&x,&sigmaA2T2,&sigmaAsqrtT,&emrt,&tmp1,&tmp2,&tmp3,&tmp4,&d1,&d2,&price);
    ddx_mkl_pricer(1,&flags,&d2,&emrt,&ddx_price);
    double reference_pricer_value1 = p.compute_call_price(x);
    double reference_pricer_value2 = p.compute_call_price(x+eps);
    /*
    std::cout << "d2 = " << d2 << std::endl;
    std::cout << "emrt = " << emrt << std::endl;
    std::cout << "ddx_price = " << ddx_price << std::endl;
    std::cout << "(reference_pricer_value2 - reference_pricer_value1)/eps  = "
              << (reference_pricer_value2 - reference_pricer_value1)/eps << std::endl;
    std::cout.flush();
    */
    REQUIRE( abs((reference_pricer_value2 - reference_pricer_value1)/eps -ddx_price )< 1.0e-4 );


    // test short call price
    flags = 2;
    mkl_pricer(1,&flags,&s,&x,&sigmaA2T2,&sigmaAsqrtT,&emrt,&tmp1,&tmp2,&tmp3,&tmp4,&d1,&d2,&price);
    ddx_mkl_pricer(1,&flags,&d2,&emrt,&ddx_price);
    reference_pricer_value1 = p.compute_call_price(x);
    reference_pricer_value2 = p.compute_call_price(x+eps);
    /*
    std::cout << "d2 = " << d2 << std::endl;
    std::cout << "emrt = " << emrt << std::endl;
    std::cout << "ddx_price = " << ddx_price << std::endl;
    std::cout << "(reference_pricer_value2 - reference_pricer_value1)/eps  = "
              << (reference_pricer_value2 - reference_pricer_value1)/eps << std::endl;
    std::cout.flush();
    */
    REQUIRE( abs((reference_pricer_value2 - reference_pricer_value1)/eps  + ddx_price )< 1.0e-4 );


    // test long put price
    flags = 1;
    mkl_pricer(1,&flags,&s,&x,&sigmaA2T2,&sigmaAsqrtT,&emrt,&tmp1,&tmp2,&tmp3,&tmp4,&d1,&d2,&price);
    ddx_mkl_pricer(1,&flags,&d2,&emrt,&ddx_price);
    reference_pricer_value1 = p.compute_put_price(x);
    reference_pricer_value2 = p.compute_put_price(x+eps);

    /*
    std::cout << "d2 = " << d2 << std::endl;
    std::cout << "emrt = " << emrt << std::endl;
    std::cout << "ddx_price = " << ddx_price << std::endl;
    std::cout << "(reference_pricer_value2 - reference_pricer_value1)/eps  = "
              << (reference_pricer_value2 - reference_pricer_value1)/eps << std::endl;
    std::cout.flush();
    */
    REQUIRE( abs((reference_pricer_value2 - reference_pricer_value1)/eps  - ddx_price )< 1.0e-4 );

    // test short put price
    flags = 3;
    mkl_pricer(1,&flags,&s,&x,&sigmaA2T2,&sigmaAsqrtT,&emrt,&tmp1,&tmp2,&tmp3,&tmp4,&d1,&d2,&price);
    ddx_mkl_pricer(1,&flags,&d2,&emrt,&ddx_price);
    reference_pricer_value1 = p.compute_put_price(x);
    reference_pricer_value2 = p.compute_put_price(x+eps);
    /*
    std::cout << "d2 = " << d2 << std::endl;
    std::cout << "emrt = " << emrt << std::endl;
    std::cout << "ddx_price = " << ddx_price << std::endl;
    std::cout << "(reference_pricer_value2 - reference_pricer_value1)/eps  = "
              << (reference_pricer_value2 - reference_pricer_value1)/eps << std::endl;
    std::cout.flush();
    */
    REQUIRE( abs((reference_pricer_value2 - reference_pricer_value1)/eps  + ddx_price )< 1.0e-4 );



}