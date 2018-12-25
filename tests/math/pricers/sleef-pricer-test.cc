//
// Created by peter on 12/25/18.
//

#include <catch2/catch.hpp>
#include <iostream>
#include <src/math/pricers/sleef_pricer.h>
#include <tests/math/pricers/Pricer.h>


#define DECLARE_AND_DEFINE(type, x, y) \
   type x[64] __attribute__((aligned(ALIGN_TO))); \
   for(UINT64 i=0;i<ALIGN_TO;i++) x[i]=y;

TEST_CASE("pricer-class equals sleef-pricer (long call)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, put_call, 1)
    DECLARE_AND_DEFINE(FLOAT, long_short, 1)


    Pricer p;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();

    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);



    // test long call price
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    FLOAT reference_pricer_value = p.compute_call_price(x[0]);
    REQUIRE(abs(reference_pricer_value - price[0]) < 1.0e-8);

}

TEST_CASE("pricer-class equals sleef-pricer (short call)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, put_call, 1)
    DECLARE_AND_DEFINE(FLOAT, long_short, -1)


    Pricer p;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();

    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);



    // test long call price
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    FLOAT reference_pricer_value = p.compute_call_price(x[0]);
    REQUIRE(abs(reference_pricer_value + price[0]) < 1.0e-8);

}

TEST_CASE("pricer-class equals sleef-pricer (long put)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, put_call, -1)
    DECLARE_AND_DEFINE(FLOAT, long_short, 1)


    Pricer p;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();

    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);



    // test long call price
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    FLOAT reference_pricer_value = p.compute_put_price(x[0]);
    REQUIRE(abs(reference_pricer_value - price[0]) < 1.0e-8);

}

TEST_CASE("pricer-class equals sleef-pricer (short put)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, put_call, -1)
    DECLARE_AND_DEFINE(FLOAT, long_short, -1)


    Pricer p;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();

    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);



    // test long call price
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    FLOAT reference_pricer_value = p.compute_put_price(x[0]);

    REQUIRE(abs(reference_pricer_value + price[0]) < 1.0e-8);

}