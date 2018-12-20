//
// Created by peter on 12/5/18.
//

#include <catch2/catch.hpp>
#include <iostream>
#include <src/math/pricers/mkl_pricer.h>
#include <tests/math/pricers/Pricer.h>


#define DECLARE_AND_DEFINE(type, x, y) \
   type x[64] __attribute__((aligned(64))); \
   for(UINT64 i=0;i<64;i++) x[i]=y;

TEST_CASE("pricer-class equals mkl-pricer (long call)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, tmp1, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp2, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp3, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp4, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp5, 0.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(UINT64, flags, 0)


    Pricer p;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_mkl_pricer();

    prepare_mkl_pricer(64, s, sigma, t, tau, r, tmp1, tmp2, tmp3, tmp4, tmp5,
                       sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);

    // test long call price
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    FLOAT reference_pricer_value = p.compute_call_price(x[0]);
    REQUIRE(abs(reference_pricer_value - price[0]) < 1.0e-8);


}


TEST_CASE("pricer-class equals mkl-pricer (short call)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, tmp1, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp2, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp3, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp4, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp5, 0.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(UINT64, flags, 2)


    Pricer p;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_mkl_pricer();

    prepare_mkl_pricer(64, s, sigma, t, tau, r, tmp1, tmp2, tmp3, tmp4, tmp5,
                       sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);

    // test long call price
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    FLOAT reference_pricer_value = p.compute_call_price(x[0]);
    REQUIRE(abs(reference_pricer_value + price[0]) < 1.0e-8);


}


TEST_CASE("pricer-class equals mkl-pricer (long put)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, tmp1, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp2, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp3, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp4, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp5, 0.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(UINT64, flags, 1)


    Pricer p;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_mkl_pricer();

    prepare_mkl_pricer(64, s, sigma, t, tau, r, tmp1, tmp2, tmp3, tmp4, tmp5,
                       sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);

    // test long call price
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    FLOAT reference_pricer_value = p.compute_put_price(x[0]);
    REQUIRE(abs(reference_pricer_value - price[0]) < 1.0e-8);

}


TEST_CASE("pricer-class equals mkl-pricer (short put)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, tmp1, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp2, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp3, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp4, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp5, 0.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(UINT64, flags, 3)


    Pricer p;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_mkl_pricer();

    prepare_mkl_pricer(64, s, sigma, t, tau, r, tmp1, tmp2, tmp3, tmp4, tmp5,
                       sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);

    // test long call price
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    FLOAT reference_pricer_value = p.compute_put_price(x[0]);
    REQUIRE(abs(reference_pricer_value + price[0]) < 1.0e-8);

}



TEST_CASE("pricer-class equals ddx-mkl-pricer (long call)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, tmp1, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp2, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp3, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp4, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp5, 0.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(UINT64, flags, 0)

    const FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_mkl_pricer();
    prepare_mkl_pricer(64, s, sigma, t, tau, r, tmp1, tmp2, tmp3, tmp4, tmp5,
                       sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    ddx_mkl_pricer(64, flags, d2, emrt, ddx_price);
    FLOAT reference_pricer_value1 = p.compute_call_price(x[0]);
    FLOAT reference_pricer_value2 = p.compute_call_price(x[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - ddx_price[0]) < 1.0e-4);

}


TEST_CASE("pricer-class equals ddx-mkl-pricer (short call)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, tmp1, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp2, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp3, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp4, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp5, 0.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(UINT64, flags, 2)

    const FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_mkl_pricer();
    prepare_mkl_pricer(64, s, sigma, t, tau, r, tmp1, tmp2, tmp3, tmp4, tmp5,
                       sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    ddx_mkl_pricer(64, flags, d2, emrt, ddx_price);
    FLOAT reference_pricer_value1 = p.compute_call_price(x[0]);
    FLOAT reference_pricer_value2 = p.compute_call_price(x[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + ddx_price[0]) < 1.0e-4);

}

TEST_CASE("pricer-class equals ddx-mkl-pricer (long put)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, tmp1, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp2, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp3, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp4, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp5, 0.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(UINT64, flags, 1)

    const FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_mkl_pricer();
    prepare_mkl_pricer(64, s, sigma, t, tau, r, tmp1, tmp2, tmp3, tmp4, tmp5,
                       sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    ddx_mkl_pricer(64, flags, d2, emrt, ddx_price);
    FLOAT reference_pricer_value1 = p.compute_put_price(x[0]);
    FLOAT reference_pricer_value2 = p.compute_put_price(x[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - ddx_price[0]) < 1.0e-4);


}

TEST_CASE("pricer-class equals ddx-mkl-pricer (short put)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, tmp1, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp2, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp3, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp4, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp5, 0.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(UINT64, flags, 3)

    const FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_mkl_pricer();
    prepare_mkl_pricer(64, s, sigma, t, tau, r, tmp1, tmp2, tmp3, tmp4, tmp5,
                       sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    ddx_mkl_pricer(64, flags, d2, emrt, ddx_price);

    FLOAT reference_pricer_value1 = p.compute_put_price(x[0]);
    FLOAT reference_pricer_value2 = p.compute_put_price(x[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + ddx_price[0]) < 1.0e-4);

}


TEST_CASE("d2dx2_pricer equals ddx-mkl-pricer diff-quot (long call)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, tmp1, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp2, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp3, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp4, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp5, 0.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(UINT64, flags, 0)
    FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_mkl_pricer();
    prepare_mkl_pricer(64, s, sigma, t, tau, r, tmp1, tmp2, tmp3, tmp4, tmp5,
                       sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    d2dx2_mkl_pricer(64, flags, s, x, d2dx2_prep, sigmaA2T2, tmp1, tmp2, d2dx2_price);

    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value1, 0.);
    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value2, 0.);

    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    ddx_mkl_pricer(64, flags, d2, emrt, reference_pricer_value1);
    for (UINT64 i = 0; i ^ 64; ++i) x[i] += eps;
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    ddx_mkl_pricer(64, flags, d2, emrt, reference_pricer_value2);


    REQUIRE(abs((reference_pricer_value2[0] - reference_pricer_value1[0]) / eps - d2dx2_price[0]) < 1.0e-4);

}


TEST_CASE("d2dx2_pricer equals ddx-mkl-pricer diff-quot (short call)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, tmp1, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp2, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp3, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp4, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp5, 0.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(UINT64, flags, 2)
    FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_mkl_pricer();
    prepare_mkl_pricer(64, s, sigma, t, tau, r, tmp1, tmp2, tmp3, tmp4, tmp5,
                       sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    d2dx2_mkl_pricer(64, flags, s, x, d2dx2_prep, sigmaA2T2, tmp1, tmp2, d2dx2_price);

    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value1, 0.);
    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value2, 0.);

    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    ddx_mkl_pricer(64, flags, d2, emrt, reference_pricer_value1);
    for (UINT64 i = 0; i ^ 64; ++i) x[i] += eps;
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    ddx_mkl_pricer(64, flags, d2, emrt, reference_pricer_value2);


    REQUIRE(abs((reference_pricer_value2[0] - reference_pricer_value1[0]) / eps - d2dx2_price[0]) < 1.0e-4);

}

TEST_CASE("d2dx2_pricer equals ddx-mkl-pricer diff-quot (long put)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, tmp1, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp2, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp3, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp4, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp5, 0.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(UINT64, flags, 1)
    FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_mkl_pricer();
    prepare_mkl_pricer(64, s, sigma, t, tau, r, tmp1, tmp2, tmp3, tmp4, tmp5,
                       sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    d2dx2_mkl_pricer(64, flags, s, x, d2dx2_prep, sigmaA2T2, tmp1, tmp2, d2dx2_price);

    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value1, 0.);
    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value2, 0.);

    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    ddx_mkl_pricer(64, flags, d2, emrt, reference_pricer_value1);
    for (UINT64 i = 0; i ^ 64; ++i) x[i] += eps;
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    ddx_mkl_pricer(64, flags, d2, emrt, reference_pricer_value2);


    REQUIRE(abs((reference_pricer_value2[0] - reference_pricer_value1[0]) / eps - d2dx2_price[0]) < 1.0e-4);

}


TEST_CASE("d2dx2_pricer equals ddx-mkl-pricer diff-quot (short put)", "[pricer]") {

    DECLARE_AND_DEFINE(FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(FLOAT, s, 70.)
    DECLARE_AND_DEFINE(FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(FLOAT, x, 72.)

    DECLARE_AND_DEFINE(FLOAT, tmp1, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp2, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp3, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp4, 0.)
    DECLARE_AND_DEFINE(FLOAT, tmp5, 0.)

    DECLARE_AND_DEFINE(FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(UINT64, flags, 3)
    FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_mkl_pricer();
    prepare_mkl_pricer(64, s, sigma, t, tau, r, tmp1, tmp2, tmp3, tmp4, tmp5,
                       sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    d2dx2_mkl_pricer(64, flags, s, x, d2dx2_prep, sigmaA2T2, tmp1, tmp2, d2dx2_price);

    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value1, 0.);
    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value2, 0.);

    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    ddx_mkl_pricer(64, flags, d2, emrt, reference_pricer_value1);
    for (UINT64 i = 0; i ^ 64; ++i) x[i] += eps;
    mkl_pricer(64, flags, s, x, sigmaA2T2, sigmaAsqrtT, emrt, tmp1, tmp2, tmp3, tmp4, d1, d2, price);
    ddx_mkl_pricer(64, flags, d2, emrt, reference_pricer_value2);


    REQUIRE(abs((reference_pricer_value2[0] - reference_pricer_value1[0]) / eps - d2dx2_price[0]) < 1.0e-4);

}

