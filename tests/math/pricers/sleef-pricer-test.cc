/*
 *
 * (c) 2019, by Peter Lennartz  // peter.lennartz@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <catch2/catch.hpp>
#include <math/pricers/sleef_pricer.h>
#include <math/pricers/Pricer.h>
#include <pricer.h>


#define DECLARE_AND_DEFINE(n, type, x, y) \
   type x[n] __attribute__((aligned(ALIGN_TO))); \
   for(UINT64 i=0;i<n;i++) x[i]=y;

TEST_CASE("pricer-class equals sleef-pricer (long call)", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, 1)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, 1)


    Pricer p;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();

    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);



    // test long call price
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    FLOAT reference_pricer_value = p.compute_call_price(x[0]);
    REQUIRE(abs(reference_pricer_value - price[0]) < 1.0e-8);

}

TEST_CASE("pricer-class equals sleef-pricer (short call)", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, 1)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, -1)


    Pricer p;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();

    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);



    // test long call price
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    FLOAT reference_pricer_value = p.compute_call_price(x[0]);
    REQUIRE(abs(reference_pricer_value + price[0]) < 1.0e-8);

}

TEST_CASE("pricer-class equals sleef-pricer (long put)", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, -1)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, 1)


    Pricer p;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();

    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);



    // test long call price
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    FLOAT reference_pricer_value = p.compute_put_price(x[0]);
    REQUIRE(abs(reference_pricer_value - price[0]) < 1.0e-8);

}

TEST_CASE("pricer-class equals sleef-pricer (short put)", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, -1)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, -1)


    Pricer p;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();

    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);



    // test long call price
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    FLOAT reference_pricer_value = p.compute_put_price(x[0]);

    REQUIRE(abs(reference_pricer_value + price[0]) < 1.0e-8);

}


TEST_CASE("pricer-class equals ddx-tw-pricer (long call)", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, ddx_price, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, 1)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, 1)

    const FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();
    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_tw_pricer(64, long_short, put_call, d2, emrt, ddx_price);
    FLOAT reference_pricer_value1 = p.compute_call_price(x[0]);
    FLOAT reference_pricer_value2 = p.compute_call_price(x[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - ddx_price[0]) < 1.0e-4);

}


TEST_CASE("pricer-class equals ddx-tw-pricer (short call)", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, ddx_price, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, -1)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, 1)

    const FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();
    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_tw_pricer(64, long_short, put_call, d2, emrt, ddx_price);
    FLOAT reference_pricer_value1 = p.compute_call_price(x[0]);
    FLOAT reference_pricer_value2 = p.compute_call_price(x[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + ddx_price[0]) < 1.0e-4);

}


TEST_CASE("pricer-class equals ddx-tw-pricer (long put)", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, ddx_price, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, 1)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, -1)

    const FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();
    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_tw_pricer(64, long_short, put_call, d2, emrt, ddx_price);
    FLOAT reference_pricer_value1 = p.compute_put_price(x[0]);
    FLOAT reference_pricer_value2 = p.compute_put_price(x[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - ddx_price[0]) < 1.0e-4);

}


TEST_CASE("pricer-class equals ddx-tw-pricer (short put)", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, ddx_price, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, -1)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, -1)

    const FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();
    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_tw_pricer(64, long_short, put_call, d2, emrt, ddx_price);
    FLOAT reference_pricer_value1 = p.compute_put_price(x[0]);
    FLOAT reference_pricer_value2 = p.compute_put_price(x[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + ddx_price[0]) < 1.0e-4);

}


TEST_CASE("d2dx2_pricer equals ddx-tw-pricer diff-quot (long call)", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, 1)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, 1)

    FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();
    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    d2dx2_tw_pricer(64, long_short, s, x, d2dx2_prep, sigmaA2T2, d2dx2_price);

    DECLARE_AND_DEFINE(64,FLOAT, reference_pricer_value1, 0.);
    DECLARE_AND_DEFINE(64,FLOAT, reference_pricer_value2, 0.);

    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_tw_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value1);
    for (UINT64 i = 0; i ^ 64; ++i) x[i] += eps;
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_tw_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value2);


    REQUIRE(abs((reference_pricer_value2[0] - reference_pricer_value1[0]) / eps - d2dx2_price[0]) < 1.0e-4);

}

TEST_CASE("d2dx2_pricer equals ddx-tw-pricer diff-quot (long put)", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, 1)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, -1)

    FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();
    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    d2dx2_tw_pricer(64, long_short, s, x, d2dx2_prep, sigmaA2T2, d2dx2_price);

    DECLARE_AND_DEFINE(64,FLOAT, reference_pricer_value1, 0.);
    DECLARE_AND_DEFINE(64,FLOAT, reference_pricer_value2, 0.);

    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_tw_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value1);
    for (UINT64 i = 0; i ^ 64; ++i) x[i] += eps;
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_tw_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value2);


    REQUIRE(abs((reference_pricer_value2[0] - reference_pricer_value1[0]) / eps - d2dx2_price[0]) < 1.0e-4);

}


TEST_CASE("d2dx2_pricer equals ddx-tw-pricer diff-quot (short call)", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, -1)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, 1)

    FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();
    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    d2dx2_tw_pricer(64, long_short, s, x, d2dx2_prep, sigmaA2T2, d2dx2_price);

    DECLARE_AND_DEFINE(64,FLOAT, reference_pricer_value1, 0.);
    DECLARE_AND_DEFINE(64,FLOAT, reference_pricer_value2, 0.);

    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_tw_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value1);
    for (UINT64 i = 0; i ^ 64; ++i) x[i] += eps;
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_tw_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value2);


    REQUIRE(abs((reference_pricer_value2[0] - reference_pricer_value1[0]) / eps - d2dx2_price[0]) < 1.0e-4);

}

TEST_CASE("used architecture", "[pricer]") {
    WARN("THIS WARNING IS SUPPOSED TO APPEAR !!! " <<
    "It is here to let the user know what kind of instruction-set is being used. " <<
    "Used instruction set is: " << getUsedInstructionSet());
    CHECK(true);
}



/*
void compute_tw_strikes_from_premiums(
        UINT64 n,
        Real_Ptr long_short_,        // 1 == long option // -1 == short option
        Real_Ptr put_call_,          // -1 == put // 1 == call
        Real_Ptr s_,
        Real_Ptr sigmaA2T2_,
        Real_Ptr sigmaAsqrtT_,
        Real_Ptr emrt_,
        Int32_Ptr to_structure,
        Real_Ptr offsets_,
        Real_Ptr prices_,
        UINT64 m,
        Real_Ptr premiums_,
        Real_Ptr instrument_prices_,
        Real_Ptr instrument_pricesl_,
        Real_Ptr instrument_pricesh_,
        Real_Ptr xl_,
        Real_Ptr xh_,
        Real_Ptr x_) {

} */

TEST_CASE("compute_tw_tsrikes_from_premiums() --- (1)", "[pricer]") {
    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(64,int32_t, to_structure, 0)
    DECLARE_AND_DEFINE(64,FLOAT, offsets, 0)
    DECLARE_AND_DEFINE(64,FLOAT, prices, 0)
    DECLARE_AND_DEFINE(64,FLOAT, x_tmp, 0)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0)

    DECLARE_AND_DEFINE(64,FLOAT, premiums, 5.)
    DECLARE_AND_DEFINE(64,FLOAT, instrument_prices,0)
    DECLARE_AND_DEFINE(64,FLOAT, instrument_pricesl,0)
    DECLARE_AND_DEFINE(64,FLOAT, instrument_pricesh,0)

    DECLARE_AND_DEFINE(64,FLOAT, xl,1)
    DECLARE_AND_DEFINE(64,FLOAT, xh,1000)

    
    DECLARE_AND_DEFINE(64,FLOAT, long_short, 1)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, 1)

    init_tw_pricer();
    prepare_tw_pricer(64, s, sigma, t, tau, r,
                      sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);

    tw_pricer(64, long_short, put_call, s, xl, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, instrument_pricesl);
    tw_pricer(64, long_short, put_call, s, xh, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, instrument_pricesh);


    for(int32_t i = 0; i < 64; ++i) {
        to_structure[i] = i;
    }

    compute_tw_strikes_from_premiums(
            64,long_short, put_call, s, sigmaA2T2, sigmaAsqrtT, emrt, to_structure, offsets, prices, x_tmp,
            64,premiums,instrument_prices, instrument_pricesl, instrument_pricesh, x, xl, xh);

    tw_pricer(64, long_short, put_call, s, xl, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, instrument_pricesl);
    tw_pricer(64, long_short, put_call, s, xh, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, instrument_pricesh);

    REQUIRE(abs(instrument_pricesl[0] - premiums[0]) < 1.0e-4);
    REQUIRE(abs(instrument_pricesh[0] - premiums[0]) < 1.0e-4);

}

TEST_CASE("d2dx2_pricer equals ddx-tw-pricer diff-quot (short put)", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, -1)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, -1)

    FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();
    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    d2dx2_tw_pricer(64, long_short, s, x, d2dx2_prep, sigmaA2T2, d2dx2_price);

    DECLARE_AND_DEFINE(64,FLOAT, reference_pricer_value1, 0.);
    DECLARE_AND_DEFINE(64,FLOAT, reference_pricer_value2, 0.);

    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_tw_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value1);
    for (UINT64 i = 0; i ^ 64; ++i) x[i] += eps;
    tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_tw_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value2);


    REQUIRE(abs((reference_pricer_value2[0] - reference_pricer_value1[0]) / eps - d2dx2_price[0]) < 1.0e-4);

}


TEST_CASE("full-tw-pricer-test", "[pricer]") {

    DECLARE_AND_DEFINE(64,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(64,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(64,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(64,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(64,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(64,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(64,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, ddx_price2, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(64,FLOAT, price, 0.)
    DECLARE_AND_DEFINE(64,FLOAT, put_call, 1)
    DECLARE_AND_DEFINE(64,FLOAT, long_short, 1)


    Pricer p;
    double eps = 1.0e-10;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_tw_pricer();
    prepare_tw_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    SECTION("long-call-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = 1;
            long_short[i] = 1;
        }
        // test long call price
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);
        FLOAT reference_pricer_value = p.compute_call_price(x[0]);
        REQUIRE(abs(reference_pricer_value - price[0]) < 1.0e-8);
    }

    SECTION("short-call-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = 1;
            long_short[i] = -1;
        }
        // test long call price
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);
        FLOAT reference_pricer_value = p.compute_call_price(x[0]);
        REQUIRE(abs(reference_pricer_value + price[0]) < 1.0e-8);
    }

    SECTION("long-put-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = -1;
            long_short[i] = 1;
        }
        // test long call price
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);
        FLOAT reference_pricer_value = p.compute_put_price(x[0]);
        REQUIRE(abs(reference_pricer_value - price[0]) < 1.0e-8);
    }

    SECTION("short-put-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = -1;
            long_short[i] = -1;
        }
        // test long call price
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);
        FLOAT reference_pricer_value = p.compute_put_price(x[0]);
        REQUIRE(abs(reference_pricer_value + price[0]) < 1.0e-8);
    }

    SECTION("long-call-derivative-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = 1;
            long_short[i] = 1;
        }
        // test long call price
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);
        FLOAT reference_pricer_value1 = p.compute_call_price(x[0]);
        FLOAT reference_pricer_value2 = p.compute_call_price(x[0] + eps);

        REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - ddx_price[0]) < 1.0e-4);
    }

    SECTION("short-call-derivative-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = 1;
            long_short[i] = -1;
        }
        // test long call price
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);
        FLOAT reference_pricer_value1 = p.compute_call_price(x[0]);
        FLOAT reference_pricer_value2 = p.compute_call_price(x[0] + eps);

        REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + ddx_price[0]) < 1.0e-4);
    }


    SECTION("long-put-derivative-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = -1;
            long_short[i] = 1;
        }
        // test long call price
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);
        FLOAT reference_pricer_value1 = p.compute_put_price(x[0]);
        FLOAT reference_pricer_value2 = p.compute_put_price(x[0] + eps);

        REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - ddx_price[0]) < 1.0e-4);
    }

    SECTION("short-put-derivative-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = -1;
            long_short[i] = -1;
        }
        // test long call price
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);
        FLOAT reference_pricer_value1 = p.compute_put_price(x[0]);
        FLOAT reference_pricer_value2 = p.compute_put_price(x[0] + eps);

        REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + ddx_price[0]) < 1.0e-4);
    }

    SECTION("long-call-second-derivative-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = 1;
            long_short[i] = 1;
        }
        // test long call price
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);

        for (UINT64 i = 0; i < 64; i++) x[i] += eps;
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price2,
                          d2dx2_price);
        for (UINT64 i = 0; i < 64; i++) x[i] -= eps;

        REQUIRE(abs((ddx_price2[0] - ddx_price[0]) / eps - d2dx2_price[0]) < 1.0e-4);
    }

    SECTION("short-call-second-derivative-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = 1;
            long_short[i] = -1;
        }
        // test long call price
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);

        for (UINT64 i = 0; i < 64; i++) x[i] += eps;
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price2,
                          d2dx2_price);
        for (UINT64 i = 0; i < 64; i++) x[i] -= eps;

        REQUIRE(abs((ddx_price2[0] - ddx_price[0]) / eps - d2dx2_price[0]) < 1.0e-4);
    }

    SECTION("long-put-second-derivative-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = -1;
            long_short[i] = 1;
        }
        // test long call price
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);

        for (UINT64 i = 0; i < 64; i++) x[i] += eps;
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price2,
                          d2dx2_price);
        for (UINT64 i = 0; i < 64; i++) x[i] -= eps;

        REQUIRE(abs((ddx_price2[0] - ddx_price[0]) / eps - d2dx2_price[0]) < 1.0e-4);
    }

    SECTION("short-put-second-derivative-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = -1;
            long_short[i] = -1;
        }
        // test long call price
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);

        for (UINT64 i = 0; i < 64; i++) x[i] += eps;
        full_tw_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price2,
                          d2dx2_price);
        for (UINT64 i = 0; i < 64; i++) x[i] -= eps;

        REQUIRE(abs((ddx_price2[0] - ddx_price[0]) / eps - d2dx2_price[0]) < 1.0e-4);
    }
}

