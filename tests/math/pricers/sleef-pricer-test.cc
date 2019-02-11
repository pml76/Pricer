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
#include <src/math/pricers/sleef_pricer.h>
#include <tests/math/pricers/Pricer.h>


#define DECLARE_AND_DEFINE(type, x, y) \
   (type x[64] __attribute__((aligned(ALIGN_TO))); \
   for(UINT64 i=0;i<ALIGN_TO;i++) x[i]=y);

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


TEST_CASE("pricer-class equals ddx-sleef-pricer (long call)", "[pricer]") {

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
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, long_short, 1)
    DECLARE_AND_DEFINE(FLOAT, put_call, 1)

    const FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();
    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_sleef_pricer(64, long_short, put_call, d2, emrt, ddx_price);
    FLOAT reference_pricer_value1 = p.compute_call_price(x[0]);
    FLOAT reference_pricer_value2 = p.compute_call_price(x[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - ddx_price[0]) < 1.0e-4);

}


TEST_CASE("pricer-class equals ddx-sleef-pricer (short call)", "[pricer]") {

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
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, long_short, -1)
    DECLARE_AND_DEFINE(FLOAT, put_call, 1)

    const FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();
    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_sleef_pricer(64, long_short, put_call, d2, emrt, ddx_price);
    FLOAT reference_pricer_value1 = p.compute_call_price(x[0]);
    FLOAT reference_pricer_value2 = p.compute_call_price(x[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + ddx_price[0]) < 1.0e-4);

}


TEST_CASE("pricer-class equals ddx-sleef-pricer (long put)", "[pricer]") {

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
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, long_short, 1)
    DECLARE_AND_DEFINE(FLOAT, put_call, -1)

    const FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();
    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_sleef_pricer(64, long_short, put_call, d2, emrt, ddx_price);
    FLOAT reference_pricer_value1 = p.compute_put_price(x[0]);
    FLOAT reference_pricer_value2 = p.compute_put_price(x[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - ddx_price[0]) < 1.0e-4);

}


TEST_CASE("pricer-class equals ddx-sleef-pricer (short put)", "[pricer]") {

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
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, long_short, -1)
    DECLARE_AND_DEFINE(FLOAT, put_call, -1)

    const FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();
    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_sleef_pricer(64, long_short, put_call, d2, emrt, ddx_price);
    FLOAT reference_pricer_value1 = p.compute_put_price(x[0]);
    FLOAT reference_pricer_value2 = p.compute_put_price(x[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + ddx_price[0]) < 1.0e-4);

}


TEST_CASE("d2dx2_pricer equals ddx-sleef-pricer diff-quot (long call)", "[pricer]") {

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
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, long_short, 1)
    DECLARE_AND_DEFINE(FLOAT, put_call, 1)

    FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();
    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    d2dx2_sleef_pricer(64, long_short, s, x, d2dx2_prep, sigmaA2T2, d2dx2_price);

    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value1, 0.);
    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value2, 0.);

    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_sleef_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value1);
    for (UINT64 i = 0; i ^ 64; ++i) x[i] += eps;
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_sleef_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value2);


    REQUIRE(abs((reference_pricer_value2[0] - reference_pricer_value1[0]) / eps - d2dx2_price[0]) < 1.0e-4);

}

TEST_CASE("d2dx2_pricer equals ddx-sleef-pricer diff-quot (long put)", "[pricer]") {

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
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, long_short, 1)
    DECLARE_AND_DEFINE(FLOAT, put_call, -1)

    FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();
    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    d2dx2_sleef_pricer(64, long_short, s, x, d2dx2_prep, sigmaA2T2, d2dx2_price);

    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value1, 0.);
    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value2, 0.);

    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_sleef_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value1);
    for (UINT64 i = 0; i ^ 64; ++i) x[i] += eps;
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_sleef_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value2);


    REQUIRE(abs((reference_pricer_value2[0] - reference_pricer_value1[0]) / eps - d2dx2_price[0]) < 1.0e-4);

}


TEST_CASE("d2dx2_pricer equals ddx-sleef-pricer diff-quot (short call)", "[pricer]") {

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
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, long_short, -1)
    DECLARE_AND_DEFINE(FLOAT, put_call, 1)

    FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();
    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    d2dx2_sleef_pricer(64, long_short, s, x, d2dx2_prep, sigmaA2T2, d2dx2_price);

    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value1, 0.);
    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value2, 0.);

    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_sleef_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value1);
    for (UINT64 i = 0; i ^ 64; ++i) x[i] += eps;
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_sleef_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value2);


    REQUIRE(abs((reference_pricer_value2[0] - reference_pricer_value1[0]) / eps - d2dx2_price[0]) < 1.0e-4);

}

TEST_CASE("d2dx2_pricer equals ddx-sleef-pricer diff-quot (short put)", "[pricer]") {

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
    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, d1, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, long_short, -1)
    DECLARE_AND_DEFINE(FLOAT, put_call, -1)

    FLOAT eps = 1.0e-10;

    Pricer p;

    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();
    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    // test long call price
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    d2dx2_sleef_pricer(64, long_short, s, x, d2dx2_prep, sigmaA2T2, d2dx2_price);

    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value1, 0.);
    DECLARE_AND_DEFINE(FLOAT, reference_pricer_value2, 0.);

    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_sleef_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value1);
    for (UINT64 i = 0; i ^ 64; ++i) x[i] += eps;
    sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price);
    ddx_sleef_pricer(64, long_short, put_call, d2, emrt, reference_pricer_value2);


    REQUIRE(abs((reference_pricer_value2[0] - reference_pricer_value1[0]) / eps - d2dx2_price[0]) < 1.0e-4);

}


TEST_CASE("full-sleef-pricer-test", "[pricer]") {

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

    DECLARE_AND_DEFINE(FLOAT, ddx_price, 0.)
    DECLARE_AND_DEFINE(FLOAT, ddx_price2, 0.)
    DECLARE_AND_DEFINE(FLOAT, d2dx2_price, 0.)

    DECLARE_AND_DEFINE(FLOAT, price, 0.)
    DECLARE_AND_DEFINE(FLOAT, put_call, 1)
    DECLARE_AND_DEFINE(FLOAT, long_short, 1)


    Pricer p;
    double eps = 1.0e-10;


    p.set_market_data(sigma[0], t[0], tau[0], r[0], s[0]);

    init_sleef_pricer();
    prepare_sleef_pricer(64, s, sigma, t, tau, r,
                         sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);


    SECTION("long-call-test") {
        for (UINT64 i = 0; i < 64; ++i) {
            put_call[i] = 1;
            long_short[i] = 1;
        }
        // test long call price
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
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
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
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
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
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
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
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
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
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
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
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
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
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
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
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
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);

        for (UINT64 i = 0; i < 64; i++) x[i] += eps;
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price2,
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
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);

        for (UINT64 i = 0; i < 64; i++) x[i] += eps;
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price2,
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
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);

        for (UINT64 i = 0; i < 64; i++) x[i] += eps;
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price2,
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
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price,
                          d2dx2_price);

        for (UINT64 i = 0; i < 64; i++) x[i] += eps;
        full_sleef_pricer(64, long_short, put_call, s, x, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep, price, ddx_price2,
                          d2dx2_price);
        for (UINT64 i = 0; i < 64; i++) x[i] -= eps;

        REQUIRE(abs((ddx_price2[0] - ddx_price[0]) / eps - d2dx2_price[0]) < 1.0e-4);
    }
}

