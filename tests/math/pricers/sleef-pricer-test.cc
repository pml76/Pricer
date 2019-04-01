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
#include <pricer-dispatcher.h>
#include <memory/context.h>

#define DECLARE_AND_DEFINE(n, x, y) \
   for(UINT64 i=0;i<n;i++) x[i]=y;

TEST_CASE("pricer-class equals sleef-pricer (long call)", "[pricer]") {

    Pricer::pricer_context context(PRICER_FLAG_TW_PRICER, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), 1)


    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
            context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);



    // test long call price
    tw_pricer(context);
    FLOAT reference_pricer_value = p.compute_call_price(context.get_x()[0]);
    REQUIRE(abs(reference_pricer_value - context.get_prices()[0]) < 1.0e-8);

}

TEST_CASE("pricer-class equals sleef-pricer (short call)", "[pricer]") {

    Pricer::pricer_context context(PRICER_FLAG_TW_PRICER, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)


    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);



    // test long call price
    tw_pricer(context);
    FLOAT reference_pricer_value = p.compute_call_price(context.get_x()[0]);
    REQUIRE(abs(reference_pricer_value + context.get_prices()[0]) < 1.0e-8);

}

TEST_CASE("pricer-class equals sleef-pricer (long put)", "[pricer]") {

    Pricer::pricer_context context(PRICER_FLAG_TW_PRICER, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), 1)


    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);



    // test long call price
    tw_pricer(context);
    FLOAT reference_pricer_value = p.compute_put_price(context.get_x()[0]);
    REQUIRE(abs(reference_pricer_value - context.get_prices()[0]) < 1.0e-8);

}

TEST_CASE("pricer-class equals sleef-pricer (short put)", "[pricer]") {

    Pricer::pricer_context context(PRICER_FLAG_TW_PRICER, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)


    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);



    // test long call price
    tw_pricer(context);
    FLOAT reference_pricer_value = p.compute_put_price(context.get_x()[0]);

    REQUIRE(abs(reference_pricer_value + context.get_prices()[0]) < 1.0e-8);

}


TEST_CASE("pricer-class equals ddx-tw-pricer (long call)", "[pricer]") {

    Pricer::pricer_context context(PRICER_FLAG_TW_COMPUTE_DDX, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_ddx_price(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), 1)

    const FLOAT eps = 1.0e-10;
    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);


    // test long call price
    tw_pricer(context);
    ddx_tw_pricer(context);
    FLOAT reference_pricer_value1 = p.compute_call_price(context.get_x()[0]);
    FLOAT reference_pricer_value2 = p.compute_call_price(context.get_x()[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[0]) < 1.0e-4);

}


TEST_CASE("pricer-class equals ddx-tw-pricer (short call)", "[pricer]") {

    Pricer::pricer_context context(PRICER_FLAG_TW_COMPUTE_DDX, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_ddx_price(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)

    const FLOAT eps = 1.0e-10;
    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);


    // test long call price
    tw_pricer(context);
    ddx_tw_pricer(context);
    FLOAT reference_pricer_value1 = p.compute_call_price(context.get_x()[0]);
    FLOAT reference_pricer_value2 = p.compute_call_price(context.get_x()[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[0]) < 1.0e-4);

}


TEST_CASE("pricer-class equals ddx-tw-pricer (long put)", "[pricer]") {

    Pricer::pricer_context context(PRICER_FLAG_TW_COMPUTE_DDX, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_ddx_price(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), 1)

    const FLOAT eps = 1.0e-10;
    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);


    // test long call price
    tw_pricer(context);
    ddx_tw_pricer(context);
    FLOAT reference_pricer_value1 = p.compute_put_price(context.get_x()[0]);
    FLOAT reference_pricer_value2 = p.compute_put_price(context.get_x()[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[0]) < 1.0e-4);

}


TEST_CASE("pricer-class equals ddx-tw-pricer (short put)", "[pricer]") {

    Pricer::pricer_context context(PRICER_FLAG_TW_COMPUTE_DDX, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_ddx_price(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)

    const FLOAT eps = 1.0e-10;
    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);


    // test long call price
    tw_pricer(context);
    ddx_tw_pricer(context);
    FLOAT reference_pricer_value1 = p.compute_put_price(context.get_x()[0]);
    FLOAT reference_pricer_value2 = p.compute_put_price(context.get_x()[0] + eps);

    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[0]) < 1.0e-4);

}


TEST_CASE("d2dx2_pricer equals ddx-tw-pricer diff-quot (long call)", "[pricer]") {


    Pricer::pricer_context context(PRICER_FLAG_TW_COMPUTE_D2DX2, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_ddx_price(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), 1)

    const FLOAT eps = 1.0e-4;
    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);


    /****************/

    // test long call price
    tw_pricer(context);
    d2dx2_tw_pricer(context);

    FLOAT reference_pricer_value1;
    FLOAT reference_pricer_value2;

    tw_pricer(context);
    ddx_tw_pricer(context);
    reference_pricer_value1 = context.get_ddx_price()[0];

    for (UINT64 i = 0; i ^ 64; ++i) context.get_x()[i] += eps;
    tw_pricer(context);
    ddx_tw_pricer(context);
    reference_pricer_value2 = context.get_ddx_price()[0];


    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2()[0]) < eps);

}

TEST_CASE("d2dx2_pricer equals ddx-tw-pricer diff-quot (long put)", "[pricer]") {


    Pricer::pricer_context context(PRICER_FLAG_TW_COMPUTE_D2DX2, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_ddx_price(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), 1)

    const FLOAT eps = 1.0e-4;
    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);


    /****************/

    // test long call price
    tw_pricer(context);
    d2dx2_tw_pricer(context);

    FLOAT reference_pricer_value1;
    FLOAT reference_pricer_value2;

    tw_pricer(context);
    ddx_tw_pricer(context);
    reference_pricer_value1 = context.get_ddx_price()[0];

    for (UINT64 i = 0; i ^ 64; ++i) context.get_x()[i] += eps;
    tw_pricer(context);
    ddx_tw_pricer(context);
    reference_pricer_value2 = context.get_ddx_price()[0];


    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2()[0]) < eps);
}


TEST_CASE("d2dx2_pricer equals ddx-tw-pricer diff-quot (short call)", "[pricer]") {

    Pricer::pricer_context context(PRICER_FLAG_TW_COMPUTE_D2DX2, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_ddx_price(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)

    const FLOAT eps = 1.0e-4;
    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);


    /****************/

    // test long call price
    tw_pricer(context);
    d2dx2_tw_pricer(context);

    FLOAT reference_pricer_value1;
    FLOAT reference_pricer_value2;

    tw_pricer(context);
    ddx_tw_pricer(context);
    reference_pricer_value1 = context.get_ddx_price()[0];

    for (UINT64 i = 0; i ^ 64; ++i) context.get_x()[i] += eps;
    tw_pricer(context);
    ddx_tw_pricer(context);
    reference_pricer_value2 = context.get_ddx_price()[0];


    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2()[0]) < eps);

}

TEST_CASE("used architecture", "[pricer]") {
    WARN("THIS WARNING IS SUPPOSED TO APPEAR !!! " <<
    "It is here to let the user know what kind of instruction-set is being used. " <<
    "Used instruction set is: " << getUsedInstructionSet());
    CHECK(true);
}



TEST_CASE("compute_tw_tsrikes_from_premiums() --- (1)", "[pricer]") {

    Pricer::pricer_context context(PRICER_FLAG_TW_PRICER, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)


    DECLARE_AND_DEFINE(ALIGN_TO, context.get_to_structure(), 0)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_offsets(), 0)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_premiums(), 5.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_instrument_prices(),0)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_instrument_pricesl(),0)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_instrument_pricesh(),0)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_xl_(),1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_xh_(),1000)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x_(),500)

    

    init_tw_pricer();
    prepare_tw_pricer(context);



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

TEST_CASE("compute_tw_tsrikes_from_premiums() --- (2)", "[pricer]") {
    DECLARE_AND_DEFINE(128,FLOAT, r, 0.01)
    DECLARE_AND_DEFINE(128,FLOAT, s, 70.)
    DECLARE_AND_DEFINE(128,FLOAT, t, 1.2)
    DECLARE_AND_DEFINE(128,FLOAT, tau, 1. / 12.)
    DECLARE_AND_DEFINE(128,FLOAT, sigma, 0.3)
    DECLARE_AND_DEFINE(128,FLOAT, x, 72.)

    DECLARE_AND_DEFINE(128,FLOAT, sigmaA, 0.)
    DECLARE_AND_DEFINE(128,FLOAT, sigmaA2T2, 0.)
    DECLARE_AND_DEFINE(128,FLOAT, sigmaAsqrtT, 0.)
    DECLARE_AND_DEFINE(128,FLOAT, emrt, 0.)
    DECLARE_AND_DEFINE(128,FLOAT, d2dx2_prep, 0.)

    DECLARE_AND_DEFINE(128,int32_t, to_structure, 0)
    DECLARE_AND_DEFINE(128,FLOAT, offsets, 0)
    DECLARE_AND_DEFINE(128,FLOAT, prices, 0)
    DECLARE_AND_DEFINE(128,FLOAT, x_tmp, 0)

    DECLARE_AND_DEFINE(128,FLOAT, d1, 0)
    DECLARE_AND_DEFINE(128,FLOAT, d2, 0)

    DECLARE_AND_DEFINE(64,FLOAT, premiums, 5.)
    DECLARE_AND_DEFINE(64,FLOAT, instrument_prices,0)
    DECLARE_AND_DEFINE(64,FLOAT, instrument_pricesl,0)
    DECLARE_AND_DEFINE(64,FLOAT, instrument_pricesh,0)

    DECLARE_AND_DEFINE(64,FLOAT, xl,1)
    DECLARE_AND_DEFINE(64,FLOAT, xh,1000)


    DECLARE_AND_DEFINE(128,FLOAT, long_short, 1)
    DECLARE_AND_DEFINE(128,FLOAT, put_call, 1)

    init_tw_pricer();
    prepare_tw_pricer(128, s, sigma, t, tau, r,
                      sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep);

    tw_pricer(128, long_short, put_call, s, xl, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, instrument_pricesl);
    tw_pricer(128, long_short, put_call, s, xh, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, instrument_pricesh);


    for(int32_t i = 0; i < 64; ++i) {
        to_structure[2*i]   = i;
        to_structure[2*i+1] = i;

        put_call[2*i]       = 1;
        put_call[2*i+1]     = -1;

        long_short[2*i]     = 1;
        long_short[2*i+1]   = -1;


    }

    compute_tw_strikes_from_premiums(
            128,long_short, put_call, s, sigmaA2T2, sigmaAsqrtT, emrt, to_structure, offsets, prices, x_tmp,
            64,premiums,instrument_prices, instrument_pricesl, instrument_pricesh, x, xl, xh);

    tw_pricer(128, long_short, put_call, s, xl, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, instrument_pricesl);
    tw_pricer(128, long_short, put_call, s, xh, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, instrument_pricesh);

    REQUIRE(abs(instrument_pricesl[0] - premiums[0]) < 1.0e-4);
    REQUIRE(abs(instrument_pricesh[0] - premiums[0]) < 1.0e-4);

}

TEST_CASE("d2dx2_pricer equals ddx-tw-pricer diff-quot (short put)", "[pricer]") {


    Pricer::pricer_context context(PRICER_FLAG_TW_COMPUTE_D2DX2, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_ddx_price(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)

    const FLOAT eps = 1.0e-4;
    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);


    /****************/

    // test long call price
    tw_pricer(context);
    d2dx2_tw_pricer(context);

    FLOAT reference_pricer_value1;
    FLOAT reference_pricer_value2;

    tw_pricer(context);
    ddx_tw_pricer(context);
    reference_pricer_value1 = context.get_ddx_price()[0];

    for (UINT64 i = 0; i ^ 64; ++i) context.get_x()[i] += eps;
    tw_pricer(context);
    ddx_tw_pricer(context);
    reference_pricer_value2 = context.get_ddx_price()[0];


    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2()[0]) < eps);
}


TEST_CASE("full-tw-pricer-test", "[pricer]") {


    Pricer::pricer_context context(PRICER_FLAG_TW_COMPUTE_D2DX2, ALIGN_TO,0);

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_r(), 0.01)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_s(), 70.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_t(), 1.2)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_tau(), 1. / 12.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigma(), 0.3)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_x(), 72.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_ddx_price(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_emrt(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d1(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_d2(), 0.)

    DECLARE_AND_DEFINE(ALIGN_TO, context.get_prices(), 0.)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
    DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)

    const FLOAT eps = 1.0e-4;
    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);


    /****************/


    SECTION("long-call-test") {
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), 1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value = p.compute_call_price(context.get_x()[0]);
        REQUIRE(abs(reference_pricer_value - context.get_prices()[0]) < 1.0e-8);
    }

    SECTION("short-call-test") {
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value = p.compute_call_price(context.get_x()[0]);
        REQUIRE(abs(reference_pricer_value + context.get_prices()[0]) < 1.0e-8);
    }

    SECTION("long-put-test") {
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), 1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value = p.compute_put_price(context.get_x()[0]);
        REQUIRE(abs(reference_pricer_value - context.get_prices()[0]) < 1.0e-8);
    }

    SECTION("short-put-test") {
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value = p.compute_put_price(context.get_x()[0]);
        REQUIRE(abs(reference_pricer_value + context.get_prices()[0]) < 1.0e-8);
    }

    SECTION("long-call-derivative-test") {
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), 1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value1 = p.compute_call_price(context.get_x()[0]);
        FLOAT reference_pricer_value2 = p.compute_call_price(context.get_x()[0] + eps);

        REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[0]) < 1.0e-4);
    }

    SECTION("short-call-derivative-test") {
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value1 = p.compute_call_price(context.get_x()[0]);
        FLOAT reference_pricer_value2 = p.compute_call_price(context.get_x()[0] + eps);

        REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[0]) < 1.0e-4);
    }


    SECTION("long-put-derivative-test") {
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), 1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value1 = p.compute_put_price(context.get_x()[0]);
        FLOAT reference_pricer_value2 = p.compute_put_price(context.get_x()[0] + eps);

        REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[0]) < 1.0e-4);
    }

    SECTION("short-put-derivative-test") {
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value1 = p.compute_put_price(context.get_x()[0]);
        FLOAT reference_pricer_value2 = p.compute_put_price(context.get_x()[0] + eps);

        REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[0]) < 1.0e-4);
    }

    SECTION("long-call-second-derivative-test") {
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), 1)
        // test long call price
        full_tw_pricer(context);

        FLOAT ddx_price2 = context.get_ddx_price()[0];

        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] += eps;
        full_tw_pricer(context);
        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] -= eps;

        REQUIRE(abs((ddx_price2 - context.get_ddx_price()[0]) / eps - context.get_d2dx2()[0]) < 1.0e-4);
    }

    SECTION("short-call-second-derivative-test") {
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), 1)
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)
        // test long call price
        full_tw_pricer(context);

        FLOAT ddx_price2 = context.get_ddx_price()[0];

        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] += eps;
        full_tw_pricer(context);
        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] -= eps;

        REQUIRE(abs((ddx_price2 - context.get_ddx_price()[0]) / eps - context.get_d2dx2()[0]) < 1.0e-4);
    }

    SECTION("long-put-second-derivative-test") {
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), 1)
        // test long call price
        full_tw_pricer(context);

        FLOAT ddx_price2 = context.get_ddx_price()[0];

        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] += eps;
        full_tw_pricer(context);
        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] -= eps;

        REQUIRE(abs((ddx_price2 - context.get_ddx_price()[0]) / eps - context.get_d2dx2()[0]) < 1.0e-4);
    }

    SECTION("short-put-second-derivative-test") {
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_put_call(), -1)
        DECLARE_AND_DEFINE(ALIGN_TO, context.get_long_short(), -1)
        // test long call price
        full_tw_pricer(context);

        FLOAT ddx_price2 = context.get_ddx_price()[0];

        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] += eps;
        full_tw_pricer(context);
        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] -= eps;

        REQUIRE(abs((ddx_price2 - context.get_ddx_price()[0]) / eps - context.get_d2dx2()[0]) < 1.0e-4);
    }
}

