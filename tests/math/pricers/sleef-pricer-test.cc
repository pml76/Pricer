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

#define SET_EQUAL_TO(n, x, y) \
   for(UINT64 i=0;i<n;i++) x[i]=y;
   
   
SCENARIO("Pricer::pricer_context implementation") {

    GIVEN("a context of a size")


    GIVEN("a context for simple price computations") {
        Pricer::pricer_context context(PRICER_FLAG_TW_PRICER, ALIGN_TO, 0);
        
        THEN("only those entries should have values not equal to zero that are used") {
            REQUIRE(context.get_s()           != nullptr);
            REQUIRE(context.get_sigma()       != nullptr);
            REQUIRE(context.get_t()           != nullptr);
            REQUIRE(context.get_tau()         != nullptr);
            REQUIRE(context.get_r()           != nullptr);
            REQUIRE(context.get_sigmaA()      != nullptr);
            REQUIRE(context.get_sigmaA2T2()   != nullptr);
            REQUIRE(context.get_sigmaAsqrtT() != nullptr);
            REQUIRE(context.get_emrt()        != nullptr);
            REQUIRE(context.get_long_short()  != nullptr);
            REQUIRE(context.get_put_call()    != nullptr);
            REQUIRE(context.get_d2dx2_prep()  != nullptr);
            REQUIRE(context.get_d1()          != nullptr);
            REQUIRE(context.get_d2()          != nullptr);
            REQUIRE(context.get_prices()      != nullptr);
            REQUIRE(context.get_x()           != nullptr);

            REQUIRE(context.get_to_structure()       == nullptr);
            REQUIRE(context.get_offsets()            == nullptr);
            REQUIRE(context.get_premiums()           == nullptr);
            REQUIRE(context.get_instrument_prices()  == nullptr);
            REQUIRE(context.get_instrument_pricesl() == nullptr);
            REQUIRE(context.get_instrument_pricesh() == nullptr);
            REQUIRE(context.get_x_()                 == nullptr);
            REQUIRE(context.get_xl_()                == nullptr);
            REQUIRE(context.get_xh_()                == nullptr);
            REQUIRE(context.get_ddx_price()          == nullptr);
            REQUIRE(context.get_d2dx2()              == nullptr);
        }

        THEN("the pointers are ALIGN_TO-bytes aligned (needed for usage of SIMD instructions)") {
            REQUIRE((((unsigned long long int) context.get_s())           & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigma())       & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_t())           & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_tau())         & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_r())           & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaA())      & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaA2T2())   & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaAsqrtT()) & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_emrt())        & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_long_short())  & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_put_call())    & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d2dx2_prep())  & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d1())          & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d2())          & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_prices())      & (ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_x())           & (ALIGN_TO-1)) == 0);
        }


    }

}
   

SCENARIO("Pricer computations are correct") {

    GIVEN("An instance of a Pricer::pricer_context class and a correspondingly initialized Vortex::Pricer object") {

        Pricer::pricer_context context(PRICER_FLAG_TW_PRICER, ALIGN_TO, 0);

        SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
        SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
        SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
        SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
        SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
        SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)

        SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
        SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
        SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
        SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
        SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

        SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
        SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

        SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)


        Vortex::Pricer p;


        p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                          context.get_r()[0], context.get_s()[0]);

        init_tw_pricer();

        prepare_tw_pricer(context);


        WHEN("we price a long call") {
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)
            // test long call price
            tw_pricer(context);
            FLOAT reference_pricer_value = p.compute_call_price(context.get_x()[0]);

            THEN("the computed prices of Vortex::Pricer and tw_pricer() are equal") {
                REQUIRE(abs(reference_pricer_value - context.get_prices()[0]) < 1.0e-8);
            }
        }

        WHEN("we price a short call") {
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)

            // compute the short call price
            tw_pricer(context);
            FLOAT reference_pricer_value = p.compute_call_price(context.get_x()[0]);
            THEN("the computed prices of the Vortex::Pricer and tw_pricer() are equal") {
                REQUIRE(abs(reference_pricer_value + context.get_prices()[0]) < 1.0e-8);
            }
        }

        WHEN("we price a long put") {
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)

            // test long call price
            tw_pricer(context);
            FLOAT reference_pricer_value = p.compute_put_price(context.get_x()[0]);

            THEN("the computed prices of the Vortex::Pricer and tw_pricer() are equal ") {
                REQUIRE(abs(reference_pricer_value - context.get_prices()[0]) < 1.0e-8);
            }
        }

        WHEN("we price a short put"){
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)

            tw_pricer(context);
            FLOAT reference_pricer_value = p.compute_put_price(context.get_x()[0]);

            THEN("the computed prices of Vortex::Pricer and tw_pricer() are equal") {
                REQUIRE(abs(reference_pricer_value + context.get_prices()[0]) < 1.0e-8);

            }
        }
    }
}


TEST_CASE("pricer-class equals sleef-pricer (short put)", "[pricer]") {

    Pricer::pricer_context context(PRICER_FLAG_TW_PRICER, ALIGN_TO,0);

    SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
    SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
    SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
    SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
    SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)

    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
    SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)


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

    SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
    SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
    SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
    SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
    SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)
    SET_EQUAL_TO(ALIGN_TO, context.get_ddx_price(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
    SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)

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

    SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
    SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
    SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
    SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
    SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)
    SET_EQUAL_TO(ALIGN_TO, context.get_ddx_price(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
    SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)

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

    SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
    SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
    SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
    SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
    SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)
    SET_EQUAL_TO(ALIGN_TO, context.get_ddx_price(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
    SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)

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

    SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
    SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
    SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
    SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
    SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)
    SET_EQUAL_TO(ALIGN_TO, context.get_ddx_price(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
    SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)

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

    SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
    SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
    SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
    SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
    SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)
    SET_EQUAL_TO(ALIGN_TO, context.get_ddx_price(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
    SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)

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

    SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
    SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
    SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
    SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
    SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)
    SET_EQUAL_TO(ALIGN_TO, context.get_ddx_price(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
    SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)

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

    SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
    SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
    SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
    SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
    SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)
    SET_EQUAL_TO(ALIGN_TO, context.get_ddx_price(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
    SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)

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



TEST_CASE("compute_tw_strikes_from_premiums() --- (1)", "[pricer]") {

    Pricer::pricer_context context(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, ALIGN_TO,ALIGN_TO);

    SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
    SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
    SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
    SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
    SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)

    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
    SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)


    SET_EQUAL_TO(ALIGN_TO, context.get_to_structure(), 0)
    SET_EQUAL_TO(ALIGN_TO, context.get_offsets(), 0)
    SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0)

    SET_EQUAL_TO(ALIGN_TO, context.get_premiums(), 5.)
    SET_EQUAL_TO(ALIGN_TO, context.get_instrument_prices(),0)
    SET_EQUAL_TO(ALIGN_TO, context.get_instrument_pricesl(),0)
    SET_EQUAL_TO(ALIGN_TO, context.get_instrument_pricesh(),0)

    SET_EQUAL_TO(ALIGN_TO, context.get_xl_(),1)
    SET_EQUAL_TO(ALIGN_TO, context.get_xh_(),1000)
    SET_EQUAL_TO(ALIGN_TO, context.get_x_(),500)

    

    init_tw_pricer();
    prepare_tw_pricer(context);

    for(int32_t i = 0; i < 64; ++i) {
        context.get_to_structure()[i] = i;
    }

    compute_tw_strikes_from_premiums(context);

    REQUIRE(abs(context.get_instrument_pricesl()[0] - context.get_premiums()[0]) < 1.0e-4);
    // REQUIRE(abs(context.get_instrument_pricesh()[0] - context.get_premiums()[0]) < 1.0e-4);

}


TEST_CASE("d2dx2_pricer equals ddx-tw-pricer diff-quot (short put)", "[pricer]") {


    Pricer::pricer_context context(PRICER_FLAG_TW_COMPUTE_D2DX2, ALIGN_TO,0);

    SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
    SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
    SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
    SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
    SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)
    SET_EQUAL_TO(ALIGN_TO, context.get_ddx_price(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
    SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)

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

    SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
    SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
    SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
    SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
    SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)
    SET_EQUAL_TO(ALIGN_TO, context.get_ddx_price(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

    SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)
    SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
    SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)

    const FLOAT eps = 1.0e-4;
    Vortex::Pricer p;


    p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                      context.get_r()[0], context.get_s()[0]);

    init_tw_pricer();

    prepare_tw_pricer(context);


    /****************/


    SECTION("long-call-test") {
        SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
        SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value = p.compute_call_price(context.get_x()[0]);
        REQUIRE(abs(reference_pricer_value - context.get_prices()[0]) < 1.0e-8);
    }

    SECTION("short-call-test") {
        SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
        SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value = p.compute_call_price(context.get_x()[0]);
        REQUIRE(abs(reference_pricer_value + context.get_prices()[0]) < 1.0e-8);
    }

    SECTION("long-put-test") {
        SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
        SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value = p.compute_put_price(context.get_x()[0]);
        REQUIRE(abs(reference_pricer_value - context.get_prices()[0]) < 1.0e-8);
    }

    SECTION("short-put-test") {
        SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
        SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value = p.compute_put_price(context.get_x()[0]);
        REQUIRE(abs(reference_pricer_value + context.get_prices()[0]) < 1.0e-8);
    }

    SECTION("long-call-derivative-test") {
        SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
        SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value1 = p.compute_call_price(context.get_x()[0]);
        FLOAT reference_pricer_value2 = p.compute_call_price(context.get_x()[0] + eps);

        FLOAT diff_quot = (reference_pricer_value2 - reference_pricer_value1) / eps;

        REQUIRE(abs(diff_quot - context.get_ddx_price()[0]) < 1.0e-4);
    }

    SECTION("short-call-derivative-test") {
        SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
        SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value1 = p.compute_call_price(context.get_x()[0]);
        FLOAT reference_pricer_value2 = p.compute_call_price(context.get_x()[0] + eps);

        REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[0]) < 1.0e-4);
    }


    SECTION("long-put-derivative-test") {
        SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
        SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value1 = p.compute_put_price(context.get_x()[0]);
        FLOAT reference_pricer_value2 = p.compute_put_price(context.get_x()[0] + eps);

        REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[0]) < 1.0e-4);
    }

    SECTION("short-put-derivative-test") {
        SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
        SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)
        // test long call price
        full_tw_pricer(context);
        FLOAT reference_pricer_value1 = p.compute_put_price(context.get_x()[0]);
        FLOAT reference_pricer_value2 = p.compute_put_price(context.get_x()[0] + eps);

        REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[0]) < 1.0e-4);
    }

    SECTION("long-call-second-derivative-test") {
        SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
        SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)
        // test long call price
        full_tw_pricer(context);

        FLOAT ddx_price = context.get_ddx_price()[0];
        FLOAT d2dx2 = context.get_d2dx2()[0];

        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] += eps;
        full_tw_pricer(context);
        FLOAT ddx_price2 = context.get_ddx_price()[0];
        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] -= eps;

        FLOAT diff_quot = (ddx_price2 - ddx_price) / eps;

        REQUIRE(abs(diff_quot - d2dx2) < 1.0e-4);
    }

    SECTION("short-call-second-derivative-test") {
        SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
        SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)
        // test long call price
        full_tw_pricer(context);

        FLOAT ddx_price = context.get_ddx_price()[0];
        FLOAT d2dx2 = context.get_d2dx2()[0];

        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] += eps;
        full_tw_pricer(context);
        FLOAT ddx_price2 = context.get_ddx_price()[0];
        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] -= eps;

        FLOAT diff_quot = (ddx_price2 - ddx_price) / eps;

        REQUIRE(abs(diff_quot - d2dx2) < 1.0e-4);    }

    SECTION("long-put-second-derivative-test") {
        SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
        SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)
        // test long call price
        full_tw_pricer(context);

        FLOAT ddx_price = context.get_ddx_price()[0];
        FLOAT d2dx2 = context.get_d2dx2()[0];

        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] += eps;
        full_tw_pricer(context);
        FLOAT ddx_price2 = context.get_ddx_price()[0];
        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] -= eps;

        FLOAT diff_quot = (ddx_price2 - ddx_price) / eps;

        REQUIRE(abs(diff_quot - d2dx2) < 1.0e-4);    }

    SECTION("short-put-second-derivative-test") {
        SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
        SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)
        // test long call price
        full_tw_pricer(context);

        FLOAT ddx_price = context.get_ddx_price()[0];
        FLOAT d2dx2 = context.get_d2dx2()[0];

        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] += eps;
        full_tw_pricer(context);
        FLOAT ddx_price2 = context.get_ddx_price()[0];
        for (UINT64 i = 0; i < ALIGN_TO; i++) context.get_x()[i] -= eps;

        FLOAT diff_quot = (ddx_price2 - ddx_price) / eps;

        REQUIRE(abs(diff_quot - d2dx2) < 1.0e-4);
    }
}

