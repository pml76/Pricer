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
#include <include/sleef_pricer.h>
#include <math/pricers/Pricer.h>
#include <pricer-dispatcher.h>
#include <include/context.h>

#define SET_EQUAL_TO(n, x, y) \
   for(uint64_t i=0;i<n;i++) x[i]=y;
   
   
SCENARIO("Pricer::pricer_context implementation") {

    GIVEN("a request for a context of size < ALIGN_TO") {
        Pricer::pricer_context context(ALIGN_TO/4);

        THEN("we get a context of size equal to ALIGN_TO") {
            REQUIRE(context.get_n_max() == ALIGN_TO);
        }
    }

    GIVEN("a request for a context of size x where x % ALIGN_TO != 0") {
        Pricer::pricer_context context(129);

        THEN("we get a context of size y where y > x and y % ALIGN_TO == 0") {
            REQUIRE(context.get_n_max() > 129);
            REQUIRE((context.get_n_max() % ALIGN_TO) == 0);
        }
    }

    GIVEN("a context for simple price computations") {
        Pricer::pricer_context context( ALIGN_TO );
        
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

            /*
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
             */
        }

        THEN("the pointers are ALIGN_TO-bytes aligned (needed for usage of SIMD instructions)") {
            REQUIRE((((unsigned long long int) context.get_s())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigma())       & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_t())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_tau())         & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_r())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaA())      & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaA2T2())   & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaAsqrtT()) & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_emrt())        & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_long_short())  & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_put_call())    & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d2dx2_prep())  & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d1())          & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d2())          & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_prices())      & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_x())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
        }


    }

}


SCENARIO("Pricer::ddx_pricer_context implementation") {

    GIVEN("a request for a context of size < ALIGN_TO") {
        Pricer::ddx_pricer_context context(ALIGN_TO/4);

        THEN("we get a context of size equal to ALIGN_TO") {
            REQUIRE(context.get_n_max() == ALIGN_TO);
        }
    }

    GIVEN("a request for a context of size x where x % ALIGN_TO != 0") {
        Pricer::ddx_pricer_context context(129);

        THEN("we get a context of size y where y > x and y % ALIGN_TO == 0") {
            REQUIRE(context.get_n_max() > 129);
            REQUIRE((context.get_n_max() % ALIGN_TO) == 0);
        }
    }


    GIVEN("a context for simple price and first derivative computations (wrt. strike)") {
        Pricer::ddx_pricer_context context( ALIGN_TO );

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
            REQUIRE(context.get_ddx_price()   != nullptr);

            /*
            REQUIRE(context.get_to_structure()       == nullptr);
            REQUIRE(context.get_offsets()            == nullptr);
            REQUIRE(context.get_premiums()           == nullptr);
            REQUIRE(context.get_instrument_prices()  == nullptr);
            REQUIRE(context.get_instrument_pricesl() == nullptr);
            REQUIRE(context.get_instrument_pricesh() == nullptr);
            REQUIRE(context.get_x_()                 == nullptr);
            REQUIRE(context.get_xl_()                == nullptr);
            REQUIRE(context.get_xh_()                == nullptr);
            REQUIRE(context.get_d2dx2()              == nullptr);
             */
        }

        THEN("the pointers are ALIGN_TO-bytes aligned (needed for usage of SIMD instructions)") {
            REQUIRE((((unsigned long long int) context.get_s())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigma())       & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_t())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_tau())         & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_r())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaA())      & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaA2T2())   & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaAsqrtT()) & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_emrt())        & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_long_short())  & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_put_call())    & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d2dx2_prep())  & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d1())          & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d2())          & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_prices())      & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_x())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_ddx_price())   & ((unsigned long long int) ALIGN_TO-1)) == 0);
        }


    }

}


SCENARIO("Pricer::d2dx2_pricer_context implementation") {

    GIVEN("a request for a context of size < ALIGN_TO") {
        Pricer::d2dx2_pricer_context context(ALIGN_TO/4);

        THEN("we get a context of size equal to ALIGN_TO") {
            REQUIRE(context.get_n_max() == ALIGN_TO);
        }
    }

    GIVEN("a request for a context of size x where x % ALIGN_TO != 0") {
        Pricer::d2dx2_pricer_context context(129);

        THEN("we get a context of size y where y > x and y % ALIGN_TO == 0") {
            REQUIRE(context.get_n_max() > 129);
            REQUIRE((context.get_n_max() % ALIGN_TO) == 0);
        }
    }

    GIVEN("a context for simple price and first derivative computations (wrt. strike)") {
        Pricer::d2dx2_pricer_context context( ALIGN_TO );

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
            REQUIRE(context.get_ddx_price()   != nullptr);
            REQUIRE(context.get_d2dx2_price() != nullptr);

            /*
            REQUIRE(context.get_to_structure()       == nullptr);
            REQUIRE(context.get_offsets()            == nullptr);
            REQUIRE(context.get_premiums()           == nullptr);
            REQUIRE(context.get_instrument_prices()  == nullptr);
            REQUIRE(context.get_instrument_pricesl() == nullptr);
            REQUIRE(context.get_instrument_pricesh() == nullptr);
            REQUIRE(context.get_x_()                 == nullptr);
            REQUIRE(context.get_xl_()                == nullptr);
            REQUIRE(context.get_xh_()                == nullptr);
             */
        }

        THEN("the pointers are ALIGN_TO-bytes aligned (needed for usage of SIMD instructions)") {
            REQUIRE((((unsigned long long int) context.get_s())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigma())       & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_t())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_tau())         & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_r())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaA())      & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaA2T2())   & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaAsqrtT()) & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_emrt())        & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_long_short())  & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_put_call())    & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d2dx2_prep())  & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d1())          & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d2())          & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_prices())      & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_x())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_ddx_price())   & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d2dx2_price()) & ((unsigned long long int) ALIGN_TO-1)) == 0);
        }


    }

}



SCENARIO("Pricer::compute_prices_of_instruments_context implementation") {

    GIVEN("a request for a context of size < ALIGN_TO") {
        Pricer::compute_prices_of_instruments_context context(ALIGN_TO/4, ALIGN_TO/4);

        THEN("we get a context of size equal to ALIGN_TO") {
            REQUIRE(context.get_n_max() == ALIGN_TO);
        }

        THEN("we get a context of size equal to ALIGN_TO") {
            REQUIRE(context.get_m_max() == ALIGN_TO);
        }
    }

    GIVEN("a request for a context of size x where x % ALIGN_TO != 0") {
        Pricer::compute_prices_of_instruments_context context(129, 129);

        THEN("we get a context of size y where y > x and y % ALIGN_TO == 0") {
            REQUIRE(context.get_n_max() > 129);
            REQUIRE((context.get_n_max() % ALIGN_TO) == 0);
        }

        THEN("we get a context of size y where y > x and y % ALIGN_TO == 0") {
            REQUIRE(context.get_m_max() > 129);
            REQUIRE((context.get_m_max() % ALIGN_TO) == 0);
        }
    }

    GIVEN("a context for simple price and first derivative computations (wrt. strike)") {
        Pricer::compute_prices_of_instruments_context context( ALIGN_TO, ALIGN_TO);

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
            REQUIRE(context.get_to_structure()!= nullptr);
            REQUIRE(context.get_instrument_prices() != nullptr);

        }

        THEN("the pointers are ALIGN_TO-bytes aligned (needed for usage of SIMD instructions)") {
            REQUIRE((((unsigned long long int) context.get_s())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigma())       & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_t())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_tau())         & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_r())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaA())      & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaA2T2())   & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_sigmaAsqrtT()) & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_emrt())        & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_long_short())  & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_put_call())    & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d2dx2_prep())  & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d1())          & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_d2())          & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_prices())      & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_x())           & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_to_structure())      & ((unsigned long long int) ALIGN_TO-1)) == 0);
            REQUIRE((((unsigned long long int) context.get_instrument_prices()) & ((unsigned long long int) ALIGN_TO-1)) == 0);
        }


    }

}


SCENARIO("Compute prices of instruments") {

    GIVEN("An instance of Pricer::compute_prices_of_instruments_context with exactly one leg per instrument") {
        Pricer::compute_prices_of_instruments_context context(100,100);

        SET_EQUAL_TO(100, context.get_r(), 0.01)
        SET_EQUAL_TO(100, context.get_s(), 70.)
        SET_EQUAL_TO(100, context.get_t(), 1.2)
        SET_EQUAL_TO(100, context.get_tau(), 1. / 12.)
        SET_EQUAL_TO(100, context.get_sigma(), 0.3)
        SET_EQUAL_TO(100, context.get_x(), 72.)

        SET_EQUAL_TO(100, context.get_sigmaA(), 0.)
        SET_EQUAL_TO(100, context.get_sigmaA2T2(), 0.)
        SET_EQUAL_TO(100, context.get_sigmaAsqrtT(), 0.)
        SET_EQUAL_TO(100, context.get_emrt(), 0.)
        SET_EQUAL_TO(100, context.get_d2dx2_prep(), 0.)

        SET_EQUAL_TO(100, context.get_d1(), 0.)
        SET_EQUAL_TO(100, context.get_d2(), 0.)

        SET_EQUAL_TO(100, context.get_prices(), 0.)

        for(uint64_t i = 0; i < 100; ++i) {
            context.get_to_structure()[i] = i;
        }

        init_tw_pricer();

        WHEN("we compute the instrument prices, the prices have to be equal to the prices of individal legs") {
            SET_EQUAL_TO(100, context.get_put_call(), 1)
            SET_EQUAL_TO(100, context.get_long_short(), 1)

            prepare_tw_pricer(context);
            tw_pricer(context);
            double tmp_prices[100];
            for(uint64_t i = 0; i < 100; ++i) {
                tmp_prices[i] = context.get_prices()[i];
            }
            compute_tw_prices_of_instruments(context);

            for(uint64_t i = 0; i < 100; ++i) {
                REQUIRE(tmp_prices[i] == context.get_instrument_prices()[i]);
            }
        }

    }
}

SCENARIO("Pricer computations are correct") {

    GIVEN("An instance of a Pricer::pricer_context class and a correspondingly initialized Vortex::Pricer object") {

        Pricer::pricer_context context( ALIGN_TO );

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

            THEN("all computed prices of tw_pricer() are equal") {
                for(uint64_t i = 0; i < ALIGN_TO; ++i) {
                    REQUIRE(abs(reference_pricer_value - context.get_prices()[i]) < 1.0e-8);
                }
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

            THEN("all computed prices of tw_pricer() are equal") {
                for(uint64_t i = 0; i < ALIGN_TO; ++i) {
                    REQUIRE(abs(reference_pricer_value + context.get_prices()[0]) < 1.0e-8);
                }
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
            THEN("all computed prices of tw_pricer() are equal") {
                for(uint64_t i = 0; i < ALIGN_TO; ++i) {
                    REQUIRE(abs(reference_pricer_value - context.get_prices()[0]) < 1.0e-8);
                }
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
            THEN("all computed prices of tw_pricer() are equal") {
                for(uint64_t i = 0; i < ALIGN_TO; ++i) {
                    REQUIRE(abs(reference_pricer_value + context.get_prices()[0]) < 1.0e-8);
                }
            }
        }

        WHEN("we price a put and a call at the money") {
            SET_EQUAL_TO(ALIGN_TO, context.get_x(), context.get_s()[0])
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)

            tw_pricer(context);
            double long_call_price = context.get_prices()[0];
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
            tw_pricer(context);
            double long_put_price = context.get_prices()[0];

            THEN("the put-call-parity holds") {
                REQUIRE(abs(long_call_price - long_put_price) < 1.0e-4);
            }

        }

        WHEN("we price a put and a call 5$ above the money") {
            SET_EQUAL_TO(ALIGN_TO, context.get_x(), context.get_s()[0]+5)
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)

            tw_pricer(context);
            double long_call_price = context.get_prices()[0];
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
            tw_pricer(context);
            double long_put_price = context.get_prices()[0];

            INFO("long_call_price = " << long_call_price);
            INFO("long_put_price  = " << long_put_price);
            INFO("long_call_price - long_put_price = " << long_call_price - long_put_price);
            INFO("5*exp(-context.get_r()[0]*context.get_t()[0]) = " << 5*exp(-context.get_r()[0]*context.get_t()[0]));

            THEN("the put-call-parity holds") {
                REQUIRE(abs(long_call_price - long_put_price+5*exp(-context.get_r()[0]*context.get_t()[0])) < 1.0e-4);
            }

        }

    }

    GIVEN("An instance of the Pricer::pricer_context class of size 129") {

        Pricer::pricer_context context( 129 );

        SET_EQUAL_TO(129, context.get_r(), 0.05)
        SET_EQUAL_TO(129, context.get_s(), 81.)
        SET_EQUAL_TO(129, context.get_t(), 1.2)
        SET_EQUAL_TO(129, context.get_tau(), 1. / 12.)
        SET_EQUAL_TO(129, context.get_sigma(), 0.2)
        SET_EQUAL_TO(129, context.get_x(), 72.)

        SET_EQUAL_TO(129, context.get_sigmaA(), 0.)
        SET_EQUAL_TO(129, context.get_sigmaA2T2(), 0.)
        SET_EQUAL_TO(129, context.get_sigmaAsqrtT(), 0.)
        SET_EQUAL_TO(129, context.get_emrt(), 0.)
        SET_EQUAL_TO(129, context.get_d2dx2_prep(), 0.)

        SET_EQUAL_TO(129, context.get_d1(), 0.)
        SET_EQUAL_TO(129, context.get_d2(), 0.)

        SET_EQUAL_TO(129, context.get_prices(), 0.)


        Vortex::Pricer p;


        p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                          context.get_r()[0], context.get_s()[0]);

        init_tw_pricer();

        prepare_tw_pricer(context);


        WHEN("we price a long call") {
            SET_EQUAL_TO(129, context.get_put_call(), 1)
            SET_EQUAL_TO(129, context.get_long_short(), 1)
            // test long call price
            tw_pricer(context);
            FLOAT reference_pricer_value = p.compute_call_price(context.get_x()[0]);

            THEN("the computed prices of Vortex::Pricer and tw_pricer() are equal") {
                REQUIRE(abs(reference_pricer_value - context.get_prices()[0]) < 1.0e-8);
            }

            THEN("all computed prices of tw_pricer() are equal") {
                for(uint64_t i = 0; i < 129; ++i) {
                    REQUIRE(abs(reference_pricer_value - context.get_prices()[i]) < 1.0e-8);
                }
            }

        }

        WHEN("we price a short call") {
            SET_EQUAL_TO(129, context.get_put_call(), 1)
            SET_EQUAL_TO(129, context.get_long_short(), -1)

            // compute the short call price
            tw_pricer(context);
            FLOAT reference_pricer_value = p.compute_call_price(context.get_x()[0]);
            THEN("the computed prices of the Vortex::Pricer and tw_pricer() are equal") {
                REQUIRE(abs(reference_pricer_value + context.get_prices()[0]) < 1.0e-8);
            }

            THEN("all computed prices of tw_pricer() are equal") {
                for(uint64_t i = 0; i < 129; ++i) {
                    REQUIRE(abs(reference_pricer_value + context.get_prices()[i]) < 1.0e-8);
                }
            }

        }

        WHEN("we price a long put") {
            SET_EQUAL_TO(129, context.get_put_call(), -1)
            SET_EQUAL_TO(129, context.get_long_short(), 1)

            // test long call price
            tw_pricer(context);
            FLOAT reference_pricer_value = p.compute_put_price(context.get_x()[0]);

            THEN("the computed prices of the Vortex::Pricer and tw_pricer() are equal ") {
                REQUIRE(abs(reference_pricer_value - context.get_prices()[0]) < 1.0e-8);
            }
            THEN("all computed prices of tw_pricer() are equal") {
                for(uint64_t i = 0; i < 129; ++i) {
                    REQUIRE(abs(reference_pricer_value - context.get_prices()[i]) < 1.0e-8);
                }
            }

        }

        WHEN("we price a short put"){
            SET_EQUAL_TO(129, context.get_put_call(), -1)
            SET_EQUAL_TO(129, context.get_long_short(), -1)

            tw_pricer(context);
            FLOAT reference_pricer_value = p.compute_put_price(context.get_x()[0]);

            THEN("the computed prices of Vortex::Pricer and tw_pricer() are equal") {
                REQUIRE(abs(reference_pricer_value + context.get_prices()[0]) < 1.0e-8);

            }
            THEN("all computed prices of tw_pricer() are equal") {
                for(uint64_t i = 0; i < 129; ++i) {
                    REQUIRE(abs(reference_pricer_value + context.get_prices()[i]) < 1.0e-8);
                }
            }

        }

    }


}


SCENARIO("First derivative computations of pricing function with respect to strike are correct") {

    GIVEN("an instance of a Pricer::ddx_pricer_context class and a correspondingly initialized Vortex::Pricer object") {

        Pricer::ddx_pricer_context context(ALIGN_TO);
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
        WHEN("we compute the first derivative with respect to the strike of the price-function for a long call") {
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)
            // test long call price
            tw_pricer(context);

            prepare_tw_pricer(context);

            // test long call price
            tw_pricer(context);
            ddx_tw_pricer(context);
            FLOAT reference_pricer_value1 = p.compute_call_price(context.get_x()[0]);
            FLOAT reference_pricer_value2 = p.compute_call_price(context.get_x()[0] + eps);

            THEN("the computed derivative has to be equal to the quotient of differences of the original prices.") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[0]) <
                        1.0e-4);
            }

            THEN("the computed first derivatives are all equal.") {
                for(uint64_t i = 0; i < ALIGN_TO; ++i) {
                    REQUIRE(abs(
                            (reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[i]) <
                            1.0e-4);
                }
            }
        }

        WHEN("we compute the first derivative with respect to the strike of the price-function for a short call") {
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)
            // test long call price
            tw_pricer(context);

            prepare_tw_pricer(context);

            // test long call price
            tw_pricer(context);
            ddx_tw_pricer(context);
            FLOAT reference_pricer_value1 = p.compute_call_price(context.get_x()[0]);
            FLOAT reference_pricer_value2 = p.compute_call_price(context.get_x()[0] + eps);

            THEN("the computed derivative has to be equal to the quotient of differences of the original prices.") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[0])
                        < 1.0e-4);
            }
            THEN("the computed first derivatives are all equal.") {
                for(uint64_t i = 0; i < ALIGN_TO; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[i])
                            < 1.0e-4);
                }
            }
        }


        WHEN("we compute the first derivative with respect to the strike of the price-function for a long put") {
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)
            // test long call price
            tw_pricer(context);

            prepare_tw_pricer(context);

            // test long call price
            tw_pricer(context);
            ddx_tw_pricer(context);
            FLOAT reference_pricer_value1 = p.compute_put_price(context.get_x()[0]);
            FLOAT reference_pricer_value2 = p.compute_put_price(context.get_x()[0] + eps);

            THEN("the computed derivative has to be equal to the quotient of differences of the original prices.") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[0])
                    < 1.0e-4);
            }
            THEN("the computed first derivatives are all equal.") {
                for(uint64_t i = 0; i < ALIGN_TO; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[i])
                            < 1.0e-4);
                }
            }
        }

        WHEN("we compute the first derivative with respect to the strike of the price-function for a short put") {
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)
            // test long call price
            tw_pricer(context);

            prepare_tw_pricer(context);

            // test long call price
            tw_pricer(context);
            ddx_tw_pricer(context);
            FLOAT reference_pricer_value1 = p.compute_put_price(context.get_x()[0]);
            FLOAT reference_pricer_value2 = p.compute_put_price(context.get_x()[0] + eps);

            THEN("the computed derivative has to be equal to the quotient of differences of the original prices.") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[0])
                        < 1.0e-4);
            }
            THEN("the computed first derivatives are all equal.") {
                for (uint64_t i = 0; i < ALIGN_TO; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[i])
                            < 1.0e-4);
                }
            }
        }

    }

    GIVEN("an instance of a Pricer::ddx_pricer_context of size 129") {

        Pricer::ddx_pricer_context context(129);
        SET_EQUAL_TO(129, context.get_r(), 0.01)
        SET_EQUAL_TO(129, context.get_s(), 70.)
        SET_EQUAL_TO(129, context.get_t(), 1.2)
        SET_EQUAL_TO(129, context.get_tau(), 1. / 12.)
        SET_EQUAL_TO(129, context.get_sigma(), 0.3)
        SET_EQUAL_TO(129, context.get_x(), 72.)
        SET_EQUAL_TO(129, context.get_ddx_price(), 0.)

        SET_EQUAL_TO(129, context.get_sigmaA(), 0.)
        SET_EQUAL_TO(129, context.get_sigmaA2T2(), 0.)
        SET_EQUAL_TO(129, context.get_sigmaAsqrtT(), 0.)
        SET_EQUAL_TO(129, context.get_emrt(), 0.)
        SET_EQUAL_TO(129, context.get_d2dx2_prep(), 0.)

        SET_EQUAL_TO(129, context.get_d1(), 0.)
        SET_EQUAL_TO(129, context.get_d2(), 0.)

        SET_EQUAL_TO(129, context.get_prices(), 0.)
        SET_EQUAL_TO(129, context.get_put_call(), 1)
        SET_EQUAL_TO(129, context.get_long_short(), 1)

        const FLOAT eps = 1.0e-10;
        Vortex::Pricer p;


        p.set_market_data(context.get_sigma()[0], context.get_t()[0], context.get_tau()[0],
                          context.get_r()[0], context.get_s()[0]);

        init_tw_pricer();
        WHEN("we compute the first derivative with respect to the strike of the price-function for a long call") {
            SET_EQUAL_TO(129, context.get_put_call(), 1)
            SET_EQUAL_TO(129, context.get_long_short(), 1)
            // test long call price
            tw_pricer(context);

            prepare_tw_pricer(context);

            // test long call price
            tw_pricer(context);
            ddx_tw_pricer(context);
            FLOAT reference_pricer_value1 = p.compute_call_price(context.get_x()[0]);
            FLOAT reference_pricer_value2 = p.compute_call_price(context.get_x()[0] + eps);

            THEN("the computed derivative has to be equal to the quotient of differences of the original prices.") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[0]) <
                        1.0e-4);
            }

            THEN("the computed first derivatives are all equal.") {
                for(uint64_t i = 0; i < 129; ++i) {
                    REQUIRE(abs(
                            (reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[i]) <
                            1.0e-4);
                }
            }
        }

        WHEN("we compute the first derivative with respect to the strike of the price-function for a short call") {
            SET_EQUAL_TO(129, context.get_put_call(), 1)
            SET_EQUAL_TO(129, context.get_long_short(), -1)
            // test long call price
            tw_pricer(context);

            prepare_tw_pricer(context);

            // test long call price
            tw_pricer(context);
            ddx_tw_pricer(context);
            FLOAT reference_pricer_value1 = p.compute_call_price(context.get_x()[0]);
            FLOAT reference_pricer_value2 = p.compute_call_price(context.get_x()[0] + eps);

            THEN("the computed derivative has to be equal to the quotient of differences of the original prices.") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[0])
                        < 1.0e-4);
            }
            THEN("the computed first derivatives are all equal.") {
                for(uint64_t i = 0; i < 129; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[i])
                            < 1.0e-4);
                }
            }
        }


        WHEN("we compute the first derivative with respect to the strike of the price-function for a long put") {
            SET_EQUAL_TO(129, context.get_put_call(), -1)
            SET_EQUAL_TO(129, context.get_long_short(), 1)
            // test long call price
            tw_pricer(context);

            prepare_tw_pricer(context);

            // test long call price
            tw_pricer(context);
            ddx_tw_pricer(context);
            FLOAT reference_pricer_value1 = p.compute_put_price(context.get_x()[0]);
            FLOAT reference_pricer_value2 = p.compute_put_price(context.get_x()[0] + eps);

            THEN("the computed derivative has to be equal to the quotient of differences of the original prices.") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[0])
                        < 1.0e-4);
            }
            THEN("the computed first derivatives are all equal.") {
                for(uint64_t i = 0; i < 129; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_ddx_price()[i])
                            < 1.0e-4);
                }
            }
        }

        WHEN("we compute the first derivative with respect to the strike of the price-function for a short put") {
            SET_EQUAL_TO(129, context.get_put_call(), -1)
            SET_EQUAL_TO(129, context.get_long_short(), -1)
            // test long call price
            tw_pricer(context);

            prepare_tw_pricer(context);

            // test long call price
            tw_pricer(context);
            ddx_tw_pricer(context);
            FLOAT reference_pricer_value1 = p.compute_put_price(context.get_x()[0]);
            FLOAT reference_pricer_value2 = p.compute_put_price(context.get_x()[0] + eps);

            THEN("the computed derivative has to be equal to the quotient of differences of the original prices.") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[0])
                        < 1.0e-4);
            }
            THEN("the computed first derivatives are all equal.") {
                for (uint64_t i = 0; i < 129; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps + context.get_ddx_price()[i])
                            < 1.0e-4);
                }
            }
        }

    }

}


SCENARIO("Second derivative computations of pricing function with respect to strike are correct") {

    GIVEN("an initialized instance of a Pricer::d2dx2_pricer_context class") {

        Pricer::d2dx2_pricer_context context(ALIGN_TO);

        SET_EQUAL_TO(ALIGN_TO, context.get_r(), 0.01)
        SET_EQUAL_TO(ALIGN_TO, context.get_s(), 70.)
        SET_EQUAL_TO(ALIGN_TO, context.get_t(), 1.2)
        SET_EQUAL_TO(ALIGN_TO, context.get_tau(), 1. / 12.)
        SET_EQUAL_TO(ALIGN_TO, context.get_sigma(), 0.3)
        SET_EQUAL_TO(ALIGN_TO, context.get_x(), 72.)
        SET_EQUAL_TO(ALIGN_TO, context.get_ddx_price(), 0.)
        SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_price(), 0.)

        SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA(), 0.)
        SET_EQUAL_TO(ALIGN_TO, context.get_sigmaA2T2(), 0.)
        SET_EQUAL_TO(ALIGN_TO, context.get_sigmaAsqrtT(), 0.)
        SET_EQUAL_TO(ALIGN_TO, context.get_emrt(), 0.)
        SET_EQUAL_TO(ALIGN_TO, context.get_d2dx2_prep(), 0.)

        SET_EQUAL_TO(ALIGN_TO, context.get_d1(), 0.)
        SET_EQUAL_TO(ALIGN_TO, context.get_d2(), 0.)

        SET_EQUAL_TO(ALIGN_TO, context.get_prices(), 0.)

        const FLOAT eps = 1.0e-4;

        WHEN("we compute the second derivative, f''(x), of the price function with respect to the strike of a long call"
             "as well as the first derivatives f'(x) and f'(x+eps) for 0<eps<<1") {
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)

            prepare_tw_pricer(context);
            tw_pricer(context);
            d2dx2_tw_pricer(context);

            FLOAT reference_pricer_value1;
            FLOAT reference_pricer_value2;

            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value1 = context.get_ddx_price()[0];

            for (uint64_t i = 0; i < ALIGN_TO; ++i) context.get_x()[i] += eps;
            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value2 = context.get_ddx_price()[0];

            THEN("the value, f''(x), of the second derivative has to be equal to (f'(x+eps)-f'(x))/eps ") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[0])
                    < eps);
            }
            THEN("the computed first derivatives are all equal.") {
                for (uint64_t i = 0; i < ALIGN_TO; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[i])
                            < eps);
                }
            }

        }

        WHEN("we compute the second derivative, f''(x), of the price function with respect to the strike of a short call"
             "as well as the first derivatives f'(x) and f'(x+eps) for 0<eps<<1") {
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), 1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)

            prepare_tw_pricer(context);
            tw_pricer(context);
            d2dx2_tw_pricer(context);

            FLOAT reference_pricer_value1;
            FLOAT reference_pricer_value2;

            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value1 = context.get_ddx_price()[0];

            for (uint64_t i = 0; i < ALIGN_TO; ++i) context.get_x()[i] += eps;
            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value2 = context.get_ddx_price()[0];

            THEN("the value, f''(x), of the second derivative has to be equal to (f'(x+eps)-f'(x))/eps ") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[0])
                        < eps);
            }
            THEN("the computed first derivatives are all equal.") {
                for (uint64_t i = 0; i < ALIGN_TO; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[i])
                            < eps);
                }
            }

        }

        WHEN("we compute the second derivative, f''(x), of the price function with respect to the strike of a long put"
             "as well as the first derivatives f'(x) and f'(x+eps) for 0<eps<<1") {
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), 1)

            prepare_tw_pricer(context);
            tw_pricer(context);
            d2dx2_tw_pricer(context);

            FLOAT reference_pricer_value1;
            FLOAT reference_pricer_value2;

            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value1 = context.get_ddx_price()[0];

            for (uint64_t i = 0; i < ALIGN_TO; ++i) context.get_x()[i] += eps;
            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value2 = context.get_ddx_price()[0];

            THEN("the value, f''(x), of the second derivative has to be equal to (f'(x+eps)-f'(x))/eps ") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[0])
                        < eps);
            }
            THEN("the computed first derivatives are all equal.") {
                for (uint64_t i = 0; i < ALIGN_TO; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[i])
                            < eps);
                }
            }

        }

        WHEN("we compute the second derivative, f''(x), of the price function with respect to the strike of a short put"
             "as well as the first derivatives f'(x) and f'(x+eps) for 0<eps<<1") {
            SET_EQUAL_TO(ALIGN_TO, context.get_put_call(), -1)
            SET_EQUAL_TO(ALIGN_TO, context.get_long_short(), -1)

            prepare_tw_pricer(context);
            tw_pricer(context);
            d2dx2_tw_pricer(context);

            FLOAT reference_pricer_value1;
            FLOAT reference_pricer_value2;

            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value1 = context.get_ddx_price()[0];

            for (uint64_t i = 0; i < ALIGN_TO; ++i) context.get_x()[i] += eps;
            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value2 = context.get_ddx_price()[0];

            THEN("the value, f''(x), of the second derivative has to be equal to (f'(x+eps)-f'(x))/eps ") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[0])
                        < eps);
            }
            THEN("the computed first derivatives are all equal.") {
                for (uint64_t i = 0; i < ALIGN_TO; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[i])
                            < eps);
                }
            }

        }

    }

    GIVEN("an initialized instance of a Pricer::d2dx2_pricer_context of size 129") {

        Pricer::d2dx2_pricer_context context(129);

        SET_EQUAL_TO(129, context.get_r(), 0.01)
        SET_EQUAL_TO(129, context.get_s(), 70.)
        SET_EQUAL_TO(129, context.get_t(), 1.2)
        SET_EQUAL_TO(129, context.get_tau(), 1. / 12.)
        SET_EQUAL_TO(129, context.get_sigma(), 0.3)
        SET_EQUAL_TO(129, context.get_x(), 72.)
        SET_EQUAL_TO(129, context.get_ddx_price(), 0.)
        SET_EQUAL_TO(129, context.get_d2dx2_price(), 0.)

        SET_EQUAL_TO(129, context.get_sigmaA(), 0.)
        SET_EQUAL_TO(129, context.get_sigmaA2T2(), 0.)
        SET_EQUAL_TO(129, context.get_sigmaAsqrtT(), 0.)
        SET_EQUAL_TO(129, context.get_emrt(), 0.)
        SET_EQUAL_TO(129, context.get_d2dx2_prep(), 0.)

        SET_EQUAL_TO(129, context.get_d1(), 0.)
        SET_EQUAL_TO(129, context.get_d2(), 0.)

        SET_EQUAL_TO(129, context.get_prices(), 0.)

        const FLOAT eps = 1.0e-4;

        WHEN("we compute the second derivative, f''(x), of the price function with respect to the strike of a long call"
             "as well as the first derivatives f'(x) and f'(x+eps) for 0<eps<<1") {
            SET_EQUAL_TO(129, context.get_put_call(), 1)
            SET_EQUAL_TO(129, context.get_long_short(), 1)

            prepare_tw_pricer(context);
            tw_pricer(context);
            d2dx2_tw_pricer(context);

            FLOAT reference_pricer_value1;
            FLOAT reference_pricer_value2;

            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value1 = context.get_ddx_price()[0];

            for (uint64_t i = 0; i < 129; ++i) context.get_x()[i] += eps;
            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value2 = context.get_ddx_price()[0];

            THEN("the value, f''(x), of the second derivative has to be equal to (f'(x+eps)-f'(x))/eps ") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[0])
                        < eps);
            }
            THEN("the computed first derivatives are all equal.") {
                for (uint64_t i = 0; i < 129; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[i])
                            < eps);
                }
            }

        }

        WHEN("we compute the second derivative, f''(x), of the price function with respect to the strike of a short call"
             "as well as the first derivatives f'(x) and f'(x+eps) for 0<eps<<1") {
            SET_EQUAL_TO(129, context.get_put_call(), 1)
            SET_EQUAL_TO(129, context.get_long_short(), -1)

            prepare_tw_pricer(context);
            tw_pricer(context);
            d2dx2_tw_pricer(context);

            FLOAT reference_pricer_value1;
            FLOAT reference_pricer_value2;

            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value1 = context.get_ddx_price()[0];

            for (uint64_t i = 0; i < 129; ++i) context.get_x()[i] += eps;
            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value2 = context.get_ddx_price()[0];

            THEN("the value, f''(x), of the second derivative has to be equal to (f'(x+eps)-f'(x))/eps ") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[0])
                        < eps);
            }
            THEN("the computed first derivatives are all equal.") {
                for (uint64_t i = 0; i < 129; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[i])
                            < eps);
                }
            }

        }

        WHEN("we compute the second derivative, f''(x), of the price function with respect to the strike of a long put"
             "as well as the first derivatives f'(x) and f'(x+eps) for 0<eps<<1") {
            SET_EQUAL_TO(129, context.get_put_call(), -1)
            SET_EQUAL_TO(129, context.get_long_short(), 1)

            prepare_tw_pricer(context);
            tw_pricer(context);
            d2dx2_tw_pricer(context);

            FLOAT reference_pricer_value1;
            FLOAT reference_pricer_value2;

            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value1 = context.get_ddx_price()[0];

            for (uint64_t i = 0; i < 129; ++i) context.get_x()[i] += eps;
            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value2 = context.get_ddx_price()[0];

            THEN("the value, f''(x), of the second derivative has to be equal to (f'(x+eps)-f'(x))/eps ") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[0])
                        < eps);
            }
            THEN("the computed first derivatives are all equal.") {
                for (uint64_t i = 0; i < 129; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[i])
                            < eps);
                }
            }

        }

        WHEN("we compute the second derivative, f''(x), of the price function with respect to the strike of a short put"
             "as well as the first derivatives f'(x) and f'(x+eps) for 0<eps<<1") {
            SET_EQUAL_TO(129, context.get_put_call(), -1)
            SET_EQUAL_TO(129, context.get_long_short(), -1)

            prepare_tw_pricer(context);
            tw_pricer(context);
            d2dx2_tw_pricer(context);

            FLOAT reference_pricer_value1;
            FLOAT reference_pricer_value2;

            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value1 = context.get_ddx_price()[0];

            for (uint64_t i = 0; i < 129; ++i) context.get_x()[i] += eps;
            tw_pricer(context);
            ddx_tw_pricer(context);
            reference_pricer_value2 = context.get_ddx_price()[0];

            THEN("the value, f''(x), of the second derivative has to be equal to (f'(x+eps)-f'(x))/eps ") {
                REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[0])
                        < eps);
            }
            THEN("the computed first derivatives are all equal.") {
                for (uint64_t i = 0; i < 129; ++i) {
                    REQUIRE(abs((reference_pricer_value2 - reference_pricer_value1) / eps - context.get_d2dx2_price()[i])
                            < eps);
                }
            }

        }

    }

}


TEST_CASE("used architecture", "[pricer]") {
    WARN("THIS WARNING IS SUPPOSED TO APPEAR !!! " <<
    "It is here to let the user know what kind of instruction-set is being used. " <<
    "Used instruction set is: " << getUsedInstructionSet());
    CHECK(true);
}


SCENARIO("Computation of upper and lower bound works") {
    GIVEN("a compute_insrtument_form_premiums_context of length (ALIGN_TO, ALIGN_TO)") {

        Pricer::compute_instrument_strikes_from_premiums_context context(ALIGN_TO, ALIGN_TO);

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
        SET_EQUAL_TO(ALIGN_TO, context.get_instrument_prices(), 0)
        SET_EQUAL_TO(ALIGN_TO, context.get_instrument_pricesl(), 0)
        SET_EQUAL_TO(ALIGN_TO, context.get_instrument_pricesh(), 0)

        SET_EQUAL_TO(ALIGN_TO, context.get_xl_(), 1)
        SET_EQUAL_TO(ALIGN_TO, context.get_xh_(), 1000)
        SET_EQUAL_TO(ALIGN_TO, context.get_x_(), 500)


        init_tw_pricer();

        WHEN("we compute the strikes of long calls with premium equal to 5") {
            prepare_tw_pricer(context);

            for (int32_t i = 0; i < ALIGN_TO; ++i) {
                context.get_to_structure()[i] = i;
            }

            compute_tw_upper_and_lower_bounds(context);

            THEN("The upper and lower bounds are priced correct") {

            }


        }
    }
}

SCENARIO("Computations of strikes from premiums works") {

    GIVEN("a compute_insrtument_form_premiums_context of length (ALIGN_TO, ALIGN_TO)") {

        Pricer::compute_instrument_strikes_from_premiums_context context(ALIGN_TO, ALIGN_TO);

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
        SET_EQUAL_TO(ALIGN_TO, context.get_instrument_prices(), 0)
        SET_EQUAL_TO(ALIGN_TO, context.get_instrument_pricesl(), 0)
        SET_EQUAL_TO(ALIGN_TO, context.get_instrument_pricesh(), 0)

        SET_EQUAL_TO(ALIGN_TO, context.get_xl_(), 1)
        SET_EQUAL_TO(ALIGN_TO, context.get_xh_(), 1000)
        SET_EQUAL_TO(ALIGN_TO, context.get_x_(), 500)


        init_tw_pricer();

        WHEN("we compute the strikes of long calls with premium equal to 5") {
            prepare_tw_pricer(context);

            for (int32_t i = 0; i < ALIGN_TO; ++i) {
                context.get_to_structure()[i] = i;
            }

            compute_tw_strikes_from_premiums(context);

            THEN("the resulting call cost a premium of 5") {
                REQUIRE(abs(context.get_instrument_prices()[0] - context.get_premiums()[0]) < 1.0e-4);
            }

            THEN("the resulting calls evaluate to 5") {
                compute_tw_prices_of_instruments(context);
                REQUIRE(abs(context.get_instrument_prices()[0] - context.get_premiums()[0]) < 1.0e-4);
            }
        }
    }


}



/*

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

*/