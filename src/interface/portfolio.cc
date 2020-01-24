//
// Created by peter on 1/22/20.
//

#include <pricer-dispatcher.h>
#include "portfolio.h"

void Portfolio::generate_portfolio(double xl, double xh) {

    ASSUME_ALIGNED(Real_Ptr ,compute_from_prices_context.get_xl_())
    ASSUME_ALIGNED(Real_Ptr ,compute_from_prices_context.get_xh_())
    ASSUME_ALIGNED(Real_Ptr ,compute_from_prices_context.get_premiums())
    ASSUME_ALIGNED(Real_Ptr ,compute_from_prices_context.get_instrument_prices())

    std::vector<uint64_t> instrument_map;

    // generate the trades
    generator.generate_instruments(market_data, &compute_from_prices_context,
            &instruments_context, instrument_map);

    // shrink n_max and m_max to the closest values
    // since we're not going to insert any more items in these
    // structures.
    compute_from_prices_context.shrink_m_max_to_actual();
    compute_from_prices_context.shrink_n_max_to_actual();
    instruments_context.shrink_m_max_to_actual();
    instruments_context.shrink_n_max_to_actual();

    // set xl and xh to their values
#ifdef NDEBUG
    #pragma omp parallel
#endif
    {

#ifdef NDEBUG
    #pragma omp for simd
#endif
        for (uint64_t i = 0; i < compute_from_prices_context.get_m_act(); ++i) {
            compute_from_prices_context.get_xl_()[i] = xl;
            compute_from_prices_context.get_xh_()[i] = xh;
        }

    }


    // compute prices in instrument_context ....
    compute_tw_prices_of_instruments(instruments_context);

    // ... and modify premiums in compute_from_premiums accordingly
#ifdef NDEBUG
#pragma omp parallel
#endif
    {

#ifdef NDEBUG
#pragma omp for simd
#endif
        for (uint64_t i = 0; i < instruments_context.get_m_act(); i++) {
            compute_from_prices_context.get_premiums()[instrument_map[i]]
                    -= instruments_context.get_instrument_prices()[i];
        }

    }


    // now we're ready to let the magic work....
    compute_tw_strikes_from_premiums(compute_from_prices_context);


    // restore the premiums, since it looks a little odd.
    //
    // strictly speaking, this is not necessary However, we do it anyway
    // since the user expects us to do so.
    //
#ifdef NDEBUG
#pragma omp parallel
#endif
    {

#ifdef NDEBUG
#pragma omp for simd
#endif
        for (uint64_t i = 0; i < instruments_context.get_m_act(); i++) {
            compute_from_prices_context.get_premiums()[instrument_map[i]]
                    += instruments_context.get_instrument_prices()[i];
        }

    }

}