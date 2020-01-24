//
// Created by peter on 1/22/20.
//

#include <pricer-dispatcher.h>
#include <stdexcept>
#include <map>
#include "portfolio.h"

void Portfolio::do_compute(double xl, double xh) {
    std::vector<uint64_t> instrument_map;

    // generate the trades
    generator.generate_instruments(market_data, &compute_from_prices_context,
            &instruments_context, instrument_map);

    // set xl and xh to their values
    for( uint64_t i = 0; i < compute_from_prices_context.get_m_act(); ++i) {
        compute_from_prices_context.get_xl_()[i] = xl;
        compute_from_prices_context.get_xh_()[i] = xh;
    }


    // compute prices in instrument_context ....
    compute_tw_prices_of_instruments(instruments_context);

    // ... and modify premiums in compute_from_premiums accordingly

    for(uint64_t i = 0; i < instruments_context.get_m_act(); i++) {
        compute_from_prices_context.get_premiums()[instrument_map[i]]
            -= instruments_context.get_instrument_prices()[i];
    }



    compute_tw_strikes_from_premiums(compute_from_prices_context);



}