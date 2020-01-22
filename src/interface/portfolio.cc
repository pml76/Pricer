//
// Created by peter on 1/22/20.
//

#include <pricer-dispatcher.h>
#include <stdexcept>
#include "portfolio.h"

void Portfolio::do_compute() {
    std::vector<uint64_t> instrument_map;

    // generate the trades
    generator.generate_instruments(market_data, &compute_from_prices_context,
            &instruments_context, instrument_map);

    // compute prices in instrument_context ....
    compute_tw_prices_of_instruments(instruments_context);

    // ... and modify premiums in compute_from_premiums accordingly
    for(uint64_t i = 0; i < instruments_context.get_m_act(); i++) {
        compute_from_prices_context.get_premiums()[instrument_map[i]]
            += instruments_context.get_instrument_prices()[i];
    }

    // now, we need to check if there is a root in the interval [xl, xh]....
    compute_tw_upper_and_lower_bounds(compute_from_prices_context);

    // ...and sort those entries out where no root exists.



}