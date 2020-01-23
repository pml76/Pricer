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

    // compute prices in instrument_context ....
    compute_tw_prices_of_instruments(instruments_context);

    // ... and modify premiums in compute_from_premiums accordingly
    for(uint64_t i = 0; i < instruments_context.get_m_act(); i++) {
        compute_from_prices_context.get_premiums()[instrument_map[i]]
            += instruments_context.get_instrument_prices()[i];
    }

    // set xl and xh to their values
    for( uint64_t i = 0; i < compute_from_prices_context.get_m_act(); ++i) {
        compute_from_prices_context.get_xl_()[i] = xl;
        compute_from_prices_context.get_xh_()[i] = xh;
    }

    // now, we need to check if there is a root in the interval [xl, xh]....
    compute_tw_upper_and_lower_bounds(compute_from_prices_context);

    // ...and sort those entries out where no root exists.
    uint64_t no_unbreaketed;
    uint64_t no_breaketed;
    check_if_roots_are_breaketed( compute_from_prices_context, &no_breaketed, &no_unbreaketed);

    Pricer::compute_instrument_strikes_from_premiums_context c(compute_from_prices_context.get_n_act(), no_breaketed);
    std::map<uint64_t, uint64_t> id_map;
    for(uint64_t k = 0, j = 0, i = 0; i < compute_from_prices_context.get_n_act(); ++i) {
        if(compute_from_prices_context.get_instrument_pricesh()[i] * compute_from_prices_context.get_instrument_pricesl()[i] < 0.) {
            c.overwrite_leg(j++, &compute_from_prices_context, i);
            if( id_map.find(compute_from_prices_context.get_to_structure()[i]) == id_map.end()) {
                id_map[compute_from_prices_context.get_to_structure()[i]] = k++;
            }
        }
    }
    for(auto &i : id_map) {
        c.overwrite_structure(i.second, &compute_from_prices_context, i.first);
    }



}