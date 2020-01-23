//
// Created by peter on 1/17/20.
//

#ifndef PRICER_IO_CSV_H
#define PRICER_IO_CSV_H

#include <iostream>
#include "market_data.h"
#include "context.h"
#include <erfa.h>

MarketData *read_market_data_from_csv(std::istream &file, char delimiter = ',');
bool checkIfAllDataWasReadSuccessful(const MarketData *m);

template<typename S>
static void output_header(S &o, const char delimiter) {
    o << "tradedate" << delimiter;
    o << "tradedate_year" << delimiter;
    o << "tradedate_month" << delimiter;
    o << "tradedate_day" << delimiter;
    o << "tradedate_jd" << delimiter;

    o << "structure_id" << delimiter;
    o << "premium_per_structure" << delimiter;
    o << "underlying" << delimiter;
    o << "volatility" << delimiter;
    o << "interest_rate" << delimiter;
    o << "long_short" << delimiter;
    o << "put_call" << delimiter;
    o << "price" << delimiter;
    o << "strike" << delimiter;

    o << "period_begin" << delimiter;
    o << "period_begin_year" << delimiter;
    o << "period_begin_month" << delimiter;
    o << "period_begin_day" << delimiter;
    o << "period_begin_jd" << delimiter;

    o << "period_end" << delimiter;
    o << "period_end_year" << delimiter;
    o << "period_end_month" << delimiter;
    o << "period_end_day" << delimiter;
    o << "period_end_jd";

    o << std::endl;
}

template<typename S>
void write_csv(Pricer::compute_instrument_strikes_from_premiums_context &context, S &o,
             const char delimiter = ',') {

    output_header(o, delimiter);

    for(uint64_t i = 0; i < context.get_n_act(); i++) {
        int y, m, d;
        double dd;

        eraJd2cal(context.get_tradedate()[i],0,&y,&m,&d,&dd);
        o << y <<"-" << ((m<10)?"0":"") << m << "-" << ((d<10)?"0":"") << d << delimiter;
        o << y << delimiter;
        o << m << delimiter;
        o << d << delimiter;
        o << context.get_tradedate()[i] << delimiter;

        o << context.get_to_structure()[i] << delimiter;
        o << context.get_premiums()[context.get_to_structure()[i]] << delimiter;
        o << context.get_s()[i] << delimiter;
        o << context.get_sigma()[i] << delimiter;
        o << context.get_r()[i] << delimiter;
        o << ((context.get_long_short()[i]==1.)?"Long":((context.get_long_short()[i]==-1.)?"Short":"N/A")) << delimiter;
        o << ((context.get_put_call()[i]==1.)?"Call":((context.get_put_call()[i]==-1.)?"Put":"N/A")) << delimiter;
        o << context.get_prices()[i] << delimiter;

        o << context.get_x()[i] + context.get_offsets()[i] << delimiter;

        eraJd2cal(context.get_tradedate()[i] + context.get_t()[i]-context.get_tau()[i]+1,0,&y,&m,&d,&dd);
        o << y <<"-" << ((m<10)?"0":"") << m << "-" << ((d<10)?"0":"") << d << delimiter;
        o << y << delimiter;
        o << m << delimiter;
        o << d << delimiter;
        o << context.get_tradedate()[i] + context.get_t()[i] - context.get_tau()[i] +1<< delimiter;

        eraJd2cal(context.get_tradedate()[i] + context.get_t()[i],0,&y,&m,&d,&dd);
        o << y <<"-" << ((m<10)?"0":"") << m << "-" << ((d<10)?"0":"") << d << delimiter;
        o << y << delimiter;
        o << m << delimiter;
        o << d << delimiter;
        o << context.get_tradedate()[i] + context.get_t()[i];

        o << std::endl;
    }

}


template<typename S>
void write_csv(Pricer::compute_prices_of_instruments_context &context, S &o,
        const char delimiter = ',') {

    output_header(o, delimiter);

    for(uint64_t i = 0; i < context.get_n_act(); i++) {
        int y, m, d;
        double dd;

        eraJd2cal(context.get_tradedate()[i],0,&y,&m,&d,&dd);
        o << y <<"-" << ((m<10)?"0":"") << m << "-" << ((d<10)?"0":"") << d << delimiter;
        o << y << delimiter;
        o << m << delimiter;
        o << d << delimiter;
        o << context.get_tradedate()[i] << delimiter;

        o << context.get_to_structure()[i] << delimiter;
        o << "NaN" << delimiter; // context.get_premiums()[context.get_to_structure()[i]] << delimiter;
        o << context.get_s()[i] << delimiter;
        o << context.get_sigma()[i] << delimiter;
        o << context.get_r()[i] << delimiter;
        o << ((context.get_long_short()[i]==1.)?"Long":((context.get_long_short()[i]==-1.)?"Short":"N/A")) << delimiter;
        o << ((context.get_put_call()[i]==1.)?"Call":((context.get_put_call()[i]==-1.)?"Put":"N/A")) << delimiter;
        o << context.get_prices()[i] << delimiter;

        o << context.get_x()[i] + context.get_offsets()[i] << delimiter;

        eraJd2cal(context.get_tradedate()[i] + context.get_t()[i]-context.get_tau()[i]+1,0,&y,&m,&d,&dd);
        o << y <<"-" << ((m<10)?"0":"") << m << "-" << ((d<10)?"0":"") << d << delimiter;
        o << y << delimiter;
        o << m << delimiter;
        o << d << delimiter;
        o << context.get_tradedate()[i] + context.get_t()[i] - context.get_tau()[i] +1<< delimiter;

        eraJd2cal(context.get_tradedate()[i] + context.get_t()[i],0,&y,&m,&d,&dd);
        o << y <<"-" << ((m<10)?"0":"") << m << "-" << ((d<10)?"0":"") << d << delimiter;
        o << y << delimiter;
        o << m << delimiter;
        o << d << delimiter;
        o << context.get_tradedate()[i] + context.get_t()[i];

        o << std::endl;
    }

}



#endif //PRICER_IO_CSV_H
