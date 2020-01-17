//
// Created by peter on 1/17/20.
//

#ifndef PRICER_MARKET_DATA_H
#define PRICER_MARKET_DATA_H


#include <c++/9.2.0/cstdint>
#include "memory/allocate.h"

typedef double jd_t;

class MarketData {
    double * __restrict__ market_prices;
    double * __restrict__ volatilities;
    double * __restrict__ interest_rates;

    jd_t * __restrict__ trade_dates;
    jd_t * __restrict__ begin_periods;
    jd_t * __restrict__ end_periods;

    uint64_t * __restrict__ term_months;

    uint64_t number_of_entries;

    uint64_t number_of_alloceted_entries;

    void extend_allocation() {
        Pricer::Private::reallocate_memory(2*number_of_alloceted_entries, number_of_alloceted_entries, &market_prices);
        Pricer::Private::reallocate_memory(2*number_of_alloceted_entries, number_of_alloceted_entries, &volatilities);
        Pricer::Private::reallocate_memory(2*number_of_alloceted_entries, number_of_alloceted_entries, &interest_rates);

        Pricer::Private::reallocate_memory(2*number_of_alloceted_entries, number_of_alloceted_entries, &trade_dates);
        Pricer::Private::reallocate_memory(2*number_of_alloceted_entries, number_of_alloceted_entries, &begin_periods);
        Pricer::Private::reallocate_memory(2*number_of_alloceted_entries, number_of_alloceted_entries, &end_periods);

        Pricer::Private::reallocate_memory(2*number_of_alloceted_entries, number_of_alloceted_entries, &term_months);

        number_of_alloceted_entries *= 2;
    }

public:

    MarketData(uint64_t number_of_entries_p = 1024*16) {
        auto n = Pricer::Private::allocate_memory<double>(number_of_entries_p, &market_prices);
        Pricer::Private::allocate_memory<double>(number_of_entries_p, &volatilities);
        Pricer::Private::allocate_memory<double>(number_of_entries_p, &interest_rates);

        Pricer::Private::allocate_memory<jd_t>(number_of_entries_p, &trade_dates);
        Pricer::Private::allocate_memory<jd_t>(number_of_entries_p, &begin_periods);
        Pricer::Private::allocate_memory<jd_t>(number_of_entries_p, &end_periods);

        Pricer::Private::allocate_memory<uint64_t>(number_of_entries_p, &term_months);

        number_of_entries = 0;
        number_of_alloceted_entries = n;
    }

    ~MarketData() {
        Pricer::Private::deallocate_memory<double>(market_prices);
        Pricer::Private::deallocate_memory<double>(volatilities);
        Pricer::Private::deallocate_memory<double>(interest_rates);

        Pricer::Private::deallocate_memory<jd_t>(trade_dates);
        Pricer::Private::deallocate_memory<jd_t>(begin_periods);
        Pricer::Private::deallocate_memory<jd_t>(end_periods);

        Pricer::Private::deallocate_memory<uint64_t>(term_months);
    }

    void add_entry( double market_price_p, double volatility_p, double interest_rate_p,
            jd_t trade_date_p, jd_t begin_period_p, jd_t end_period_p, uint64_t term_month_p) {


        if( number_of_entries >= number_of_alloceted_entries) {
            extend_allocation();
        }

        market_prices[number_of_entries]  = market_price_p;
        volatilities[number_of_entries]   = volatility_p;
        interest_rates[number_of_entries] = interest_rate_p;
        trade_dates[number_of_entries]    = trade_date_p;
        begin_periods[number_of_entries]  = begin_period_p;
        end_periods[number_of_entries++]  = end_period_p;

    }

    double *getMarketPrices() const {
        return market_prices;
    }

    double *getVolatilities() const {
        return volatilities;
    }

    double *getInterestRates() const {
        return interest_rates;
    }


    jd_t *getTradeDates() const {
        return trade_dates;
    }

    jd_t *getBeginPeriods() const {
        return begin_periods;
    }

    jd_t *getEndPeriods() const {
        return end_periods;
    }

    uint64_t *getTermMonths() const {
        return term_months;
    }

    uint64_t getNumberOfEntries() const {
        return number_of_entries;
    }

};


#endif //PRICER_MARKET_DATA_H
