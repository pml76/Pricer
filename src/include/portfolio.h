//
// Created by peter on 1/22/20.
//

#ifndef PRICER_PORTFOLIO_H
#define PRICER_PORTFOLIO_H

#include <vector>
#include "context.h"
#include "market_data.h"

class PortfolioGeneratorInterface {
public:

    virtual void generate_instruments(
            const MarketData &market_data,
            Pricer::compute_instrument_strikes_from_premiums_context *c1,
            Pricer::compute_prices_of_instruments_context *c2,
            std::vector<uint64_t> &instrument_map
            ) const = 0;
};

class Portfolio {
    Pricer::compute_instrument_strikes_from_premiums_context compute_from_prices_context;
    Pricer::compute_prices_of_instruments_context instruments_context;

    const MarketData &market_data;
    const PortfolioGeneratorInterface &generator;

public:
    Portfolio(const MarketData &market_data_p, const PortfolioGeneratorInterface &generator_p) :
        market_data(market_data_p), generator(generator_p) {}

    void do_compute(double xl, double xh);
};

#endif //PRICER_PORTFOLIO_H
