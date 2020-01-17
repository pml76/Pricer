//
// Created by peter on 1/17/20.
//
#include <sstream>

#include <catch2/catch.hpp>
#include <io_csv.h>
#include <erfa.h>


SCENARIO("Import from CSV") {

    GIVEN("A good CSV ") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate,trade_date,begin_period,end_period" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31" << std::endl;
        std::istringstream f(s.str());

        THEN("the routine should not throw an exception") {
            MarketData *m;
            REQUIRE_NOTHROW(m = read_market_data_from_csv(f,','));
            delete m;
        }

        GIVEN("The corresponding MarketData object") {
            MarketData *m;
            m = read_market_data_from_csv(f,',');

            THEN("We have exactly one set of market data read") {
                REQUIRE(m->getNumberOfEntries() == 1);
            }

            THEN("The market_price equals 70.5") {
                REQUIRE((m->getMarketPrices()[0] - 70.5) < 1.0e-12);
            }

            THEN("The volatility equals 0.34") {
                REQUIRE((m->getVolatilities()[0] - 0.34) < 1.0e-12);
            }

            THEN("The interest_rate equals 0.01") {
                REQUIRE((m->getInterestRates()[0] - 0.01) < 1.0e-12);
            }

            THEN("The trade_date equals 2012-01-01") {
                double d1 = m->getTradeDates()[0];
                double dd;
                int y, m, d;
                eraJd2cal(d1, 0, &y, &m, &d, &dd);
                bool b = (y == 2012) && (m == 1) && (d ==1);
                REQUIRE( b == true );
            }

            THEN("The begin_period equals 2012-12-01") {
                double d1 = m->getBeginPeriods()[0];
                double dd;
                int y, m, d;
                eraJd2cal(d1, 0, &y, &m, &d, &dd);
                bool b = (y == 2012) && (m == 12) && (d ==1);
                REQUIRE( b == true );
            }

            THEN("The end_period equals 2012-12-31") {
                double d1 = m->getEndPeriods()[0];
                double dd;
                int y, m, d;
                eraJd2cal(d1, 0, &y, &m, &d, &dd);
                bool b = (y == 2012) && (m == 12) && (d ==31);
                REQUIRE( b == true );
            }

            THEN("The term_month equals 11") {
                auto t = m->getTermMonths()[0];
                REQUIRE( t == 11 );
            }

            REQUIRE(checkIfAllDataWasReadSuccessful(m) == true);

            delete m;
        }

    }

    GIVEN("A CSV without a 'market_price' - field") {
        std::stringstream s;
        s << "market_price1,volatility,interest_rate,trade_date,begin_period,end_period" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }

    GIVEN("A CSV without a 'volatility' - field") {
        std::stringstream s;
        s << "market_price,volatility1,interest_rate,trade_date,begin_period,end_period" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }

    GIVEN("A CSV without a 'interest_rate' - field") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate1,trade_date,begin_period,end_period" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }

    GIVEN("A CSV without a 'trade_date' - field") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate,trade_date1,begin_period,end_period" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }

    GIVEN("A CSV without a 'begin_period' - field") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate,trade_date,begin_period1,end_period" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }

    GIVEN("A CSV without a 'end_period' - field") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate,trade_date,begin_period,end_period1" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }

    GIVEN("A CSV with invalid numbers") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate,trade_date,begin_period,end_period" << std::endl;
        s << "x70.5,x0.34,v0.01,v2012-01-01,2012-v12-01,2012-12-v31" << std::endl;
        std::istringstream f(s.str());

        MarketData *m = read_market_data_from_csv(f, ',');

        THEN("the result is always NAN") {
            REQUIRE(std::isnan(m->getMarketPrices()[0])==true);
            REQUIRE(std::isnan(m->getVolatilities()[0])==true);
            REQUIRE(std::isnan(m->getInterestRates()[0])==true);
            REQUIRE(std::isnan(m->getTradeDates()[0])==true);
            REQUIRE(std::isnan(m->getBeginPeriods()[0])==true);
            REQUIRE(std::isnan(m->getEndPeriods()[0])==true);

            REQUIRE(checkIfAllDataWasReadSuccessful(m) == false);

        }

        delete m;
    }
}