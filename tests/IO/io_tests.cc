//
// Created by peter on 1/17/20.
//
#include <sstream>

#include <catch2/catch.hpp>
#include <io_csv.h>
#include <erfa.h>


SCENARIO("Export to CSV") {
    GIVEN("Some Trades") {

        Pricer::compute_instrument_strikes_from_premiums_context context(1,1);
        uint64_t structure = context.add_structure(72, 5);
        double d1, d2, b1,b2,e1,e2;
        eraCal2jd(2012,1,1,&d1,&d2);
        eraCal2jd(2012,6,1,&b1,&b2);
        eraCal2jd(2012,6,30,&e1,&e2);
        context.add_leg(d1+d2, 60, 70, 0.03, 0.01, (e1+e2)-(d1+d2), (e1+e2)-(b1+b2)+1, -1., -1., 0., structure);
        context.add_leg(d1+d2, 60, 70, 0.03, 0.01, (e1+e2)-(d1+d2), (e1+e2)-(b1+b2)+1, 1., -1., 0., structure);
        std::stringstream s;


        const char *const header = "tradedate,tradedate_year,tradedate_month,tradedate_day,tradedate_jd,structure_id,premium_per_structure,underlying,volatility,interest_rate,long_short,put_call,price,strike,period_begin,period_begin_year,period_begin_month,period_begin_day,period_begin_jd,period_end,period_end_year,period_end_month,period_end_day,period_end_jd";

        THEN("the output is correct when context is a Pricer::compute_instrument_strikes_from_premiums_context") {
            write_csv(context, s);

            std::string line;
            std::getline(s,line);
            REQUIRE(line == header);

            std::getline(s, line);
            const char *const put_line = "2012-01-01,2012,1,1,2.45593e+06,1,5,70,0.03,0.01,Short,Put,0,60,2012-06-01,2012,6,1,2.45608e+06,2012-06-30,2012,6,30,2.45611e+06";
            REQUIRE(line == put_line);

            std::getline(s, line);
            const char *const call_line = "2012-01-01,2012,1,1,2.45593e+06,1,5,70,0.03,0.01,Short,Call,0,60,2012-06-01,2012,6,1,2.45608e+06,2012-06-30,2012,6,30,2.45611e+06";
            REQUIRE(line == call_line);

        }

        THEN("the output is correct when the context is a Pricer::compute_prices_of_instruments_context ") {
            write_csv(static_cast<Pricer::compute_prices_of_instruments_context &>(context), s);

            std::string line;
            std::getline(s,line);
            REQUIRE(line == header);

            std::getline(s, line);
            const char *const put_line = "2012-01-01,2012,1,1,2.45593e+06,1,NaN,70,0.03,0.01,Short,Put,0,60,2012-06-01,2012,6,1,2.45608e+06,2012-06-30,2012,6,30,2.45611e+06";
            REQUIRE(line == put_line);

            std::getline(s, line);
            const char *const call_line = "2012-01-01,2012,1,1,2.45593e+06,1,NaN,70,0.03,0.01,Short,Call,0,60,2012-06-01,2012,6,1,2.45608e+06,2012-06-30,2012,6,30,2.45611e+06";
            REQUIRE(line == call_line);

        }


    }
}

SCENARIO("Import from CSV") {

    GIVEN("A good CSV ") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate,trade_date,begin_period,end_period,premium" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31,3.14" << std::endl;
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

            THEN("The premium equals 3.14") {
                REQUIRE((m->getPremiums()[0] - 3.14) < 1.0e-12);
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
        s << "market_price1,volatility,interest_rate,trade_date,begin_period,end_period,premium" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31,3.14" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }

    GIVEN("A CSV without a 'volatility' - field") {
        std::stringstream s;
        s << "market_price,volatility1,interest_rate,trade_date,begin_period,end_period,premium" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31,3.14" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }

    GIVEN("A CSV without a 'interest_rate' - field") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate1,trade_date,begin_period,end_period,premium" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31,3.14" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }

    GIVEN("A CSV without a 'trade_date' - field") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate,trade_date1,begin_period,end_period,premium" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31,3.14" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }

    GIVEN("A CSV without a 'begin_period' - field") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate,trade_date,begin_period1,end_period,premium" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31,3.14" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }

    GIVEN("A CSV without a 'end_period' - field") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate,trade_date,begin_period,end_period1,premium" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31,3.14" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }

    GIVEN("A CSV without a 'premium' - field") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate,trade_date,begin_period,end_period,premium1" << std::endl;
        s << "70.5,0.34,0.01,2012-01-01,2012-12-1,2012-12-31,3.14" << std::endl;
        std::istringstream f(s.str());

        THEN("throw runtime-error") {
            REQUIRE_THROWS_AS(read_market_data_from_csv(f,','), std::runtime_error);
        }
    }


    GIVEN("A CSV with invalid numbers") {
        std::stringstream s;
        s << "market_price,volatility,interest_rate,trade_date,begin_period,end_period,premium" << std::endl;
        s << "x70.5,x0.34,v0.01,v2012-01-01,2012-v12-01,2012-12-v31,v3.14" << std::endl;
        std::istringstream f(s.str());

        MarketData *m = read_market_data_from_csv(f, ',');

        THEN("the result is always NAN") {
            REQUIRE(std::isnan(m->getMarketPrices()[0])==true);
            REQUIRE(std::isnan(m->getVolatilities()[0])==true);
            REQUIRE(std::isnan(m->getInterestRates()[0])==true);
            REQUIRE(std::isnan(m->getTradeDates()[0])==true);
            REQUIRE(std::isnan(m->getBeginPeriods()[0])==true);
            REQUIRE(std::isnan(m->getEndPeriods()[0])==true);
            REQUIRE(std::isnan(m->getPremiums()[0])==true);

            REQUIRE(checkIfAllDataWasReadSuccessful(m) == false);

        }

        delete m;
    }
}