//
// Created by peter on 1/17/20.
//

#include <c++/9.2.0/fstream>
#include "io_csv.h"
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <erfa.h>

/**
 *  Returns true iff the all market-dada in m was read successful.
 *
 * @param m
 * @return
 */
bool checkIfAllDataWasReadSuccessful(const MarketData *m) {
    uint64_t n = m->getNumberOfEntries();

    for(uint64_t i = 0; i < n; i++ ) {
        if(std::isnan(m->getMarketPrices()[i])) return false;
        if(std::isnan(m->getVolatilities()[i])) return false;
        if(std::isnan(m->getInterestRates()[i])) return false;
        if(std::isnan(m->getTradeDates()[i])) return false;
        if(std::isnan(m->getBeginPeriods()[i])) return false;
        if(std::isnan(m->getEndPeriods()[i])) return false;
        if(std::isnan(m->getPremiums()[i])) return false;
    }

    return true;
}


double getDoubleValue(unsigned int idx, const std::vector<std::string> &v);
jd_t getJD(unsigned int idx, const std::vector<std::string> &v);

unsigned int getIndex(const std::string name, std::map<std::string, unsigned int> &column_map) {
    if(column_map.find(name) == column_map.end()) {
        throw std::runtime_error("No volatility column in input found.");
    }
    return column_map[name];
}

uint64_t getMonthsBetween( jd_t t1, jd_t t2) {
    double dd;
    int y, m, d;
    uint64_t res;

    eraJd2cal(t1, 0, &y, &m, &d, &dd);
    res = 12*y+m;

    eraJd2cal(t2, 0, &y, &m, &d, &dd);
    res -= (12*y+m);

    return res;


}

/**
 * 
 * @param file a csv-file with the following columns:
 *              double market_price  --- price @ which the market trades
 *              double volatility    --- volatility of market_price
 *              double interest_rate --- the interest_rate for the end of the tenor of market_price
 *              double premiums      --- the premiums to be spent
 *              
 *              std::string trade_date   --- date of this trade in YYYY-MM-DD format
 *              std::string begin_period --- begin of settlement-period in YYYY—MM-DD format (i.e., the begin of the 
 *                                              period for which market_price is traded)
 *              std::string end_period   --- end of settlement-period in YYYY-MM-DD format.
 *                                              
 * @param delimiter. The delimiter of the csv-file.
 * @return a pointer to a MarketData object containing the read data. The object can contain NANs in case some data was unreadible.
 */
MarketData *read_market_data_from_csv(std::istream &file, char delimiter) {
    MarketData *m = new MarketData;

    // get the column-names
    std::string line;
    
    std::getline(file, line);
    std::map<std::string, unsigned int> column_map;
    std::string e;

    {
        std::istringstream f(line);
        unsigned int i = 0;
        while(std::getline(f, e, delimiter)) {
            column_map[e] = i++;
        }
    }
    
    unsigned int market_price_idx  = getIndex("market_price", column_map);
    unsigned int volatility_idx    = getIndex( "volatility", column_map);
    unsigned int interest_rate_idx = getIndex( "interest_rate", column_map);
    unsigned int premiums_idx      = getIndex( "premium", column_map);

    unsigned int trade_date_idx    = getIndex( "trade_date", column_map);
    unsigned int begin_period_idx  = getIndex( "begin_period", column_map);
    unsigned int end_period_idx     = getIndex( "end_period", column_map);

    
    while(std::getline(file, line)) {
        std::istringstream f(line);
        std::vector<std::string> v;
        
        while(std::getline(f, e, delimiter)) {
            v.push_back(e);
        }
        
        
        double market_price  = getDoubleValue(market_price_idx, v);
        double volatility    = getDoubleValue(volatility_idx, v);
        double interest_rate = getDoubleValue(interest_rate_idx, v);
        double premium       = getDoubleValue(premiums_idx, v);

        jd_t trade_date      = getJD(trade_date_idx, v);
        jd_t begin_period    = getJD(begin_period_idx, v);
        jd_t end_period      = getJD(end_period_idx, v);

        uint64_t term_month  = getMonthsBetween(begin_period, trade_date);
        
        m->add_entry(market_price, volatility, interest_rate, premium, trade_date, begin_period, end_period, term_month);
    }

    return m;

}

double getDoubleValue(unsigned int idx, const std::vector<std::string> &v) {
    double value;
    {
        std::string s;
        if (v.size() <= idx) {
            s = "NAN";
        } else {
            s = v[idx];
        }

        const char *begin_string = s.c_str();
        char *end_string = nullptr;
        value = strtod(begin_string, &end_string);
        if (begin_string == end_string) {
            value = strtod("NAN", nullptr);
        }
    }
    return value;
}


jd_t getJD(unsigned int idx, const std::vector<std::string> &v) {
    double d1, d2;
    {
        std::string s;
        if (v.size() <= idx) {
            return std::strtod( "NAN", nullptr );
        } else {
            s = v[idx];
        }

        const char *begin_string = s.c_str();
        char *end_string = nullptr;
        int y = strtoul(begin_string, &end_string, 10);
        if (begin_string == end_string) {
            return strtod("NAN", nullptr);
        }

        end_string = nullptr;
        begin_string += 5;
        int m = strtoul(begin_string, &end_string, 10);
        if (begin_string == end_string) {
            return strtod("NAN", nullptr);
        }

        end_string = nullptr;
        begin_string += 3;
        int d = strtoul(begin_string, &end_string, 10);
        if (begin_string == end_string) {
            return strtod("NAN", nullptr);
        }

        eraCal2jd(y, m, d, &d1, &d2);
        return d1+d2;

    }

}