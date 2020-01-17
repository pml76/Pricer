//
// Created by peter on 1/17/20.
//

#ifndef PRICER_IO_CSV_H
#define PRICER_IO_CSV_H

#include "market_data.h"


MarketData *read_market_data_from_csv(std::istream &file, char delimiter = ',');
bool checkIfAllDataWasReadSuccessful(const MarketData *m);

#endif //PRICER_IO_CSV_H
