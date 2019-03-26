//
// Created by VortexUser on 3/26/19.
//

#ifndef PRICER_PRICER_H
#define PRICER_PRICER_H

#include <pricer-dispatcher.h>

namespace Pricer {


    class pricer {

        enum state_t : uint64_t {
            initialized         = 0,
            prepared            = 1,
            price_computed      = 2
        };

    public:
        pricer( class pricer_context *context) : context(context){
            init_tw_pricer();
            state = initialized;
        }

        ~pricer();

        void compute_price( class pricer_context *context);


    };
}

#endif //PRICER_PRICER_H
