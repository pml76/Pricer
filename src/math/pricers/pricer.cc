//
// Created by VortexUser on 3/26/19.
//

#include <memory/context.h>
#include "pricer.h"

namespace Pricer {

    void pricer::prepare_pricer(class pricer_context *context) {
        prepare_tw_pricer(
                context->get_n(),
                context->get_s(),
                context->get_sigma(),
                context->get_t(),
                context->get_tau(),
                context->get_r(),
                context->get_sigmaA(),
                context->get_sigmaA2T2(),
                context->get_sigmaAsqrtT(),
                context->get_emrt(),
                context->get_d2dx2_prep());
    }

    void pricer::compute_price(class pricer_context *context) {
        tw_pricer(
                context->get_n(),
                context->get_long_short(),
                context->get_put_call(),
                context->get_s(),
                context->get_x(),
                context->get_sigmaA2T2(),
                context->get_sigmaAsqrtT(),
                context->get_emrt(),
                context->get_d1(),
                context->get_d2dx2_prep(),
                context->get_price() );

    }

}