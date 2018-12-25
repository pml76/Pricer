//
// Created by peter on 12/21/18.
//

#ifndef PRICER_SLEEF_PRICER_H
#define PRICER_SLEEF_PRICER_H

#include <src/math/pricers/pricer-base.h>

#ifdef __cplusplus
#define restrict
extern "C" {
#endif

void init_sleef_pricer();

void prepare_sleef_pricer(
        UINT64 n,
        Real_Ptr s_,                 /// [in] future price
        Real_Ptr sigma_,             /// [in] vola
        Real_Ptr t_,                 /// [in] time to maturity
        Real_Ptr tau_,               /// [in] time of avg. period
        Real_Ptr r_,                 /// [in] interest rate
        Real_Ptr sigmaA_,           /// [out] adjusted vola
        Real_Ptr sigmaA2T2_,        /// [out] sigmaA^2t/2
        Real_Ptr sigmaAsqrtT_,      /// [out] sigmaA*sqrt(t)
        Real_Ptr emrt_,             /// [out] exp(-rt)/2
        Real_Ptr d2dx2_prep_);


void sleef_pricer(
        UINT64 n,
        Real_Ptr long_short_,     // 1 == long option // -1 == short option
        Real_Ptr put_call_,       // -1 == put // 1 == call
        Real_Ptr s_,                /// [in] stock price
        Real_Ptr x_,                /// [in] strike
        Real_Ptr sigmaA2T2_,        /// [in] sigmaA^2t/2
        Real_Ptr sigmaAsqrtT_,      /// [in] sigmaA*sqrt(t)
        Real_Ptr emrt_,
        Real_Ptr d1_,               /// [out] d1
        Real_Ptr d2_,
        Real_Ptr price_);

#ifdef __cplusplus
};
#endif

#endif //PRICER_SLEEF_PRICER_H
