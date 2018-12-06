//
// Created by peter on 12/3/18.
//

#ifndef PRICER_MKL_PRICER_H
#define PRICER_MKL_PRICER_H

#define FLOAT double

#include <mkl.h>

#ifdef __cplusplus
#define restrict
extern "C" {
#endif

void init_mkl_pricer();

void prepare_mkl_pricer(
        MKL_INT64 n,
        FLOAT * restrict s,                 /// [in] stock price
        FLOAT * restrict sigma,             /// [in] vola
        FLOAT * restrict t,                 /// [in] time to maturity
        FLOAT * restrict tau,               /// [in] time of avg. period
        FLOAT * restrict r,                 /// [in] interest rate
        FLOAT * restrict tmp1,
        FLOAT * restrict tmp2,
        FLOAT * restrict tmp3,
        FLOAT * restrict tmp4,
        FLOAT * restrict tmp5,
        FLOAT * restrict sigmaA,           /// [out] adjusted vola
        FLOAT * restrict sigmaA2T2,        /// [out] sigmaA^2t/2
        FLOAT * restrict sigmaAsqrtT,      /// [out] sigmaA*sqrt(t)
        FLOAT * restrict emrt);


void mkl_pricer(
        MKL_INT64 n,
        MKL_INT64 * restrict flags,
        FLOAT * restrict s,                /// [in] stock price
        FLOAT * restrict x,                /// [in] strike
        FLOAT * restrict sigmaA2T2,        /// [in] sigmaA^2t/2
        FLOAT * restrict sigmaAsqrtT,      /// [in] sigmaA*sqrt(t)
        FLOAT * restrict emrt,
        FLOAT * restrict tmp1,
        FLOAT * restrict tmp2,
        FLOAT * restrict tmp3,
        FLOAT * restrict tmp4,
        FLOAT * restrict d1,               /// [out] d1
        FLOAT * restrict d2,
        FLOAT * restrict price);

void ddx_mkl_pricer(
        MKL_INT64 n,
        MKL_INT64 *restrict flags,
        FLOAT *restrict d2,
        FLOAT *restrict emrt,
        FLOAT *restrict ddx_price);

#ifdef __cplusplus
};
#endif

#endif //PRICER_MKL_PRICER_H
