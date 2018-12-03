//
// Created by peter on 12/3/18.
//

#ifndef PRICER_MKL_PRICER_H
#define PRICER_MKL_PRICER_H

#define FLOAT double

void init_mkl_pricer();

void prepare_mkl_pricer(
        MKL_INT64 n,
        FLOAT *__attribute__((aligned(32))) restrict s,                 /// [in] stock price
        FLOAT *__attribute__((aligned(32))) restrict x,                 /// [in] strike
        FLOAT *__attribute__((aligned(32))) restrict sigma,             /// [in] vola
        FLOAT *__attribute__((aligned(32))) restrict t,                 /// [in] time to maturity
        FLOAT *__attribute__((aligned(32))) restrict tau,               /// [in] time of avg. period
        FLOAT *__attribute__((aligned(32))) restrict r,                 /// [in] interest rate
        FLOAT *__attribute__((aligned(32))) restrict tmp1,
        FLOAT *__attribute__((aligned(32))) restrict tmp2,
        FLOAT *__attribute__((aligned(32))) restrict tmp3,
        FLOAT *__attribute__((aligned(32))) restrict tmp4,
        FLOAT *__attribute__((aligned(32))) restrict tmp5,
        FLOAT *__attribute__((aligned(32))) restrict sigmaA,           /// [out] adjusted vola
        FLOAT *__attribute__((aligned(32))) restrict sigmaA2T2,        /// [out] sigmaA^2t/2
        FLOAT *__attribute__((aligned(32))) restrict sigmaAsqrtT,      /// [out] sigmaA*sqrt(t)
        FLOAT *__attribute__((aligned(32))) restrict emrt);


void mkl_pricer(
        MKL_INT64 n,
        MKL_INT64 *__attribute__((aligned(32))) restrict flags,
        FLOAT *__attribute__((aligned(32))) restrict s,                /// [in] stock price
        FLOAT *__attribute__((aligned(32))) restrict x,                /// [in] strike
        FLOAT *__attribute__((aligned(32))) restrict sigmaA2T2,        /// [in] sigmaA^2t/2
        FLOAT *__attribute__((aligned(32))) restrict sigmaAsqrtT,      /// [in] sigmaA*sqrt(t)
        FLOAT *__attribute__((aligned(32))) restrict emrt,
        FLOAT *__attribute__((aligned(32))) restrict tmp1,
        FLOAT *__attribute__((aligned(32))) restrict tmp2,
        FLOAT *__attribute__((aligned(32))) restrict d1,               /// [out] d1
        FLOAT *__attribute__((aligned(32))) restrict d2,
        FLOAT *__attribute__((aligned(32))) restrict price);

#endif //PRICER_MKL_PRICER_H
