


#pragma GCC optimize("O3", "unroll-loops", "no-omit-frame-pointer", "inline") //Optimization flags
#pragma GCC option("arch=native", "tune=native", "no-zero-upper") //Enable AVX

#include <mkl.h>
#include <math.h>
#include <src/math/mkl_pricer.h>


static FLOAT ln_of_2, msqrt2;


static inline void __attribute__((always_inline)) compute_d_values(
        MKL_INT64 n,
        FLOAT *__attribute__((aligned(32))) restrict s,                /// [in] stock price
        FLOAT *__attribute__((aligned(32))) restrict x,                /// [in] strike
        FLOAT *__attribute__((aligned(32))) restrict sigmaA2T2,        /// [in] sigmaA^2t/2
        FLOAT *__attribute__((aligned(32))) restrict sigmaAsqrtT,      /// [in] sigmaA*sqrt(t)
        FLOAT *__attribute__((aligned(32))) restrict tmp1,
        FLOAT *__attribute__((aligned(32))) restrict tmp2,
        FLOAT *__attribute__((aligned(32))) restrict d1,               /// [out] d1
        FLOAT *__attribute__((aligned(32))) restrict d2) {            /// [out] d2


    for (MKL_INT64 i = 0; i < n; ++i) {
        tmp1[i] = s[i] / x[i];
    }

    vdLn(n, tmp1, tmp2);

    for (MKL_INT64 i = 0; i < n; ++i) {
        d1[i] = (tmp2[i] + sigmaA2T2[i]) / (sigmaAsqrtT[i] * msqrt2);
        d2[i] = d1[i] - sigmaAsqrtT[i] / msqrt2;
    }

}


void init_mkl_pricer() {
    ln_of_2 = log(2.);
    msqrt2 = -sqrt(2);
}

/**
 * @brief prepares the pricing routine.
 *
 * Has to be called before the pricer is called. All pointers must point to arrays of size n.
 *
 * @param n         number of items in lists.
 * @param s         stock price
 * @param x         strike
 * @param sigma     vola
 * @param t         time to maturity
 * @param tau       time of averaging period
 * @param r         interest rate
 *
 *
 *
 */
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
        FLOAT *__attribute__((aligned(32))) restrict emrt) {          /// [out] exp(-rt)/2

    FLOAT tt1, tt2;
    for (MKL_INT64 i = 0; i < n; ++i) {
        tt1 = sigma[i] * sigma[i];
        tt2 = tt1 + ln_of_2;
        tmp1[i] = tt2 * t[i];
        tmp2[i] = tt2 * (t[i] - tau[i]);
        tmp5[i] = tt1 * tt1 * tau[i] * tau[i];
    }

    vdExp(n, tmp1, tmp3);
    vdExp(n, tmp2, tmp4);

    FLOAT tt3;
    for (MKL_INT64 i = 0; i < n; ++i) {
        tt3 = tmp4[i] * (1. - sigma[i] * sigma[i] * tau[i]);
        tmp2[i] = (tmp3[i] - tt3) / tmp5[i];
    }

    vdLn(n, tmp2, tmp1);
    for (MKL_INT64 i = 0; i < n; ++i) {
        tmp1[i] /= t[i];
    }

    vdSqrt(n, tmp1, sigmaA);
    vdSqrt(n, t, tmp1);


    for (MKL_INT64 i = 0; i < n; ++i) {
        sigmaA2T2[i] = sigmaA[i] * sigmaA[i] * t[i] / 2.;
        sigmaAsqrtT[i] = sigmaA[i] * tmp1[i];
        tmp2[i] = -r[i] * t[i] - ln_of_2;
    }

    vdExp(n, tmp2, emrt);

}

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
        FLOAT *__attribute__((aligned(32))) restrict price) {            /// [out] d2

    compute_d_values(n, s, x, sigmaA2T2,
                     sigmaAsqrtT,
                     tmp1,
                     tmp2,
                     d1,
                     d2);

    for (MKL_INT64 i = 0; i < n; ++i) {
        if ((flags[i] & 1) != 0) {
            d1[i] = -d1[i];
            d2[i] = -d2[i];
        }
    }


    vdErfc(n, d1, tmp1);
    vdErfc(n, d2, tmp2);

    double sx1, sx2;
    for (MKL_INT64 i = 0; i < n; ++i) {
        if ((flags[i] & 1) != 0) {
            sx1 = x[i];
            sx2 = s[i];
        } else {
            sx1 = s[i];
            sx2 = x[i];
        }
        price[i] = sx1 * tmp1[i] - sx2 * tmp2[i];
        price[i] *= emrt[i];

        if ((flags[i] & 2) != 0) {
            price[i] = -price[i];
        }
    }

}