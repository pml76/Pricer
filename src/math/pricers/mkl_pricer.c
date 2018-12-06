


// #pragma GCC optimize("O3", "unroll-loops", "no-omit-frame-pointer", "inline") //Optimization flags
// #pragma GCC option("arch=native", "tune=native", "no-zero-upper", "funsafe-math") //Enable AVX

#include <mkl.h>
#include <math.h>
#include <src/math/pricers/mkl_pricer.h>
#include <stdio.h>


static FLOAT ln_of_2, msqrt2;


// __attribute__((aligned(32)))

/**
 * @brief computes the d-vealues of the black scholes formula.
 *
 * On return, the values stored in d1 and d2 are equal to
 *
 *    d1 = (log(s/x)+sigmaA2T2) / (sigmaAsqrtT * (-sqrt(2)))
 *    d2 = d1 - sigmaAsqrtT / -sqrt(2)
 *
 * @note All arrays given to this function are assumed to be of length n
 *
 * @param n [in] Size of the arrays
 * @param s [in] Array of stock prices
 * @param x [in] Array of strikes
 * @param sigmaA2T2 [in] Array computed by prepare_mkl_pricer()
 * @param sigmaAsqrtT [in[ Array computed by prepare_mkl_pricer()
 * @param tmp1
 * @param tmp2
 * @param d1 [out]
 * @param d2 [out]
 */
static inline void __attribute__((always_inline)) compute_d_values(
        MKL_INT64 n,
        FLOAT *restrict s,                /// [in] stock price
        FLOAT *restrict x,                /// [in] strike
        FLOAT *restrict sigmaA2T2,        /// [in] sigmaA^2t/2
        FLOAT *restrict sigmaAsqrtT,      /// [in] sigmaA*sqrt(t)
        FLOAT *restrict tmp1,
        FLOAT *restrict tmp2,
        FLOAT *restrict d1,               /// [out] d1
        FLOAT *restrict d2) {            /// [out] d2


    for (MKL_INT64 i = 0; i < n; ++i) {
        tmp1[i] = s[i] / x[i];
    }

    vdLn(n, tmp1, tmp2);

    for (MKL_INT64 i = 0; i < n; ++i) {
        d1[i] = (tmp2[i] + sigmaA2T2[i]) / (sigmaAsqrtT[i] * msqrt2);
        d2[i] = d1[i] - sigmaAsqrtT[i] / msqrt2;
    }

}

/**
 * @brief Set some constants before computing.
 *
 * @note: This function has to be called before any of the other function of this file can be called
 */
void init_mkl_pricer() {
    ln_of_2 = log(2.);
    msqrt2 = -sqrt(2);
}

/**
 * @brief Prepare the computation of the option prices.
 *
 * This routine computes some values that are used during the computation of the of the option prices that are not
 * dependent on the strike. Upon return the following will be true for each entry in the arrays
 *
 * sigmaA      = sqrt(log(M)/t)  where M = ( 2exp(sigma^2t)-2exp(sigma^2(t-tau))(1+sigma^2 tau) ) / (sigma^4 tau^2)
 * sigmaA2T2   = sigmaA^2*t/2
 * sigmaAsqrtT = sigmaA*sqrt(t)
 * emrt        = exp(-r*t)/2
 *
 * @note: All arrays given to this function are assumed to be of length n
 *
 * @param n [in] Size of the arrays
 * @param s [in] Price of stock
 * @param sigma [in] Vola of stock-price
 * @param t [in] Time to maturity in years
 * @param tau [in] Time of averaging period in years
 * @param r [in] interest rate
 * @param tmp1
 * @param tmp2
 * @param tmp3
 * @param tmp4
 * @param tmp5
 * @param sigmaA [out]
 * @param sigmaA2T2 [out]
 * @param sigmaAsqrtT [out]
 * @param emrt [out]
 */
void prepare_mkl_pricer(
        MKL_INT64 n,
        FLOAT *restrict s,                 /// [in] future price
        FLOAT *restrict sigma,             /// [in] vola
        FLOAT *restrict t,                 /// [in] time to maturity
        FLOAT *restrict tau,               /// [in] time of avg. period
        FLOAT *restrict r,                 /// [in] interest rate
        FLOAT *restrict tmp1,
        FLOAT *restrict tmp2,
        FLOAT *restrict tmp3,
        FLOAT *restrict tmp4,
        FLOAT *restrict tmp5,
        FLOAT *restrict sigmaA,           /// [out] adjusted vola
        FLOAT *restrict sigmaA2T2,        /// [out] sigmaA^2t/2
        FLOAT *restrict sigmaAsqrtT,      /// [out] sigmaA*sqrt(t)
        FLOAT *restrict emrt) {          /// [out] exp(-rt)/2

    FLOAT tt1;
    for (MKL_INT64 i = 0; i < n; ++i) {
        tt1 = sigma[i] * sigma[i];
        tmp1[i] = tt1 * t[i] + ln_of_2;
        tmp2[i] = tt1 * (t[i] - tau[i]) + ln_of_2;
        tmp5[i] = tt1 * tt1 * tau[i] * tau[i];
    }

    vdExp(n, tmp1, tmp3);
    vdExp(n, tmp2, tmp4);


    FLOAT tt3;
    for (MKL_INT64 i = 0; i < n; ++i) {
        tt3 = tmp4[i] * (1. + sigma[i] * sigma[i] * tau[i]);
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

/**
 * @brief Compute the prices of the asian options.
 *
 * This routine computes the prices of asian options. What kind of options is encoded in the
 * flags.
 *
 *  put     when flags & 1 != 0
 *  call    when flags & 1 == 0
 *
 *  short   when flags & 2 != 0
 *  long    when flags & 1 == 0
 *
 *  @note All arrays given to the routine are assumed to be of length n.
 *
 * @param n [in] Size of the arrays
 * @param flags [in]
 * @param s [in] stock prices
 * @param x [in] strikes
 * @param sigmaA2T2 [in] values computed by prepare_mkl_pricer()
 * @param sigmaAsqrtT [in] values computed by prepare_mkl_pricer()
 * @param emrt [in] values_computed by prepare_mkl_pricer()
 * @param tmp1
 * @param tmp2
 * @param d1 [out] computed by compute_d_values()
 * @param d2 [out] computed by compute_d_velues()
 * @param price [out]
 */

void mkl_pricer(
        MKL_INT64 n,
        MKL_INT64 *restrict flags,
        FLOAT *restrict s,                /// [in] stock price
        FLOAT *restrict x,                /// [in] strike
        FLOAT *restrict sigmaA2T2,        /// [in] sigmaA^2t/2
        FLOAT *restrict sigmaAsqrtT,      /// [in] sigmaA*sqrt(t)
        FLOAT *restrict emrt,
        FLOAT *restrict tmp1,
        FLOAT *restrict tmp2,
        FLOAT *restrict tmp3,
        FLOAT *restrict tmp4,
        FLOAT *restrict d1,               /// [out] d1
        FLOAT *restrict d2,
        FLOAT *restrict price) {            /// [out] d2

    compute_d_values(n, s, x, sigmaA2T2,
                     sigmaAsqrtT,
                     tmp1,
                     tmp2,
                     d1,
                     d2);

    for (MKL_INT64 i = 0; i < n; ++i) {
        if ((flags[i] & 1) != 0) {
            tmp3[i] = -d1[i];
            tmp4[i] = -d2[i];
        } else {
            tmp3[i] = d1[i];
            tmp4[i] = d2[i];
        }
    }


    vdErfc(n, tmp3, tmp1);
    vdErfc(n, tmp4, tmp2);

    double sx1, sx2;
    double w1, w2;
    for (MKL_INT64 i = 0; i < n; ++i) {
        if ((flags[i] & 1) != 0) {
            sx1 = x[i];
            sx2 = s[i];
            w1 = tmp2[i];
            w2 = tmp1[i];
        } else {
            sx1 = s[i];
            sx2 = x[i];
            w1 = tmp1[i];
            w2 = tmp2[i];
        }
        price[i] = sx1 * w1 - sx2 * w2;
        price[i] *= emrt[i];

        if ((flags[i] & 2) != 0) {
            price[i] = -price[i];
        }
    }

}

/**
 * @brief Computes the first partial derivative of the price according to the strike.
 *
 * @note All arrays are assumed to be of length n.
 *
 * @param n [in] Size of the arrays
 * @param flags [in] Flags, see mkl_pricer()
 * @param d2  [in] Values computed by mkl_pricer()
 * @param emrt  [in] Values computed by prepare_mkl_pricer()
 * @param ddx_price [out]
 */
void ddx_mkl_pricer(
        MKL_INT64 n,
        MKL_INT64 *restrict flags,
        FLOAT *restrict d2,
        FLOAT *restrict emrt,
        FLOAT *restrict ddx_price) {            /// [out] d2


    vdErfc(n, d2, ddx_price);

    double d;
    for (MKL_INT64 i = 0; i < n; ++i) {

        if ((flags[i] & 1) != 0) {
            d = 2.;
        } else {
            d = 0;
        }

        ddx_price[i] = d - ddx_price[i];
        ddx_price[i] *= emrt[i];

        if ((flags[i] & 2) != 0) {
            ddx_price[i] = -ddx_price[i];
        }
    }

}




