/*
 *
 * (c) 2019, by Peter Lennartz  // peter.lennartz@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */



// #pragma GCC optimize("O3", "unroll-loops", "no-omit-frame-pointer", "inline") //Optimization flags
// #pragma GCC option("arch=native", "tune=native", "no-zero-upper", "funsafe-math") //Enable AVX

#include <mkl.h>
#include <math.h>
#include <math/pricers/mkl_pricer.h>


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
        UINT64 n,
        Real_Ptr s,                /// [in] stock price
        Real_Ptr x,                /// [in] strike
        Real_Ptr sigmaA2T2,        /// [in] sigmaA^2t/2
        Real_Ptr sigmaAsqrtT,      /// [in] sigmaA*sqrt(t)
        Real_Ptr tmp1,
        Real_Ptr tmp2,
        Real_Ptr d1,               /// [out] d1
        Real_Ptr d2) {            /// [out] d2


    ASSUME(n % 64 == 0)

    // for optimization assume that every
    // pointer handed to this function is aligned.
    ASSUME_ALIGNED(s)
    ASSUME_ALIGNED(x)
    ASSUME_ALIGNED(sigmaA2T2)
    ASSUME_ALIGNED(sigmaAsqrtT)
    ASSUME_ALIGNED(tmp1)
    ASSUME_ALIGNED(tmp2)
    ASSUME_ALIGNED(d1)
    ASSUME_ALIGNED(d2)


    for (UINT64 i = 0; i < n; ++i) {
        tmp1[i] = s[i] / x[i];
    }

    vdLn(n, tmp1, tmp2);

    for (UINT64 i = 0; i < n; ++i) {
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
        UINT64 n,
        Real_Ptr s,                 /// [in] future price
        Real_Ptr sigma,             /// [in] vola
        Real_Ptr t,                 /// [in] time to maturity
        Real_Ptr tau,               /// [in] time of avg. period
        Real_Ptr r,                 /// [in] interest rate
        Real_Ptr tmp1,
        Real_Ptr tmp2,
        Real_Ptr tmp3,
        Real_Ptr tmp4,
        Real_Ptr tmp5,
        Real_Ptr sigmaA,           /// [out] adjusted vola
        Real_Ptr sigmaA2T2,        /// [out] sigmaA^2t/2
        Real_Ptr sigmaAsqrtT,      /// [out] sigmaA*sqrt(t)
        Real_Ptr emrt,             /// [out] exp(-rt)/2
        Real_Ptr d2dx2_prep) {

    ASSUME(n % 64 == 0)

    ASSUME_ALIGNED(s)
    ASSUME_ALIGNED(sigma)
    ASSUME_ALIGNED(t)
    ASSUME_ALIGNED(tau)
    ASSUME_ALIGNED(r)
    ASSUME_ALIGNED(tmp1)
    ASSUME_ALIGNED(tmp2)
    ASSUME_ALIGNED(tmp3)
    ASSUME_ALIGNED(tmp4)
    ASSUME_ALIGNED(tmp5)
    ASSUME_ALIGNED(sigmaA)
    ASSUME_ALIGNED(sigmaA2T2)
    ASSUME_ALIGNED(sigmaAsqrtT)
    ASSUME_ALIGNED(emrt)
    ASSUME_ALIGNED(d2dx2_prep)

    FLOAT tt1;
    for (UINT64 i = 0; i < n; ++i) {
        tt1 = sigma[i] * sigma[i];
        tmp1[i] = tt1 * t[i] + ln_of_2;
        tmp2[i] = tt1 * (t[i] - tau[i]) + ln_of_2;
        tmp5[i] = tt1 * tt1 * tau[i] * tau[i];
    }

    vdExp(n, tmp1, tmp3);
    vdExp(n, tmp2, tmp4);


    FLOAT tt3;
    for (UINT64 i = 0; i < n; ++i) {
        tt3 = tmp4[i] * (1. + sigma[i] * sigma[i] * tau[i]);
        tmp2[i] = (tmp3[i] - tt3) / tmp5[i];
    }

    vdLn(n, tmp2, tmp1);
    for (UINT64 i = 0; i < n; ++i) {
        tmp1[i] /= t[i];
    }

    vdSqrt(n, tmp1, sigmaA);
    vdSqrt(n, t, tmp1);


    for (UINT64 i = 0; i < n; ++i) {
        sigmaA2T2[i] = sigmaA[i] * sigmaA[i] * t[i] / 2.;
        sigmaAsqrtT[i] = sigmaA[i] * tmp1[i];
        tmp2[i] = -r[i] * t[i] - ln_of_2;
        tmp3[i] = -sigmaA2T2[i] / 4;
    }

    vdExp(n, tmp2, emrt);
    vdExp(n, tmp3, tmp2);

    for (UINT64 i = 0; i < n; ++i) {
        tmp1[i] = 2. * s[i] / (2. * M_PI * sigmaA2T2[i]);
    }

    vdSqrt(n, tmp1, d2dx2_prep);

    for (UINT64 i = 0; i < n; ++i) {
        d2dx2_prep[i] *= tmp2[i] * emrt[i];
    }


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
        UINT64 n,
        Uint64_Ptr flags,
        Real_Ptr s,                /// [in] stock price
        Real_Ptr x,                /// [in] strike
        Real_Ptr sigmaA2T2,        /// [in] sigmaA^2t/2
        Real_Ptr sigmaAsqrtT,      /// [in] sigmaA*sqrt(t)
        Real_Ptr emrt,
        Real_Ptr tmp1,
        Real_Ptr tmp2,
        Real_Ptr tmp3,
        Real_Ptr tmp4,
        Real_Ptr d1,               /// [out] d1
        Real_Ptr d2,
        Real_Ptr price) {            /// [out] d2

    ASSUME(n % 64 == 0)

    ASSUME_ALIGNED(flags)
    ASSUME_ALIGNED(x)
    ASSUME_ALIGNED(s)
    ASSUME_ALIGNED(sigmaA2T2)
    ASSUME_ALIGNED(sigmaAsqrtT)
    ASSUME_ALIGNED(emrt)
    ASSUME_ALIGNED(tmp1)
    ASSUME_ALIGNED(tmp2)
    ASSUME_ALIGNED(tmp3)
    ASSUME_ALIGNED(tmp4)
    ASSUME_ALIGNED(d1)
    ASSUME_ALIGNED(d2)
    ASSUME_ALIGNED(price)

    compute_d_values(n, s, x, sigmaA2T2,
                     sigmaAsqrtT,
                     tmp1,
                     tmp2,
                     d1,
                     d2);

    for (UINT64 i = 0; i < n; ++i) {
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
    for (UINT64 i = 0; i < n; ++i) {
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
        UINT64 n,
        Uint64_Ptr flags,
        Real_Ptr d2,
        Real_Ptr emrt,
        Real_Ptr ddx_price) {            /// [out] d2

    ASSUME(n % 64 == 0)

    ASSUME_ALIGNED(flags)
    ASSUME_ALIGNED(d2)
    ASSUME_ALIGNED(emrt)
    ASSUME_ALIGNED(ddx_price)


    vdErfc(n, d2, ddx_price);

    double d;
    for (UINT64 i = 0; i < n; ++i) {

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


void d2dx2_mkl_pricer(
        UINT64 n,
        Uint64_Ptr flags,
        Real_Ptr s,
        Real_Ptr x,
        Real_Ptr d2dx2_prep,
        Real_Ptr sigmaA2T2,
        Real_Ptr tmp1,
        Real_Ptr tmp2,
        Real_Ptr d2dx2
) {

    ASSUME(n % 64 == 0)

    ASSUME_ALIGNED(flags)
    ASSUME_ALIGNED(s)
    ASSUME_ALIGNED(x)
    ASSUME_ALIGNED(d2dx2_prep)
    ASSUME_ALIGNED(sigmaA2T2)
    ASSUME_ALIGNED(tmp1)
    ASSUME_ALIGNED(tmp2)
    ASSUME_ALIGNED(d2dx2)

    for (UINT64 i = 0; i < n; ++i) {
        tmp1[i] = s[i] / x[i];
    }

    vdLn(n, tmp1, d2dx2);

    for (UINT64 i = 0; i < n; ++i) {
        tmp1[i] = -d2dx2[i] * d2dx2[i] / (4. * sigmaA2T2[i]);
    }

    vdExp(n, tmp1, d2dx2);

    for (UINT64 i = 0; i < n; ++i) {
        tmp1[i] = x[i] * x[i] * x[i];
    }

    vdSqrt(n, tmp1, tmp2);

    for (UINT64 i = 0; i < n; ++i) {
        d2dx2[i] = d2dx2[i] * d2dx2_prep[i] / tmp2[i];

        if ((flags[i] & 2) != 0) {
            d2dx2[i] = -d2dx2[i];
        }
    }

}




