//
// Created by peter on 1/25/20.
//

#include "my_pricer.h"



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

// #include <my.h>
#include <math.h>
#include <mkl.h>
#include <pricer-base.h>

static FLOAT ln_of_2, msqrt2;

static const double ln2 = log(2);


#ifndef __INTEL_COMPILER
#define CONST const
#else
#define CONST
#endif

#define MLA mla
#define C2V(x) (x)

#define POLY2(x, c1, c0) MLA(x, C2V(c1), C2V(c0))
#define POLY3(x, x2, c2, c1, c0) MLA(x2, C2V(c2), MLA(x, C2V(c1), C2V(c0)))
#define POLY4(x, x2, c3, c2, c1, c0) MLA(x2, MLA(x, C2V(c3), C2V(c2)), MLA(x, C2V(c1), C2V(c0)))
#define POLY5(x, x2, x4, c4, c3, c2, c1, c0) MLA(x4, C2V(c4), POLY4(x, x2, c3, c2, c1, c0))
#define POLY6(x, x2, x4, c5, c4, c3, c2, c1, c0) MLA(x4, POLY2(x, c5, c4), POLY4(x, x2, c3, c2, c1, c0))
#define POLY7(x, x2, x4, c6, c5, c4, c3, c2, c1, c0) MLA(x4, POLY3(x, x2, c6, c5, c4), POLY4(x, x2, c3, c2, c1, c0))
#define POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0) MLA(x4, POLY4(x, x2, c7, c6, c5, c4), POLY4(x, x2, c3, c2, c1, c0))
#define POLY9(x, x2, x4, x8, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, C2V(c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY10(x, x2, x4, x8, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY2(x, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY11(x, x2, x4, x8, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY3(x, x2, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY12(x, x2, x4, x8, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY4(x, x2, cb, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY13(x, x2, x4, x8, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY5(x, x2, x4, cc, cb, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY14(x, x2, x4, x8, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY6(x, x2, x4, cd, cc, cb, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY15(x, x2, x4, x8, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY7(x, x2, x4, ce, cd, cc, cb, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY8(x, x2, x4, cf, ce, cd, cc, cb, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY17(x, x2, x4, x8, x16, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, C2V(d0), POLY16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY18(x, x2, x4, x8, x16, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, POLY2(x, d1, d0), POLY16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY19(x, x2, x4, x8, x16, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, POLY3(x, x2, d2, d1, d0), POLY16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))

static INLINE CONST double mla(double x, double y, double z) { return x * y + z; }

static INLINE CONST int64_t doubleToRawLongBits(double d) {
    union {
        double f;
        int64_t i;
    } tmp;
    tmp.f = d;
    return tmp.i;
}

static INLINE CONST double longBitsToDouble(int64_t i) {
    union {
        double f;
        int64_t i;
    } tmp;
    tmp.i = i;
    return tmp.f;
}


static INLINE CONST int ilogb2k(double d) {
    return ((doubleToRawLongBits(d) >> 52) & 0x7ff) - 0x3ff;
}

static INLINE CONST double ldexp3k(double d, int e) { // very fast, no denormal
    return longBitsToDouble(doubleToRawLongBits(d) + (((int64_t)e) << 52));
}


static CONST double my_ln(double d) {
    double x, x2, t, m;
    int e;

    e = ilogb2k(d * (1.0/0.75));
    m = ldexp3k(d, -e);

    x = (m-1) / (m+1);
    x2 = x * x;

    double x4 = x2 * x2, x8 = x4 * x4;

    t = POLY7(x2, x4, x8,
              0.153487338491425068243146,
              0.152519917006351951593857,
              0.181863266251982985677316,
              0.222221366518767365905163,
              0.285714294746548025383248,
              0.399999999950799600689777,
              0.6666666666667778740063);

    x = x * 2 + 0.693147180559945286226764 * e + x * x2 * t;

    return x;
}


/*
inline
static double my_ln(double x) {

    static const double a1 = 0.9999964239, a2 = -0.4998741238, a3 = 0.3317990258, a4 = -0.2407338084,
        a5 = 0.1676540711, a6 = -0.0953293897, a7 = 0.0360884937, a8 = -0.0064535442;

    double x0 = 0.;
    while (x > 2.) {
        x /= 2.;
        x0 += ln2;
    }

    while (x < 1.) {
        x *= 2.;
        x0 -= ln2;
    }

    x = x - 1.;

    return ((((((a8*x + a7)*x + a6)*a5 + a4)*x + a3)*x + a2)*x+a1)*x;
}

*/


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
 * @param sigmaA2T2 [in] Array computed by prepare_my_pricer()
 * @param sigmaAsqrtT [in[ Array computed by prepare_my_pricer()
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
    ASSUME_ALIGNED(Real_Ptr ,s)
    ASSUME_ALIGNED(Real_Ptr ,x)
    ASSUME_ALIGNED(Real_Ptr ,sigmaA2T2)
    ASSUME_ALIGNED(Real_Ptr ,sigmaAsqrtT)
    ASSUME_ALIGNED(Real_Ptr ,tmp1)
    ASSUME_ALIGNED(Real_Ptr ,tmp2)
    ASSUME_ALIGNED(Real_Ptr ,d1)
    ASSUME_ALIGNED(Real_Ptr ,d2)


    for (UINT64 i = 0; i < n; ++i) {
        tmp1[i] = s[i] / x[i];
        tmp2[i] = my_ln(tmp1[i]);
    }

//    vdLn(n, tmp1, tmp2);

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
void init_my_pricer() {
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
void prepare_my_pricer(
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

    ASSUME_ALIGNED(Real_Ptr ,s)
    ASSUME_ALIGNED(Real_Ptr ,sigma)
    ASSUME_ALIGNED(Real_Ptr ,t)
    ASSUME_ALIGNED(Real_Ptr ,tau)
    ASSUME_ALIGNED(Real_Ptr ,r)
    ASSUME_ALIGNED(Real_Ptr ,tmp1)
    ASSUME_ALIGNED(Real_Ptr ,tmp2)
    ASSUME_ALIGNED(Real_Ptr ,tmp3)
    ASSUME_ALIGNED(Real_Ptr ,tmp4)
    ASSUME_ALIGNED(Real_Ptr ,tmp5)
    ASSUME_ALIGNED(Real_Ptr ,sigmaA)
    ASSUME_ALIGNED(Real_Ptr ,sigmaA2T2)
    ASSUME_ALIGNED(Real_Ptr ,sigmaAsqrtT)
    ASSUME_ALIGNED(Real_Ptr ,emrt)
    ASSUME_ALIGNED(Real_Ptr ,d2dx2_prep)

    FLOAT tt1;
    for (UINT64 i = 0; i < n; ++i) {
        tt1 = sigma[i] * sigma[i];
        tmp1[i] = tt1 * t[i] + ln_of_2;
        tmp2[i] = tt1 * (t[i] - tau[i]) + ln_of_2;
        tmp5[i] = tt1 * tt1 * tau[i] * tau[i];
    }


    for( uint64_t i = 0; i < n; ++i) {
        tmp3[i] = exp(tmp1[i]);
        tmp4[i] = exp(tmp2[i]);
    }
//    vdExp(n, tmp1, tmp3);
//    vdExp(n, tmp2, tmp4);


    FLOAT tt3;
    for (UINT64 i = 0; i < n; ++i) {
        tt3 = tmp4[i] * (1. + sigma[i] * sigma[i] * tau[i]);
        tmp2[i] = (tmp3[i] - tt3) / tmp5[i];
    }

//    vdLn(n, tmp2, tmp1);

    for (UINT64 i = 0; i < n; ++i) {
        tmp1[i] = log(tmp2[i]);
        tmp1[i] /= t[i];
    }

    for(uint64_t i = 0; i < n; ++i) {
        sigmaA[i] = sqrt(tmp1[i]);
        tmp1[i] = sqrt(t[i]);
    }
    //vdSqrt(n, tmp1, sigmaA);
    //vdSqrt(n, t, tmp1);


    for (UINT64 i = 0; i < n; ++i) {
        sigmaA2T2[i] = sigmaA[i] * sigmaA[i] * t[i] / 2.;
        sigmaAsqrtT[i] = sigmaA[i] * tmp1[i];
        tmp2[i] = -r[i] * t[i] - ln_of_2;
       // tmp3[i] = -sigmaA2T2[i] / 4;
    }

    for(uint64_t i = 0; i < n; i++) {
        emrt[i] = exp(tmp2[i]);

    }

    //vdExp(n, tmp2, emrt);

    // vdExp(n, tmp3, tmp2);

    /*
    for (UINT64 i = 0; i < n; ++i) {
        tmp1[i] = 2. * s[i] / (2. * M_PI * sigmaA2T2[i]);
    }

    vdSqrt(n, tmp1, d2dx2_prep);

    for (UINT64 i = 0; i < n; ++i) {
        d2dx2_prep[i] *= tmp2[i] * emrt[i];
    }
*/

}

inline
static double my_erfc(double x) {
    union {
        double d;
        uint64_t i;
    } d1, d2;

    static const double  a1 = 0.0705230784, a2 = 0.0422820123, a3 = 0.0092705272,
        a4 = 0.0001520143, a5 = 0.0002765672, a6 = 0.0000430638;
    double x2 = x*x;
    d1.d = x;
    d1.i &= 0x7FFFFFFF;

    d2.d = 1./((1.+d1.d*a1+x2*a2)+x2*(d1.d*a3 + x2*(a4 + d1.d*a5 +x2*a6)));
    d2.d *= d2.d;
    d2.d *= d2.d;
    d2.d *= d2.d;
    d2.d *= d2.d;

    return d2.d;

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
 * @param sigmaA2T2 [in] values computed by prepare_my_pricer()
 * @param sigmaAsqrtT [in] values computed by prepare_my_pricer()
 * @param emrt [in] values_computed by prepare_my_pricer()
 * @param tmp1
 * @param tmp2
 * @param d1 [out] computed by compute_d_values()
 * @param d2 [out] computed by compute_d_velues()
 * @param price [out]
 */

void my_pricer(
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

    ASSUME_ALIGNED(Uint64_Ptr,flags)
    ASSUME_ALIGNED(Real_Ptr ,x)
    ASSUME_ALIGNED(Real_Ptr ,s)
    ASSUME_ALIGNED(Real_Ptr ,sigmaA2T2)
    ASSUME_ALIGNED(Real_Ptr ,sigmaAsqrtT)
    ASSUME_ALIGNED(Real_Ptr ,emrt)
    ASSUME_ALIGNED(Real_Ptr ,tmp1)
    ASSUME_ALIGNED(Real_Ptr ,tmp2)
    ASSUME_ALIGNED(Real_Ptr ,tmp3)
    ASSUME_ALIGNED(Real_Ptr ,tmp4)
    ASSUME_ALIGNED(Real_Ptr ,d1)
    ASSUME_ALIGNED(Real_Ptr ,d2)
    ASSUME_ALIGNED(Real_Ptr ,price)

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


    for(uint64_t i = 0; i < n; i++) {
        tmp1[i] = my_erfc(tmp3[i]);
        tmp2[i] = my_erfc(tmp4[i]);
        
        // vdErfc(n, tmp3, tmp1);
        // vdErfc(n, tmp4, tmp2);
    }

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
