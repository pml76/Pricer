//
// Created by peter on 11/20/18.
//

/****************************************************/

#include <float.h>
#include <math.h>
#include <errno.h>


double ieee754_exp(double);

typedef unsigned char __u_char;
typedef unsigned short int __u_short;
typedef unsigned int __u_int;
typedef unsigned long int __u_long;


typedef signed char __int8_t;
typedef unsigned char __uint8_t;
typedef signed short int __int16_t;
typedef unsigned short int __uint16_t;
typedef signed int __int32_t;
typedef unsigned int __uint32_t;

typedef signed long int __int64_t;
typedef unsigned long int __uint64_t;

typedef __int8_t int8_t;
typedef __int16_t int16_t;
typedef __int32_t int32_t;
typedef __int64_t int64_t;

typedef __uint8_t uint8_t;
typedef __uint16_t uint16_t;
typedef __uint32_t uint32_t;
typedef __uint64_t uint64_t;


typedef union {
    double value;
    struct {
        uint32_t lsw;
        uint32_t msw;
    } parts;
    uint64_t word;
} ieee_double_shape_type;

/************************************************************/

static const double
        tiny = 1e-300,
        half = 5.00000000000000000000e-01,
        one = 1.00000000000000000000e+00,
        two = 2.00000000000000000000e+00,

        erx = 8.45062911510467529297e-01,


        efx = 1.28379167095512586316e-01,
        pp[] = {1.28379167095512558561e-01,
                -3.25042107247001499370e-01,
                -2.84817495755985104766e-02,
                -5.77027029648944159157e-03,
                -2.37630166566501626084e-05},
        qq[] = {0.0, 3.97917223959155352819e-01,
                6.50222499887672944485e-02,
                5.08130628187576562776e-03,
                1.32494738004321644526e-04,
                -3.96022827877536812320e-06},


        pa[] = {-2.36211856075265944077e-03,
                4.14856118683748331666e-01,
                -3.72207876035701323847e-01,
                3.18346619901161753674e-01,
                -1.10894694282396677476e-01,
                3.54783043256182359371e-02,
                -2.16637559486879084300e-03},
        qa[] = {0.0, 1.06420880400844228286e-01,
                5.40397917702171048937e-01,
                7.18286544141962662868e-02,
                1.26171219808761642112e-01,
                1.36370839120290507362e-02,
                1.19844998467991074170e-02},


        ra[] = {-9.86494403484714822705e-03,
                -6.93858572707181764372e-01,
                -1.05586262253232909814e+01,
                -6.23753324503260060396e+01,
                -1.62396669462573470355e+02,
                -1.84605092906711035994e+02,
                -8.12874355063065934246e+01,
                -9.81432934416914548592e+00},
        sa[] = {0.0, 1.96512716674392571292e+01,
                1.37657754143519042600e+02,
                4.34565877475229228821e+02,
                6.45387271733267880336e+02,
                4.29008140027567833386e+02,
                1.08635005541779435134e+02,
                6.57024977031928170135e+00,
                -6.04244152148580987438e-02},


        rb[] = {-9.86494292470009928597e-03,
                -7.99283237680523006574e-01,
                -1.77579549177547519889e+01,
                -1.60636384855821916062e+02,
                -6.37566443368389627722e+02,
                -1.02509513161107724954e+03,
                -4.83519191608651397019e+02},
        sb[] = {0.0, 3.03380607434824582924e+01,
                3.25792512996573918826e+02,
                1.53672958608443695994e+03,
                3.19985821950859553908e+03,
                2.55305040643316442583e+03,
                4.74528541206955367215e+02,
                -2.24409524465858183362e+01};

double glibc_erf(double x) {
    int32_t hx, ix, i;
    double R, S, P, Q, s, y, z, r;
    do {
        ieee_double_shape_type gh_u;
        gh_u.value = (x);
        (hx) = gh_u.parts.msw;
    }
    while (0);
    ix = hx & 0x7fffffff;
    if (ix >= 0x7ff00000) {
        i = ((uint32_t) hx >> 31) << 1;
        return (double) (1 - i) + one / x;
    }

    if (ix < 0x3feb0000) {
        double r1, r2, s1, s2, s3, z2, z4;
        if (ix < 0x3e300000) {
            if (ix < 0x00800000) {

                double ret = 0.0625 * (16.0 * x + (16.0 * efx) * x);
                do {
                    __typeof(ret) force_underflow_tmp = (ret);
                    if (_Generic (((force_underflow_tmp)), float: (__typeof(force_underflow_tmp)) __builtin_fabsf(
                            force_underflow_tmp), _Float32: (__typeof(force_underflow_tmp)) __builtin_fabsf(
                            force_underflow_tmp), default: (__typeof(force_underflow_tmp)) __builtin_fabs(
                            force_underflow_tmp), long double: (__typeof(force_underflow_tmp)) __builtin_fabsl(
                            force_underflow_tmp), _Float64x: (__typeof(force_underflow_tmp)) __builtin_fabsl(
                            force_underflow_tmp), _Float128: (__typeof(force_underflow_tmp)) __builtin_fabsf128(
                            force_underflow_tmp)) <
                        _Generic (((force_underflow_tmp)), float: (__typeof(force_underflow_tmp)) 1.17549435082228750796873653722224568e-38F, _Float32: (__typeof(force_underflow_tmp)) 1.17549435082228750796873653722224568e-38F, default: (__typeof(force_underflow_tmp)) ((double) 2.22507385850720138309023271733240406e-308L), long double: (__typeof(force_underflow_tmp)) 3.36210314311209350626267781732175260e-4932L, _Float64x: (__typeof(force_underflow_tmp)) 3.36210314311209350626267781732175260e-4932L, _Float128: (__typeof(force_underflow_tmp)) 3.36210314311209350626267781732175260e-4932F128)) {
                        __typeof(force_underflow_tmp) force_underflow_tmp2 = force_underflow_tmp * force_underflow_tmp;
                        do {
                            if (sizeof(force_underflow_tmp2) <= sizeof(double) ||
                                __builtin_types_compatible_p (__typeof(force_underflow_tmp2), _Float128))
                                    __asm __volatile ("" : : "x" (force_underflow_tmp2));
                            else __asm __volatile ("" : : "f" (force_underflow_tmp2));
                        }
                        while (0);
                    }
                }
                while (0);
                return ret;
            }
            return x + efx * x;
        }
        z = x * x;
        r1 = pp[0] + z * pp[1];
        z2 = z * z;
        r2 = pp[2] + z * pp[3];
        z4 = z2 * z2;
        s1 = one + z * qq[1];
        s2 = qq[2] + z * qq[3];
        s3 = qq[4] + z * qq[5];
        r = r1 + z2 * r2 + z4 * pp[4];
        s = s1 + z2 * s2 + z4 * s3;
        y = r / s;
        return x + x * y;
    }
    if (ix < 0x3ff40000) {
        double s2, s4, s6, P1, P2, P3, P4, Q1, Q2, Q3, Q4;
        s = fabs(x) - one;
        P1 = pa[0] + s * pa[1];
        s2 = s * s;
        Q1 = one + s * qa[1];
        s4 = s2 * s2;
        P2 = pa[2] + s * pa[3];
        s6 = s4 * s2;
        Q2 = qa[2] + s * qa[3];
        P3 = pa[4] + s * pa[5];
        Q3 = qa[4] + s * qa[5];
        P4 = pa[6];
        Q4 = qa[6];
        P = P1 + s2 * P2 + s4 * P3 + s6 * P4;
        Q = Q1 + s2 * Q2 + s4 * Q3 + s6 * Q4;
        if (hx >= 0)
            return erx + P / Q;
        else
            return -erx - P / Q;
    }
    if (ix >= 0x40180000) {
        if (hx >= 0)
            return one - tiny;
        else
            return tiny - one;
    }
    x = fabs(x);
    s = one / (x * x);
    if (ix < 0x4006DB6E) {
        double R1, R2, R3, R4, S1, S2, S3, S4, s2, s4, s6, s8;
        R1 = ra[0] + s * ra[1];
        s2 = s * s;
        S1 = one + s * sa[1];
        s4 = s2 * s2;
        R2 = ra[2] + s * ra[3];
        s6 = s4 * s2;
        S2 = sa[2] + s * sa[3];
        s8 = s4 * s4;
        R3 = ra[4] + s * ra[5];
        S3 = sa[4] + s * sa[5];
        R4 = ra[6] + s * ra[7];
        S4 = sa[6] + s * sa[7];
        R = R1 + s2 * R2 + s4 * R3 + s6 * R4;
        S = S1 + s2 * S2 + s4 * S3 + s6 * S4 + s8 * sa[8];
    } else {
        double R1, R2, R3, S1, S2, S3, S4, s2, s4, s6;
        R1 = rb[0] + s * rb[1];
        s2 = s * s;
        S1 = one + s * sb[1];
        s4 = s2 * s2;
        R2 = rb[2] + s * rb[3];
        s6 = s4 * s2;
        S2 = sb[2] + s * sb[3];
        R3 = rb[4] + s * rb[5];
        S3 = sb[4] + s * sb[5];
        S4 = sb[6] + s * sb[7];
        R = R1 + s2 * R2 + s4 * R3 + s6 * rb[6];
        S = S1 + s2 * S2 + s4 * S3 + s6 * S4;
    }
    z = x;
    do {
        ieee_double_shape_type sl_u;
        sl_u.value = (z);
        sl_u.parts.lsw = (0);
        (z) = sl_u.value;
    }
    while (0);
    r = ieee754_exp(-z * z - 0.5625) *
        ieee754_exp((z - x) * (z + x) + R / S);
    if (hx >= 0)
        return one - r / x;
    else
        return r / x - one;
}

double glibc_erfc(double x) {
    int32_t hx, ix;
    double R, S, P, Q, s, y, z, r;
    do {
        ieee_double_shape_type gh_u;
        gh_u.value = (x);
        (hx) = gh_u.parts.msw;
    }
    while (0);
    ix = hx & 0x7fffffff;
    if (ix >= 0x7ff00000) {
        double ret = (double) (((uint32_t) hx >> 31) << 1) + one / x;
        if (0 && ret == 0.0)
            return 0.0;
        return ret;
    }

    if (ix < 0x3feb0000) {
        double r1, r2, s1, s2, s3, z2, z4;
        if (ix < 0x3c700000)
            return one - x;
        z = x * x;
        r1 = pp[0] + z * pp[1];
        z2 = z * z;
        r2 = pp[2] + z * pp[3];
        z4 = z2 * z2;
        s1 = one + z * qq[1];
        s2 = qq[2] + z * qq[3];
        s3 = qq[4] + z * qq[5];
        r = r1 + z2 * r2 + z4 * pp[4];
        s = s1 + z2 * s2 + z4 * s3;
        y = r / s;
        if (hx < 0x3fd00000) {
            return one - (x + x * y);
        } else {
            r = x * y;
            r += (x - half);
            return half - r;
        }
    }
    if (ix < 0x3ff40000) {
        double s2, s4, s6, P1, P2, P3, P4, Q1, Q2, Q3, Q4;
        s = fabs(x) - one;
        P1 = pa[0] + s * pa[1];
        s2 = s * s;
        Q1 = one + s * qa[1];
        s4 = s2 * s2;
        P2 = pa[2] + s * pa[3];
        s6 = s4 * s2;
        Q2 = qa[2] + s * qa[3];
        P3 = pa[4] + s * pa[5];
        Q3 = qa[4] + s * qa[5];
        P4 = pa[6];
        Q4 = qa[6];
        P = P1 + s2 * P2 + s4 * P3 + s6 * P4;
        Q = Q1 + s2 * Q2 + s4 * Q3 + s6 * Q4;
        if (hx >= 0) {
            z = one - erx;
            return z - P / Q;
        } else {
            z = erx + P / Q;
            return one + z;
        }
    }
    if (ix < 0x403c0000) {
        x = fabs(x);
        s = one / (x * x);
        if (ix < 0x4006DB6D) {
            double R1, R2, R3, R4, S1, S2, S3, S4, s2, s4, s6, s8;
            R1 = ra[0] + s * ra[1];
            s2 = s * s;
            S1 = one + s * sa[1];
            s4 = s2 * s2;
            R2 = ra[2] + s * ra[3];
            s6 = s4 * s2;
            S2 = sa[2] + s * sa[3];
            s8 = s4 * s4;
            R3 = ra[4] + s * ra[5];
            S3 = sa[4] + s * sa[5];
            R4 = ra[6] + s * ra[7];
            S4 = sa[6] + s * sa[7];
            R = R1 + s2 * R2 + s4 * R3 + s6 * R4;
            S = S1 + s2 * S2 + s4 * S3 + s6 * S4 + s8 * sa[8];
        } else {
            double R1, R2, R3, S1, S2, S3, S4, s2, s4, s6;
            if (hx < 0 && ix >= 0x40180000)
                return two - tiny;
            R1 = rb[0] + s * rb[1];
            s2 = s * s;
            S1 = one + s * sb[1];
            s4 = s2 * s2;
            R2 = rb[2] + s * rb[3];
            s6 = s4 * s2;
            S2 = sb[2] + s * sb[3];
            R3 = rb[4] + s * rb[5];
            S3 = sb[4] + s * sb[5];
            S4 = sb[6] + s * sb[7];
            R = R1 + s2 * R2 + s4 * R3 + s6 * rb[6];
            S = S1 + s2 * S2 + s4 * S3 + s6 * S4;
        }
        z = x;
        do {
            ieee_double_shape_type sl_u;
            sl_u.value = (z);
            sl_u.parts.lsw = (0);
            (z) = sl_u.value;
        }
        while (0);
        r = ieee754_exp(-z * z - 0.5625) *
            ieee754_exp((z - x) * (z + x) + R / S);
        if (hx > 0) {
            double ret = (r / x);
            if (ret == 0)
                (errno = (
                        34
                ));
            return ret;
        } else
            return two - r / x;
    } else {
        if (hx > 0) {
            (errno = (
                    34
            ));
            return tiny * tiny;
        } else
            return two - tiny;
    }
}






