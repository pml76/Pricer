//
// Created by peter on 12/21/18.
//



#include "misc.h"

extern const double rempitabdp[];

#define __SLEEFSIMDDP_C__

#if (defined(_MSC_VER))
#pragma fp_contract (off)
#endif

#ifdef ENABLE_SSE2
#define CONFIG 2
#include "helpersse2.h"
#ifdef DORENAME
#ifdef ENABLE_GNUABI
#include "renamesse2_gnuabi.h"
#else
#include "renamesse2.h"
#endif
#endif
#endif

#ifdef ENABLE_SSE4
#define CONFIG 4
#include "helpersse2.h"
#ifdef DORENAME
#include "renamesse4.h"
#endif
#endif

#ifdef ENABLE_AVX
#define CONFIG 1

#include "helperavx.h"

#ifdef DORENAME
#ifdef ENABLE_GNUABI
#include "renameavx_gnuabi.h"
#else
#include "renameavx.h"
#endif
#endif
#endif

#ifdef ENABLE_FMA4
#define CONFIG 4
#include "helperavx.h"
#ifdef DORENAME
#ifdef ENABLE_GNUABI
#include "renamefma4_gnuabi.h"
#else
#include "renamefma4.h"
#endif
#endif
#endif

#ifdef ENABLE_AVX2
#define CONFIG 1
#include "helperavx2.h"
#ifdef DORENAME
#ifdef ENABLE_GNUABI
#include "renameavx2_gnuabi.h"
#else
#include "renameavx2.h"
#endif
#endif
#endif

#ifdef ENABLE_AVX2128
#define CONFIG 1
#include "helperavx2_128.h"
#ifdef DORENAME
#include "renameavx2128.h"
#endif
#endif

#ifdef ENABLE_AVX512F
#define CONFIG 1
#include "helperavx512f.h"
#ifdef DORENAME
#ifdef ENABLE_GNUABI
#include "renameavx512f_gnuabi.h"
#else
#include "renameavx512f.h"
#endif
#endif
#endif

#ifdef ENABLE_ADVSIMD
#define CONFIG 1
#include "helperadvsimd.h"
#ifdef DORENAME
#ifdef ENABLE_GNUABI
#include "renameadvsimd_gnuabi.h"
#else
#include "renameadvsimd.h"
#endif
#endif
#endif

#ifdef ENABLE_VSX
#define CONFIG 1
#include "helperpower_128.h"
#ifdef DORENAME
#include "renamevsx.h"
#endif
#endif

//

#ifdef ENABLE_VECEXT
#define CONFIG 1
#include "helpervecext.h"
#ifdef DORENAME
#include "renamevecext.h"
#endif
#endif

#ifdef ENABLE_PUREC
#define CONFIG 1
#include "helperpurec.h"
#ifdef DORENAME
#include "renamepurec.h"
#endif
#endif

#ifdef ENABLE_SVE
#define CONFIG 1
#include "helpersve.h"
#ifdef DORENAME
#ifdef ENABLE_GNUABI
#include "renamesve_gnuabi.h"
#else
#include "renamesve.h"
#endif /* ENABLE_GNUABI */
#endif /* DORENAME */
#endif /* ENABLE_SVE */

//

#include <src/math/pricers/pricer-base.h>
#include <tgmath.h>

static vdouble ln_of_2, msqrt2;
vdouble one, two, four, pi2, one_over_two, zero;


vdouble xexp(vdouble d);

vdouble xlog(vdouble d);

vdouble xsqrt(vdouble d);

vdouble xerfc_u15(vdouble a);


void init_sleef_pricer() {
    ln_of_2 = vcast_vd_d(log(2.));
    msqrt2 = vcast_vd_d(-sqrt(2));
    one = vcast_vd_d(1.);
    two = vcast_vd_d(2.);
    four = vcast_vd_d(4.);
    pi2 = vmul_vd_vd_vd(two, vcast_vd_d(M_PI));
    one_over_two = vdiv_vd_vd_vd(one, two);
    zero = vcast_vd_d(0.);

}


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
        Real_Ptr d2dx2_prep_) {

    ASSUME(n % 64 == 0)

    ASSUME_ALIGNED(s_)
    ASSUME_ALIGNED(sigma_)
    ASSUME_ALIGNED(t_)
    ASSUME_ALIGNED(tau_)
    ASSUME_ALIGNED(r_)
    ASSUME_ALIGNED(sigmaA_)
    ASSUME_ALIGNED(sigmaA2T2_)
    ASSUME_ALIGNED(sigmaAsqrtT_)
    ASSUME_ALIGNED(emrt_)
    ASSUME_ALIGNED(d2dx2_prep_)

    vdouble tmp1, tmp2, tmp3, tmp4, tmp5;
    vdouble tt1, tt3;
    vdouble s, sigma, t, tau, r, sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep;

    for (UINT64 i = 0; i < n; i += sizeof(vdouble) / sizeof(double)) {

        s = vload_vd_p(&s_[i]);
        sigma = vload_vd_p(&sigma_[i]);
        t = vload_vd_p(&t_[i]);
        tau = vload_vd_p(&tau_[i]);
        r = vload_vd_p(&r_[i]);

        tt1 = vmul_vd_vd_vd(sigma, sigma);
        tmp1 = vadd_vd_vd_vd(ln_of_2, vmul_vd_vd_vd(tt1, t));
        tmp2 = vadd_vd_vd_vd(ln_of_2, vmul_vd_vd_vd(tt1, vsub_vd_vd_vd(t, tau)));
        tmp5 = vmul_vd_vd_vd(vmul_vd_vd_vd(tt1, tt1), vmul_vd_vd_vd(tau, tau));

        tmp3 = xexp(tmp1);
        tmp4 = xexp(tmp2);

        tt3 = vmul_vd_vd_vd(tmp4, vadd_vd_vd_vd(one, vmul_vd_vd_vd(sigma, vmul_vd_vd_vd(sigma, tau))));
        tmp2 = vdiv_vd_vd_vd(vsub_vd_vd_vd(tmp3, tt3), tmp5);

        sigmaA = xsqrt(vdiv_vd_vd_vd(xlog(tmp2), t));
        tmp1 = xsqrt(t);

        sigmaA2T2 = vmul_vd_vd_vd(vmul_vd_vd_vd(sigmaA, sigmaA), vdiv_vd_vd_vd(t, two));
        emrt = xexp(vsub_vd_vd_vd(vneg_vd_vd(vmul_vd_vd_vd(r, t)), ln_of_2));
        sigmaAsqrtT = vmul_vd_vd_vd(sigmaA, tmp1);
        tmp2 = xexp(vneg_vd_vd(vdiv_vd_vd_vd(sigmaA2T2, four)));
        tmp1 = vdiv_vd_vd_vd(vmul_vd_vd_vd(two, s), vmul_vd_vd_vd(pi2, sigmaA2T2));
        d2dx2_prep = vmul_vd_vd_vd(xsqrt(tmp1), vmul_vd_vd_vd(tmp2, emrt));

        vstore_v_p_vd(&sigmaA_[i], sigmaA);
        vstore_v_p_vd(&emrt_[i], emrt);
        vstore_v_p_vd(&sigmaA2T2_[i], sigmaA2T2);
        vstore_v_p_vd(&sigmaAsqrtT_[i], sigmaAsqrtT);
        vstore_v_p_vd(&d2dx2_prep_[i], d2dx2_prep);

    }

}

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
        Real_Ptr price_) {            /// [out] d2

    ASSUME(n % 64 == 0)

    ASSUME_ALIGNED(x_)
    ASSUME_ALIGNED(s_)
    ASSUME_ALIGNED(sigmaA2T2_)
    ASSUME_ALIGNED(sigmaAsqrtT_)
    ASSUME_ALIGNED(emrt_)
    ASSUME_ALIGNED(d1_)
    ASSUME_ALIGNED(d2_)
    ASSUME_ALIGNED(price_)
    ASSUME_ALIGNED(long_short_)
    ASSUME_ALIGNED(put_call_)

    vdouble tmp1, tmp2, tmp3, tmp4, x, s, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price;
    vdouble long_short, put_call;

    for (UINT64 i = 0; i < n; i += sizeof(vdouble) / sizeof(double)) {
        s = vload_vd_p(&s_[i]);
        x = vload_vd_p(&x_[i]);
        sigmaA2T2 = vload_vd_p(&sigmaA2T2_[i]);
        sigmaAsqrtT = vload_vd_p(&sigmaAsqrtT_[i]);
        emrt = vload_vd_p(&emrt_[i]);
        long_short = vload_vd_p(&long_short_[i]);
        put_call = vload_vd_p(&put_call_[i]);


        tmp2 = xlog(vdiv_vd_vd_vd(s, x));
        d1 = vdiv_vd_vd_vd(vadd_vd_vd_vd(tmp2, sigmaA2T2), vmul_vd_vd_vd(sigmaAsqrtT, msqrt2));
        d2 = vsub_vd_vd_vd(d1, vdiv_vd_vd_vd(sigmaAsqrtT, msqrt2));

        tmp3 = vmul_vd_vd_vd(put_call, d1);
        tmp4 = vmul_vd_vd_vd(put_call, d2);

        tmp1 = xerfc_u15(tmp3);
        tmp2 = xerfc_u15(tmp4);

        price = vsub_vd_vd_vd(
                vmul_vd_vd_vd(put_call, vmul_vd_vd_vd(s, tmp1)),
                vmul_vd_vd_vd(put_call, vmul_vd_vd_vd(x, tmp2)));

        price = vmul_vd_vd_vd(vmul_vd_vd_vd(emrt, long_short), price);

        vstore_v_p_vd(&d1_[i], d1);
        vstore_v_p_vd(&d2_[i], d2);
        vstore_v_p_vd(&price_[i], price);
    }


}