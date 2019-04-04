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

#include <memory/context.h>
#include <iostream>

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

#include <math/pricers/pricer-base.h>
#include <pricer-renamer.h>
#include <tgmath.h>

#include <omp.h>

static vdouble ln_of_2, msqrt2;
vdouble one, two, four, pi2, mone, zero;


extern "C" {

    vdouble xexp(vdouble d);

    vdouble xlog(vdouble d);

    vdouble xsqrt(vdouble d);

    vdouble xerfc_u15(vdouble a);
}

void init_tw_pricer() {
    ln_of_2 = vcast_vd_d(log(2.));
    msqrt2 = vcast_vd_d(-sqrt(2));
    one = vcast_vd_d(1.);
    two = vcast_vd_d(2.);
    four = vcast_vd_d(4.);
    pi2 = vmul_vd_vd_vd(two, vcast_vd_d(M_PI));
    mone = vcast_vd_d(-1.);
    zero = vcast_vd_d(0.);

}


void prepare_tw_pricer( Pricer::pricer_context &context ) {

/*
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
*/

    ASSUME(context.get_n_max() % 64 == 0)

    ASSUME_ALIGNED(Real_Ptr ,context.get_s())
    ASSUME_ALIGNED(Real_Ptr ,context.get_sigma())
    ASSUME_ALIGNED(Real_Ptr ,context.get_t())
    ASSUME_ALIGNED(Real_Ptr ,context.get_tau())
    ASSUME_ALIGNED(Real_Ptr ,context.get_r())
    ASSUME_ALIGNED(Real_Ptr ,context.get_sigmaA())
    ASSUME_ALIGNED(Real_Ptr ,context.get_sigmaA2T2())
    ASSUME_ALIGNED(Real_Ptr ,context.get_sigmaAsqrtT())
    ASSUME_ALIGNED(Real_Ptr ,context.get_emrt())
    ASSUME_ALIGNED(Real_Ptr ,context.get_d2dx2_prep())



#pragma omp parallel
    {
        vdouble tmp1, tmp2, tmp3, tmp4, tmp5;
        vdouble tt1, tt3;
        vdouble s, sigma, t, tau, r, sigmaA, sigmaA2T2, sigmaAsqrtT, emrt, d2dx2_prep;

        uint64_t n2 = context.get_n_max() / (64 / sizeof(double));
        uint64_t tid = omp_get_thread_num();
        uint64_t num_threads = omp_get_num_threads();
        uint64_t begin = ((tid * n2) / num_threads) * (64 / sizeof(double));
        uint64_t end = (((tid + 1) * n2) / num_threads) * (64 / sizeof(double));

        for (uint64_t i = begin; i < end; i += sizeof(vdouble) / sizeof(double)) {

            s = vload_vd_p(&context.get_s()[i]);
            sigma = vload_vd_p(&context.get_sigma()[i]);
            t = vload_vd_p(&context.get_t()[i]);
            tau = vload_vd_p(&context.get_tau()[i]);
            r = vload_vd_p(&context.get_r()[i]);

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

            vstore_v_p_vd(&context.get_sigmaA()[i], sigmaA);
            vstore_v_p_vd(&context.get_emrt()[i], emrt);
            vstore_v_p_vd(&context.get_sigmaA2T2()[i], sigmaA2T2);
            vstore_v_p_vd(&context.get_sigmaAsqrtT()[i], sigmaAsqrtT);
            vstore_v_p_vd(&context.get_d2dx2_prep()[i], d2dx2_prep);

        }
    }

}

void tw_pricer( Pricer::pricer_context &context ) {

/*        UINT64 n,
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
*/

    if( (context.get_flags() & PRICER_FLAG_TW_PRICER ) != PRICER_FLAG_TW_PRICER ) {
        std::cerr << "tw_pricer() called with improper context." << std::endl;
        exit(-1);
    }


    ASSUME(context.get_n_max() % 64 == 0)

    ASSUME_ALIGNED(Real_Ptr ,context.get_x())
    ASSUME_ALIGNED(Real_Ptr ,context.get_s())
    ASSUME_ALIGNED(Real_Ptr ,context.get_sigmaA2T2())
    ASSUME_ALIGNED(Real_Ptr ,context.get_sigmaAsqrtT())
    ASSUME_ALIGNED(Real_Ptr ,context.get_emrt())
    ASSUME_ALIGNED(Real_Ptr ,context.get_d1())
    ASSUME_ALIGNED(Real_Ptr ,context.get_d2())
    ASSUME_ALIGNED(Real_Ptr ,context.get_prices())
    ASSUME_ALIGNED(Real_Ptr ,context.get_long_short())
    ASSUME_ALIGNED(Real_Ptr ,context.get_put_call())

// #pragma omp parallel
    {
        vdouble tmp1, tmp2, tmp3, tmp4, x, s, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price;
        vdouble long_short, put_call;

        uint64_t n2 = context.get_n_max() / (64 / sizeof(double));
        uint64_t tid = omp_get_thread_num();
        uint64_t num_threads = omp_get_num_threads();
        uint64_t begin = ((tid * n2) / num_threads) * (64 / sizeof(double));
        uint64_t end = (((tid + 1) * n2) / num_threads) * (64 / sizeof(double));

        for (uint64_t i = begin; i < end; i += sizeof(vdouble) / sizeof(double)) {

            s = vload_vd_p(&context.get_s()[i]);
            x = vload_vd_p(&context.get_x()[i]);
            sigmaA2T2 = vload_vd_p(&context.get_sigmaA2T2()[i]);
            sigmaAsqrtT = vload_vd_p(&context.get_sigmaAsqrtT()[i]);
            emrt = vload_vd_p(&context.get_emrt()[i]);
            long_short = vload_vd_p(&context.get_long_short()[i]);
            put_call = vload_vd_p(&context.get_put_call()[i]);


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

            vstore_v_p_vd(&context.get_d1()[i], d1);
            vstore_v_p_vd(&context.get_d2()[i], d2);
            vstore_v_p_vd(&context.get_prices()[i], price);
        }
    }


}


void ddx_tw_pricer( Pricer::pricer_context &context ) {

/*        UINT64 n,
        Real_Ptr long_short_,     // 1 == long option // -1 == short option
        Real_Ptr put_call_,       // -1 == put // 1 == call
        Real_Ptr d2_,
        Real_Ptr emrt_,
        Real_Ptr ddx_price_) {
  */

    if((context.get_flags() & PRICER_FLAG_TW_COMPUTE_DDX) != PRICER_FLAG_TW_COMPUTE_DDX ) {
        std::cerr << "ddx_tw_pricer() called with improper context." << std::endl;
        exit(-1);
    }

    ASSUME(context.get_n_max() % 64 == 0)

    ASSUME_ALIGNED(Real_Ptr ,context.get_long_short())
    ASSUME_ALIGNED(Real_Ptr ,context.get_put_call())
    ASSUME_ALIGNED(Real_Ptr ,context.get_d2())
    ASSUME_ALIGNED(Real_Ptr ,context.get_emrt())
    ASSUME_ALIGNED(Real_Ptr ,context.get_ddx_price())

#pragma omp parallel
    {
        vdouble long_short, put_call, d2, emrt, ddx_price;

        uint64_t n2 = context.get_n_max() / (64 / sizeof(double));
        uint64_t tid = omp_get_thread_num();
        uint64_t num_threads = omp_get_num_threads();
        uint64_t begin = ((tid * n2) / num_threads) * (64 / sizeof(double));
        uint64_t end = (((tid + 1) * n2) / num_threads) * (64 / sizeof(double));

        for (uint64_t i = begin; i < end; i += sizeof(vdouble) / sizeof(double)) {

            d2 = vload_vd_p(&context.get_d2()[i]);
            put_call = vload_vd_p(&context.get_put_call()[i]);
            long_short = vload_vd_p(&context.get_long_short()[i]);
            emrt = vload_vd_p(&context.get_emrt()[i]);

            ddx_price = xerfc_u15(d2);
            ddx_price = vsub_vd_vd_vd(vadd_vd_vd_vd(one, vmul_vd_vd_vd(put_call, mone)), ddx_price);
            ddx_price = vmul_vd_vd_vd(vmul_vd_vd_vd(emrt, long_short), ddx_price);

            vstore_v_p_vd(&context.get_ddx_price()[i], ddx_price);
        }
    }

}


void d2dx2_tw_pricer( Pricer::pricer_context &context ) {

/*        UINT64 n,
        Real_Ptr long_short_,
        Real_Ptr s_,
        Real_Ptr x_,
        Real_Ptr d2dx2_prep_,
        Real_Ptr sigmaA2T2_,
        Real_Ptr d2dx2_
) {

 */
    if((context.get_flags() & PRICER_FLAG_TW_COMPUTE_D2DX2) != PRICER_FLAG_TW_COMPUTE_D2DX2) {
        std::cerr << "d2dx2_tw_pricer() called with improper context." << std::endl;
        exit(-1);
    }


    ASSUME(context.get_n_max() % 64 == 0)

    ASSUME_ALIGNED(Real_Ptr ,context.get_long_short())
    ASSUME_ALIGNED(Real_Ptr ,context.get_s())
    ASSUME_ALIGNED(Real_Ptr ,context.get_x())
    ASSUME_ALIGNED(Real_Ptr ,context.get_d2dx2_prep())
    ASSUME_ALIGNED(Real_Ptr ,context.get_sigmaA2T2())
    ASSUME_ALIGNED(Real_Ptr ,context.get_d2dx2())

#pragma omp parallel
    {
        vdouble long_short, s, x, d2dx2_prep, sigmaA2T2, d2dx2;

        uint64_t n2 = context.get_n_max() / (64 / sizeof(double));
        uint64_t tid = omp_get_thread_num();
        uint64_t num_threads = omp_get_num_threads();
        uint64_t begin = ((tid * n2) / num_threads) * (64 / sizeof(double));
        uint64_t end = (((tid + 1) * n2) / num_threads) * (64 / sizeof(double));

        for (uint64_t i = begin; i < end; i += sizeof(vdouble) / sizeof(double)) {

            long_short = vload_vd_p(&context.get_long_short()[i]);
            s = vload_vd_p(&context.get_s()[i]);
            x = vload_vd_p(&context.get_x()[i]);
            d2dx2_prep = vload_vd_p(&context.get_d2dx2_prep()[i]);
            sigmaA2T2 = vload_vd_p(&context.get_sigmaA2T2()[i]);

            d2dx2 = xlog(vdiv_vd_vd_vd(s, x));
            d2dx2 = xexp(vneg_vd_vd(vdiv_vd_vd_vd(vmul_vd_vd_vd(d2dx2, d2dx2), vmul_vd_vd_vd(four, sigmaA2T2))));
            d2dx2 = vdiv_vd_vd_vd(vmul_vd_vd_vd(d2dx2, d2dx2_prep), xsqrt(vmul_vd_vd_vd(vmul_vd_vd_vd(x, x), x)));
            d2dx2 = vmul_vd_vd_vd(d2dx2, long_short);

            vstore_v_p_vd(&context.get_d2dx2()[i], d2dx2);
        }
    }

}

void full_tw_pricer( Pricer::pricer_context &context ) {

/*        UINT64 n,
        Real_Ptr long_short_,     // 1 == long option // -1 == short option
        Real_Ptr put_call_,       // -1 == put // 1 == call
        Real_Ptr s_,                /// [in] stock price
        Real_Ptr x_,                /// [in] strike
        Real_Ptr sigmaA2T2_,        /// [in] sigmaA^2t/2
        Real_Ptr sigmaAsqrtT_,      /// [in] sigmaA*sqrt(t)
        Real_Ptr emrt_,
        Real_Ptr d2dx2_prep_,
        Real_Ptr price_,
        Real_Ptr ddx_price_,
        Real_Ptr d2dx2_) {
*/

    if((context.get_flags() & (PRICER_FLAG_TW_PRICER | PRICER_FLAG_TW_COMPUTE_DDX | PRICER_FLAG_TW_COMPUTE_D2DX2) ) !=
            (PRICER_FLAG_TW_PRICER | PRICER_FLAG_TW_COMPUTE_DDX | PRICER_FLAG_TW_COMPUTE_D2DX2) ) {
        std::cerr << "full_tw_pricer() called with improper context." << std::endl;
        exit(-1);
    }

    ASSUME(context.get_n_max() % 64 == 0)

    ASSUME_ALIGNED(Real_Ptr ,context.get_x())
    ASSUME_ALIGNED(Real_Ptr ,context.get_s())
    ASSUME_ALIGNED(Real_Ptr ,context.get_sigmaA2T2())
    ASSUME_ALIGNED(Real_Ptr ,context.get_sigmaAsqrtT())
    ASSUME_ALIGNED(Real_Ptr ,context.get_emrt())
    ASSUME_ALIGNED(Real_Ptr ,context.get_prices())
    ASSUME_ALIGNED(Real_Ptr ,context.get_long_short())
    ASSUME_ALIGNED(Real_Ptr ,context.get_put_call())
    ASSUME_ALIGNED(Real_Ptr ,context.get_ddx_price())
    ASSUME_ALIGNED(Real_Ptr ,context.get_d2dx2())
    ASSUME_ALIGNED(Real_Ptr ,context.get_d2dx2_prep())

#pragma omp parallel
    {

        vdouble tmp1, tmp2, tmp3, tmp4, x, s, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price;
        vdouble long_short, put_call;
        vdouble ddx_price;
        vdouble d2dx2, d2dx2_prep;

        uint64_t n2 = context.get_n_max() / (64 / sizeof(double));
        uint64_t tid = omp_get_thread_num();
        uint64_t num_threads = omp_get_num_threads();
        uint64_t begin = ((tid * n2) / num_threads) * (64 / sizeof(double));
        uint64_t end = (((tid + 1) * n2) / num_threads) * (64 / sizeof(double));

        for (uint64_t i = begin; i < end; i += sizeof(vdouble) / sizeof(double)) {

            s = vload_vd_p(&context.get_s()[i]);
            x = vload_vd_p(&context.get_x()[i]);
            sigmaA2T2 = vload_vd_p(&context.get_sigmaA2T2()[i]);
            sigmaAsqrtT = vload_vd_p(&context.get_sigmaAsqrtT()[i]);
            emrt = vload_vd_p(&context.get_emrt()[i]);
            long_short = vload_vd_p(&context.get_long_short()[i]);
            put_call = vload_vd_p(&context.get_put_call()[i]);
            d2dx2_prep = vload_vd_p(&context.get_d2dx2_prep()[i]);

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

            ddx_price = xerfc_u15(d2);
            ddx_price = vsub_vd_vd_vd(vadd_vd_vd_vd(one, vmul_vd_vd_vd(put_call, mone)), ddx_price);
            ddx_price = vmul_vd_vd_vd(vmul_vd_vd_vd(emrt, long_short), ddx_price);


            d2dx2 = xlog(vdiv_vd_vd_vd(s, x));
            d2dx2 = xexp(vneg_vd_vd(vdiv_vd_vd_vd(vmul_vd_vd_vd(d2dx2, d2dx2), vmul_vd_vd_vd(four, sigmaA2T2))));
            d2dx2 = vdiv_vd_vd_vd(vmul_vd_vd_vd(d2dx2, d2dx2_prep), xsqrt(vmul_vd_vd_vd(vmul_vd_vd_vd(x, x), x)));
            d2dx2 = vmul_vd_vd_vd(d2dx2, long_short);

            vstore_v_p_vd(&context.get_d2dx2()[i], d2dx2);
            vstore_v_p_vd(&context.get_ddx_price()[i], ddx_price);
            vstore_v_p_vd(&context.get_prices()[i], price);

        }
    }

}


void compute_tw_strikes_from_premiums( Pricer::pricer_context &context ) {

    if( (context.get_flags() & PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES)
                != PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES ) {
        std::cerr << "compute_tw_strices_from_premiums() called with wrong initialized context" << std::endl;
        exit(-1);
    }

/*
        UINT64 n,
        Real_Ptr long_short_,        // 1 == long option // -1 == short option
        Real_Ptr put_call_,          // -1 == put // 1 == call
        Real_Ptr s_,
        Real_Ptr sigmaA2T2_,
        Real_Ptr sigmaAsqrtT_,
        Real_Ptr emrt_,
        Int32_Ptr to_structure,
        Real_Ptr offsets_,
        Real_Ptr prices_,
        Real_Ptr x_tmp_,
        uint64_t m,
        Real_Ptr premiums_,
        Real_Ptr instrument_prices_,
        Real_Ptr instrument_pricesl_,
        Real_Ptr instrument_pricesh_,
        Real_Ptr x_,
        Real_Ptr xl_,
        Real_Ptr xh_) {
*/

    ASSUME(context.get_n_max() % 64 == 0)
    ASSUME(context.get_m_max() % 64 == 0)

    ASSUME_ALIGNED(Real_Ptr ,context.get_long_short())
    ASSUME_ALIGNED(Real_Ptr ,context.get_put_call())
    ASSUME_ALIGNED(Real_Ptr ,context.get_s())
    ASSUME_ALIGNED(Real_Ptr ,context.get_sigmaA2T2())
    ASSUME_ALIGNED(Real_Ptr ,context.get_sigmaAsqrtT())
    ASSUME_ALIGNED(Real_Ptr ,context.get_emrt())
    ASSUME_ALIGNED(Int32_Ptr,context.get_to_structure())
    ASSUME_ALIGNED(Real_Ptr ,context.get_offsets())
    ASSUME_ALIGNED(Real_Ptr ,context.get_prices())
    ASSUME_ALIGNED(Real_Ptr ,context.get_x())

    ASSUME_ALIGNED(Real_Ptr ,context.get_premiums())
    ASSUME_ALIGNED(Real_Ptr ,context.get_instrument_prices())
    ASSUME_ALIGNED(Real_Ptr ,context.get_instrument_pricesh())
    ASSUME_ALIGNED(Real_Ptr ,context.get_instrument_pricesl())
    ASSUME_ALIGNED(Real_Ptr ,context.get_x_())
    ASSUME_ALIGNED(Real_Ptr ,context.get_xh_())
    ASSUME_ALIGNED(Real_Ptr ,context.get_xl_())

    vdouble err;
    double buffer[sizeof(vdouble)/ sizeof(double)] __attribute__((aligned(ALIGN_TO)));

    // fill the instrument_priceh and instrument_pricel arrays with the prices evaluated at xh_ and xl_,
    // respectively.
#pragma omp parallel
    {
        vdouble tmp, tmp1, tmp2, tmp3, tmp4, x, s, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price;
        vdouble long_short, put_call;


        uint64_t m2 =context.get_m_max() / (64 / sizeof(double));
        uint64_t tid = omp_get_thread_num();
        uint64_t num_threads = omp_get_num_threads();
        uint64_t m_begin = ((tid * m2) / num_threads) * (64 / sizeof(double));
        uint64_t m_end = (((tid + 1) * m2) / num_threads) * (64 / sizeof(double));
        uint64_t n2 = context.get_n_max() / (64 / sizeof(double));
        uint64_t n_begin = ((tid * n2) / num_threads) * (64 / sizeof(double));
        uint64_t n_end = (((tid + 1) * n2) / num_threads) * (64 / sizeof(double));

        ///
        /// set the variables to zero where we are going to aggregate the numbers
        ///
        for (uint64_t i = m_begin; i < m_end; i += sizeof(vdouble) / sizeof(double)) {
            tmp = vneg_vd_vd(vload_vd_p(&context.get_premiums()[i]));
            vstore_v_p_vd(&context.get_instrument_pricesh()[i], tmp);
            vstore_v_p_vd(&context.get_instrument_pricesl()[i], tmp);
        }

#pragma omp barrier

        ///
        /// aggregate the prices and its first two derivatives
        ///
        double d;

#pragma omp for schedule(static)
        for (uint64_t i = 0; i < context.get_n_max(); ++i) {
#pragma omp atomic update
            context.get_instrument_prices()[context.get_to_structure()[i]] += context.get_prices()[i];
        }

#pragma omp barrier

        for (uint64_t i = n_begin; i < n_end; i += sizeof(vdouble) / sizeof(double)) {


            tmp = vload_vd_p(&context.get_offsets()[i]);

            s = vload_vd_p(&context.get_s()[i]);
            x = vadd_vd_vd_vd(vload_vd_p(&context.get_xh_()[i]),tmp);
            sigmaA2T2 = vload_vd_p(&context.get_sigmaA2T2()[i]);
            sigmaAsqrtT = vload_vd_p(&context.get_sigmaAsqrtT()[i]);
            emrt = vload_vd_p(&context.get_emrt()[i]);
            long_short = vload_vd_p(&context.get_long_short()[i]);
            put_call = vload_vd_p(&context.get_put_call()[i]);


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

            vstore_v_p_vd(&context.get_prices()[i], price);



            x = vadd_vd_vd_vd(vload_vd_p(&context.get_xl_()[i]),tmp);

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

            vstore_v_p_vd(&context.get_instrument_pricesl()[i], price);

        }


    }


    /*
     *  Main loop starts here.
     */

    do {
        err = zero;
#pragma omp parallel
        {

            vdouble err1 = zero;
            vdouble tmp, tmp1, tmp2, tmp3, tmp4, x, xl, xh, s, sigmaA2T2, sigmaAsqrtT, emrt, d1, d2, price;
            vdouble long_short, put_call;
            vopmask op;
            vdouble pricel, priceh;

            uint64_t m2 =context.get_m_max() / (64 / sizeof(double));
            uint64_t tid = omp_get_thread_num();
            uint64_t num_threads = omp_get_num_threads();
            uint64_t m_begin = ((tid * m2) / num_threads) * (64 / sizeof(double));
            uint64_t m_end = (((tid + 1) * m2) / num_threads) * (64 / sizeof(double));
            uint64_t n2 = context.get_n_max() / (64 / sizeof(double));
            uint64_t n_begin = ((tid * n2) / num_threads) * (64 / sizeof(double));
            uint64_t n_end = (((tid + 1) * n2) / num_threads) * (64 / sizeof(double));


            for (uint64_t i = n_begin; i < n_end; i += sizeof(vdouble) / sizeof(double)) {

                xh = vadd_vd_vd_vd( vgather_vd_p_vi(context.get_xh_(),vloadu_vi_p(&context.get_to_structure()[i])), vload_vd_p(&context.get_offsets()[i]));
                xl = vadd_vd_vd_vd( vgather_vd_p_vi(context.get_xl_(),vloadu_vi_p(&context.get_to_structure()[i])), vload_vd_p(&context.get_offsets()[i]));
                x = vdiv_vd_vd_vd(vadd_vd_vd_vd(xl,xh),two);

                s = vload_vd_p(&context.get_s()[i]);

                sigmaA2T2 = vload_vd_p(&context.get_sigmaA2T2()[i]);
                sigmaAsqrtT = vload_vd_p(&context.get_sigmaAsqrtT()[i]);
                emrt = vload_vd_p(&context.get_emrt()[i]);
                long_short = vload_vd_p(&context.get_long_short()[i]);
                put_call = vload_vd_p(&context.get_put_call()[i]);

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

                vstore_v_p_vd(&context.get_prices()[i], price);
                vstore_v_p_vd(&context.get_x()[i],x);
            }




            ///
            /// set the variables to zero where we are going to aggregate the numbers
            ///
            for (uint64_t i = m_begin; i < m_end; i += sizeof(vdouble) / sizeof(double)) {
                tmp = vload_vd_p(&context.get_premiums()[i]);
                vstore_v_p_vd(&context.get_instrument_prices()[i], vneg_vd_vd(tmp));
            }

#pragma omp barrier

            ///
            /// aggregate the prices and its first two derivatives
            ///
            double d;

#pragma omp for schedule(static)
            for (uint64_t i = 0; i < context.get_n_max(); ++i) {

                d = context.get_x()[i]-context.get_offsets()[i];

#pragma omp atomic update
                context.get_instrument_prices()[context.get_to_structure()[i]] += context.get_prices()[i];

#pragma omp atomic write
                context.get_x_()[context.get_to_structure()[i]] = d;

            }

#pragma omp barrier

            for (uint64_t i = m_begin; i < m_end; i += sizeof(vdouble) / sizeof(double)) {

                price  = vload_vd_p(&context.get_instrument_prices()[i]);
                pricel = vload_vd_p(&context.get_instrument_pricesl()[i]);
                priceh = vload_vd_p(&context.get_instrument_pricesh()[i]);

                x  = vload_vd_p(&context.get_x_()[i]);
                xl = vload_vd_p(&context.get_xl_()[i]);
                xh = vload_vd_p(&context.get_xh_()[i]);

                op     = vgt_vo_vd_vd(vmul_vd_vd_vd(price, pricel), zero);

                pricel = vsel_vd_vo_vd_vd(op, price, pricel);
                priceh = vsel_vd_vo_vd_vd(op, priceh, price);

                xl     = vsel_vd_vo_vd_vd(op, x, xl);
                xh     = vsel_vd_vo_vd_vd(op, xh, x);

                vstore_v_p_vd(&context.get_instrument_pricesh()[i], priceh);
                vstore_v_p_vd(&context.get_instrument_pricesl()[i], pricel);
                vstore_v_p_vd(&context.get_xl_()[i], xl);
                vstore_v_p_vd(&context.get_xh_()[i], xh);

                tmp = vsub_vd_vd_vd(vmax_vd_vd_vd(xh,xl),vmin_vd_vd_vd(xh,xl));
                err1 = vadd_vd_vd_vd(err1,vmul_vd_vd_vd(tmp,tmp));

            }

#pragma omp critical
            {
                err = vadd_vd_vd_vd(err, err1);
            }
        }

        vstore_v_p_vd(buffer, err);
        for(uint64_t i = 1; i < sizeof(vdouble)/ sizeof(double); ++i) {
            buffer[0] += buffer[i];
        }


    } while(sqrt(buffer[0]) > 1.0e-5);


}