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


#ifndef PRICER_MKL_PRICER_H
#define PRICER_MKL_PRICER_H


#include <mkl.h>

#include <math/pricers/pricer-base.h>

#ifdef __cplusplus
#define restrict
extern "C" {
#endif

void init_mkl_pricer();

void prepare_mkl_pricer(
        UINT64 n,
        Real_Ptr s,                 /// [in] stock price
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
        Real_Ptr emrt,
        Real_Ptr d2dx2_prep);


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
        Real_Ptr price);

void ddx_mkl_pricer(
        UINT64 n,
        Uint64_Ptr flags,
        Real_Ptr d2,
        Real_Ptr emrt,
        Real_Ptr ddx_price);


void d2dx2_mkl_pricer(
        UINT64 n,
        Uint64_Ptr flags,
        Real_Ptr s,
        Real_Ptr x,
        Real_Ptr d2dx2_prep,
        Real_Ptr sigmaA2T2,
        Real_Ptr tmp1,
        Real_Ptr tmp2,
        Real_Ptr d2dx2);

#ifdef __cplusplus
};
#endif

#endif //PRICER_MKL_PRICER_H
