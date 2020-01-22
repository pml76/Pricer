/*
 *  (c) by Peter Lennartz
 * 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <context.h>
#include <cstdint>
#include <cstring>

#define MEM_REALLOC(type,var,old_n, new_n, p, v) {get_ ## var ## _max() = Private::reallocate_memory<type>(new_n,old_n,&p); for (uint64_t n = old_n; n < get_ ## var ## _max(); ++n){p[n] = v;}}
#define MEM_ALLOC(type,n,p,v) {get_ ## n ## _max() = Private::allocate_memory<type>(n, &p); if (get_ ## n ## _max() == 0) {return;} for(uint64_t i = 0; i < get_ ## n ## _max(); ++i) {p[i]=v;}}
#define MEM_DEALLOC(type, p) {Private::deallocate_memory<type>(p);};

namespace Pricer {

    void pricer_context::realloc_mem(uint64_t n_p) {
        uint64_t n_old = get_n_max();

        MEM_REALLOC(FLOAT, n, n_old, n_p, m__tradedate, 0.);
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__s, 70.);
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__sigma,0.3)
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__t,1.)
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__tau,1./12.)
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__r,0.01)
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__sigmaA,0.3)
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__sigmaA2T2,0.3*0.3)
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__sigmaAsqrtT,0.3)
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__emrt,1.)
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__long_short,1.)
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__put_call,1.)

        MEM_REALLOC(FLOAT, n, n_old, n_p, m__d1, 0.)
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__d2, 0.)
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__prices, 0.)
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__x, 72.)

    }


    void pricer_context::alloc_mem(uint64_t n) {

        MEM_ALLOC(FLOAT, n, m__tradedate, 0.)
        MEM_ALLOC(FLOAT, n, m__s, 70.)
        MEM_ALLOC(FLOAT, n, m__sigma, 0.3)
        MEM_ALLOC(FLOAT, n, m__t, 1.)
        MEM_ALLOC(FLOAT, n, m__tau, 1./12.)
        MEM_ALLOC(FLOAT, n, m__r,0.01)
        MEM_ALLOC(FLOAT, n, m__sigmaA,0.3)
        MEM_ALLOC(FLOAT, n, m__sigmaA2T2,0.3*0.3)
        MEM_ALLOC(FLOAT, n, m__sigmaAsqrtT, 0.3)
        MEM_ALLOC(FLOAT, n, m__emrt, 1.)
        MEM_ALLOC(FLOAT, n, m__long_short, 1.)
        MEM_ALLOC(FLOAT, n, m__put_call, 1.)
//        MEM_ALLOC(FLOAT, n, m__d2dx2_prep)

        MEM_ALLOC(FLOAT, n, m__d1, 0.)
        MEM_ALLOC(FLOAT, n, m__d2, 0.)
        MEM_ALLOC(FLOAT, n, m__prices, 0.)
        MEM_ALLOC(FLOAT, n, m__x, 72.)

    }

    void pricer_context::dealloc_mem() {
        MEM_DEALLOC(FLOAT, m__tradedate)
        MEM_DEALLOC(FLOAT, m__s)
        MEM_DEALLOC(FLOAT, m__sigma)
        MEM_DEALLOC(FLOAT, m__t)
        MEM_DEALLOC(FLOAT, m__tau)
        MEM_DEALLOC(FLOAT, m__r)
        MEM_DEALLOC(FLOAT, m__sigmaA)
        MEM_DEALLOC(FLOAT, m__sigmaA2T2)
        MEM_DEALLOC(FLOAT, m__sigmaAsqrtT)
        MEM_DEALLOC(FLOAT, m__emrt)
        MEM_DEALLOC(FLOAT, m__long_short)
        MEM_DEALLOC(FLOAT, m__put_call)
//        MEM_DEALLOC(FLOAT, m__d2dx2_prep)

        MEM_DEALLOC(FLOAT, m__d1)
        MEM_DEALLOC(FLOAT, m__d2)
        MEM_DEALLOC(FLOAT, m__prices)
        MEM_DEALLOC(FLOAT, m__x)

    }



    void ddx_pricer_context::alloc_mem(uint64_t n) {

        MEM_ALLOC(FLOAT, n, m__ddx_price, 0.)

    }

    void ddx_pricer_context::dealloc_mem() {
        MEM_DEALLOC(FLOAT, m__ddx_price)
    }

    void ddx_pricer_context::realloc_mem(uint64_t n_p) {
        uint64_t n_old = get_n_max();
        pricer_context::realloc_mem(n_p);
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__ddx_price, 0.);
    }
    



    void d2dx2_pricer_context::alloc_mem(uint64_t n) {

        MEM_ALLOC(FLOAT, n, m__d2dx2_price, 0.)
    }

    void d2dx2_pricer_context::dealloc_mem() {
        MEM_DEALLOC(FLOAT, m__d2dx2_price)
    }

    void d2dx2_pricer_context::realloc_mem(uint64_t n_p) {
        uint64_t n_old = get_n_max();
        
        ddx_pricer_context::realloc_mem(n_p);
        MEM_REALLOC(FLOAT, n, n_old, n_p, m__d2dx2_price, 0.);

    }





    void compute_prices_of_instruments_context::alloc_mem(uint64_t n, uint64_t m) {
        m++;
        MEM_ALLOC(FLOAT, m, m__instrument_prices, 0.)
        MEM_ALLOC(FLOAT, m, m__x_, 72.)
        MEM_ALLOC(int32_t , n, m__to_structure, get_m_max())
        MEM_ALLOC(FLOAT, n, m__offsets, 0.)
    }

    void compute_prices_of_instruments_context::dealloc_mem() {
        MEM_DEALLOC(int32_t , m__to_structure)
        MEM_DEALLOC(FLOAT, m__instrument_prices)
        MEM_DEALLOC(FLOAT, m__offsets)
        MEM_DEALLOC(FLOAT, m__x_)
    }


    void compute_prices_of_instruments_context::realloc_mem(uint64_t n_p, uint64_t m_p) {

        if( get_m_max() != m_p) {
            m_p++;
            uint64_t m_old = get_m_max();
            MEM_REALLOC(FLOAT, m, m_old, m_p, m__instrument_prices, 0.)
            MEM_REALLOC(FLOAT, m, m_old, m_p, m__x_, 72.)

            // update the unused fields of m__to_structure[].
            for(uint64_t i = get_n_act(); i < get_n_max(); i++) {
                m__to_structure[i] = get_m_max() - 1;
            }
        }

        if( get_n_max() != n_p) {
            uint64_t n_old = get_n_max();
            pricer_context::realloc_mem(n_p);
            MEM_REALLOC(int32_t, n, n_old, n_p, m__to_structure, get_m_max()-1)
            MEM_REALLOC(FLOAT, n, n_old, n_p, m__offsets, 0.)
        }

    }

    void compute_instrument_strikes_from_premiums_context::alloc_mem(uint64_t n, uint64_t m) {


        MEM_ALLOC(FLOAT, m, m__premiums, 0.)
        MEM_ALLOC(FLOAT, m, m__instrument_pricesl, 0.)
        MEM_ALLOC(FLOAT, m, m__instrument_pricesh, 0.)

        MEM_ALLOC(FLOAT, m, m__xl_, 0.)
        MEM_ALLOC(FLOAT, m, m__xh_, 0.)
    }

    void compute_instrument_strikes_from_premiums_context::realloc_mem(uint64_t n_p, uint64_t m_p) {
        uint64_t m_old = get_m_max();
        compute_prices_of_instruments_context::realloc_mem(n_p, m_p);

        if( get_m_max() != m_p) {
            m_p++;
            MEM_REALLOC(FLOAT, m, m_old, m_p, m__premiums, 0.)
            MEM_REALLOC(FLOAT, m, m_old, m_p, m__instrument_pricesl, 0.)
            MEM_REALLOC(FLOAT, m, m_old, m_p, m__instrument_pricesh, 0.)

            MEM_REALLOC(FLOAT, m, m_old, m_p, m__xl_, 0.)
            MEM_REALLOC(FLOAT, m, m_old, m_p, m__xh_, 0.)
        }


    }


    void compute_instrument_strikes_from_premiums_context::dealloc_mem() {

        MEM_DEALLOC(FLOAT, m__premiums)
        MEM_DEALLOC(FLOAT, m__instrument_pricesl)
        MEM_DEALLOC(FLOAT, m__instrument_pricesh)
        MEM_DEALLOC(FLOAT, m__xl_)
        MEM_DEALLOC(FLOAT, m__xh_)
    }




}