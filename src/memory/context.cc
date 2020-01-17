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


#define MEM_ALLOC(type,n,p) {get_ ## n ## _max() = Private::allocate_memory<type>(n, &p); if (get_ ## n ## _max() == 0) return;}
#define MEM_DEALLOC(type, p) {Private::deallocate_memory<type>(p);};
#define MEM_INIT(x,i,v) {x[i]=v;}

namespace Pricer {
    void pricer_context::alloc_mem(uint64_t n) {


        MEM_ALLOC(FLOAT, n, m__s)
        MEM_ALLOC(FLOAT, n, m__sigma)
        MEM_ALLOC(FLOAT, n, m__t)
        MEM_ALLOC(FLOAT, n, m__tau)
        MEM_ALLOC(FLOAT, n, m__r)
        MEM_ALLOC(FLOAT, n, m__sigmaA)
        MEM_ALLOC(FLOAT, n, m__sigmaA2T2)
        MEM_ALLOC(FLOAT, n, m__sigmaAsqrtT)
        MEM_ALLOC(FLOAT, n, m__emrt)
        MEM_ALLOC(FLOAT, n, m__long_short)
        MEM_ALLOC(FLOAT, n, m__put_call)
        MEM_ALLOC(FLOAT, n, m__d2dx2_prep)

        MEM_ALLOC(FLOAT, n, m__d1)
        MEM_ALLOC(FLOAT, n, m__d2)
        MEM_ALLOC(FLOAT, n, m__prices)
        MEM_ALLOC(FLOAT, n, m__x)

        init_memory(0, m__n_max);
    }

    void ddx_pricer_context::alloc_mem(uint64_t n) {

        MEM_ALLOC(FLOAT, n, m__ddx_price)

        init_memory(0, get_n_max());
    }

    void d2dx2_pricer_context::alloc_mem(uint64_t n) {

        MEM_ALLOC(FLOAT, n, m__d2dx2_price)

        init_memory(0, get_n_max());

    }

    void compute_prices_of_instruments_context::alloc_mem(uint64_t n, uint64_t m) {
        MEM_ALLOC(int32_t , n, m__to_structure)
        MEM_ALLOC(FLOAT, m, m__instrument_prices)
        MEM_ALLOC(FLOAT, n, m__offsets)
        MEM_ALLOC(FLOAT, m, m__x_)

        init_memory(0, get_n_max(), 0, get_m_max());
    }

    void compute_instrument_strikes_from_premiums_context::alloc_mem(uint64_t n, uint64_t m) {


        MEM_ALLOC(FLOAT, m, m__premiums)
        MEM_ALLOC(FLOAT, m, m__instrument_pricesl)
        MEM_ALLOC(FLOAT, m, m__instrument_pricesh)

        MEM_ALLOC(FLOAT, m, m__xl_)
        MEM_ALLOC(FLOAT, m, m__xh_)

        init_memory(0, get_n_max(), 0, get_m_max());
    }


    void compute_instrument_strikes_from_premiums_context::dealloc_mem() {

        MEM_DEALLOC(FLOAT, m__premiums)
        MEM_DEALLOC(FLOAT, m__instrument_pricesl)
        MEM_DEALLOC(FLOAT, m__instrument_pricesh)
        MEM_DEALLOC(FLOAT, m__xl_)
        MEM_DEALLOC(FLOAT, m__xh_)
    }

    void ddx_pricer_context::dealloc_mem() {
        MEM_DEALLOC(FLOAT, m__ddx_price)
    }

    void d2dx2_pricer_context::dealloc_mem() {
        MEM_DEALLOC(FLOAT, m__d2dx2_price)
    }

    void compute_prices_of_instruments_context::dealloc_mem() {
        MEM_DEALLOC(int32_t , m__to_structure)
        MEM_DEALLOC(FLOAT, m__instrument_prices)
        MEM_DEALLOC(FLOAT, m__offsets)
        MEM_DEALLOC(FLOAT, m__x_)
    }


    void pricer_context::dealloc_mem() {
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
        MEM_DEALLOC(FLOAT, m__d2dx2_prep)

        MEM_DEALLOC(FLOAT, m__d1)
        MEM_DEALLOC(FLOAT, m__d2)
        MEM_DEALLOC(FLOAT, m__prices)
        MEM_DEALLOC(FLOAT, m__x)

    }


    void ddx_pricer_context::init_memory(uint64_t n1, uint64_t n2) {
        pricer_context::init_memory(n1, n2);

        for (uint64_t n = n1; n < n2; ++n) {
            MEM_INIT(m__ddx_price, n, 0.)
        }
    }

    void d2dx2_pricer_context::init_memory(uint64_t n1, uint64_t n2) {
        ddx_pricer_context::init_memory(n1,n2);

        for (uint64_t n = n1; n < n2; ++n) {
            MEM_INIT(m__d2dx2_price, n, 0.)
        }
    }

    void compute_prices_of_instruments_context::init_memory(uint64_t n1, uint64_t n2, uint64_t m1, uint64_t m2) {
        pricer_context::init_memory(n1, n2);

        for(uint64_t i = n1; i < n2; ++i) {
            MEM_INIT(m__to_structure, i, m2-1)
            MEM_INIT(m__offsets, i, 0.)
        }

        for(uint64_t i = m1; i < m2; ++i) {
            MEM_INIT(m__instrument_prices, i, 0.)
            MEM_INIT(m__x_, i, 0.)
        }

    }

    void compute_instrument_strikes_from_premiums_context::init_memory(uint64_t n1, uint64_t n2, uint64_t m1, uint64_t m2) {

        compute_prices_of_instruments_context::init_memory(n1, n2, m1, m2);
        d2dx2_pricer_context::init_memory(n1, n2);


        for(uint64_t i = m1; i < m2; ++i) {
            MEM_INIT(m__premiums, i, 0.)
            MEM_INIT(m__instrument_pricesl, i, 0.)
            MEM_INIT(m__instrument_pricesh, i, 0.)
            MEM_INIT(m__xl_, i, 0.)
            MEM_INIT(m__xh_, i, 0.)
        }
    }

    void pricer_context::init_memory(uint64_t n1, uint64_t n2) {

        for(uint64_t n = n1; n < n2; ++n) {
            MEM_INIT(m__s,n,70.)
            MEM_INIT(m__sigma,n,0.3)
            MEM_INIT(m__t,n,1.)
            MEM_INIT(m__tau,n,1./12.)
            MEM_INIT(m__r,n,0.01)
            MEM_INIT(m__sigmaA,n,0.3)
            MEM_INIT(m__sigmaA2T2,n,0.3*0.3)
            MEM_INIT(m__sigmaAsqrtT,n,0.3)
            MEM_INIT(m__emrt,n,1.)
            MEM_INIT(m__long_short,n,1.)
            MEM_INIT(m__put_call,n,1.)
            MEM_INIT(m__d2dx2_prep,n,0)
        }

    }

}