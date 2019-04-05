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

#include <memory/context.h>
#include <cstdint>
#include <cstring>


#define MEM_ALLOC(type,n,p) {get_ ## n ## _max() = Private::allocate_memory<type>(n, &p); if (get_ ## n ## _max() == 0) return;}
#define MEM_REALLOC(type,type2, p) {type2 tmp; m__n_max=Private::allocate_memory<type>(n,&tmp); if(m__n_max==0); std::memcpy(tmp,p,m__n); p=tmp;}
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

        /*
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, int32_t, n, m__to_structure)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, n, m__offsets)

        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__premiums)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__instrument_prices)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__instrument_pricesl)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__instrument_pricesh)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__x_)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__xl_)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__xh_)

        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_DDX, FLOAT, n, m__ddx_price)

        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_D2DX2, FLOAT, n, m__d2dx2)
*/
        init_memory(0, m__n_max);
    }

    void ddx_pricer_context::alloc_mem(uint64_t n) {

        MEM_ALLOC(FLOAT, n, m__ddx_price)
    }

    void d2dx2_pricer_context::alloc_mem(uint64_t n) {

        MEM_ALLOC(FLOAT, n, m__d2dx2_price)

    }

    void ddx_pricer_context::dealloc_mem() {
        MEM_DEALLOC(FLOAT, m__ddx_price)
    }

    void d2dx2_pricer_context::dealloc_mem() {
        MEM_DEALLOC(FLOAT, m__d2dx2_price)
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

        /*
        COND_MEM_DEALLOC(PRICER_FLAG_TW_COMPUTE_DDX, FLOAT, m__ddx_price)

        COND_MEM_DEALLOC(PRICER_FLAG_TW_COMPUTE_D2DX2, FLOAT, m__d2dx2)
        */
    }

    void ddx_pricer_context::realloc_mem(uint64_t n) {

        pricer_context::realloc_mem(n);

        uint64_t n_max_old = get_n_max();

        MEM_REALLOC(FLOAT, Real_Ptr, m__ddx_price)

        init_memory( n_max_old, get_n_max());

    }

    void d2dx2_pricer_context::realloc_mem(uint64_t n) {

        ddx_pricer_context::realloc_mem(n);

        uint64_t n_max_old = get_n_max();

        MEM_REALLOC(FLOAT, Real_Ptr, m__d2dx2_price)

        init_memory( n_max_old, get_n_max());

    }

    void pricer_context::realloc_mem(uint64_t n) {
        uint64_t n_max_old = m__n_max;
     //   uint64_t m_max_old = m__m_max;

        MEM_REALLOC(FLOAT, Real_Ptr , m__s)
        MEM_REALLOC(FLOAT, Real_Ptr , m__sigma)
        MEM_REALLOC(FLOAT, Real_Ptr , m__t)
        MEM_REALLOC(FLOAT, Real_Ptr , m__tau)
        MEM_REALLOC(FLOAT, Real_Ptr , m__r)
        MEM_REALLOC(FLOAT, Real_Ptr , m__sigmaA)
        MEM_REALLOC(FLOAT, Real_Ptr , m__sigmaA2T2)
        MEM_REALLOC(FLOAT, Real_Ptr , m__sigmaAsqrtT)
        MEM_REALLOC(FLOAT, Real_Ptr , m__emrt)
        MEM_REALLOC(FLOAT, Real_Ptr , m__long_short)
        MEM_REALLOC(FLOAT, Real_Ptr , m__put_call)
        MEM_REALLOC(FLOAT, Real_Ptr , m__d2dx2_prep)

        MEM_REALLOC(FLOAT, Real_Ptr , m__d1)
        MEM_REALLOC(FLOAT, Real_Ptr , m__d2)
        MEM_REALLOC(FLOAT, Real_Ptr , m__prices)
        MEM_REALLOC(FLOAT, Real_Ptr , m__x)

        /*
        COND_MEM_REALLOC(PRICER_FLAG_TW_COMPUTE_DDX, FLOAT, Real_Ptr, m__ddx_price)

        COND_MEM_REALLOC(PRICER_FLAG_TW_COMPUTE_D2DX2, FLOAT, Real_Ptr, m__d2dx2)
*/
        init_memory(n_max_old, m__n_max);
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

            /*
            COND_MEM_INIT(PRICER_FLAG_TW_PRICER, m__d1,n,1)
            COND_MEM_INIT(PRICER_FLAG_TW_PRICER, m__d2,n,1)
            COND_MEM_INIT(PRICER_FLAG_TW_PRICER, m__prices,n,1)
            COND_MEM_INIT(PRICER_FLAG_TW_PRICER, m__x,n,72.)

            COND_MEM_INIT(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, m__to_structure,n,get_m_max()-1)
            COND_MEM_INIT(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, m__offsets,n,0.)

            COND_MEM_INIT(PRICER_FLAG_TW_COMPUTE_DDX, m__ddx_price, n, 0.)

            COND_MEM_INIT(PRICER_FLAG_TW_COMPUTE_D2DX2, m__d2dx2, n, 0.)
             */
        }

        /*
        for(uint64_t m = m1; m < m2; ++m) {
           COND_MEM_INIT(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, m__premiums, m, 2.)
           COND_MEM_INIT(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, m__instrument_prices, m, 0.)
           COND_MEM_INIT(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, m__instrument_pricesl, m, 0.)
           COND_MEM_INIT(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, m__instrument_pricesh, m, 0.)
           COND_MEM_INIT(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, m__x_, m, 0.)
           COND_MEM_INIT(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, m__xl_, m, 0.)
           COND_MEM_INIT(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, m__xh_, m, 0.)
        }
    */
    }

}