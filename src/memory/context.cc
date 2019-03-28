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


#define MEM_ALLOC(type,n,p) {m__ ## n ## _max = Private::allocate_memory<type>(n * sizeof(type), &p); if (m__ ## n ## _max == 0) return;}
#define COND_MEM_ALLOC(flag,type,n,p) {if((m__flags & flag) == flag) MEM_ALLOC(type,n,p) else p = nullptr;}

#define MEM_REALLOC(type,type2, p) {type2 tmp; m__n_max=Private::allocate_memory<type>(n*sizeof(type),&tmp); if(m__n_max==0); std::memcpy(tmp,p,m__n); p=tmp;}
#define COND_MEM_REALLOC(flag, type, type2, p) {if((m__flags & flag) == flag) MEM_REALLOC(type,type2,p)}

#define MEM_DEALLOC(type, p) {Private::deallocate_memory<type>(p);};
#define COND_MEM_DEALLOC(flags, type, p) {if((m__flags & flags) == flags) MEM_DEALLOC(type, p)}

namespace Pricer {
    void pricer_context::alloc_mem(uint64_t n, uint64_t m) {
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

        COND_MEM_ALLOC(PRICER_FLAG_TW_PRICER, FLOAT, n, m__d1)
        COND_MEM_ALLOC(PRICER_FLAG_TW_PRICER, FLOAT, n, m__d2)
        COND_MEM_ALLOC(PRICER_FLAG_TW_PRICER, FLOAT, n, m__prices)
        COND_MEM_ALLOC(PRICER_FLAG_TW_PRICER, FLOAT, n, m__x)

        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, int32_t, n, m__to_structure)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, n, m__offsets)



        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__premiums)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__instrument_prices)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__instrument_pricesl)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__instrument_pricesh)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__x_)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__xl_)
        COND_MEM_ALLOC(PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES, FLOAT, m, m__xh_)

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

        COND_MEM_DEALLOC(PRICER_FLAG_TW_PRICER, FLOAT, m__d1)
        COND_MEM_DEALLOC(PRICER_FLAG_TW_PRICER, FLOAT, m__d2)
        COND_MEM_DEALLOC(PRICER_FLAG_TW_PRICER, FLOAT, m__prices)
        COND_MEM_DEALLOC(PRICER_FLAG_TW_PRICER, FLOAT, m__x)

    }

    void pricer_context::realloc_mem(uint64_t n) {
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

        COND_MEM_REALLOC(PRICER_FLAG_TW_PRICER, FLOAT, Real_Ptr , m__d1)
        COND_MEM_REALLOC(PRICER_FLAG_TW_PRICER, FLOAT, Real_Ptr , m__d2)
        COND_MEM_REALLOC(PRICER_FLAG_TW_PRICER, FLOAT, Real_Ptr , m__prices)
        COND_MEM_REALLOC(PRICER_FLAG_TW_PRICER, FLOAT, Real_Ptr , m__x)

    }



}