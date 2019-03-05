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

#include <memory.h>
#include <cstdint>

#define MEM_ALLOC(type,p) {m__n_max = Private::allocate_memory<FLOAT>(n * sizeof(type), &p); if (m__n_max == 0) return;}
#define COND_MEM_ALLOC(flag,type,p) {if(m__flags & flag) MEM_ALLOC(type,p)}

namespace Pricer {
    void pricer_context::alloc_mem(uint64_t n) {
        MEM_ALLOC(FLOAT, m__s)
        MEM_ALLOC(FLOAT, m__sigma)
        MEM_ALLOC(FLOAT, m__t)
        MEM_ALLOC(FLOAT, m__tau)
        MEM_ALLOC(FLOAT, m__r)
        MEM_ALLOC(FLOAT, m__sigmaA)
        MEM_ALLOC(FLOAT, m__sigmaA2T2)
        MEM_ALLOC(FLOAT, m__sigmaAsqrtT)
        MEM_ALLOC(FLOAT, m__emrt)
        MEM_ALLOC(FLOAT, m__long_short)
        MEM_ALLOC(FLOAT, m__put_call)
        MEM_ALLOC(FLOAT, m__d2dx2_prep)

        COND_MEM_ALLOC(PRICER_FLAG_TW_PRICER, FLOAT, m__d1)
        COND_MEM_ALLOC(PRICER_FLAG_TW_PRICER, FLOAT, m__d2)
        COND_MEM_ALLOC(PRICER_FLAG_TW_PRICER, FLOAT, m__price)

    }
}