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

#ifndef PRICER_MEMORY_H
#define PRICER_MEMORY_H

#include <cstdint>
#include <cstdlib>
#include <math/pricers/pricer-base.h>

namespace Pricer {

    namespace Private {

        template <typename T> inline
        uint64_t allocate_memory( int64_t length, T * __restrict__* ptr ) {

            size_t size = (length / ALIGN_TO + 1) * ALIGN_TO;
            *ptr = static_cast<T *>(aligned_alloc(ALIGN_TO, size));

            if (*ptr)
                return size;
            else
                return 0;
        }

        template <typename T> inline
        void deallocate_memory(T *ptr) {
            if(ptr) {
                free(ptr);
            }
            ptr = nullptr;
        }
    }

#define DEFINE_VARIABLE(type,x) private:type m__ ## x; public:inline type& INLINE get_ ## x() {return m__ ## x;}

#define PRICER_FLAG_TW_PRICER                           1
#define PRICER_FLAG_TW_COMPUTE_STRIKES_OF_MICROHEDGES   2
    class pricer_context {

    public:
        pricer_context(uint64_t flags, uint64_t n_max, uint64_t m_max = 0) :
            m__flags(flags),
            m__n(0), m__m(0) {
                alloc_mem(n_max, m_max);
        }

        ~pricer_context() {
            dealloc_mem();
        };

        void realloc_mem(uint64_t n);

    private:

        void alloc_mem(uint64_t n, uint64_t m);
        void dealloc_mem();
        void init_memory(uint64_t n1, uint64_t n2, uint64_t m1, uint64_t m2);

    public:
        uint64_t get_n_ALIGNED() {
            if((get_n()%ALIGN_TO) == 0) {
                return get_n();
            }
            return get_n() - (get_n() % ALIGN_TO) + ALIGN_TO;
        }

        uint64_t get_m_ALIGNED() {
            if((get_m()%ALIGN_TO) == 0) {
                return get_m();
            }
            return get_m() - (get_m() % ALIGN_TO) + ALIGN_TO;
        }


        DEFINE_VARIABLE(uint64_t,flags)
        DEFINE_VARIABLE(uint64_t, n_max)
        DEFINE_VARIABLE(uint64_t, n)

        DEFINE_VARIABLE(Real_Ptr, prices)
        DEFINE_VARIABLE(Real_Ptr, d1)
        DEFINE_VARIABLE(Real_Ptr, d2)
        DEFINE_VARIABLE(Real_Ptr, long_short)
        DEFINE_VARIABLE(Real_Ptr, put_call)
        DEFINE_VARIABLE(Real_Ptr, x)
        DEFINE_VARIABLE(Real_Ptr, s)
        DEFINE_VARIABLE(Real_Ptr, sigma)
        DEFINE_VARIABLE(Real_Ptr, t)
        DEFINE_VARIABLE(Real_Ptr, tau)
        DEFINE_VARIABLE(Real_Ptr, r)
        DEFINE_VARIABLE(Real_Ptr, sigmaA)
        DEFINE_VARIABLE(Real_Ptr, sigmaA2T2)
        DEFINE_VARIABLE(Real_Ptr, sigmaAsqrtT)
        DEFINE_VARIABLE(Real_Ptr, emrt)
        DEFINE_VARIABLE(Real_Ptr, d2dx2_prep)
        DEFINE_VARIABLE(Real_Ptr, ddx_price)
        DEFINE_VARIABLE(Real_Ptr, d2dx2)

        DEFINE_VARIABLE(Int32_Ptr, to_structure)
        DEFINE_VARIABLE(Real_Ptr, offsets)

        DEFINE_VARIABLE(uint64_t, m)
        DEFINE_VARIABLE(uint64_t, m_max)
        DEFINE_VARIABLE(Real_Ptr, premiums)
        DEFINE_VARIABLE(Real_Ptr, instrument_prices)
        DEFINE_VARIABLE(Real_Ptr, instrument_pricesl)
        DEFINE_VARIABLE(Real_Ptr, instrument_pricesh)
        DEFINE_VARIABLE(Real_Ptr, x_)
        DEFINE_VARIABLE(Real_Ptr, xl_)
        DEFINE_VARIABLE(Real_Ptr, xh_)


    };

}


#endif //PRICER_MEMORY_H
