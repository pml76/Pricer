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
            free(ptr);
        }
    }

#define DEFINE_VARIABLE(type,x) private:type m__ ## x; public:type get_ ## x() const{return m__ ## x;}

#define PRICER_FLAG_TW_PRICER 1
    class pricer_context {

    public:
        pricer_context(uint64_t flags, uint64_t n_max) :
            m__flags(flags),
            m__n(0) {
                alloc_mem(n_max);
        }

        ~pricer_context() {
            dealloc_mem();
        };

        void realloc_mem(uint64_t n);

    private:

        void alloc_mem(uint64_t n);
        void dealloc_mem();



        DEFINE_VARIABLE(uint64_t,flags)
        DEFINE_VARIABLE(uint64_t, n_max)
        DEFINE_VARIABLE(uint64_t, n)

        DEFINE_VARIABLE(Real_Ptr, price)
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
    };

}


#endif //PRICER_MEMORY_H
