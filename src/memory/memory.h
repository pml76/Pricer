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

        inline template <typename T>
        uint64_t allocate_memory( int64_t length, T **ptr ) {

            size_t size = (length / ALIGN_TO + 1) * ALIGN_TO;
            *ptr = static_cast<T *>(aligned_alloc(ALIGN_TO, size));

            if (*ptr)
                return size;
            else
                return 0;
        }

        inline template <typename T>
        void deallocate_memory(T *ptr) {
            free(ptr);
        }
    }

    class pricer_context {
    public:
        pricer_context(uint64_t flags, uint64_t n_max, uint64_t inc = 1024 * 1024 * ALIGN_TO);

        void add_data(uint64_t n, Real_Ptr s, Real_Ptr sigma, Real_Ptr t, Real_Ptr tau, Real_Ptr r);

    private:

#define PRICER_FLAG_TW_PRICER 1

        void alloc_mem(uint64_t n);


        uint64_t m__flags;
        uint64_t m__n_max;
        uint64_t m__n;

        Real_Ptr  m__price;
        Real_Ptr m__d1;
        Real_Ptr m__d2;
        Real_Ptr m__long_short;
        Real_Ptr m__put_call;
        Real_Ptr m__s;
        Real_Ptr m__sigma;
        Real_Ptr m__t;
        Real_Ptr m__tau;
        Real_Ptr m__r;
        Real_Ptr m__sigmaA;
        Real_Ptr m__sigmaA2T2;
        Real_Ptr m__sigmaAsqrtT;
        Real_Ptr m__emrt;
        Real_Ptr m__d2dx2_prep;
    };

}


#endif //PRICER_MEMORY_H
