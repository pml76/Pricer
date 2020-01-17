//
// Created by peter on 1/17/20.
//

#ifndef PRICER_ALLOCATE_H
#define PRICER_ALLOCATE_H

#include <cstdint>
#include <stdlib.h>
#include <cstring>

#include <include/pricer-base.h>

namespace Pricer {

    namespace Private {
        template<typename T>
        inline
        uint64_t allocate_memory(uint64_t length, T *__restrict__ *ptr) {

            size_t size;
            if ((length % ALIGN_TO) == 0) {
                size = length;
            } else {
                size = length - (length % ALIGN_TO) + ALIGN_TO;
            }
#ifdef __WIN64
            *ptr = static_cast<T *>(_aligned_malloc(size * sizeof(T), ALIGN_TO));
#else
            *ptr = static_cast<T *>(aligned_alloc(ALIGN_TO, size * sizeof(T)));
#endif
            if (*ptr)
                return size;
            else
                return 0;
        }

        template<typename T>
        inline
        void deallocate_memory(T *ptr) {
            if (ptr) {
                free(ptr);
            }
            ptr = nullptr;
        }

        template<typename T>
        inline
        uint64_t reallocate_memory(uint64_t new_length, uint64_t old_length, T *__restrict__ *ptr) {
            T *__restrict__ p;
            auto n = allocate_memory<T>(new_length, &p);
            if (n > 0) {
                memcpy(p, *ptr, old_length * sizeof(T));
                deallocate_memory(*ptr);
                *ptr = p;
            }
            return n;

        }


    }

}


#endif //PRICER_ALLOCATE_H
