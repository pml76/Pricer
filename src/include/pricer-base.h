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


#ifndef PRICER_PRICER_BASE_H
#define PRICER_PRICER_BASE_H

#include <stdint.h>


#define ALIGN_TO 64


#ifdef __cplusplus
#define restrict
extern "C" {
#endif

#ifdef __GNUC__

#ifdef INLINE
#undef INLINE
#endif
#define INLINE __attribute__((always_inline))

#ifdef ASSUME_ALIGNED
#undef ASSUME_ALIGNED
#endif
#define ASSUME_ALIGNED(type,x) x=(type) __builtin_assume_aligned((void *)x, ALIGN_TO);

#ifdef ASSUME
#undef ASSUME
#endif
#define ASSUME(cond) if(!(cond)) __builtin_unreachable();

#ifdef LIKELY
#undef LIKELY
#endif
#define LIKELY(x)   __builtin_expect((x), 1)

#ifdef UNLIKELY
#undef UNLIKELY
#endif
#define UNLIKELY(x) __builtin_expect((x), 0)

#ifdef PREFETCH
#undef PREFETCH
#endif
#define PREFETCH(a, b, c) __builtin_prefetch(a,b,c)
#endif

#ifdef __clang__

#ifdef INLINE
#undef INLINE
#endif
#define INLINE __attribute__((always_inline))


#ifdef ASSUME_ALIGNED
#undef ASSUME_ALIGNED
#endif
#define ASSUME_ALIGNED(type,x) x=(type) __builtin_assume_aligned(x,ALIGN_TO);

#ifdef ASSUME
#undef ASSUME
#endif
#define ASSUME(cond) __builtin_assume(cond);

#ifdef LIKELY
#undef LIKELY
#endif
#define LIKELY(x)   __builtin_expect((x), 1)

#ifdef UNLIKELY
#undef UNLIKELY
#endif
#define UNLIKELY(x) __builtin_expect((x), 0)
#endif

#ifdef __INTEL_COMPILER

#ifdef INLINE
#undef INLINE
#endif
#define INLINE __attribute__((always_inline))


#ifdef ASSUME_ALIGNED
#undef ASSUME_ALIGNED
#endif
#define ASSUME_ALIGNED(type,x) __assume_aligned(x,64);


#ifdef ASSUME
#undef ASSUME
#endif
#define ASSUME(cond) __assume(cond);

#ifdef LIKELY
#undef LIKELY
#endif
#define LIKELY(x)   __builtin_expect((x), 1)


#ifdef UNLIKELY
#undef UNLIKELY
#endif
#define UNLIKELY(x) __builtin_expect((x), 0)

#endif

#define DECLARE_AND_DEFINE_ARRAY(type, size, x, y) \
   type x[64] __attribute__((aligned(ALIGN_TO))); \
   for(UINT64 i=0;i<ALIGN_TO;i++) x[i]=y;


typedef double FLOAT;
typedef unsigned long long UINT64;

typedef FLOAT *__restrict__ Real_Ptr;
typedef uint64_t *__restrict__ Uint64_Ptr;
typedef int32_t *__restrict__ Int32_Ptr;

#ifdef __cplusplus
};
#endif

#endif //PRICER_PRICER_BASE_H
