//
// Created by peter on 12/22/18.
//

#ifndef PRICER_PRICER_BASE_H
#define PRICER_PRICER_BASE_H


#define ALIGN_TO 64


#ifdef __cplusplus
#define restrict
extern "C" {
#endif

#ifdef __GNUC__
#ifdef ASSUME_ALIGNED
#undef ASSUME_ALIGNED
#endif
#define ASSUME_ALIGNED(x) x=__builtin_assume_aligned(x,ALIGN_TO);

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

#ifdef ASSUME_ALIGNED
#undef ASSUME_ALIGNED
#endif
#define ASSUME_ALIGNED(x) x=__builtin_assume_aligned(x,64);

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

#ifdef ASSUME_ALIGNED
#undef ASSUME_ALIGNED
#endif
#define ASSUME_ALIGNED(x) __assume_aligned(x,64);


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
typedef UINT64 *__restrict__ Uint64_Ptr;
typedef int32_t *__restrict Int32_Ptr;

#ifdef __cplusplus
};
#endif

#endif //PRICER_PRICER_BASE_H
