//
// Created by peter on 12/19/18.
//

#include <src/math/include/erfc.h>


#ifdef __GNUC__
#define BUILTIN_ASSUME_ALIGNED(x) x=__builtin_assume_aligned(x,64);
#define ASSUME(cond) if(!(cond)) __builtin_unreachable();
#define LIKELY(x)   __builtin_expect((x), 1)
#define UNLIKELY(x) __builtin_expect((x), 0)
#endif

#ifdef __clang__
#define BUILTIN_ASSUME_ALIGNED(x) x=__builtin_assume_aligned(x,64);
#define ASSUME(cond) __builtin_assume(cond);
#define LIKELY(x)   __builtin_expect((x), 1)
#define UNLIKELY(x) __builtin_expect((x), 0)
#endif

#ifdef __INTEL_COMPILER
#define BUILTIN_ASSUME_ALIGNED(x) __assume_aligned(x,64);
#define ASSUME(cond) __assume(cond);
#define LIKELY(x)   __builtin_expect((x), 1)
#define UNLIKELY(x) __builtin_expect((x), 0)
#endif

#define DECLARE_AND_DEFINE(type, x, y) \
   type x[64] __attribute__((aligned(64))); \
   for(UINT64 i=0;i<64;i++) x[i]=y;

typedef double FLOAT;
typedef unsigned long long UINT64;

typedef FLOAT *__restrict__ Real_Ptr;
typedef UINT64 *__restrict__ Uint64_Ptr;


void ERfc_Test(UINT64 n, Real_Ptr a, Real_Ptr b) {


    BUILTIN_ASSUME_ALIGNED(b)
    BUILTIN_ASSUME_ALIGNED(a)

    ASSUME(n % 64 == 0)

    for (UINT64 i = 0; i < 64; ++i) {
        b[i] = xerfc_u15(a[i]);
    }

}