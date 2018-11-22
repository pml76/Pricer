//
// Created by peter on 11/20/18.
//

#ifndef PRICER_SQRT_WRAPPER_H
#define PRICER_SQRT_WRAPPER_H


#define SQRT_VERSION 3


#if SQRT_VERSION == 1

#ifdef __cplusplus
extern "C"
#endif
double ieee754_sqrt(double x);

#endif

#if SQRT_VERSION == 2


__inline __attribute__ ((__always_inline__)) double
ieee754_sqrt (double d)
{
  double res;



  asm ("sqrtsd %1, %0" : "=x" (res) : "xm" (d));

  return res;
}

#endif


#if SQRT_VERSION = 3

#define ieee754_sqrt(x) sqrt(x)

#endif

#endif //PRICER_SQRT_WRAPPER_H
