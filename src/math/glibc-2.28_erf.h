//
// Created by peter on 11/20/18.
//

#ifndef PRICER_ERF_WRAPPER_H
#define PRICER_ERF_WRAPPER_H


#ifdef __cplusplus
extern "C"
#endif
double glibc_erf(double x);

#ifdef __cplusplus
extern "C"
#endif
double glibc_erfc(double x);


#endif //PRICER_ERF_WRAPPER_H
