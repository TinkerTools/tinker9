#ifndef TINKER_GPU_ACC_MATHFUNC_H_
#define TINKER_GPU_ACC_MATHFUNC_H_

#include "decl_mathfunc.h"

#pragma acc routine(abs) seq

#pragma acc routine(sqrt) seq
#pragma acc routine(exp) seq
#pragma acc routine(floor) seq
#pragma acc routine(fabs) seq
#pragma acc routine(pow) seq
#pragma acc routine(cos) seq
#pragma acc routine(sin) seq
#pragma acc routine(erf) seq
#pragma acc routine(erfc) seq
#pragma acc routine(fmin) seq
#pragma acc routine(fmax) seq

#pragma acc routine(sqrtf) seq
#pragma acc routine(expf) seq
#pragma acc routine(floorf) seq
#pragma acc routine(fabsf) seq
#pragma acc routine(powf) seq
#pragma acc routine(cosf) seq
#pragma acc routine(sinf) seq
#pragma acc routine(erff) seq
#pragma acc routine(erfcf) seq
#pragma acc routine(fminf) seq
#pragma acc routine(fmaxf) seq

#endif
