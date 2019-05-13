#ifndef TINKER_GPU_ACC_MATHFUNC_H_
#define TINKER_GPU_ACC_MATHFUNC_H_

#include "util/macro.h"
#include <cmath>

#pragma acc routine(sqrt) seq
#pragma acc routine(exp) seq
#pragma acc routine(floor) seq
#pragma acc routine(fabs) seq
#pragma acc routine(pow) seq

#pragma acc routine(sqrtf) seq
#pragma acc routine(expf) seq
#pragma acc routine(floorf) seq
#pragma acc routine(fabsf) seq
#pragma acc routine(powf) seq

#ifdef TINKER_GPU_DOUBLE
#  define REAL_SQRT sqrt
#  define REAL_EXP exp
#  define REAL_FLOOR floor
#  define REAL_ABS fabs
#  define REAL_POW pow
#  define REAL_RECIP(x) (1.0 / (x))
#endif

#ifdef TINKER_GPU_SINGLE
#  define REAL_SQRT sqrtf
#  define REAL_EXP expf
#  define REAL_FLOOR floorf
#  define REAL_ABS fabsf
#  define REAL_POW powf
#  define REAL_RECIP(x) (1.0f / (x))
#endif

#define REAL_SQ(x) ((x) * (x))
#define REAL_CUBE(x) ((x) * (x) * (x))

#endif
