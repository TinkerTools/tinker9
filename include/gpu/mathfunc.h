#ifndef TINKER_GPU_MATHFUNC_H_
#define TINKER_GPU_MATHFUNC_H_

#include "util/macro.h"
#include <cmath>

#pragma acc routine(sqrt) seq
#pragma acc routine(exp) seq
#pragma acc routine(floor) seq
#pragma acc routine(fabs) seq

#pragma acc routine(sqrtf) seq
#pragma acc routine(expf) seq
#pragma acc routine(floorf) seq
#pragma acc routine(fabsf) seq

#ifdef TINKER_GPU_DOUBLE
#  define REAL_SQRT sqrt
#  define REAL_EXP exp
#  define REAL_FLOOR floor
#  define REAL_ABS fabs
#endif

#ifdef TINKER_GPU_SINGLE
#  define REAL_SQRT sqrtf
#  define REAL_EXP expf
#  define REAL_FLOOR floorf
#  define REAL_ABS fabsf
#endif

#define REAL_SQ(x) ((x) * (x))

#endif
