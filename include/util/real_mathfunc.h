#ifndef TINKER_UTIL_REAL_MATHFUNC_H_
#define TINKER_UTIL_REAL_MATHFUNC_H_

#include "macro.h"
#include <cmath>

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
