#ifndef TINKER_GPU_MATHFUNC_H
#  define TINKER_GPU_MATHFUNC_H_

#  include "util/macro.h"
#  include <cmath>

#  pragma acc routine(sqrt) seq
#  pragma acc routine(sqrtf) seq
#  pragma acc routine(exp) seq
#  pragma acc routine(expf) seq

#  ifdef TINKER_GPU_DOUBLE
#    define REAL_SQRT sqrt
#    define REAL_EXP exp
#  endif

#  ifdef TINKER_GPU_SINGLE
#    define REAL_SQRT sqrtf
#    define REAL_EXP expf
#  endif

#endif
