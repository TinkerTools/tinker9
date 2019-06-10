#ifndef TINKER_GPU_ACC_MATHFUNC_H_
#define TINKER_GPU_ACC_MATHFUNC_H_

#include "util/real_mathfunc.h"

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

#pragma acc routine(abs) seq

#define TINKER_ACC_PARALLEL_PRAGMA_STRING_(...)                                \
  TINKER_STR_VA_ARGS_IMPL(acc parallel loop independent deviceptr(__VA_ARGS__))
#define TINKER_ACC_PARALLEL(ITER, TOTAL, OPER, ...)                            \
  {                                                                            \
    _Pragma(TINKER_ACC_PARALLEL_PRAGMA_STRING_(                                \
        __VA_ARGS__)) for (int ITER = 0; ITER < TOTAL; ++ITER) {               \
      OPER;                                                                    \
    }                                                                          \
  }

TINKER_NAMESPACE_BEGIN
// ans = a dot b
void dotprod(float* cpu_ans, const float* gpu_a, const float* gpu_b, int cpu_n);
void dotprod(double* cpu_ans, const double* gpu_a, const double* gpu_b,
             int cpu_n);
TINKER_NAMESPACE_END

#endif
