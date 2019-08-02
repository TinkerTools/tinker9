#ifndef TINKER_UTIL_RT_CUDART_H_
#define TINKER_UTIL_RT_CUDART_H_

#include "util_macro.h"
#include <cuda_runtime.h>
#include <cufft.h>

TINKER_NAMESPACE_BEGIN
typedef cufftHandle fft_plan_t;
TINKER_NAMESPACE_END

#endif
