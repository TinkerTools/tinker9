#ifndef TINKER_UTIL_CUDART_H_
#define TINKER_UTIL_CUDART_H_

#ifndef TINKER_HOSTONLY
#  include <cuda_runtime.h>
#  include <cufft.h>

#  include "util_macro.h"
TINKER_NAMESPACE_BEGIN
typedef cufftHandle fft_plan_t;
TINKER_NAMESPACE_END

#endif

#endif
