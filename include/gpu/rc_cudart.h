#ifndef TINKER_SRC_GPU_RC_CUDART_H_
#define TINKER_SRC_GPU_RC_CUDART_H_

#ifndef TINKER_HOSTONLY
#  include <cuda_runtime.h>
#  include <cufft.h>

#  include "util_macro.h"
TINKER_NAMESPACE_BEGIN
namespace gpu {
typedef cufftHandle fft_plan_t;
}
TINKER_NAMESPACE_END

#endif

#endif
