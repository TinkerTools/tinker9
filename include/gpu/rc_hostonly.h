#ifndef TINKER_SRC_GPU_RC_HOSTONLY_H_
#define TINKER_SRC_GPU_RC_HOSTONLY_H_

#ifdef TINKER_HOSTONLY

#  include "util/macro.h"
#  include <cstdlib>
#  include <fftw3.h>

TINKER_NAMESPACE_BEGIN
namespace gpu {
enum hostonly_cuda_enum__ {
  cudaSuccess = 0,
  cudaErrorInvalidHostPointer,

  cudaMemcpyDeviceToDevice,
  cudaMemcpyDeviceToHost,
  cudaMemcpyHostToDevice
};

typedef hostonly_cuda_enum__ cudaError_t;
const char* cudaGetErrorName(cudaError_t error);
const char* cudaGetErrorString(cudaError_t error);

cudaError_t cudaFree(void* ptr);
cudaError_t cudaMalloc(void** devPtr, size_t size);
template <class T>
cudaError_t cudaMalloc(T** devPtr, size_t size) {
  return cudaMalloc(reinterpret_cast<void**>(devPtr), size);
}
cudaError_t cudaMemcpy(void* dst, const void* src, size_t count,
                       hostonly_cuda_enum__ kind);
cudaError_t cudaMemset(void* devPtr, int value, size_t count);

typedef struct hostonly_fftw_plans_st__ {
#  if defined(TINKER_GPU_SINGLE)
  fftwf_plan
#  elif defined(TINKER_GPU_DOUBLE)
  fftw_plan
#  else
  static_assert(false, "");
#  endif
      planf, /// fft front plan
      planb; /// fft back plan
} fft_plan_t;
}
TINKER_NAMESPACE_END

#endif

#endif
