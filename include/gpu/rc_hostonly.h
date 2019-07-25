#ifndef TINKER_SRC_GPU_RC_HOSTONLY_H_
#define TINKER_SRC_GPU_RC_HOSTONLY_H_

#ifdef TINKER_HOSTONLY

#  include "util_macro.h"
#  include <cstdlib>
#  include <fftw3.h>

TINKER_NAMESPACE_BEGIN
namespace gpu {
enum hostonly_cuda_enum_ {
  cudaSuccess = 0,
  cudaErrorInvalidHostPointer,

  cudaMemcpyDeviceToDevice,
  cudaMemcpyDeviceToHost,
  cudaMemcpyHostToDevice
};

typedef hostonly_cuda_enum_ cudaError_t;
const char* cudaGetErrorName(cudaError_t error);
const char* cudaGetErrorString(cudaError_t error);

typedef hostonly_cuda_enum_* cudaStream_t;
cudaError_t cudaStreamCreate(cudaStream_t*);
cudaError_t cudaStreamDestroy(cudaStream_t);
cudaError_t cudaStreamSynchronize(cudaStream_t);

cudaError_t cudaFree(void* ptr);
cudaError_t cudaMalloc(void** devPtr, size_t size);
template <class T>
cudaError_t cudaMalloc(T** devPtr, size_t size) {
  return cudaMalloc(reinterpret_cast<void**>(devPtr), size);
}
cudaError_t cudaMemcpy(void* dst, const void* src, size_t count,
                       hostonly_cuda_enum_ kind);
cudaError_t cudaMemcpyAsync(void* dst, const void* src, size_t count,
                            hostonly_cuda_enum_ kind, cudaStream_t stream);
cudaError_t cudaMemset(void* devPtr, int value, size_t count);

typedef struct hostonly_fftw_plans_st_ {
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
