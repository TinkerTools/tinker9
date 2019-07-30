#ifdef TINKER_HOSTONLY

#  include "util_hostonly.h"
#  include <cstring>

TINKER_NAMESPACE_BEGIN
const char* cudaGetErrorName(cudaError_t error) {
  const char* s1 = "HostOnlySuccess";
  const char* s2 = "HostOnlyGenericError";
  if (error == cudaSuccess)
    return s1;
  else
    return s2;
}

const char* cudaGetErrorString(cudaError_t error) {
  const char* s1 = "[HostOnly] Success";
  const char* s2 = "[HostOnly] Opps, something was wrong";
  if (error == cudaSuccess)
    return s1;
  else
    return s2;
}

cudaError_t cudaStreamCreate(cudaStream_t* s) {
  *s = nullptr;
  return cudaSuccess;
}

cudaError_t cudaStreamDestroy(cudaStream_t) { return cudaSuccess; }

cudaError_t cudaStreamSynchronize(cudaStream_t) { return cudaSuccess; }

cudaError_t cudaFree(void* ptr) {
  std::free(ptr);
  return cudaSuccess;
}

cudaError_t cudaMalloc(void** devPtr, size_t size) {
  *devPtr = std::malloc(size);
  return cudaSuccess;
}

cudaError_t cudaMemcpy(void* dst, const void* src, size_t count,
                       hostonly_cuda_enum_) {
  std::memcpy(dst, src, count);
  return cudaSuccess;
}

cudaError_t cudaMemcpyAsync(void* dst, const void* src, size_t count,
                            hostonly_cuda_enum_, cudaStream_t) {
  std::memcpy(dst, src, count);
  return cudaSuccess;
}

cudaError_t cudaMemset(void* devPtr, int value, size_t count) {
  std::memset(devPtr, value, count);
  return cudaSuccess;
}
TINKER_NAMESPACE_END

#endif
