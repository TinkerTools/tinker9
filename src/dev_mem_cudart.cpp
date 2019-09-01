#include "dev_mem.h"

#ifdef TINKER_CUDART
#  include "error.h"
#  include <cuda_runtime.h>

TINKER_NAMESPACE_BEGIN
void DeviceMemory::copyin_bytes(void* dst, const void* src, size_t nbytes) {
  check_rt(cudaMemcpy(dst, src, nbytes, cudaMemcpyHostToDevice));
}

void DeviceMemory::copyout_bytes(void* dst, const void* src, size_t nbytes) {
  check_rt(cudaMemcpy(dst, src, nbytes, cudaMemcpyDeviceToHost));
}

void DeviceMemory::copy_bytes(void* dst, const void* src, size_t nbytes) {
  check_rt(cudaMemcpy(dst, src, nbytes, cudaMemcpyDeviceToDevice));
}

void DeviceMemory::zero_bytes(void* dst, size_t nbytes) {
  check_rt(cudaMemset(dst, 0, nbytes));
}

void DeviceMemory::deallocate_bytes(void* ptr) { check_rt(cudaFree(ptr)); }

void DeviceMemory::allocate_bytes(void** pptr, size_t nbytes) {
  *pptr = nullptr;
  check_rt(cudaMalloc(pptr, nbytes));
}
TINKER_NAMESPACE_END
#endif
