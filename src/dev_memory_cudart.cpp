#include "macro.h"

#if TINKER_CUDART
#   include "error.h"
#   include <cstring>
#   include <cuda_runtime.h>

TINKER_NAMESPACE_BEGIN
void device_memory_copyin_bytes(void* dst, const void* src, size_t nbytes)
{
   check_rt(cudaMemcpy(dst, src, nbytes, cudaMemcpyHostToDevice));
}

void device_memory_copyout_bytes(void* dst, const void* src, size_t nbytes)
{
   check_rt(cudaMemcpy(dst, src, nbytes, cudaMemcpyDeviceToHost));
}

void device_memory_copy_bytes(void* dst, const void* src, size_t nbytes)
{
   check_rt(cudaMemcpy(dst, src, nbytes, cudaMemcpyDeviceToDevice));
}

void device_memory_zero_bytes(void* dst, size_t nbytes)
{
   check_rt(cudaMemset(dst, 0, nbytes));
}

void device_memory_deallocate_bytes(void* ptr)
{
   check_rt(cudaFree(ptr));
}

void device_memory_allocate_bytes(void** pptr, size_t nbytes)
{
   *pptr = nullptr;
   check_rt(cudaMalloc(pptr, nbytes));
}
TINKER_NAMESPACE_END
#endif
