#include "dev_array.h"
#include "error.h"
#include <cstring>
#include <cuda_runtime.h>


TINKER_NAMESPACE_BEGIN
void device_memory_copyin_bytes(void* dst, const void* src, size_t nbytes,
                                void* stream)
{
   auto s = reinterpret_cast<cudaStream_t>(stream);
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyHostToDevice, s));
   check_rt(cudaStreamSynchronize(s));
}


void device_memory_copyout_bytes(void* dst, const void* src, size_t nbytes,
                                 void* stream)
{
   auto s = reinterpret_cast<cudaStream_t>(stream);
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToHost, s));
   check_rt(cudaStreamSynchronize(s));
}


void device_memory_copy_bytes(void* dst, const void* src, size_t nbytes,
                              void* stream)
{
   auto s = reinterpret_cast<cudaStream_t>(stream);
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToDevice, s));
   check_rt(cudaStreamSynchronize(s));
}


void device_memory_zero_bytes_async(void* dst, size_t nbytes, void* stream)
{
   check_rt(cudaMemsetAsync(dst, 0, nbytes, (cudaStream_t)stream));
}


void device_memory_zero_bytes(void* dst, size_t nbytes)
{
   device_memory_zero_bytes_async(dst, nbytes);
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
