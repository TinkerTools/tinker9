#include "dev_array.h"
#include "cudalib.h"
#include "error.h"
#include <cstring>
#include <cuda_runtime.h>


TINKER_NAMESPACE_BEGIN
void device_memory_copyin_bytes(void* dst, const void* src, size_t nbytes,
                                int sync)
{
   cudaStream_t s;
   if (sync)
      s = nullptr;
   else
      s = nonblk;
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyHostToDevice, s));
   if (sync)
      check_rt(cudaStreamSynchronize(nullptr));
}


void device_memory_copyout_bytes(void* dst, const void* src, size_t nbytes,
                                 DMFlag flag)
{
   cudaStream_t s;
   if (flag & DMFlag::DEFAULT_Q)
      s = nullptr;
   else
      s = nonblk;
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToHost, s));
   if (flag & DMFlag::WAIT)
      check_rt(cudaStreamSynchronize(s));
}


void device_memory_copy_bytes(void* dst, const void* src, size_t nbytes,
                              DMFlag flag)
{
   cudaStream_t s;
   if (flag & DMFlag::DEFAULT_Q)
      s = nullptr;
   else
      s = nonblk;
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToDevice, s));
   if (flag & DMFlag::WAIT)
      check_rt(cudaStreamSynchronize(s));
}


void device_memory_zero_bytes(void* dst, size_t nbytes, DMFlag flag)
{
   cudaStream_t s;
   if (flag & DMFlag::DEFAULT_Q)
      s = nullptr;
   else
      s = nonblk;
   check_rt(cudaMemsetAsync(dst, 0, nbytes, s));
   if (flag & DMFlag::WAIT)
      check_rt(cudaStreamSynchronize(s));
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
