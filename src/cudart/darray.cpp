#include "tool/darray.h"
#include "glob.accasync.h"
#include "tool/cudalib.h"
#include "tool/error.h"
#include <cstring>
#include <cuda_runtime.h>


namespace tinker {
void wait_for(int queue)
{
   cudaStream_t s = queue == sync_queue ? nullptr : nonblk;
   check_rt(cudaStreamSynchronize(s));
}

void device_memory_copyin_bytes(void* dst, const void* src, size_t nbytes,
                                LPFlag flag)
{
   cudaStream_t s;
   if (flag & LPFlag::DEFAULT_Q)
      s = nullptr;
   else
      s = nonblk;
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyHostToDevice, s));
   if (flag & LPFlag::WAIT)
      check_rt(cudaStreamSynchronize(s));
}

void device_memory_copyout_bytes_async(void* dst, const void* src,
                                       size_t nbytes, int queue)
{
   cudaStream_t s = queue == sync_queue ? nullptr : nonblk;
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToHost, s));
}

void device_memory_copy_bytes(void* dst, const void* src, size_t nbytes,
                              LPFlag flag)
{
   cudaStream_t s;
   if (flag & LPFlag::DEFAULT_Q)
      s = nullptr;
   else
      s = nonblk;
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToDevice, s));
   if (flag & LPFlag::WAIT)
      check_rt(cudaStreamSynchronize(s));
}


void device_memory_zero_bytes(void* dst, size_t nbytes, LPFlag flag)
{
   if (dst == nullptr)
      return;


   cudaStream_t s;
   if (flag & LPFlag::DEFAULT_Q)
      s = nullptr;
   else
      s = nonblk;
   check_rt(cudaMemsetAsync(dst, 0, nbytes, s));
   if (flag & LPFlag::WAIT)
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
}
