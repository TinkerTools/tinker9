#include "tool/darray.h"
#include "glob.accasync.h"
#include "tool/cudalib.h"
#include "tool/error.h"
#include <cstring>
#include <cuda_runtime.h>


namespace tinker {
void wait_for(int queue)
{
   cudaStream_t st = queue == g::q1 ? g::s1 : g::s0;
   check_rt(cudaStreamSynchronize(st));
}

void device_memory_copyin_bytes_async(void* dst, const void* src, size_t nbytes,
                                      int queue)
{
   cudaStream_t st = queue == g::q1 ? g::s1 : g::s0;
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyHostToDevice, st));
}

void device_memory_copyout_bytes_async(void* dst, const void* src,
                                       size_t nbytes, int queue)
{
   cudaStream_t st = queue == g::q1 ? g::s1 : g::s0;
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToHost, st));
}

void device_memory_copy_bytes_async(void* dst, const void* src, size_t nbytes,
                                    int queue)
{
   cudaStream_t st = queue == g::q1 ? g::s1 : g::s0;
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToDevice, st));
}


void device_memory_zero_bytes_async(void* dst, size_t nbytes, int queue)
{
   if (dst == nullptr)
      return;


   cudaStream_t st = queue == g::q1 ? g::s1 : g::s0;
   check_rt(cudaMemsetAsync(dst, 0, nbytes, st));
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
