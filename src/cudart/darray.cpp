#include "tool/darray.h"
#include "accasync.h"
#include "tool/cudalib.h"
#include "tool/error.h"
#include <cstring>
#include <cuda_runtime.h>

namespace tinker {
void waitFor(int queue)
{
   cudaStream_t st = queue == g::q1 ? g::s1 : g::s0;
   check_rt(cudaStreamSynchronize(st));
}

void deviceMemoryCopyinBytesAsync(void* dst, const void* src, size_t nbytes, int queue)
{
   cudaStream_t st = queue == g::q1 ? g::s1 : g::s0;
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyHostToDevice, st));
}

void deviceMemoryCopyoutBytesAsync(void* dst, const void* src, size_t nbytes, int queue)
{
   cudaStream_t st = queue == g::q1 ? g::s1 : g::s0;
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToHost, st));
}

void deviceMemoryCopyBytesAsync(void* dst, const void* src, size_t nbytes, int queue)
{
   cudaStream_t st = queue == g::q1 ? g::s1 : g::s0;
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToDevice, st));
}

void deviceMemoryZeroBytesAsync(void* dst, size_t nbytes, int queue)
{
   if (dst == nullptr)
      return;

   cudaStream_t st = queue == g::q1 ? g::s1 : g::s0;
   check_rt(cudaMemsetAsync(dst, 0, nbytes, st));
}

void deviceMemoryDeallocate(void* ptr)
{
   check_rt(cudaFree(ptr));
}

void deviceMemoryAllocateBytes(void** pptr, size_t nbytes)
{
   *pptr = nullptr;
   check_rt(cudaMalloc(pptr, nbytes));
}
}
