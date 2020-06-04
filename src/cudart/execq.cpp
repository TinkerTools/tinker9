#include "execq.h"
#include "tool/error.h"
#include <cuda_runtime.h>


namespace tinker {
class ExecQ::Impl
{
public:
   cudaStream_t ss;
};


void ExecQ::deallocate()
{
   check_rt(cudaStreamDestroy(ptr->ss));
   delete ptr;
}


void ExecQ::allocate()
{
   ptr = new ExecQ::Impl;
   check_rt(cudaStreamCreate(&ptr->ss));
}


void ExecQ::synchronize()
{
   check_rt(cudaStreamSynchronize(ptr->ss));
}


void ExecQ::copy_bytes(void* dst, const void* src, size_t nbytes)
{
   check_rt(
      cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToDevice, ptr->ss));
}
}
