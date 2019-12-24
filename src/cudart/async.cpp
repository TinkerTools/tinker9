#include "async.h"
#include "error.h"
#include <cuda_runtime.h>


TINKER_NAMESPACE_BEGIN
void deallocate_stream(void* s)
{
   auto ss = reinterpret_cast<cudaStream_t>(s);
   check_rt(cudaStreamDestroy(ss));
}


void allocate_stream(void** s)
{
   auto ss = reinterpret_cast<cudaStream_t*>(s);
   check_rt(cudaStreamCreate(ss));
}


void synchronize_stream(void* s)
{
   auto ss = reinterpret_cast<cudaStream_t>(s);
   check_rt(cudaStreamSynchronize(ss));
}


void copy_bytes_async(void* dst, const void* src, size_t nbytes, void* s)
{
   auto ss = reinterpret_cast<cudaStream_t>(s);
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToDevice, ss));
}
TINKER_NAMESPACE_END
