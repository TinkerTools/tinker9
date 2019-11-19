#include "async.h"
#include "error.h"
#include <cuda_runtime.h>


TINKER_NAMESPACE_BEGIN
class StreamSt
{};


void* async_acc;


void deallocate_stream(Stream s)
{
   check_rt(cudaStreamDestroy(reinterpret_cast<cudaStream_t>(s)));
}


void allocate_stream(Stream* s)
{
   check_rt(cudaStreamCreate(reinterpret_cast<cudaStream_t*>(s)));
}


void synchronize_stream(Stream s)
{
   check_rt(cudaStreamSynchronize(reinterpret_cast<cudaStream_t>(s)));
}


void copy_bytes_async(void* dst, const void* src, size_t nbytes, Stream s)
{
   check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToDevice,
                            reinterpret_cast<cudaStream_t>(s)));
}
TINKER_NAMESPACE_END
