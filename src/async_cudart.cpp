#include "async.h"
#include "error.h"

#if TINKER_CUDART
TINKER_NAMESPACE_BEGIN
void deallocate_stream (Stream s)
{
   check_rt (cudaStreamDestroy (s));
}

void allocate_stream (Stream* s)
{
   check_rt (cudaStreamCreate (s));
}

void synchronize_stream (Stream s)
{
   check_rt (cudaStreamSynchronize (s));
}

void copy_bytes_async (void* dst, const void* src, size_t nbytes, Stream s)
{
   check_rt (cudaMemcpyAsync (dst, src, nbytes, cudaMemcpyDeviceToDevice, s));
}
TINKER_NAMESPACE_END
#endif
