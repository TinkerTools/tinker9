#ifndef TINKER_ASYNC_H_
#define TINKER_ASYNC_H_

#include "macro.h"
#include <cstring>

#if TINKER_CUDART
#   include <cuda_runtime.h>
TINKER_NAMESPACE_BEGIN
typedef cudaStream_t Stream;
TINKER_NAMESPACE_END
#endif

#if TINKER_HOST
TINKER_NAMESPACE_BEGIN
class StreamSt;
typedef StreamSt* Stream;
TINKER_NAMESPACE_END
#endif

TINKER_NAMESPACE_BEGIN
/// @brief
/// deallocate, allocate, and synchronize the asynchronous stream
/// @{
void deallocate_stream (Stream);
void allocate_stream (Stream*);
void synchronize_stream (Stream);
/// @}

/// @brief
/// copy between two device addresses without blocking the calling thread
void copy_bytes_async (void* dst, const void* src, size_t nbytes, Stream s);
TINKER_NAMESPACE_END

#endif
