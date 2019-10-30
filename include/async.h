#pragma once
#include "macro.h"
#include <cstring>


#if TINKER_CUDART
#   include <cuda_runtime.h>
TINKER_NAMESPACE_BEGIN
using Stream = cudaStream_t;
TINKER_NAMESPACE_END
#endif


#if TINKER_HOST
TINKER_NAMESPACE_BEGIN
class StreamSt;
typedef StreamSt* Stream;
TINKER_NAMESPACE_END
#endif


TINKER_NAMESPACE_BEGIN
/// \brief Deallocate the asynchronous stream.
void deallocate_stream(Stream);
/// \brief Allocate the asynchronous stream.
void allocate_stream(Stream*);
/// \brief Synchronize the asynchronous stream.
void synchronize_stream(Stream);


/// \brief
/// Copy between two device addresses without blocking the calling thread.
void copy_bytes_async(void* dst, const void* src, size_t nbytes, Stream s);
TINKER_NAMESPACE_END
