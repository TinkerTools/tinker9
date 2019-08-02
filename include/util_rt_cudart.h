#ifndef TINKER_UTIL_RT_CUDART_H_
#define TINKER_UTIL_RT_CUDART_H_

#include "util_macro.h"
#include <cuda_runtime.h>
#include <cufft.h>

TINKER_NAMESPACE_BEGIN
typedef cudaStream_t Stream;
void dealloc_stream(Stream);
void alloc_stream(Stream*);
void sync_stream(Stream);

enum class CopyDirection { DeviceToDevice, DeviceToHost, HostToDevice };
void copy_memory(void* dst, const void* src, size_t count, CopyDirection);
void copy_memory_async(void* dst, const void* src, size_t count, CopyDirection,
                       Stream);

typedef cufftHandle FFTPlan;
TINKER_NAMESPACE_END

#endif
