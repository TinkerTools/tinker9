#ifndef TINKER_UTIL_RT_HOSTONLY_H_
#define TINKER_UTIL_RT_HOSTONLY_H_

#include "util_macro.h"
#include <cstdlib>
#include <fftw3.h>

TINKER_NAMESPACE_BEGIN
typedef int* Stream;
void dealloc_stream(Stream);
void alloc_stream(Stream*);
void sync_stream(Stream);

enum class CopyDirection { DeviceToDevice, DeviceToHost, HostToDevice };
void copy_memory(void* dst, const void* src, size_t count, CopyDirection);
void copy_memory_async(void* dst, const void* src, size_t count, CopyDirection,
                       Stream);

struct FFTPlan {
#if defined(TINKER_SINGLE_PRECISION)
  fftwf_plan
#elif defined(TINKER_DOUBLE_PRECISION)
  fftw_plan
#else
  static_assert(false, "");
#endif
      planf, ///< fft front plan
      planb; ///< fft back plan
};
TINKER_NAMESPACE_END

#endif
