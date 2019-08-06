#ifndef TINKER_UTIL_RT_H_
#define TINKER_UTIL_RT_H_

#ifdef TINKER_HOST
#  include "rt_host.h"
#else
#  include "rt_cudart.h"
#endif

#include "gen_unit.h"

TINKER_NAMESPACE_BEGIN
typedef GenericUnit<FFTPlan, GenericUnitVersion::V0> FFTPlanUnit;

/// @brief
/// zero-out, deallocate, or allocate bytes on device
/// @{
void zero_bytes(void* ptr, size_t nbytes);
void dealloc_bytes(void* ptr);
void alloc_bytes(void** ptr, size_t nbytes);
template <class T>
void alloc_bytes(T** ptr, size_t nbytes) {
  return alloc_bytes(reinterpret_cast<void**>(ptr), nbytes);
}
/// @}

/// @brief
/// copy @c nbytes from host to device (copyin), from device to host (copyout),
/// or between two device addresses (copy);
/// will block the calling thread
/// @{
void copyin_bytes(void* dst, const void* src, size_t nbytes);
void copyout_bytes(void* dst, const void* src, size_t nbytes);
void copy_bytes(void* dst, const void* src, size_t nbytes);
/// @}

template <>
struct GenericUnitAlloc<GenericUnitVersion::V1> {
  struct Dealloc {
    void operator()(void* ptr) { dealloc_bytes(ptr); }
  };

  struct Alloc {
    void operator()(void** ptr, size_t nbytes) { alloc_bytes(ptr, nbytes); }
  };

  struct Copyin {
    void operator()(void* dst, const void* src, size_t nbytes) {
      copyin_bytes(dst, src, nbytes);
    }
  };
};

/// @brief
/// deallocate, allocate, and synchronize the asynchronous stream
/// @{
void dealloc_stream(Stream);
void alloc_stream(Stream*);
void sync_stream(Stream);
/// @}
/// @brief
/// copy between two device addresses without blocking the calling thread
void copy_bytes_async(void* dst, const void* src, size_t nbytes, Stream s);
TINKER_NAMESPACE_END

#endif
