#include "dev_mem.h"

#ifdef TINKER_HOST
#  include <cstdlib>

TINKER_NAMESPACE_BEGIN
void DeviceMemory::copyin_bytes(void* dst, const void* src, size_t nbytes) {
  std::memcpy(dst, src, nbytes);
}

void DeviceMemory::copyout_bytes(void* dst, const void* src, size_t nbytes) {
  std::memcpy(dst, src, nbytes);
}

void DeviceMemory::copy_bytes(void* dst, const void* src, size_t nbytes) {
  std::memcpy(dst, src, nbytes);
}

void DeviceMemory::zero_bytes(void* dst, size_t nbytes) {
  std::memset(dst, 0, nbytes);
}

void DeviceMemory::deallocate_bytes(void* ptr) { std::free(ptr); }

void DeviceMemory::allocate_bytes(void** pptr, size_t nbytes) {
  *pptr = std::malloc(nbytes);
}
TINKER_NAMESPACE_END
#endif
