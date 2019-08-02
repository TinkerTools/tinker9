#ifdef TINKER_HOSTONLY

#  include "util_rt_hostonly.h"
#  include <cstring>

TINKER_NAMESPACE_BEGIN
void dealloc_stream(Stream) {}

void alloc_stream(Stream* s) { *s = nullptr; }

void sync_stream(Stream) {}

void copy_memory(void* dst, const void* src, size_t count, CopyDirection) {
  std::memcpy(dst, src, count);
}

void copy_memory_async(void* dst, const void* src, size_t count, CopyDirection,
                       Stream) {
  std::memcpy(dst, src, count);
}
TINKER_NAMESPACE_END

#endif
