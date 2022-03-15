#include "tool/darray.h"
#include <cstdlib>
#include <cstring>

namespace tinker {
void wait_for(int) {}

void device_memory_copyin_bytes_async(void* dst, const void* src, size_t nbytes, int)
{
   std::memcpy(dst, src, nbytes);
}

void device_memory_copyout_bytes_async(void* dst, const void* src, size_t nbytes, int)
{
   std::memcpy(dst, src, nbytes);
}

void device_memory_copy_bytes_async(void* dst, const void* src, size_t nbytes, int)
{
   std::memcpy(dst, src, nbytes);
}

void device_memory_zero_bytes_async(void* dst, size_t nbytes, int)
{
   if (dst == nullptr)
      return;

   std::memset(dst, 0, nbytes);
}

void device_memory_deallocate_bytes(void* ptr)
{
   std::free(ptr);
}

void device_memory_allocate_bytes(void** pptr, size_t nbytes)
{
   *pptr = std::malloc(nbytes);
}
}
