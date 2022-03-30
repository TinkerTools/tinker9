#include "tool/darray.h"
#include <cstdlib>
#include <cstring>

namespace tinker {
void waitFor(int) {}

void deviceMemoryCopyinBytesAsync(void* dst, const void* src, size_t nbytes, int)
{
   std::memcpy(dst, src, nbytes);
}

void deviceMemoryCopyoutBytesAsync(void* dst, const void* src, size_t nbytes, int)
{
   std::memcpy(dst, src, nbytes);
}

void deviceMemoryCopyBytesAsync(void* dst, const void* src, size_t nbytes, int)
{
   std::memcpy(dst, src, nbytes);
}

void deviceMemoryZeroBytesAsync(void* dst, size_t nbytes, int)
{
   if (dst == nullptr)
      return;

   std::memset(dst, 0, nbytes);
}

void deviceMemoryDeallocate(void* ptr)
{
   std::free(ptr);
}

void deviceMemoryAllocateBytes(void** pptr, size_t nbytes)
{
   *pptr = std::malloc(nbytes);
}
}
