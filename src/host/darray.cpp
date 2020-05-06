#include "darray.h"
#include <cstdlib>
#include <cstring>


namespace tinker {
void device_memory_copyin_bytes(void* dst, const void* src, size_t nbytes,
                                DMFlag)
{
   std::memcpy(dst, src, nbytes);
}


void device_memory_copyout_bytes(void* dst, const void* src, size_t nbytes,
                                 DMFlag)
{
   std::memcpy(dst, src, nbytes);
}


void device_memory_copy_bytes(void* dst, const void* src, size_t nbytes, DMFlag)
{
   std::memcpy(dst, src, nbytes);
}


void device_memory_zero_bytes(void* dst, size_t nbytes, DMFlag)
{
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
