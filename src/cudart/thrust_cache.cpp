#include "thrust_cache.h"
#include "darray.h"
#include "md.h"


namespace tinker {
ThrustCache::ThrustCache()
   : nbytes(0)
   , ptr(nullptr)
{}


auto ThrustCache::allocate(ptrdiff_t numbyte) -> value_type*
{
   if (numbyte > nbytes) {
      nbytes = numbyte;
      device_memory_deallocate_bytes(ptr);
      device_memory_allocate_bytes(reinterpret_cast<void**>(&ptr), nbytes);
   }
   return ptr;
}


void ThrustCache::deallocate(value_type*, size_t)
{
   // does not do anything
   return;
}


void ThrustCache::clear()
{
   device_memory_deallocate_bytes(ptr);
   ptr = nullptr;
   nbytes = 0;
}


ThrustCache thrust_cache;


void thrust_cache_dealloc()
{
   thrust_cache.clear();
}


void thrust_cache_alloc()
{
   const size_t numbyte = 10 * n * sizeof(real);
   thrust_cache.allocate(numbyte);
}
}
