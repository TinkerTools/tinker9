#include "md/inc.h"
#include "tool/darray.h"
#include "tool/thrustcache.h"

namespace tinker {
ThrustCache::ThrustCache()
   : nbytes(0)
   , ptr(nullptr)
{}

auto ThrustCache::allocate(ptrdiff_t numbyte) -> value_type*
{
   if ((size_t)numbyte > nbytes) {
      nbytes = numbyte;
      deviceMemoryDeallocate(ptr);
      deviceMemoryAllocateBytes(reinterpret_cast<void**>(&ptr), nbytes);
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
   deviceMemoryDeallocate(ptr);
   ptr = nullptr;
   nbytes = 0;
}

ThrustCache& ThrustCache::instance()
{
   static ThrustCache thrust_cache;
   return thrust_cache;
}

void ThrustCache::deallocate()
{
   ThrustCache::instance().clear();
}

void ThrustCache::allocate()
{
   const size_t numbyte = 10 * n * sizeof(real);
   ThrustCache::instance().allocate(numbyte);
}
}
