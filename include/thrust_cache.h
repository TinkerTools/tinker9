#pragma once
#include "macro.h"
#include <cstddef>

namespace tinker {
class ThrustCache
{
public:
   using value_type = char;
   ThrustCache();
   value_type* allocate(ptrdiff_t);
   void deallocate(value_type*, size_t);
   void clear();

private:
   value_type* ptr;
   size_t nbytes;
};
extern ThrustCache thrust_cache;

void thrust_cache_alloc();
void thrust_cache_dealloc();
}
