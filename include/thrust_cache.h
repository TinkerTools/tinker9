#pragma once
#include "macro.h"
#include <cstddef>


TINKER_NAMESPACE_BEGIN
class ThrustCache
{
public:
   using value_type = char;
   ThrustCache();
   value_type* allocate(ptrdiff_t);
   void deallocate(value_type*, size_t);
   void clear();


private:
   size_t nbytes;
   value_type* ptr;
};
TINKER_EXTERN ThrustCache thrust_cache;
TINKER_NAMESPACE_END
