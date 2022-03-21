#pragma once
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

   static ThrustCache& instance();
   static void allocate();
   static void deallocate();

private:
   value_type* ptr;
   size_t nbytes;
};
}
