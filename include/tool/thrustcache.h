#pragma once
#include <cstddef>

namespace tinker {
/// \ingroup cuda_syntax
/// \brief Device memory cache for the Thrust Library.
class ThrustCache
{
public:
   // required
   using value_type = char;
   ThrustCache();
   value_type* allocate(ptrdiff_t);
   void deallocate(value_type*, size_t);
   void clear();

   // custom
   static ThrustCache& instance();
   static void allocate();
   static void deallocate();

private:
   value_type* ptr;
   size_t nbytes;
};
}
