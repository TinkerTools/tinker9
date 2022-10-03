#pragma once
#include <cstddef>

namespace tinker {
/// \ingroup platform
/// \brief Device memory cache (allocator) for the Thrust Library.
///
/// The strategy of this allocator is as follows:
/// it only expands the size of the allocated memory,
/// and does not free the memory on calling \c deallocate(value_type*, size_t),
/// until \c clear() is called.
class ThrustCache
{
public:
   using value_type = char;
   ThrustCache();
   value_type* allocate(ptrdiff_t);      ///< Expands and returns the allocated memory upon request.
   void deallocate(value_type*, size_t); ///< Does not do anything.
   void clear();                         ///< Deallocates the memory.

   static ThrustCache& instance(); ///< Returns the reference to the singleton object.
   static void allocate();         ///< Allocates the memory managed by the singleton object.
   static void deallocate();       ///< Deallocates the memory managed by the singleton object.

private:
   value_type* ptr;
   size_t nbytes;
};
}
