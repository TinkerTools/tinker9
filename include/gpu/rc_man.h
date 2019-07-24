#ifndef TINKER_SRC_GPU_RC_MAN_H_
#define TINKER_SRC_GPU_RC_MAN_H_

#include "util/cxx.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
/// resource management operations
typedef enum {
  rc_dealloc = 0x001, /// deallocate device memory
  rc_alloc = 0x002,   /// allocate device memory
  rc_copyin = 0x004,  /// update device data from host memory, or directly
                      /// initialize device data
  rc_evolve = 0x010,
} rc_t;

/**
 * Resource Management
 * 
 * To deallocate resource in reverse order of allocation, use named objects.
 * @code
 * rc_man<foo_data> foo_random_name{rc};
 * rc_man<bar_data> bar_random_name{rc};
 * @endcode
 * 
 * To deallocate resource in the same order of allocation, use unnamed objects.
 * @code
 * rc_man<foo_data>{rc};
 * rc_man<bar_data>{rc};
 * @endcode
 */
template <void (*F)(rc_t)>
class rc_man {
private:
  rc_t rc_;
  bool dealloc_() const { return rc_ & rc_dealloc; }

public:
  rc_man(rc_t rc) : rc_(rc) {
    if (!dealloc_()) {
      F(rc_);
    }
  }

  ~rc_man() {
    if (dealloc_()) {
      F(rc_);
    }
  }
};
}
TINKER_NAMESPACE_END

#endif
