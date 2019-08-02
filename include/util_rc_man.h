#ifndef TINKER_UTIL_RC_MAN_H_
#define TINKER_UTIL_RC_MAN_H_

#include "util_cxx.h"

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * wrappers of intel @c for_rtl_init_, @c for_rtl_finish_ functions;
 * of gnu @c _gfortran_set_args function; or
 * of other fortran runtime functions.
 */
/// @{
void fortran_runtime_initialize(int, char**);
void fortran_runtime_finish();
/// @}

/**
 * @brief
 * set up and clean up host and device environment
 */
/// @{
void initialize();
void finish();
/// @}

/**
 * @brief
 * resource management
 *
 * To deallocate resource in reverse order of allocation, use named objects.
 * @code
 * rc_man foo42_{foo_data, op};
 * rc_man bar42_{bar_data, op};
 * @endcode
 *
 * To deallocate resource in the same order of allocation, use unnamed objects.
 * @code
 * rc_man{foo_data, op};
 * rc_man{bar_data, op};
 * @endcode
 */
class ResourceManagement {
public:
  typedef enum {
    dealloc = 0x001, ///< deallocation
    alloc = 0x002,   ///< allocation
    init = 0x004,    ///< initialization
    evolve = 0x008   ///< evolution
  } rc_op;

private:
  void (*f_)(rc_op);
  rc_op op_;
  bool will_dealloc_() const { return op_ & dealloc; }
  bool only_dealloc_() const { return op_ == dealloc; }

public:
  ResourceManagement(void (*f)(rc_op), rc_op op)
      : f_(f)
      , op_(op) {
    if (!will_dealloc_()) {
      f_(op_);
    }
  }

  ~ResourceManagement() {
    if (only_dealloc_()) {
      f_(op_);
    }
  }
};
typedef ResourceManagement rc_man;
typedef rc_man::rc_op rc_op;
constexpr rc_op rc_dealloc = rc_man::dealloc;
constexpr rc_op rc_alloc = rc_man::alloc;
constexpr rc_op rc_init = rc_man::init;

void host_data(rc_op);
void device_data(rc_op);
TINKER_NAMESPACE_END

#endif
