#pragma once
#include "enum_op.h"


namespace tinker {
enum class rc_op
{
   DEALLOC = 0x001,
   ALLOC = 0x002,
   INIT = 0x004
};
TINKER_ENABLE_ENUM_BITMASK(rc_op);
constexpr rc_op rc_dealloc = rc_op::DEALLOC;
constexpr rc_op rc_alloc = rc_op::ALLOC;
constexpr rc_op rc_init = rc_op::INIT;


/**
 * \ingroup mem
 * \brief Resource management.
 *
 * To deallocate resource in reverse order of allocation, use named objects.
 * \code
 * rc_man foo42_{foo_data, op};
 * rc_man bar42_{bar_data, op};
 * \endcode
 *
 * To deallocate resource in the same order of allocation, use unnamed objects.
 * \code
 * rc_man{foo_data, op};
 * rc_man{bar_data, op};
 * \endcode
 */
class ResourceManagement
{
private:
   void (*f_)(rc_op);
   rc_op op_;
   bool will_dealloc_() const;
   bool only_dealloc_() const;


public:
   ResourceManagement(void (*f)(rc_op), rc_op op);
   ~ResourceManagement();
};
/// \ingroup mem
using rc_man = ResourceManagement;


/// \ingroup mem
/// \brief Set up host and device environment.
void initialize();
/// \ingroup mem
/// \brief Clean up host and device environment.
void finish();
/// \ingroup mem
/// \brief Set up and clean up host environment.
void host_data(rc_op);
/// \brief Set up and clean up device environment.
void device_data(rc_op);


/**
 * \ingroup mem
 * \brief
 * Wrappers of Intel `for_rtl_init_` function;
 * of GNU `_gfortran_set_args` function;
 * or of other Fortran runtime functions.
 */
void fortran_runtime_initialize(int, char**);
/**
 * \ingroup mem
 * \brief
 * Wrappers of Intel `for_rtl_finish_` function;
 * or of other Fortran runtime functions.
 */
void fortran_runtime_finish();
}
