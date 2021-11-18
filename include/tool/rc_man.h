#pragma once
#include "tool/enum_op.h"


namespace tinker {
enum class ResourceOperation
{
   DEALLOC = 0x001,
   ALLOC = 0x002,
   INIT = 0x004
};
TINKER_ENABLE_ENUM_BITMASK(ResourceOperation);
using rc_op = ResourceOperation;
/**
 * \ingroup rc
 * Deallocates resource.
 */
constexpr rc_op rc_dealloc = rc_op::DEALLOC;
/**
 * \ingroup rc
 * Allocates resource.
 */
constexpr rc_op rc_alloc = rc_op::ALLOC;
/**
 * \ingroup rc
 * Initializes resource.
 */
constexpr rc_op rc_init = rc_op::INIT;


/**
 * \ingroup rc
 * Resource management. Allocates resource in the object constructor and
 * deallocates resource in the object destructor.
 *
 * To deallocate resource in reverse order of allocation, use named objects.
 * \code
 * rc_man foo42{foo_data, op};
 * rc_man bar42{bar_data, op};
 * \endcode
 *
 * To deallocate resource in the same order of allocation, use unnamed objects.
 * \code
 * rc_man {foo_data, op};
 * rc_man {bar_data, op};
 * \endcode
 */
class ResourceManagement
{
private:
   void (*m_f)(rc_op);
   rc_op m_op;
   bool will_dealloc() const;
   bool only_dealloc() const;


public:
   /**
    * \param f   Function to (de)allocate and/or initialize resource.
    * \param op  Resource operation flag.
    */
   ResourceManagement(void (*f)(rc_op), rc_op op);
   ~ResourceManagement();
};
/**
 * \ingroup rc
 * Type alias.
 */
using rc_man = ResourceManagement;


/**
 * \ingroup rc
 * Sets up host and device environment.
 */
void initialize();
/**
 * \ingroup rc
 * Cleans up host and device environment.
 */
void finish();
/**
 * \ingroup rc
 * Set up and clean up device environment.
 */
void device_data(rc_op);
}
