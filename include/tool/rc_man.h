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
using RcOp = ResourceOperation;

/// \ingroup rc
/// Deallocates resource.
constexpr RcOp rc_dealloc = RcOp::DEALLOC;

/// \ingroup rc
/// Allocates resource.
constexpr RcOp rc_alloc = RcOp::ALLOC;

/// \ingroup rc
/// Initializes resource.
constexpr RcOp rc_init = RcOp::INIT;

/// \ingroup rc
/// Resource management. Allocates resource in the object constructor and
/// deallocates resource in the object destructor.
///
/// To deallocate resource in reverse order of allocation, use named objects.
/// \code
/// RcMan foo42{fooData, op};
/// RcMan bar42{barData, op};
/// \endcode
///
/// To deallocate resource in the same order of allocation, use unnamed objects.
/// \code
/// RcMan {fooData, op};
/// RcMan {barData, op};
/// \endcode
class ResourceManagement
{
private:
   void (*m_f)(RcOp);
   RcOp m_op;
   bool will_dealloc() const;
   bool only_dealloc() const;

public:
   /// \param f   Function to (de)allocate and/or initialize resource.
   /// \param op  Resource operation flag.
   ResourceManagement(void (*f)(RcOp), RcOp op);
   ~ResourceManagement();
};

/// \ingroup rc
/// Type alias.
using RcMan = ResourceManagement;

/// \ingroup rc
/// Sets up host and device environment.
void initialize();

/// \ingroup rc
/// Cleans up host and device environment.
void finish();

/// \ingroup rc
/// Set up and clean up device environment.
void device_data(RcOp);
}
