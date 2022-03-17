#pragma once
#include "macro.h"
#include <type_traits>

namespace tinker {
/// \ingroup rc
/// Direct mathematical calculation of enum class is prohibited in C++ syntax.
template <class E>
struct EnableEnumBitMask
{
   static constexpr bool value = false;
};

/// \def TINKER_ENABLE_ENUM_BITMASK
/// \ingroup rc
/// Explicitly enables mathematical calculation by casting enum class to integer.
#define TINKER_ENABLE_ENUM_BITMASK(x)                                                              \
   template <>                                                                                     \
   struct EnableEnumBitMask<x>                                                                     \
   {                                                                                               \
      static constexpr bool value = true;                                                          \
   }

template <class E>
constexpr typename std::enable_if<EnableEnumBitMask<E>::value, E>::type operator|(E lhs, E rhs)
{
   using ut = typename std::underlying_type<E>::type;
   return static_cast<E>(static_cast<ut>(lhs) | static_cast<ut>(rhs));
}

template <class E>
constexpr bool operator&(E lhs, E rhs)
{
   using ut = typename std::underlying_type<E>::type;
   return static_cast<bool>(static_cast<ut>(lhs) & static_cast<ut>(rhs));
}

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
