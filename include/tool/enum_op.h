#pragma once
#include "macro.h"
#include <type_traits>


namespace tinker {
/**
 * \ingroup rc
 * Direct mathematical calculation of enum class is prohibited in C++ syntax.
 */
template <class E>
struct EnableEnumBitMask
{
   static constexpr bool value = false;
};


/**
 * \def TINKER_ENABLE_ENUM_BITMASK
 * \ingroup rc
 * Explicitly enables mathematical calculation by casting enum class to integer.
 */
#define TINKER_ENABLE_ENUM_BITMASK(x)                                          \
   template <>                                                                 \
   struct EnableEnumBitMask<x>                                                 \
   {                                                                           \
      static constexpr bool value = true;                                      \
   }


template <class E>
constexpr typename std::enable_if<EnableEnumBitMask<E>::value, E>::type
operator|(E lhs, E rhs)
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
}
