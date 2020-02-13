#pragma once
#include "macro.h"
#include <type_traits>


TINKER_NAMESPACE_BEGIN
template <class E>
struct EnableEnumBitMask
{
   static constexpr bool value = false;
};


#define TINKER_ENABLE_ENUM_BITMASK(x)                                          \
   template <>                                                                 \
   struct EnableEnumBitMask<x>                                                 \
   {                                                                           \
      static constexpr bool value = true;                                      \
   };


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
TINKER_NAMESPACE_END
