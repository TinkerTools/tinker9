#pragma once
#include "ff/precision.h"
#include <cstddef>
#include <type_traits>

namespace tinker {
/// \ingroup acc_syntax
/// \brief Adds `value` to `buffer[offset]`.
/// \note `value` and `buffer` elements are of the same type.
#pragma acc routine seq
template <class T>
inline void atomic_add(T value, T* buffer, size_t offset = 0)
{
   #pragma acc atomic update
   buffer[offset] += value;
}

/// \ingroup acc_syntax
/// \brief Adds `value` to `buffer[offset]` via fixed-point arithmetic.
/// \tparam T  Must be of floating-point type.
#pragma acc routine seq
template <class T,
   class = typename std::enable_if<std::is_same<T, float>::value ||
      std::is_same<T, double>::value>::type>
inline void atomic_add(T value, fixed* buffer, size_t offset = 0)
{
   // float -> (signed) long long -> fixed
   #pragma acc atomic update
   buffer[offset] += static_cast<fixed>(static_cast<long long>(value * 0x100000000ull));
}

/// \ingroup acc_syntax
/// \brief Adds virial `{xx,yx,zx,yy,zy,zz}` to `buffer[offset][0 to 7]`.
#pragma acc routine seq
template <class T>
inline void atomic_add(T vxx, T vyx, T vzx, T vyy, T vzy, T vzz, T (*buffer)[8], size_t offset = 0)
{
   atomic_add(vxx, buffer[offset], 0);
   atomic_add(vyx, buffer[offset], 1);
   atomic_add(vzx, buffer[offset], 2);
   atomic_add(vyy, buffer[offset], 3);
   atomic_add(vzy, buffer[offset], 4);
   atomic_add(vzz, buffer[offset], 5);
}

/// \ingroup acc_syntax
/// \brief Adds virial `{xx,yx,zx,yy,zy,zz}` to `buffer[offset][0 to 7]` via fixed-point arithmetic.
#pragma acc routine seq
template <class T,
   class = typename std::enable_if<std::is_same<T, float>::value ||
      std::is_same<T, double>::value>::type>
inline void atomic_add(
   T vxx, T vyx, T vzx, T vyy, T vzy, T vzz, fixed (*buffer)[8], size_t offset = 0)
{
   atomic_add(vxx, buffer[offset], 0);
   atomic_add(vyx, buffer[offset], 1);
   atomic_add(vzx, buffer[offset], 2);
   atomic_add(vyy, buffer[offset], 3);
   atomic_add(vzy, buffer[offset], 4);
   atomic_add(vzz, buffer[offset], 5);
}

/// \ingroup acc_syntax
/// \brief Converts `val` of floating-point from its original type to type `G`.
#pragma acc routine seq
template <class G, class T>
inline G floatTo(T val)
{
   static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value, "");
   if CONSTEXPR (std::is_same<G, fixed>::value)
      return static_cast<G>(static_cast<long long>(val * 0x100000000ull));
   else
      return static_cast<G>(val);
}

/// \ingroup acc_syntax
/// \brief Converts `val` of fixed-point to floating-point.
#pragma acc routine seq
template <class T>
inline T fixedTo(fixed val)
{
   return static_cast<T>(static_cast<long long>(val)) / 0x100000000ull;
}

/// \ingroup acc_syntax
/// \brief Used as `eq<T1,T2>()` for two type identifiers.
template <class T, class U>
constexpr bool eq()
{
   return std::is_same<T, U>::value;
}
}
