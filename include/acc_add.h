#pragma once

#include "macro.h"
#include <cstring>
#include <type_traits>

TINKER_NAMESPACE_BEGIN
/// \brief
/// Add \c value to `buffer[offset]`.
/// \note
/// \c value and \c buffer elements are of the same type.
#pragma acc routine seq
template <class T>
void atomic_add_value (T value, T* buffer, size_t offset = 0)
{
   #pragma acc atomic update
   buffer[offset] += value;
}

/// \brief
/// Add \c value to `buffer[offset]` via fixed-point arithmetic.
/// \tparam T
/// Must be a floating point type.
#pragma acc routine seq
template <class T>
void atomic_add_value (T value, unsigned long long* buffer, size_t offset = 0)
{
   static_assert (
      std::is_same<T, float>::value || std::is_same<T, double>::value, "");
   // float -> (signed) long long (int) -> unsigned long long (int)
   #pragma acc atomic update
   buffer[offset] += static_cast<unsigned long long> (
      static_cast<long long> (value * fixed_point));
}

#pragma acc routine seq
template <class T>
void atomic_add_value (T vxx, T vyx, T vzx, T vyy, T vzy, T vzz, T* buffer,
                       size_t offset = 0)
{
   size_t offv = offset * 8;
   atomic_add_value (vxx, buffer, offv + 0);
   atomic_add_value (vyx, buffer, offv + 1);
   atomic_add_value (vzx, buffer, offv + 2);
   atomic_add_value (vyy, buffer, offv + 3);
   atomic_add_value (vzy, buffer, offv + 4);
   atomic_add_value (vzz, buffer, offv + 5);
}

#pragma acc routine seq
template <class T>
void atomic_add_value (T vxx, T vyx, T vzx, T vyy, T vzy, T vzz,
                       unsigned long long* buffer, size_t offset = 0)
{
   size_t offv = offset * 8;
   atomic_add_value (vxx, buffer, offv + 0);
   atomic_add_value (vyx, buffer, offv + 1);
   atomic_add_value (vzx, buffer, offv + 2);
   atomic_add_value (vyy, buffer, offv + 3);
   atomic_add_value (vzy, buffer, offv + 4);
   atomic_add_value (vzz, buffer, offv + 5);
}
TINKER_NAMESPACE_END
