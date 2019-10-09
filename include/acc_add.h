#pragma once


#include "macro.h"
#include <cstring>
#include <type_traits>


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup util
 * \brief
 * Add \c value to `buffer[offset]`.
 * \note
 * \c value and \c buffer elements are of the same type.
 */
#pragma acc routine seq
template <class T>
void atomic_add_value(T value, T* buffer, size_t offset = 0)
{
   #pragma acc atomic update
   buffer[offset] += value;
}


/**
 * \ingroup util
 * \brief
 * Add \c value to `buffer[offset]` via fixed-point arithmetic.
 * \tparam T
 * Must be a floating point type.
 */
#pragma acc routine seq
template <class T>
void atomic_add_value(T value, unsigned long long* buffer, size_t offset = 0)
{
   static_assert(
      std::is_same<T, float>::value || std::is_same<T, double>::value, "");
   // float -> (signed) long long (int) -> unsigned long long (int)
   #pragma acc atomic update
   buffer[offset] += static_cast<unsigned long long>(
      static_cast<long long>(value * fixed_point));
}


/**
 * \ingroup util
 * \brief
 * Add virial to the virial buffer.
 * \note
 * Virial and virial buffer are of the same type.
 */
#pragma acc routine seq
template <class T>
void atomic_add_value(T vxx, T vyx, T vzx, T vyy, T vzy, T vzz, T (*buffer)[8],
                      size_t offset = 0)
{
   atomic_add_value(vxx, buffer[offset], 0);
   atomic_add_value(vyx, buffer[offset], 1);
   atomic_add_value(vzx, buffer[offset], 2);
   atomic_add_value(vyy, buffer[offset], 3);
   atomic_add_value(vzy, buffer[offset], 4);
   atomic_add_value(vzz, buffer[offset], 5);
}


/**
 * \ingroup util
 * \brief
 * Add virial to the virial buffer via fixed-point arithmetic.
 */
#pragma acc routine seq
template <class T>
void atomic_add_value(T vxx, T vyx, T vzx, T vyy, T vzy, T vzz,
                      unsigned long long (*buffer)[8], size_t offset = 0)
{
   atomic_add_value(vxx, buffer[offset], 0);
   atomic_add_value(vyx, buffer[offset], 1);
   atomic_add_value(vzx, buffer[offset], 2);
   atomic_add_value(vyy, buffer[offset], 3);
   atomic_add_value(vzy, buffer[offset], 4);
   atomic_add_value(vzz, buffer[offset], 5);
}
TINKER_NAMESPACE_END
