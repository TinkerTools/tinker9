#pragma once
#include "mathfunc_libfunc.h"
#include <cstring>
#include <type_traits>


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup atomic
 * \brief Add `value` to `buffer[offset]`.
 * \note `value` and `buffer` elements are of the same type.
 */
#pragma acc routine seq
template <class T>
void atomic_add(T value, T* buffer, size_t offset = 0)
{
   #pragma acc atomic update
   buffer[offset] += value;
}


/**
 * \ingroup atomic
 * \brief Add `value` to `buffer[offset]` via fixed-point arithmetic.
 * \tparam T Must be a floating point type.
 */
#pragma acc routine seq
template <class T>
void atomic_add(T value, unsigned long long* buffer, size_t offset = 0)
{
   static_assert(
      std::is_same<T, float>::value || std::is_same<T, double>::value, "");
   // float -> (signed) long long (int) -> unsigned long long (int)
   #pragma acc atomic update
   buffer[offset] += static_cast<unsigned long long>(
      static_cast<long long>(value * 0x100000000ull));
}


/**
 * \ingroup atomic
 * \brief Add virial to the virial buffer.
 * \note Virial and virial buffer are of the same type.
 */
#pragma acc routine seq
template <class T>
void atomic_add(T vxx, T vyx, T vzx, T vyy, T vzy, T vzz, T (*buffer)[8],
                size_t offset = 0)
{
   atomic_add(vxx, buffer[offset], 0);
   atomic_add(vyx, buffer[offset], 1);
   atomic_add(vzx, buffer[offset], 2);
   atomic_add(vyy, buffer[offset], 3);
   atomic_add(vzy, buffer[offset], 4);
   atomic_add(vzz, buffer[offset], 5);
}


/**
 * \ingroup atomic
 * \brief Add virial to the virial buffer via fixed-point arithmetic.
 * \tparam T Must be a floating point type.
 */
#pragma acc routine seq
template <class T>
void atomic_add(T vxx, T vyx, T vzx, T vyy, T vzy, T vzz,
                unsigned long long (*buffer)[8], size_t offset = 0)
{
   atomic_add(vxx, buffer[offset], 0);
   atomic_add(vyx, buffer[offset], 1);
   atomic_add(vzx, buffer[offset], 2);
   atomic_add(vyy, buffer[offset], 3);
   atomic_add(vzy, buffer[offset], 4);
   atomic_add(vzz, buffer[offset], 5);
}
TINKER_NAMESPACE_END
