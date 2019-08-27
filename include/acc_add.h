#ifndef TINKER_ACC_ADD_H_
#define TINKER_ACC_ADD_H_

#include "macro.h"
#include <type_traits>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * add @c value to @c buffer[@c offset] atomically
 */
#pragma acc routine seq
template <class T>
inline void atomic_add_value(T value, T* buffer, int offset = 0) {
  #pragma acc atomic update
  buffer[offset] += value;
}

/**
 * @brief
 * add @c value to @c buffer[@c offset] atomically via fixed-point arithmetic
 *
 * @tparam T
 * must be a floating point type
 */
#pragma acc routine seq
template <
    class T,
    class = typename std::enable_if<std::is_floating_point<T>::value>::type>
inline void atomic_add_value(T value, unsigned long long* buffer,
                             int offset = 0) {
  #pragma acc atomic update
  buffer[offset] += static_cast<unsigned long long>(
      static_cast<long long>(value * TINKER_FIXED_POINT));
}
TINKER_NAMESPACE_END

#endif
