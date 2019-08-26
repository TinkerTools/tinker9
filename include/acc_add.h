#ifndef TINKER_ACC_ADD_H_
#define TINKER_ACC_ADD_H_

#include "macro.h"
#include <type_traits>

TINKER_NAMESPACE_BEGIN
#pragma acc routine seq
template <class T>
void atomic_add_value(T* p, T e, int offset) {
  #pragma acc atomic update
  p[offset] += e;
}

#pragma acc routine seq
template <
    class T,
    class = typename std::enable_if<std::is_floating_point<T>::value>::type>
void atomic_add_value(unsigned long long* p, T e, int offset) {
  #pragma acc atomic update
  p[offset] += static_cast<unsigned long long>(
      static_cast<long long>(e * 0x100000000ull));
}
TINKER_NAMESPACE_END

#endif
