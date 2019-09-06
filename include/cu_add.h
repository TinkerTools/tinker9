#ifndef TINKER_CU_ADD_H_
#define TINKER_CU_ADD_H_

#include "macro.h"
#include <type_traits>

#if __CUDA_ARCH__ < 600
inline __device__ double //
atomicAdd(double* address, double val) {
  unsigned long long* address_as_ull = (unsigned long long*)address;
  unsigned long long old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif

TINKER_NAMESPACE_BEGIN
template <class T>
__device__ void //
atomic_add_value(T value, T* buffer, int offset = 0) {
  atomicAdd(&buffer[offset], value);
}

template <class T>
__device__ void //
atomic_add_value(T value, unsigned long long* buffer, int offset = 0) {
  static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value,
                "");
  atomicAdd(&buffer[offset],
            static_cast<unsigned long long>(
                static_cast<long long>(value * fixed_point)));
}
TINKER_NAMESPACE_END

#endif
