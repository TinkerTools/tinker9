#pragma once
#include "macro.h"
#include <cstring>
#include <type_traits>


// https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
__device__
double atomicAdd(double* ptr, double v)
{
   unsigned long long int* ullptr = (unsigned long long int*)ptr;
   unsigned long long int old = *ullptr, assumed;
   do {
      assumed = old;
      old = atomicCAS(ullptr, assumed,
                      __double_as_longlong(v + __longlong_as_double(assumed)));
   } while (assumed != old);
   // using floating-point comparison will hang in case of NaN
   // (since NaN != NaN)
   return __longlong_as_double(old);
}
#endif


TINKER_NAMESPACE_BEGIN
template <class T>
__device__
void atomic_add_value(T value, T* buffer, size_t offset = 0)
{
   atomicAdd(&buffer[offset], value);
}


template <class T>
__device__
void atomic_add_value(T value, unsigned long long* buffer, size_t offset = 0)
{
   static_assert(
      std::is_same<T, float>::value || std::is_same<T, double>::value, "");
   // float -> (signed) long long (int) -> unsigned long long (int)
   atomicAdd(&buffer[offset],
             static_cast<unsigned long long>(
                static_cast<long long>(value * 0x100000000ull)));
}


template <class T>
__device__
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


template <class T>
__device__
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
