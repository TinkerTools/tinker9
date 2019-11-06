#pragma once
#include "energy_buffer.h"
#if !TINKER_CUDART
#   error TINKER_CUDART must be true.
#endif


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
   atomicAdd(&buffer[offset],
             static_cast<unsigned long long>(static_cast<long long>(
                value * buffer_traits<float, 1>::fixed_point)));
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
