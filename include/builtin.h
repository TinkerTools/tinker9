#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
namespace builtin {
// int clz(int);
// int clz(long long);
// int ffs(int);
// int popc(int);


// CUDA
#ifdef __CUDACC__
__device__
__host__
inline int clz(int v)
{
#   ifdef __CUDA_ARCH__
   return __clz(v);
#   else
   return __builtin_clz(v);
#   endif
}


__device__
__host__
inline int clz(long long v)
{
#   ifdef __CUDA_ARCH__
   return __clzll(v);
#   else
   return __builtin_clzll(v);
#   endif
}


__device__
__host__
inline int ffs(int v)
{
#   ifdef __CUDA_ARCH__
   return __ffs(v);
#   else
   return __builtin_ffs(v);
#   endif
}


__device__
__host__
inline int popc(int v)
{
#   ifdef __CUDA_ARCH__
   return __popc(v);
#   else
   return __builtin_popcount(v);
#   endif
}
#else
// GCC
inline int clz(int v)
{
   return __builtin_clz(v);
}


inline int clz(long long v)
{
   return __builtin_clzll(v);
}


inline int ffs(int v)
{
   return __builtin_ffs(v);
}


inline int popc(int v)
{
   return __builtin_popcount(v);
}
#endif
}
TINKER_NAMESPACE_END
