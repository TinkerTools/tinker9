#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup math
 * \brief
 * If `b < 0`, return `-abs(a)`, otherwise, return `abs(a)`.
 *
 * Reverse-engineered from the 32-bit integer version of Fortran `SIGN(A,B)`.
 * Standard C and CUDA math libraries only have float and double versions.
 */
__device__
inline int int_copysign(int a, int b)
{
   int mask = (a ^ b) >> 31;
   return (mask + a) ^ mask;
}


__device__
inline int int_copysign_if(int a, int b)
{
   int ans = abs(a);
   if (b < 0)
      ans = -ans;
   return ans;
}
TINKER_NAMESPACE_END
