#pragma once
#include "macro.h"


namespace tinker {
/**
 * \ingroup cuda_syntax
 * If `b < 0`, returns `-abs(a)`, otherwise, returns `abs(a)`.
 * Similar to the 32-bit integer version of Fortran `SIGN(A,B)`.
 * \note Standard C and CUDA math libraries only have float and double versions.
 */
__device__
inline int int_copysign_shift(int a, int b)
{
   int mask = (a ^ b) >> 31;
   return (mask + a) ^ mask;
}


/**
 * \ingroup cuda_syntax
 * If `b < 0`, returns `-abs(a)`, otherwise, returns `abs(a)`.
 */
__device__
inline int int_copysign_if(int a, int b)
{
   int ans = abs(a);
   if (b < 0)
      ans = -ans;
   return ans;
}


/**
 * \def INT_COPYSIGN
 * \ingroup cuda_syntax
 * Defines the implementation of `int copysign(int,int);` function.
 */
#define INT_COPYSIGN int_copysign_shift
}
