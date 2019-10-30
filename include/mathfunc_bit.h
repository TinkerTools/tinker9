#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
#define FFSN_CODE_(val_, n_, ffs_func_)                                        \
   unsigned int v = val_;                                                      \
   int ans = 0;                                                                \
   int i = 0;                                                                  \
   while (v && (i++ < n_)) {                                                   \
      int pos = ffs_func_(v);                                                  \
      ans += pos;                                                              \
      v >>= pos;                                                               \
   }                                                                           \
   return ans


/**
 * \ingroup math
 * \brief
 * Find the position of the n-th least significant bit set in a 32-bit integer.
 *
 * \param val A 32-bit integer.
 * \param n Ranges from 1 to 32.
 * \return
 * A value from 0 to 32.
 *    - If `val` equals 0, always returns 0.
 *    - If `n` equals 0, always returns 0.
 *    - If `n` is greater than `POPC`, which is the number of bits that are set
 *      to 1 in `val`, returns as if `n` equals `POPC`.
 */
inline int builtin_ffsn(int val, int n)
{
   FFSN_CODE_(val, n, __builtin_ffs);
}


#ifdef __CUDACC__
__device__
inline int ffsn(int val, int n)
{
   FFSN_CODE_(val, n, __ffs);
}
#endif
#undef FFSN_CODE_
TINKER_NAMESPACE_END
