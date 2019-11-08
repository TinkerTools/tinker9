#pragma once
#include "builtin.h"


TINKER_NAMESPACE_BEGIN
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
#ifdef __CUDACC__
__host__
__device__
#endif
inline int ffsn(int val, int n)
{
   int ans = 0;
   int i = 0;
   while (val && (i++ < n)) {
      int pos = builtin::ffs(val);
      ans += pos;
      val >>= pos;
   }
   return ans;
}
TINKER_NAMESPACE_END
