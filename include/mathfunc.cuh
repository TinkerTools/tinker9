#pragma once
#include "macro.h"


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
 *      to 1 in `val`, returns an undefined value.
 */
__device__
inline int ffsn_loop(int val, int n)
{
   int ans = 0;
   int i = 0;
   while (val && (i++ < n)) {
      int pos = __ffs(val);
      ans += pos;
      val >>= pos;
   }
   return ans;
}


__device__
inline int ffsn(int c1, int n)
{
   int c2 = c1 - ((c1 >> 1) & 0x55555555);
   int c4 = ((c2 >> 2) & 0x33333333) + (c2 & 0x33333333);
   int c8 = ((c4 >> 4) + c4) & 0x0f0f0f0f;
   int c16 = ((c8 >> 8) + c8);
   // int c32 = ((c16 >> 16) + c16) & 0x3f; // __popc
   int r = 0, i = n - 1, t;
   // clang-format off
   t = (c16    ) & 0x1f; if (i >= t) { r += 16; i -= t; }
   t = (c8 >> r) & 0x0f; if (i >= t) { r +=  8; i -= t; }
   t = (c4 >> r) & 0x07; if (i >= t) { r +=  4; i -= t; }
   t = (c2 >> r) & 0x03; if (i >= t) { r +=  2; i -= t; }
   t = (c1 >> r) & 0x01; if (i >= t) { r +=  1;         }
   if (n == 0)                         r = -1;
   // if (n > c32)                     r = 32;
   // clang-format on
   return r + 1;
}
TINKER_NAMESPACE_END
