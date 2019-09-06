#ifndef TINKER_CU_SWITCH_H_
#define TINKER_CU_SWITCH_H_

#include "mathfunc.h"

TINKER_NAMESPACE_BEGIN
template <int DO_DTAPER>
inline __device__ void //
switch_taper5(real rik, real cut, real off, real& __restrict__ taper,
              real& __restrict__ dtaper) {
  real _1_ab = REAL_RECIP(cut - off);
  real x = (rik - off) * _1_ab;
  real x2 = x * x;
  real x3 = x2 * x;
  taper = x3 * (6 * x2 - 15 * x + 10);
  if_constexpr(DO_DTAPER) { dtaper = 30 * REAL_SQ(x * (1 - x)) * _1_ab; }
}
TINKER_NAMESPACE_END

#endif
