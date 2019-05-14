#ifndef TINKER_GPU_ACC_SWITCH_H_
#define TINKER_GPU_ACC_SWITCH_H_

#include "acc_mathfunc.h"
#include "decl_real.h"
#include "decl_switch.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
#pragma acc routine seq
template <int DO_DTAPER>
void switch_taper5(real rik, real cut, real off, real& taper, real& dtaper) {
  // S2(x) = 6 x**5 - 15 x**4 + 10 x**3
  // S2(x): [0,1] :-> [0,1]
  // taper5: [cut,off] :-> [1,0]
  real _1_ab = REAL_RECIP(cut - off);
  real x = (rik - off) * _1_ab;
  real x2 = x * x;
  real x3 = x2 * x;
  taper = x3 * (6 * x2 - 15 * x + 10);
  if_constexpr(DO_DTAPER) { dtaper = 30 * REAL_SQ(x * (1 - x)) * _1_ab; }
}
}
TINKER_NAMESPACE_END

#endif
