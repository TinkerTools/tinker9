#ifndef TINKER_GPU_ACC_H_
#define TINKER_GPU_ACC_H_

//======================================================================
// acc math functions

#include "decl_mathfunc.h"

#pragma acc routine(abs) seq

#pragma acc routine(sqrt) seq
#pragma acc routine(exp) seq
#pragma acc routine(floor) seq
#pragma acc routine(fabs) seq
#pragma acc routine(pow) seq
#pragma acc routine(cos) seq
#pragma acc routine(sin) seq
#pragma acc routine(acos) seq
#pragma acc routine(erf) seq
#pragma acc routine(erfc) seq
#pragma acc routine(fmin) seq
#pragma acc routine(fmax) seq

#pragma acc routine(sqrtf) seq
#pragma acc routine(expf) seq
#pragma acc routine(floorf) seq
#pragma acc routine(fabsf) seq
#pragma acc routine(powf) seq
#pragma acc routine(cosf) seq
#pragma acc routine(sinf) seq
#pragma acc routine(acosf) seq
#pragma acc routine(erff) seq
#pragma acc routine(erfcf) seq
#pragma acc routine(fminf) seq
#pragma acc routine(fmaxf) seq

#include "decl_box.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
//======================================================================
// switch

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

//======================================================================
// image and imagen

/**
 * applys periodic boundary conditions to displacement (xr, yr, zr) and
 * preserves the correct signs
 */
#pragma acc routine seq
void image(real& __restrict__ xr, real& __restrict__ yr, real& __restrict__ zr,
           const box_t* __restrict__ pb);

/**
 * applys periodic boundary conditions to displacement (xr, yr, zr) but only
 * guarantee the lengths are correct
 */
#pragma acc routine seq
void imagen(real& __restrict__ xr, real& __restrict__ yr, real& __restrict__ zr,
            const box_t* __restrict__ pb);
}
TINKER_NAMESPACE_END

#endif
